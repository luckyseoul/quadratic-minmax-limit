#!/usr/bin/env bash
# Non-destructive, verified full-repository backup to Nick's large drive.
set -euo pipefail
umask 077

if [[ $# -ne 2 ]]; then
    printf 'Usage: bash milestone_backup.sh ABSOLUTE_REPOSITORY_PATH MILESTONE_LABEL\n' >&2
    exit 2
fi
source_repo=$(realpath -e -- "$1")
milestone_label=$2
if [[ ! $milestone_label =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$ ]]; then
    printf 'Milestone label must be a short, safe filename component.\n' >&2
    exit 2
fi
case "$source_repo" in
    /|/home|/home/nick|/tmp|/mnt|/mnt/storage|/mnt/storage/*)
        printf 'Refusing a broad source or a source inside the backup drive.\n' >&2
        exit 2
        ;;
esac
repo_top=$(git -C "$source_repo" rev-parse --show-toplevel)
if [[ $(realpath -e -- "$repo_top") != "$source_repo" ]]; then
    printf 'Source must be the exact Git working-tree root.\n' >&2
    exit 2
fi
git_dir=$(realpath -e -- "$(git -C "$source_repo" rev-parse --absolute-git-dir)")
git_common_dir=$(realpath -e -- "$(git -C "$source_repo" rev-parse --path-format=absolute --git-common-dir)")
for required_tool in git tar zstd sha256sum findmnt df du jq mktemp cmp; do
    command -v "$required_tool" >/dev/null
done
mount_target=$(findmnt -n -T /mnt/storage -o TARGET)
mount_options=$(findmnt -n -T /mnt/storage -o OPTIONS)
if [[ $mount_target != /mnt/storage || ,$mount_options, != *,rw,* ]]; then
    printf 'The real /mnt/storage mount must be present and writable; do not remount automatically.\n' >&2
    exit 3
fi
source_kib=$(du -sk -- "$source_repo" | awk '{print $1}')
git_dir_kib=$(du -sk -- "$git_dir" | awk '{print $1}')
git_common_kib=$(du -sk -- "$git_common_dir" | awk '{print $1}')
free_kib=$(df -Pk /mnt/storage | awk 'NR==2 {print $4}')
if (( free_kib < 2 * (source_kib + git_dir_kib + git_common_kib) + 1048576 )); then
    printf 'Insufficient free space for a full archive plus history bundle and verification.\n' >&2
    exit 3
fi
snapshot_head=$(git -C "$source_repo" rev-parse HEAD)
snapshot_branch=$(git -C "$source_repo" symbolic-ref --quiet --short HEAD || true)
snapshot_utc=$(date -u +%Y%m%dT%H%M%SZ)
repo_name=$(basename -- "$source_repo")
repo_key=$(printf '%s' "$source_repo" | sha256sum | cut -c1-12)
backup_parent="/mnt/storage/backups/codex/${repo_name}-${repo_key}"
for backup_component in /mnt/storage/backups /mnt/storage/backups/codex "$backup_parent"; do
    if [[ -L $backup_component ]]; then
        printf 'Refusing a symlinked backup destination component: %s\n' "$backup_component" >&2
        exit 3
    fi
    mkdir -p -- "$backup_component"
done
backup_parent=$(realpath -e -- "$backup_parent")
if [[ $backup_parent != /mnt/storage/backups/codex/* ||
      $(findmnt -n -T "$backup_parent" -o TARGET) != /mnt/storage ]]; then
    printf 'Backup parent escaped the verified large-drive mount.\n' >&2
    exit 3
fi
snapshot_dir=$(mktemp -d "$backup_parent/${snapshot_utc}-${snapshot_head:0:12}-${milestone_label}.XXXXXX")
if [[ $(realpath -e -- "$snapshot_dir") != "$snapshot_dir" ||
      $(findmnt -n -T "$snapshot_dir" -o TARGET) != /mnt/storage ]]; then
    printf 'Fresh destination is not on the verified large-drive mount: %s\n' "$snapshot_dir" >&2
    exit 3
fi
printf 'Backup directory: %s\n' "$snapshot_dir"
printf 'Partial artifacts are retained at that exact path if any check fails.\n'

git -C "$source_repo" status --porcelain=v1 -z > "$snapshot_dir/status.before.z"
git -C "$source_repo" show-ref > "$snapshot_dir/refs.before.txt"
git -C "$source_repo" diff --binary > "$snapshot_dir/unstaged.patch"
git -C "$source_repo" diff --cached --binary > "$snapshot_dir/staged.patch"
git -C "$source_repo" bundle create "$snapshot_dir/repository.bundle" --all HEAD
tar -I 'zstd -T2 -3' -cf "$snapshot_dir/full-worktree.tar.zst" -C "$source_repo" .
# The separate metadata archives also preserve linked/split/sparse/conflicted
# indexes, worktree metadata, configuration, hooks, reflogs and unreachable objects.
tar -I 'zstd -T2 -3' -cf "$snapshot_dir/git-common-dir.tar.zst" -C "$git_common_dir" .
tar -I 'zstd -T2 -3' -cf "$snapshot_dir/git-dir.tar.zst" -C "$git_dir" .

# Preserve the global instruction file and the exact implementation too.
if [[ -f /home/nick/.codex/AGENTS.md ]]; then
    cp -- /home/nick/.codex/AGENTS.md "$snapshot_dir/global-AGENTS.md"
fi
cp -- "${BASH_SOURCE[0]}" "$snapshot_dir/backup-tool.sh"

git -C "$source_repo" bundle verify "$snapshot_dir/repository.bundle"
zstd -t -- "$snapshot_dir/full-worktree.tar.zst"
zstd -t -- "$snapshot_dir/git-common-dir.tar.zst" "$snapshot_dir/git-dir.tar.zst"
tar -I zstd --compare -f "$snapshot_dir/full-worktree.tar.zst" -C "$source_repo"
tar -I zstd --compare -f "$snapshot_dir/git-common-dir.tar.zst" -C "$git_common_dir"
tar -I zstd --compare -f "$snapshot_dir/git-dir.tar.zst" -C "$git_dir"
git -C "$source_repo" status --porcelain=v1 -z > "$snapshot_dir/status.after.z"
git -C "$source_repo" show-ref > "$snapshot_dir/refs.after.txt"
cmp -- "$snapshot_dir/status.before.z" "$snapshot_dir/status.after.z"
cmp -- "$snapshot_dir/refs.before.txt" "$snapshot_dir/refs.after.txt"
[[ $(git -C "$source_repo" rev-parse HEAD) == "$snapshot_head" ]]

jq -n --arg source "$source_repo" --arg head "$snapshot_head" \
    --arg branch "$snapshot_branch" --arg utc "$snapshot_utc" \
    --arg label "$milestone_label" --arg destination "$snapshot_dir" \
    --arg git_dir "$git_dir" --arg git_common_dir "$git_common_dir" \
    --arg mount_source "$(findmnt -n -T /mnt/storage -o SOURCE)" \
    '{schema:"codex_milestone_backup_v1",source:$source,head:$head,branch:$branch,
      utc:$utc,label:$label,destination:$destination,mount_source:$mount_source,
      git_dir:$git_dir,git_common_dir:$git_common_dir,
      scope:"Entire working tree including ignored/untracked files; complete Git dir and common-dir archives; all-ref Git bundle and staged/unstaged patches",
      full_archive_compared_to_source:true,bundle_verified:true,source_status_and_refs_unchanged:true,
      deletions:0,complete:true}' > "$snapshot_dir/receipt.json"
(
    cd -- "$snapshot_dir"
    checksum_inputs=(*)
    sha256sum -- "${checksum_inputs[@]}" > SHA256SUMS
    sha256sum -c SHA256SUMS
)
printf 'VERIFIED_BACKUP=%s\n' "$snapshot_dir"
sha256sum -- "$snapshot_dir/SHA256SUMS" "$snapshot_dir/receipt.json"

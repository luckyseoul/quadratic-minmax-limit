"""Route-neutral proof status for the original MathOverflow 413935 question.

The production registry is deliberately empty: no reviewed completion proof
is currently registered.  Admission is an explicit source-code/review change,
not a public ``proved=True`` argument, an environment setting, or the result of
AND-ing optional route diagnostics.  Each future entry must identify the
all-orders theorem and a separate review, with both artifact bytes pinned.

Hash checks establish provenance, not mathematical correctness.  The named
theorem and review must establish the stated conclusion for the actual
sequence before an entry is admitted.  Tests may replace the private registry
boundary with clearly synthetic fixtures; those are never production proofs.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from hashlib import sha256
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
PROBLEM = "MathOverflow 413935: the actual all-orders sequence alpha_n"


class Conclusion(str, Enum):
    EXISTENCE = "existence"
    VALUE = "value"
    NONEXISTENCE = "nonexistence"


@dataclass(frozen=True)
class ReviewedCompletionProof:
    """A source-admitted theorem/review pair, not a computational proof flag.

    ``formal_proof_path``/``formal_proof_sha256`` are an optional third
    pinned artifact for a machine-checked formalization (e.g. Lean 4 +
    Mathlib) of ``theorem``, following the provenance pattern OpenAI used
    for its Navier-Stokes/Euler blowup proofs: a human-reviewed writeup
    plus a hash-pinned formal proof file, independently re-checkable by a
    second kernel (their "Comparator" tool against an alternate trusted
    kernel, not just the primary Lean elaborator). Supplying a formal proof
    here is optional and does not by itself change how an entry is
    admitted; the theorem/review pair remains the required core, and a
    stale or unhashed formal-proof artifact still fails the whole entry
    closed, exactly like a stale theorem or review artifact.

    ``formal_verification_note`` is free text recording how the formal
    artifact was checked (toolchain, axioms it depends on, whether an
    independent kernel re-verified it) -- provenance only, not a
    machine-parsed correctness flag.
    """

    proof_id: str
    conclusion: Conclusion
    theorem: str
    theorem_path: str
    theorem_sha256: str
    review_path: str
    review_sha256: str
    limit_value: str | None = None
    problem: str = PROBLEM
    formal_proof_path: str | None = None
    formal_proof_sha256: str | None = None
    formal_verification_note: str | None = None


# Add an entry only after review of a complete proof for the original problem.
# Local propositions, finite certificates, and conditional reductions do not
# belong here.  No route (Paley, amplification, or otherwise) is mandatory.
_REVIEWED_COMPLETION_PROOFS: tuple[ReviewedCompletionProof, ...] = ()


def _reviewed_completion_entries() -> tuple[ReviewedCompletionProof, ...]:
    """Auditable registry boundary; no caller-supplied acceptance parameters."""
    return _REVIEWED_COMPLETION_PROOFS


def _pinned_artifact(relative_path: str, digest: str) -> Path:
    if not isinstance(relative_path, str) or not relative_path.strip():
        raise ValueError("a nonempty repository-relative artifact path is required")
    path = Path(relative_path)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError("proof artifacts must stay inside the repository")
    if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise ValueError("proof artifacts require an exact lowercase SHA-256")
    resolved = (ROOT / path).resolve(strict=True)
    if not resolved.is_relative_to(ROOT.resolve()) or not resolved.is_file():
        raise ValueError("proof artifacts must be repository files")
    if sha256(resolved.read_bytes()).hexdigest() != digest:
        raise ValueError(f"stale proof artifact: {relative_path}")
    return resolved


def _validated_completion_entries() -> tuple[ReviewedCompletionProof, ...]:
    entries = _reviewed_completion_entries()
    if not isinstance(entries, tuple):
        raise TypeError("the reviewed completion registry must be an explicit tuple")
    identifiers: set[str] = set()
    for entry in entries:
        if not isinstance(entry, ReviewedCompletionProof):
            raise TypeError("only reviewed theorem/review records can enter the registry")
        if not isinstance(entry.proof_id, str) or not entry.proof_id.strip():
            raise ValueError("each reviewed proof needs a nonempty identifier")
        if entry.proof_id in identifiers:
            raise ValueError(f"duplicate reviewed proof identifier: {entry.proof_id}")
        identifiers.add(entry.proof_id)
        if entry.problem != PROBLEM:
            raise ValueError("a completion proof must address the actual all-orders problem")
        if not isinstance(entry.conclusion, Conclusion):
            raise ValueError("the proof conclusion must be existence, value, or nonexistence")
        if not isinstance(entry.theorem, str) or not entry.theorem.strip():
            raise ValueError("each reviewed proof needs an explicit theorem identifier")
        if entry.conclusion is Conclusion.VALUE:
            if not isinstance(entry.limit_value, str) or not entry.limit_value.strip():
                raise ValueError("an identified-value proof must state its exact limit value")
        elif entry.limit_value is not None:
            raise ValueError("only an identified-value proof may supply a limit value")
        theorem_path = _pinned_artifact(entry.theorem_path, entry.theorem_sha256)
        review_path = _pinned_artifact(entry.review_path, entry.review_sha256)
        if theorem_path == review_path:
            raise ValueError("the proof and its review must be separate artifacts")
        has_formal_path = entry.formal_proof_path is not None
        has_formal_hash = entry.formal_proof_sha256 is not None
        if has_formal_path != has_formal_hash:
            raise ValueError(
                "a formal proof artifact needs both a path and a matching SHA-256; "
                "supply both formal_proof_path and formal_proof_sha256, or neither"
            )
        if has_formal_path:
            formal_path = _pinned_artifact(entry.formal_proof_path, entry.formal_proof_sha256)
            if formal_path in (theorem_path, review_path):
                raise ValueError(
                    "the formal proof artifact must be separate from the theorem and review"
                )
        elif entry.formal_verification_note is not None:
            raise ValueError(
                "a formal_verification_note requires a pinned formal_proof_path/sha256; "
                "provenance notes may not stand in for a hashed artifact"
            )

    conclusions = {entry.conclusion for entry in entries}
    if Conclusion.NONEXISTENCE in conclusions and conclusions.intersection(
        {Conclusion.EXISTENCE, Conclusion.VALUE}
    ):
        raise ValueError("conflicting existence and nonexistence proof records")
    values = {entry.limit_value for entry in entries if entry.conclusion is Conclusion.VALUE}
    if len(values) > 1:
        raise ValueError("conflicting exact limit values in reviewed proof records")
    return entries


def original_mo_status() -> dict:
    """Report reviewed global conclusions independently of every optional route.

    An invalid/stale/conflicting registry fails closed and reports its error;
    it cannot be converted into a global conclusion by other theorem flags.
    """
    errors: list[str] = []
    try:
        entries = _validated_completion_entries()
    except (OSError, TypeError, ValueError) as exc:
        entries = ()
        errors.append(str(exc))
    conclusions = {entry.conclusion for entry in entries}
    existence = bool(conclusions.intersection({Conclusion.EXISTENCE, Conclusion.VALUE}))
    nonexistence = Conclusion.NONEXISTENCE in conclusions
    values = [entry.limit_value for entry in entries if entry.conclusion is Conclusion.VALUE]
    value = values[0] if values else None
    settled = existence or nonexistence
    if nonexistence:
        result = "NONEXISTENCE_PROVED"
    elif value is not None:
        result = "VALUE_PROVED"
    elif existence:
        result = "EXISTENCE_PROVED"
    else:
        result = "OPEN"
    return {
        "problem": PROBLEM,
        "status": result,
        "limit_status": "NONEXISTENT" if nonexistence else "CLOSED" if existence else "OPEN",
        "problem_settled": settled,
        "existence_proved": existence,
        "nonexistence_proved": nonexistence,
        "value_proved": value is not None,
        "limit_value": value,
        "registry_valid": not errors,
        "registry_errors": errors,
        "reviewed_completion_proofs": [asdict(entry) for entry in entries],
        "required_optional_routes": [],
        "rule": (
            "Only explicit reviewed all-orders proof entries determine global status. "
            "Optional route diagnostics neither certify nor veto completion. "
            "Artifact hashes verify provenance, not mathematics."
        ),
    }

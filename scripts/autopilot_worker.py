#!/usr/bin/env python3
"""Local-model worker for the research autopilot (see autopilot/README.md).

Invoked by research_autopilot.py as AUTOPILOT_WORKER_CMD. It reads the task
JSON at argv[1], runs a tool-calling agent loop against a local Ollama model,
and writes a schema-valid result to AUTOPILOT_RESULT_JSON. The model never
gets raw shell access: it can only read repository files and run a fixed,
allow-listed set of existing scripts, and it must finish by calling
submit_result with the exact fields the controller requires.
"""
from __future__ import annotations
import json, os, subprocess, sys, urllib.request

OLLAMA_HOST = os.environ.get("OLLAMA_HOST", "http://localhost:11434")
MODEL = os.environ.get("AUTOPILOT_MODEL", "qwen2.5:14b-instruct")
REPO = os.environ.get("AUTOPILOT_REPO", ".")
MAX_TURNS = int(os.environ.get("AUTOPILOT_MAX_TURNS", "20"))

REQUIRED_RESULT_KEYS = (
    "status", "claim", "implication", "evidence", "scope",
    "dependencies", "attack_targets", "next_tasks", "compute",
)
VALID_STATUS = {"success", "counterexample", "partial", "dead_end", "error"}
VALID_SCOPE = {"finite", "asymptotic", "conditional", "all-orders"}

# Existing repo scripts a worker is allowed to invoke, not arbitrary shell.
ALLOWED_COMMANDS = {
    "run_tests": ["python3", "-m", "pytest", "-q"],
    "milestone_backup": ["bash", "scripts/milestone_backup.sh"],
}

TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "read_file",
            "description": "Read a UTF-8 text file from the repository, relative to the repo root.",
            "parameters": {
                "type": "object",
                "properties": {"path": {"type": "string"}},
                "required": ["path"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "list_dir",
            "description": "List entries in a repository directory, relative to the repo root.",
            "parameters": {
                "type": "object",
                "properties": {"path": {"type": "string"}},
                "required": ["path"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_command",
            "description": (
                "Run one allow-listed repository command: "
                + ", ".join(sorted(ALLOWED_COMMANDS))
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "name": {"type": "string", "enum": sorted(ALLOWED_COMMANDS)},
                    "args": {"type": "array", "items": {"type": "string"}},
                },
                "required": ["name"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "submit_result",
            "description": "Finish the task by submitting the final result. Call this exactly once, last.",
            "parameters": {
                "type": "object",
                "properties": {
                    "status": {"type": "string", "enum": sorted(VALID_STATUS)},
                    "claim": {"type": "string"},
                    "implication": {"type": "string"},
                    "evidence": {"type": "array", "items": {"type": "string"}},
                    "scope": {"type": "string", "enum": sorted(VALID_SCOPE)},
                    "dependencies": {"type": "array", "items": {"type": "string"}},
                    "attack_targets": {"type": "array", "items": {"type": "string"}},
                    "next_tasks": {"type": "array", "items": {"type": "string"}},
                    "compute": {"type": "object"},
                },
                "required": list(REQUIRED_RESULT_KEYS),
            },
        },
    },
]


def repo_path(rel: str) -> str:
    full = os.path.realpath(os.path.join(REPO, rel))
    root = os.path.realpath(REPO)
    if not (full == root or full.startswith(root + os.sep)):
        raise ValueError("path escapes repository root")
    return full


def call_read_file(args):
    path = repo_path(args["path"])
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()[:20000]


def call_list_dir(args):
    path = repo_path(args.get("path", "."))
    return "\n".join(sorted(os.listdir(path)))


def call_run_command(args):
    name = args["name"]
    if name not in ALLOWED_COMMANDS:
        raise ValueError(f"command not allow-listed: {name}")
    cmd = ALLOWED_COMMANDS[name] + list(args.get("args", []))
    proc = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=600)
    return json.dumps({"returncode": proc.returncode, "stdout": proc.stdout[-8000:], "stderr": proc.stderr[-4000:]})


DISPATCH = {"read_file": call_read_file, "list_dir": call_list_dir, "run_command": call_run_command}


def ollama_chat(messages):
    body = json.dumps({"model": MODEL, "messages": messages, "tools": TOOLS, "stream": False}).encode()
    req = urllib.request.Request(f"{OLLAMA_HOST}/api/chat", data=body, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        return json.loads(resp.read())


def validate_result(result: dict):
    missing = [k for k in REQUIRED_RESULT_KEYS if k not in result]
    if missing:
        raise ValueError(f"submit_result missing keys: {missing}")
    if result["status"] not in VALID_STATUS:
        raise ValueError(f"invalid status: {result['status']!r}")
    if result["scope"] not in VALID_SCOPE:
        raise ValueError(f"invalid scope: {result['scope']!r}")


def main():
    task_path = sys.argv[1]
    result_path = os.environ["AUTOPILOT_RESULT_JSON"]
    task = json.loads(open(task_path).read())

    system = (
        "You are a research-autopilot worker. Do the exact objective given. "
        "Read the listed files first via read_file. Use run_command only for "
        "the allow-listed commands. Finite computation is evidence, never an "
        "all-orders proof. When finished, call submit_result exactly once "
        "with every required field filled in accurately; do not fabricate "
        "evidence paths."
    )
    messages = [
        {"role": "system", "content": system},
        {"role": "user", "content": json.dumps(task, indent=2)},
    ]

    for _ in range(MAX_TURNS):
        response = ollama_chat(messages)
        msg = response["message"]
        messages.append(msg)
        tool_calls = msg.get("tool_calls") or []
        if not tool_calls:
            continue
        for call in tool_calls:
            fn = call["function"]["name"]
            args = call["function"].get("arguments") or {}
            if isinstance(args, str):
                args = json.loads(args)
            if fn == "submit_result":
                validate_result(args)
                with open(result_path, "w") as f:
                    json.dump(args, f, indent=2)
                return
            try:
                output = DISPATCH[fn](args)
            except Exception as e:
                output = json.dumps({"error": str(e)})
            messages.append({"role": "tool", "content": str(output)})

    with open(result_path, "w") as f:
        json.dump({
            "status": "error", "claim": "worker did not submit a result within turn budget",
            "implication": "none", "evidence": [], "scope": "finite",
            "dependencies": [], "attack_targets": [], "next_tasks": [],
            "compute": {"turns": MAX_TURNS},
        }, f, indent=2)


if __name__ == "__main__":
    main()

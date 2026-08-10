from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import re
import shlex
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from glob import glob

try:
    from rich import box
    from rich.console import Console
    from rich.live import Live
    from rich.panel import Panel
    from rich.progress import BarColumn, Progress, TextColumn
    from rich.table import Table
    from rich.text import Text
    from rich.console import Group

    RICH_AVAILABLE = True
except ImportError:
    box = None
    Console = object
    Live = None
    Panel = None
    Progress = None
    Table = None
    Text = None
    Group = None
    RICH_AVAILABLE = False

STATUS_FILENAME = "job_status.json"
AUTOMATION_WORKDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUTOMATION_DEPENDENCIES = {
    "checkoutputs.py": [
        "BTVNanoCommissioning.utils.xrootdtools",
        "alive_progress",
        "rich",
    ],
    "haddoutputs.py": [
        "alive_progress",
        "rich",
        "uproot",
    ],
}
ICONS = {
    "ok": "✅",
    "warn": "⚠️",
    "error": "❌",
    "info": "ℹ️",
    "combine": "🧩",
    "submit": "🚀",
    "check": "🔎",
    "resubmit": "🔁",
    "hadd": "📦",
    "done": "🎉",
    "clock": "🕒",
}
CONDOR_EVENT_RE = re.compile(r"^(\d{3}) \((\d+)\.(\d+)\.\d+\)")
CONDOR_CLUSTER_RE = re.compile(r"(?:job|hadd)\.(?:log|out|err)_(\d+)(?:-\d+)?$")
CONDOR_SUBMIT_RE = re.compile(
    r"\b(?P<count>\d+)\s+job\(s\)\s+submitted\s+to\s+cluster\s+" r"(?P<cluster>\d+)\b",
    re.IGNORECASE,
)
AUTOMATION_MAX_SUBMIT_RETRIES = 20
AUTOMATION_SUBMIT_RETRY_SHORT_ATTEMPTS = 5
AUTOMATION_SUBMIT_RETRY_SHORT_DELAY = 30.0
AUTOMATION_SUBMIT_RETRY_LONG_DELAY = 300.0
AUTOMATION_OUTPUT_LIMIT = 1200


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def parse_condor_submit_output(output: str) -> dict | None:
    match = CONDOR_SUBMIT_RE.search(output or "")
    if not match:
        return None
    return {
        "count": int(match.group("count")),
        "cluster_id": int(match.group("cluster")),
    }


def _status_path(job_dir: str) -> str:
    return os.path.join(job_dir, STATUS_FILENAME)


def _stage_entry(stage: str, summary: str = "", detail: str = "") -> dict:
    return {
        "stage": stage,
        "summary": summary,
        "detail": detail,
        "timestamp": _now_iso(),
    }


def _default_status(job_dir: str, job_name: str = "", output_dir: str = "") -> dict:
    return {
        "job_dir": os.path.abspath(job_dir) if job_dir else "",
        "job_name": (
            job_name or os.path.basename(os.path.abspath(job_dir)) if job_dir else ""
        ),
        "output_dir": os.path.abspath(output_dir) if output_dir else "",
        "created_at": "",
        "updated_at": "",
        "submitted_jobs": None,
        "submission_cluster_id": None,
        "checked_outputs": {
            "status": "",
            "missing_outputs": None,
            "checked_jobs": None,
            "job_ids": [],
        },
        "resubmitted_jobs": {
            "count": None,
            "job_ids": [],
            "pending_reset": False,
            "cluster_id": None,
        },
        "hadd_jobs_submitted": None,
        "hadd_submission_cluster_id": None,
        "hadd_jobs_checked": {
            "status": "",
            "failed": None,
            "checked_jobs": None,
        },
        "all_done": False,
        "automation": {
            "step": "waiting",
            "last_action_at": "",
            "last_action": "",
            "last_returncode": None,
            "last_stdout": "",
            "last_stderr": "",
            "active_pid": None,
            "active_command": "",
            "warnings": [],
            "blocked_reason": "",
            "blocked_state": "",
            "processed_job_clusters": [],
            "processed_hadd_clusters": [],
            "unreconciled_submission": None,
        },
        "stages_passed": [],
    }


def _merge_status(data: dict, job_dir: str) -> dict:
    base = _default_status(job_dir)
    if not isinstance(data, dict):
        return base
    for key, value in base.items():
        if key not in data:
            data[key] = value
    if not isinstance(data.get("checked_outputs"), dict):
        data["checked_outputs"] = base["checked_outputs"]
    else:
        for key, value in base["checked_outputs"].items():
            data["checked_outputs"].setdefault(key, value)
    if not isinstance(data.get("resubmitted_jobs"), dict):
        data["resubmitted_jobs"] = base["resubmitted_jobs"]
    else:
        for key, value in base["resubmitted_jobs"].items():
            data["resubmitted_jobs"].setdefault(key, value)
    if not isinstance(data.get("hadd_jobs_checked"), dict):
        data["hadd_jobs_checked"] = base["hadd_jobs_checked"]
    else:
        for key, value in base["hadd_jobs_checked"].items():
            data["hadd_jobs_checked"].setdefault(key, value)
    if not isinstance(data.get("automation"), dict):
        data["automation"] = base["automation"]
    else:
        for key, value in base["automation"].items():
            data["automation"].setdefault(key, value)
        for key in ("warnings", "processed_job_clusters", "processed_hadd_clusters"):
            if not isinstance(data["automation"].get(key), list):
                data["automation"][key] = []
    if not isinstance(data.get("stages_passed"), list):
        data["stages_passed"] = []
    return data


def load_status(job_dir: str, create: bool = True) -> dict:
    path = _status_path(job_dir)
    if os.path.exists(path):
        try:
            with open(path) as f:
                data = json.load(f)
        except json.JSONDecodeError:
            return _default_status(job_dir)
        data = _merge_status(data, job_dir)
        if _normalize_status_metadata(job_dir, data):
            save_status(job_dir, data)
        return data
    if not create:
        data = _default_status(job_dir)
        _normalize_status_metadata(job_dir, data)
        return data
    data = _default_status(job_dir)
    _normalize_status_metadata(job_dir, data)
    save_status(job_dir, data)
    return data


def save_status(job_dir: str, data: dict) -> None:
    data["job_dir"] = os.path.abspath(job_dir)
    data["updated_at"] = _now_iso()
    os.makedirs(job_dir, exist_ok=True)
    path = _status_path(job_dir)
    tmp_path = f"{path}.tmp"
    with open(tmp_path, "w") as f:
        json.dump(data, f, indent=4, sort_keys=False)
    os.replace(tmp_path, path)


def ensure_status(job_dir: str, job_name: str = "", output_dir: str = "") -> dict:
    data = load_status(job_dir, create=True)
    if job_name and not data.get("job_name"):
        data["job_name"] = job_name
    if output_dir and not data.get("output_dir"):
        data["output_dir"] = os.path.abspath(output_dir)
    if not data.get("created_at"):
        data["created_at"] = _now_iso()
    save_status(job_dir, data)
    return data


def _normalize_status_metadata(job_dir: str, data: dict) -> bool:
    changed = False

    if job_dir:
        job_name = os.path.basename(os.path.abspath(job_dir))
        if job_name and not data.get("job_name"):
            data["job_name"] = job_name
            changed = True

    args = _arguments(job_dir) if job_dir else {}
    output_dir = data.get("output_dir") or args.get("outputDir", "")
    if output_dir:
        normalized_output_dir = os.path.abspath(str(output_dir))
        if data.get("output_dir") != normalized_output_dir:
            data["output_dir"] = normalized_output_dir
            changed = True

    jobnum_path = os.path.join(job_dir, "jobnum_list.txt") if job_dir else ""
    if jobnum_path and os.path.isfile(jobnum_path):
        try:
            with open(jobnum_path) as f:
                job_ids = [line.strip() for line in f if line.strip()]
        except OSError:
            job_ids = []
        if job_ids and data.get("submitted_jobs") is None:
            data["submitted_jobs"] = len(job_ids)
            changed = True
        checked_outputs = data.get("checked_outputs")
        if (
            isinstance(checked_outputs, dict)
            and checked_outputs.get("checked_jobs") is None
            and job_ids
        ):
            checked_outputs["checked_jobs"] = len(job_ids)
            changed = True

    return changed


def _append_stage(data: dict, stage: str, summary: str = "", detail: str = "") -> None:
    data.setdefault("stages_passed", [])
    data["stages_passed"].append(_stage_entry(stage, summary=summary, detail=detail))


def _trim_output(value: str, limit: int = AUTOMATION_OUTPUT_LIMIT) -> str:
    value = value or ""
    if len(value) <= limit:
        return value
    return value[-limit:]


def _automation_retry_delay(attempt: int, expect_condor_submit: bool) -> float:
    if not expect_condor_submit:
        return 0.0
    if attempt < AUTOMATION_SUBMIT_RETRY_SHORT_ATTEMPTS:
        return AUTOMATION_SUBMIT_RETRY_SHORT_DELAY
    return AUTOMATION_SUBMIT_RETRY_LONG_DELAY


def _automation_dependency_modules() -> list[str]:
    modules: list[str] = []
    for dependency_list in AUTOMATION_DEPENDENCIES.values():
        for module in dependency_list:
            if module not in modules:
                modules.append(module)
    return modules


def _validate_automation_environment() -> None:
    if not AUTOMATION_DEPENDENCIES:
        return

    probe = (
        "import importlib, json, sys\n"
        f"modules = {json.dumps(_automation_dependency_modules())}\n"
        "missing = []\n"
        "for module in modules:\n"
        "    try:\n"
        "        importlib.import_module(module)\n"
        "    except Exception as exc:\n"
        '        missing.append({"module": module, "error": f"{type(exc).__name__}: {exc}"})\n'
        "print(json.dumps(missing))\n"
        "sys.exit(1 if missing else 0)\n"
    )
    completed = subprocess.run(
        [sys.executable, "-c", probe],
        cwd=AUTOMATION_WORKDIR,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode == 0:
        return

    missing: list[dict[str, str]] = []
    output = (completed.stdout or "").strip()
    if output:
        try:
            missing = json.loads(output)
        except json.JSONDecodeError:
            missing = []

    if not missing:
        stderr = (completed.stderr or "").strip()
        details = stderr or "unknown import error"
    else:
        details = "\n".join(
            f"- {item.get('module', 'unknown')}: {item.get('error', 'unknown error')}"
            for item in missing
        )

    message = (
        "Automation preflight failed: one or more modules required by subjobs "
        "could not be imported from the dashboard working directory "
        f"({AUTOMATION_WORKDIR}).\n"
        f"{details}\n"
        "Source the BTVNanoCommissioning environment before starting live automation."
    )
    raise RuntimeError(message)


def _blocked_automation_action(data: dict) -> str | None:
    auto = _automation(data)
    reason = str(auto.get("blocked_reason") or "")
    last_action = str(auto.get("last_action") or "")

    if "All checked jobs failed" in reason or last_action == "checkoutputs":
        return "checkoutputs"
    if "no arrays_* root outputs were found" in reason:
        return "hadd_check"
    if last_action in {"hadd_submit", "hadd_check"}:
        return last_action
    return None


def _current_output_scan(job_dir: str, data: dict) -> tuple[int, int, list[str]] | None:
    output_dir = _output_dir(job_dir, data)
    job_list_path = os.path.join(job_dir, "jobnum_list.txt")
    if not output_dir or not os.path.isfile(job_list_path):
        return None

    try:
        with open(job_list_path) as f:
            job_ids = [line.strip() for line in f if line.strip()]
    except OSError:
        return None

    missing: list[str] = []
    for job_id in job_ids:
        outfile = os.path.join(output_dir, f"hists_{job_id}", f"hists_{job_id}.coffea")
        if not os.path.isfile(outfile):
            missing.append(job_id)
            continue
        for rootfile in glob(
            os.path.join(output_dir, f"arrays_hists_{job_id}", "*", "*", "*.root")
        ):
            try:
                if os.path.getsize(rootfile) < 10 * 1024:
                    missing.append(job_id)
                    break
            except OSError:
                missing.append(job_id)
                break

    return len(job_ids), len(missing), missing


def _resume_blocked_automation(job_dir: str, data: dict) -> bool:
    auto = _automation(data)
    if auto.get("step") != "blocked" or not auto.get("blocked_reason"):
        return False
    if auto.get("unreconciled_submission"):
        return False

    if data.get("all_done"):
        auto["step"] = "done"
        auto["blocked_reason"] = ""
        auto["blocked_state"] = ""
        auto["last_action_at"] = _now_iso()
        auto["active_pid"] = None
        auto["active_command"] = ""
        _append_stage(data, "Automation resumed", "Cleared stale blocked state", "done")
        save_status(job_dir, data)
        return True

    action = _blocked_automation_action(data)
    if action is None:
        return False

    reason = str(auto.get("blocked_reason") or "")
    if action == "checkoutputs" and "All checked jobs failed" in reason:
        current = _current_output_scan(job_dir, data)
        if current is None:
            return False
        checked_jobs, missing_jobs, _ = current
        if checked_jobs == 0 or missing_jobs >= checked_jobs:
            return False
    if action == "hadd_check" and "no arrays_* root outputs were found" in reason:
        output_dir = _output_dir(job_dir, data)
        if not _has_array_outputs(output_dir):
            return False

    auto["step"] = "waiting"
    auto["blocked_reason"] = ""
    auto["blocked_state"] = ""
    auto["last_action_at"] = _now_iso()
    auto["active_pid"] = None
    auto["active_command"] = ""
    if action == "checkoutputs":
        auto["processed_job_clusters"] = []
    elif action == "hadd_submit":
        auto["processed_hadd_clusters"] = []
        data["hadd_jobs_submitted"] = None
        data["all_done"] = False
    elif action == "hadd_check":
        auto["processed_hadd_clusters"] = []
        data["all_done"] = False

    _append_stage(
        data,
        "Automation resumed",
        "Retrying blocked automation",
        auto.get("last_action", ""),
    )
    save_status(job_dir, data)
    return True


def _resume_stale_checkoutputs(
    job_dir: str, data: dict, condor_status: dict | None
) -> bool:
    auto = _automation(data)
    if (
        auto.get("step") != "waiting"
        or str(auto.get("last_action") or "") != "checkoutputs"
    ):
        return False

    checked = data.get("checked_outputs", {})
    missing_outputs = checked.get("missing_outputs")
    if not isinstance(missing_outputs, int) or missing_outputs <= 0:
        return False

    if data.get("all_done"):
        return False

    condor_state = (condor_status or {}).get("state")
    condor_cluster = (condor_status or {}).get("cluster")
    stale_processed_cluster = (
        _cluster_processed(data, "processed_job_clusters", condor_cluster)
        and _condor_queue_empty(condor_cluster) is True
    )
    if (
        condor_state not in (None, "unsubmitted")
        and condor_cluster is not None
        and not stale_processed_cluster
    ):
        return False

    _append_stage(
        data,
        "Automation resumed",
        "Retrying stale checkoutputs with no live condor cluster",
        auto.get("last_action", ""),
    )
    save_status(job_dir, data)
    return _run_automation_command(
        job_dir,
        [sys.executable, _script_path("checkoutputs.py"), job_dir, "--condor"],
        "checkoutputs",
        expect_condor_submit=True,
    )


def resume_blocked_automation_jobs(root: str) -> int:
    resumed = 0
    for job_dir in _job_dirs(root):
        data = load_status(job_dir, create=False)
        if _resume_blocked_automation(job_dir, data):
            resumed += 1
    return resumed


def _automation(data: dict) -> dict:
    auto = data.setdefault("automation", {})
    base = _default_status(data.get("job_dir", ""))["automation"]
    for key, value in base.items():
        auto.setdefault(key, value)
    for key in ("warnings", "processed_job_clusters", "processed_hadd_clusters"):
        if not isinstance(auto.get(key), list):
            auto[key] = []
    return auto


def _automation_active_pid(data: dict) -> int | None:
    auto = _automation(data)
    pid = auto.get("active_pid")
    if isinstance(pid, int) and pid > 0:
        return pid
    try:
        pid_int = int(pid)
    except Exception:
        return None
    return pid_int if pid_int > 0 else None


def _process_alive(pid: int, expected_fragment: str = "") -> bool:
    try:
        completed = subprocess.run(
            ["ps", "-p", str(pid), "-o", "args="],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError:
        return False
    if completed.returncode != 0:
        return False
    args = (completed.stdout or "").strip()
    if not args:
        return False
    if expected_fragment and expected_fragment not in args:
        return False
    return True


def _find_matching_automation_pid(job_dir: str, script_name: str) -> int | None:
    try:
        completed = subprocess.run(
            ["ps", "-eo", "pid=,args="],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError:
        return None
    if completed.returncode != 0:
        return None
    job_dir_abs = os.path.abspath(job_dir)
    for line in (completed.stdout or "").splitlines():
        line = line.strip()
        if not line:
            continue
        pid_text, _, args = line.partition(" ")
        try:
            pid = int(pid_text)
        except ValueError:
            continue
        args = args.strip()
        if script_name in args and job_dir_abs in args:
            return pid
    return None


def _automation_subprocess_alive(data: dict, expected_fragment: str = "") -> bool:
    pid = _automation_active_pid(data)
    if pid is None:
        return False
    return _process_alive(pid, expected_fragment=expected_fragment)


def _automation_step_is_active(action: str, step: str) -> bool:
    if action == "checkoutputs":
        return step == "checking_outputs" or step.startswith("checkoutputs_attempt_")
    if action in {"hadd_submit", "hadd_check"}:
        return step in {"submitting_hadd", "checking_hadd"} or step.startswith(
            f"{action}_attempt_"
        )
    return False


def _set_automation_action(
    data: dict,
    action: str,
    *,
    step: str | None = None,
    returncode: int | None = None,
    stdout: str = "",
    stderr: str = "",
    active_pid: int | None = None,
    active_command: str | None = None,
) -> None:
    auto = _automation(data)
    if step:
        auto["step"] = step
    auto["last_action"] = action
    auto["last_action_at"] = _now_iso()
    auto["last_returncode"] = returncode
    auto["last_stdout"] = _trim_output(stdout)
    auto["last_stderr"] = _trim_output(stderr)
    auto["active_pid"] = active_pid
    if active_command is not None:
        auto["active_command"] = active_command


def _add_automation_warning(data: dict, message: str) -> None:
    auto = _automation(data)
    warnings = auto.setdefault("warnings", [])
    if message not in warnings:
        warnings.append(message)
    auto["last_action_at"] = _now_iso()


def record_automation_warning(job_dir: str, message: str) -> dict:
    data = load_status(job_dir)
    _add_automation_warning(data, message)
    save_status(job_dir, data)
    return data


def record_automation_blocked(job_dir: str, reason: str) -> dict:
    data = load_status(job_dir)
    auto = _automation(data)
    auto["step"] = "blocked"
    auto["blocked_reason"] = reason
    auto["blocked_state"] = ""
    auto["active_pid"] = None
    auto["active_command"] = ""
    _add_automation_warning(data, reason)
    _append_stage(data, "Automation blocked", reason, "")
    save_status(job_dir, data)
    return data


def record_processing_complete(job_dir: str, summary: str = "All done") -> dict:
    data = load_status(job_dir)
    data["all_done"] = True
    auto = _automation(data)
    auto["step"] = "done"
    auto["blocked_reason"] = ""
    auto["blocked_state"] = ""
    auto["active_pid"] = None
    auto["active_command"] = ""
    _append_stage(data, "All done", summary, "")
    save_status(job_dir, data)
    return data


def record_submission(
    job_dir: str,
    total_jobs: int,
    job_name: str = "",
    output_dir: str = "",
    cluster_id: int | None = None,
) -> dict:
    data = ensure_status(job_dir, job_name=job_name, output_dir=output_dir)
    data["submitted_jobs"] = int(total_jobs)
    data["submission_cluster_id"] = cluster_id
    detail = f"cluster {cluster_id}" if cluster_id is not None else ""
    _append_stage(data, "Submitted", f"{total_jobs} jobs", detail)
    save_status(job_dir, data)
    return data


def record_checked_outputs(
    job_dir: str,
    checked_jobs: int,
    missing_jobs: int,
    job_ids: list | None = None,
) -> dict:
    data = load_status(job_dir)
    data["checked_outputs"] = {
        "status": (
            "All jobs done" if missing_jobs == 0 else f"{missing_jobs} outputs missing"
        ),
        "missing_outputs": int(missing_jobs),
        "checked_jobs": int(checked_jobs),
        "job_ids": list(job_ids or []),
    }
    data["resubmitted_jobs"] = {
        "count": None,
        "job_ids": [],
        "pending_reset": False,
        "cluster_id": None,
    }
    if missing_jobs == 0:
        _append_stage(
            data, "Checked outputs", "All jobs done", f"{checked_jobs} jobs checked"
        )
    else:
        _append_stage(
            data,
            "Checked outputs",
            f"{missing_jobs} outputs missing",
            f"{checked_jobs} jobs checked",
        )
    save_status(job_dir, data)
    return data


def record_resubmitted_jobs(
    job_dir: str, job_ids: list, cluster_id: int | None = None
) -> dict:
    data = load_status(job_dir)
    count = len(job_ids)
    data["resubmitted_jobs"] = {
        "count": count,
        "job_ids": [str(job_id) for job_id in job_ids],
        "pending_reset": True,
        "cluster_id": cluster_id,
    }
    detail = f"cluster {cluster_id}" if cluster_id is not None else ""
    _append_stage(data, "Jobs resubmitted", f"{count} jobs resubmitted", detail)
    save_status(job_dir, data)
    return data


def record_hadd_submitted(
    job_dir: str, n_jobs: int, cluster_id: int | None = None
) -> dict:
    data = load_status(job_dir)
    data["hadd_jobs_submitted"] = int(n_jobs)
    data["hadd_submission_cluster_id"] = cluster_id
    detail = f"cluster {cluster_id}" if cluster_id is not None else ""
    _append_stage(data, "Hadd jobs submitted", f"{n_jobs} jobs submitted", detail)
    save_status(job_dir, data)
    return data


def record_hadd_checked(job_dir: str, checked_jobs: int, failed_jobs: int) -> dict:
    data = load_status(job_dir)
    data["hadd_jobs_checked"] = {
        "status": "All done" if failed_jobs == 0 else f"{failed_jobs} failed",
        "failed": int(failed_jobs),
        "checked_jobs": int(checked_jobs),
    }
    if failed_jobs == 0:
        data["all_done"] = True
        _append_stage(
            data,
            "Hadd jobs checked",
            "All done",
            f"{checked_jobs} jobs checked",
        )
        _append_stage(data, "All done", "All done", "")
    else:
        data["all_done"] = False
        _append_stage(
            data,
            "Hadd jobs checked",
            f"{failed_jobs} failed",
            f"{checked_jobs} jobs checked",
        )
    save_status(job_dir, data)
    return data


def find_job_dir_by_output(
    output_dir: str, search_root: str | None = None
) -> str | None:
    target = os.path.abspath(output_dir)
    module_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    roots = [search_root] if search_root else [os.getcwd(), module_root]
    roots = [r for r in roots if r]
    for root in roots:
        for job_dir in sorted(
            glob(os.path.join(root, "jobs_*")) + glob(os.path.join(root, "job_*"))
        ):
            args_path = os.path.join(job_dir, "arguments.json")
            if not os.path.isfile(args_path):
                continue
            try:
                with open(args_path) as f:
                    args = json.load(f)
            except Exception:
                continue
            if os.path.abspath(str(args.get("outputDir", ""))) == target:
                return job_dir
    return None


def _job_dirs(root: str) -> list[str]:
    dirs = []
    module_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    search_roots = [root, module_root] if root else [module_root]
    for search_root in search_roots:
        for pattern in ("jobs_*", "job_*"):
            dirs.extend(glob(os.path.join(search_root, pattern)))
    return sorted({os.path.abspath(d) for d in dirs if os.path.isdir(d)})


def _fmt_value(value, default=""):
    if value is None:
        return default
    return str(value)


def _state_icon(state: str) -> str:
    return {
        "done": ICONS["ok"],
        "combine": ICONS["combine"],
        "warn": ICONS["warn"],
        "error": ICONS["error"],
        "info": ICONS["info"],
    }.get(state, ICONS["info"])


def _state_style(state: str) -> str:
    return {
        "done": "bold green",
        "combine": "bold magenta",
        "warn": "bold yellow",
        "error": "bold red",
        "info": "bold cyan",
    }.get(state, "bold cyan")


def _job_overall_state(data: dict) -> str:
    if data.get("all_done"):
        return "done"
    if _automation(data).get("step") == "blocked":
        return "error"
    if data.get("hadd_jobs_submitted") is not None:
        return "combine"
    if data.get("checked_outputs", {}).get("missing_outputs") not in (None, 0):
        return "warn"
    if data.get("submitted_jobs") is not None:
        return "info"
    return "info"


def _job_title(data: dict) -> str:
    name = data.get("job_name") or os.path.basename(data.get("job_dir", "")) or "job"
    if str(name).startswith("jobs_"):
        name = str(name)[5:]
    return str(name)


def _job_progress_rank(
    data: dict, condor: dict | None, hadd_condor: dict | None
) -> tuple[int, int, str]:
    title = _job_title(data).lower()

    if data.get("all_done"):
        return (50, 0, title)

    auto = _automation(data)
    checked = data.get("checked_outputs", {})
    checked_status = str(checked.get("status") or "")
    missing_outputs = checked.get("missing_outputs")
    hadd_checked = data.get("hadd_jobs_checked", {})
    hadd_checked_status = str(hadd_checked.get("status") or "")
    hadd_submitted = data.get("hadd_jobs_submitted")
    submitted_jobs = data.get("submitted_jobs")
    blocked_reason = str(auto.get("blocked_reason") or "")
    blocked_action = _blocked_automation_action(data)

    if (
        data.get("submitted_jobs") is None
        and not checked_status
        and hadd_submitted is None
    ):
        return (0, 0, title)

    if checked_status == "All jobs done" or missing_outputs == 0:
        if _is_array_job(data.get("job_dir", ""), data) and hadd_submitted is None:
            return (3, 0, title)
        if (
            _is_array_job(data.get("job_dir", ""), data)
            and hadd_checked_status != "All done"
        ):
            return (4, 0, title)
        return (2, 0, title)

    if checked_status or isinstance(missing_outputs, int):
        if isinstance(missing_outputs, int) and missing_outputs > 0:
            if isinstance(submitted_jobs, int) and missing_outputs >= submitted_jobs:
                return (1, 0, title)
            return (2, 0, title)
        return (1, 0, title)

    condor_state = (condor or {}).get("state")
    if condor_state in ("running", "warn"):
        return (1, 0, title)
    if condor_state == "unsubmitted":
        return (0, 1, title)

    if blocked_reason:
        if (
            blocked_action == "checkoutputs"
            or "All checked jobs failed" in blocked_reason
        ):
            return (2, 1, title)
        if (
            blocked_action in {"hadd_submit", "hadd_check"}
            or "hadd" in blocked_reason.lower()
        ):
            return (4, 1, title)
        return (2, 1, title)

    if hadd_submitted is not None:
        if hadd_checked_status == "All done":
            return (50, 0, title)
        return (4, 0, title)

    if submitted_jobs is not None:
        return (1, 0, title)

    if hadd_condor and hadd_condor.get("state") in ("running", "warn"):
        return (4, 1, title)

    return (1, 9, title)


def _job_stage_label(stage: str) -> tuple[str, str]:
    mapping = {
        "Submitted": (ICONS["submit"], "green"),
        "Checked outputs": (ICONS["check"], "cyan"),
        "Jobs resubmitted": (ICONS["resubmit"], "yellow"),
        "Hadd jobs submitted": (ICONS["hadd"], "magenta"),
        "Hadd jobs checked": (ICONS["hadd"], "blue"),
        "Automation blocked": (ICONS["warn"], "red"),
        "All done": (ICONS["done"], "green"),
    }
    return mapping.get(stage, (ICONS["info"], "white"))


def _stage_summary(data: dict) -> str:
    pieces = []
    submitted = data.get("submitted_jobs")
    if submitted is not None:
        pieces.append(f"{ICONS['submit']} Submitted {submitted}")
    checked = data.get("checked_outputs", {})
    if checked.get("status"):
        pieces.append(f"{ICONS['check']} Checked: {checked['status']}")
    resubmitted = data.get("resubmitted_jobs", {})
    if resubmitted.get("count"):
        pieces.append(f"{ICONS['resubmit']} Resubmitted {resubmitted['count']}")
    hadd_sub = data.get("hadd_jobs_submitted")
    if hadd_sub is not None:
        pieces.append(f"{ICONS['hadd']} Hadd submitted {hadd_sub}")
    hadd_checked = data.get("hadd_jobs_checked", {})
    if hadd_checked.get("status"):
        pieces.append(f"{ICONS['hadd']} Hadd checked: {hadd_checked['status']}")
    if data.get("all_done"):
        pieces.append(f"{ICONS['done']} All done")
    auto = _automation(data)
    if auto.get("blocked_reason"):
        pieces.append(f"{ICONS['warn']} Blocked: {auto['blocked_reason']}")
    return (
        f" {ICONS['clock']} " + "  →  ".join(pieces)
        if pieces
        else f"{ICONS['clock']} No status yet"
    )


def _latest_condor_cluster(log_dir: str, prefix: str = "job") -> int | None:
    if not os.path.isdir(log_dir):
        return None
    clusters = []
    for path in glob(os.path.join(log_dir, f"{prefix}.*_*")):
        match = CONDOR_CLUSTER_RE.search(os.path.basename(path))
        if match:
            clusters.append(int(match.group(1)))
    return max(clusters) if clusters else None


def read_condor_status(
    job_dir: str,
    data: dict | None = None,
    *,
    prefix: str = "job",
    log_subdir: str = "log",
    unavailable_label: str = "unsubmitted",
) -> dict:
    log_dir = os.path.join(job_dir, log_subdir)
    cluster = _latest_condor_cluster(log_dir, prefix=prefix)
    if cluster is None:
        return {
            "state": "unsubmitted",
            "label": unavailable_label,
            "detail": "logs unavailable",
            "cluster": None,
            "done": None,
            "total": None,
        }

    log_path = os.path.join(log_dir, f"{prefix}.log_{cluster}")
    if not os.path.isfile(log_path):
        return {
            "state": "unsubmitted",
            "label": unavailable_label,
            "detail": f"{prefix}.log_{cluster} unavailable",
            "cluster": cluster,
            "done": None,
            "total": None,
        }

    submitted = set()
    terminated = set()
    held = set()
    removed = set()
    try:
        with open(log_path, errors="replace") as f:
            for line in f:
                match = CONDOR_EVENT_RE.match(line)
                if not match:
                    continue
                event, event_cluster, proc = match.groups()
                if int(event_cluster) != cluster:
                    continue
                if event == "000":
                    submitted.add(proc)
                elif event == "005":
                    terminated.add(proc)
                elif event == "012":
                    held.add(proc)
                elif event == "009":
                    removed.add(proc)
    except OSError as exc:
        return {
            "state": "error",
            "label": "status unavailable",
            "detail": str(exc),
            "cluster": cluster,
            "done": None,
            "total": None,
        }

    total = len(submitted)
    if total == 0:
        return {
            "state": "unsubmitted",
            "label": unavailable_label,
            "detail": f"{prefix}.log_{cluster} has no submit events",
            "cluster": cluster,
            "done": None,
            "total": None,
        }

    done = len(terminated)
    state = "done" if done >= total else "warn" if held or removed else "running"
    job_kind = "jobs"
    if prefix == "hadd":
        job_kind = "hadd jobs"
    elif data:
        resubmitted = data.get("resubmitted_jobs", {})
        resubmitted_count = resubmitted.get("count")
        submitted_jobs = data.get("submitted_jobs")
        if (resubmitted_count and int(resubmitted_count) == total) or (
            isinstance(submitted_jobs, int) and total < submitted_jobs
        ):
            job_kind = "resubmission jobs"
    extras = []
    if held:
        extras.append(f"{len(held)} held")
    if removed:
        extras.append(f"{len(removed)} removed")
    detail = f"cluster {cluster}"
    if extras:
        detail += f" ({', '.join(extras)})"
    return {
        "state": state,
        "label": f"{done}/{total} {job_kind} done",
        "job_kind": job_kind,
        "detail": detail,
        "cluster": cluster,
        "done": done,
        "total": total,
    }


def read_hadd_condor_status(job_dir: str) -> dict:
    return read_condor_status(
        job_dir,
        prefix="hadd",
        log_subdir=os.path.join("haddscripts", "hadd_logs"),
        unavailable_label="Unavailable",
    )


def _condor_queue_state(cluster: int | None) -> str:
    if cluster is None:
        return "unknown"
    try:
        completed = subprocess.run(
            ["condor_q", str(cluster), "-json"],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError:
        return "error"
    if completed.returncode != 0:
        return "error"
    try:
        jobs = json.loads(completed.stdout or "")
    except (json.JSONDecodeError, TypeError):
        return "error"
    if not isinstance(jobs, list):
        return "error"
    return "empty" if len(jobs) == 0 else "present"


def _condor_queue_empty(cluster: int | None) -> bool | None:
    state = _condor_queue_state(cluster)
    if state == "empty":
        return True
    if state == "present":
        return False
    return None


def _script_path(script_name: str) -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), script_name)


def _submission_status_recorded(data: dict, action: str, submission: dict) -> bool:
    expected_count = submission["count"]
    expected_cluster = submission["cluster_id"]
    if action == "checkoutputs":
        recorded = data.get("resubmitted_jobs", {})
        return (
            recorded.get("count") == expected_count
            and recorded.get("cluster_id") == expected_cluster
        )
    if action in {"hadd_submit", "hadd_check"}:
        return (
            data.get("hadd_jobs_submitted") == expected_count
            and data.get("hadd_submission_cluster_id") == expected_cluster
        )
    return False


def _record_unreconciled_submission(
    job_dir: str,
    data: dict,
    action: str,
    submission: dict,
    completed: subprocess.CompletedProcess,
) -> dict:
    cluster_id = submission["cluster_id"]
    submitted_jobs = submission["count"]
    state = "submitted_but_status_update_failed"
    reason = (
        f"Condor submitted {submitted_jobs} job(s) to cluster {cluster_id}, but "
        "the submission status was not persisted. The submission command will not "
        "be retried; reconcile this cluster before resuming automation."
    )
    auto = _automation(data)
    auto["step"] = "blocked"
    auto["blocked_reason"] = reason
    auto["blocked_state"] = state
    auto["active_pid"] = None
    auto["active_command"] = ""
    auto["unreconciled_submission"] = {
        "state": state,
        "action": action,
        "cluster_id": cluster_id,
        "submitted_jobs": submitted_jobs,
        "detected_at": _now_iso(),
    }
    _add_automation_warning(data, reason)
    _set_automation_action(
        data,
        action,
        step="blocked",
        returncode=completed.returncode,
        stdout=completed.stdout or "",
        stderr=completed.stderr or "",
        active_pid=None,
        active_command="",
    )
    _append_stage(
        data,
        "Automation blocked",
        state,
        f"{action}: cluster {cluster_id}, {submitted_jobs} jobs",
    )
    save_status(job_dir, data)
    return data


def _run_automation_command(
    job_dir: str,
    command: list[str],
    action: str,
    *,
    expect_condor_submit: bool = True,
    max_retries: int = AUTOMATION_MAX_SUBMIT_RETRIES,
) -> dict:
    attempts = max(1, max_retries if expect_condor_submit else 1)
    last_completed: subprocess.CompletedProcess | None = None
    for attempt in range(1, attempts + 1):
        proc = subprocess.Popen(
            command,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        data = load_status(job_dir)
        _set_automation_action(
            data,
            action,
            step=f"{action}_attempt_{attempt}",
            active_pid=proc.pid,
            active_command=shlex.join(command),
        )
        save_status(job_dir, data)
        stdout, stderr = proc.communicate()
        completed = subprocess.CompletedProcess(
            command, proc.returncode, stdout, stderr
        )
        last_completed = completed
        data = load_status(job_dir)
        stdout = completed.stdout or ""
        stderr = completed.stderr or ""
        _set_automation_action(
            data,
            action,
            step=f"{action}_attempt_{attempt}",
            returncode=completed.returncode,
            stdout=stdout,
            stderr=stderr,
            active_pid=None,
            active_command=shlex.join(command),
        )
        save_status(job_dir, data)

        post_data = load_status(job_dir, create=False)
        checked = post_data.get("checked_outputs", {})
        hadd_checked = post_data.get("hadd_jobs_checked", {})
        submission_output = "\n".join(part for part in (stdout, stderr) if part)
        submission = (
            parse_condor_submit_output(submission_output)
            if expect_condor_submit
            else None
        )
        if submission is not None:
            if not _submission_status_recorded(post_data, action, submission):
                return _record_unreconciled_submission(
                    job_dir,
                    post_data,
                    action,
                    submission,
                    completed,
                )
            data = load_status(job_dir)
            auto = _automation(data)
            auto["blocked_reason"] = ""
            auto["blocked_state"] = ""
            auto["unreconciled_submission"] = None
            _set_automation_action(
                data,
                action,
                step=action,
                returncode=completed.returncode,
                stdout=stdout,
                stderr=stderr,
                active_pid=None,
                active_command="",
            )
            save_status(job_dir, data)
            return data

        action_completed_without_submit = (
            (action == "checkoutputs" and checked.get("status") == "All jobs done")
            or (action == "hadd_check" and hadd_checked.get("status") == "All done")
            or (
                action == "hadd_submit"
                and (
                    post_data.get("all_done")
                    or hadd_checked.get("status") == "All done"
                )
            )
        )
        action_blocked_without_submit = action == "checkoutputs" and _all_jobs_failed(
            post_data
        )
        submit_ok = (
            (not expect_condor_submit)
            or action_completed_without_submit
            or action_blocked_without_submit
        )
        if (completed.returncode == 0 and submit_ok) or action_blocked_without_submit:
            data = load_status(job_dir)
            if action_blocked_without_submit:
                auto = _automation(data)
                auto["step"] = "blocked"
                auto["blocked_state"] = ""
                auto["blocked_reason"] = (
                    "All checked jobs failed; skipping condor resubmission because this likely indicates a code/config bug."
                )
                _add_automation_warning(data, auto["blocked_reason"])
            _set_automation_action(
                data,
                action,
                step="blocked" if action_blocked_without_submit else action,
                returncode=completed.returncode,
                stdout=stdout,
                stderr=stderr,
            )
            save_status(job_dir, data)
            return data

        if attempt < attempts:
            time.sleep(_automation_retry_delay(attempt, expect_condor_submit))

    stdout = last_completed.stdout if last_completed else ""
    stderr = last_completed.stderr if last_completed else ""
    detail = (
        f"{action} failed after {attempts} attempt(s); "
        "condor_submit confirmation was not found in stdout or stderr"
        if expect_condor_submit
        else f"{action} failed after {attempts} attempt(s)"
    )
    data = load_status(job_dir)
    auto = _automation(data)
    auto["step"] = "blocked"
    auto["blocked_state"] = ""
    auto["blocked_reason"] = detail
    _add_automation_warning(data, detail)
    _set_automation_action(
        data,
        action,
        step="blocked",
        returncode=last_completed.returncode if last_completed else None,
        stdout=stdout,
        stderr=stderr,
    )
    save_status(job_dir, data)
    return data


def _arguments(job_dir: str) -> dict:
    args_path = os.path.join(job_dir, "arguments.json")
    if not os.path.isfile(args_path):
        return {}
    try:
        with open(args_path) as f:
            return json.load(f)
    except Exception:
        return {}


def _truthy(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _output_dir(job_dir: str, data: dict) -> str:
    output_dir = data.get("output_dir") or ""
    if not output_dir:
        output_dir = str(_arguments(job_dir).get("outputDir", ""))
    return os.path.abspath(output_dir) if output_dir else ""


def _has_array_outputs(output_dir: str) -> bool:
    if not output_dir or not os.path.isdir(output_dir):
        return False
    for array_dir in glob(os.path.join(output_dir, "arrays_*")):
        if glob(os.path.join(array_dir, "*", "*", "*.root")):
            return True
    return False


def _is_array_job(job_dir: str, data: dict) -> bool:
    return _truthy(_arguments(job_dir).get("isArray", False))


def _mark_cluster_processed(data: dict, key: str, cluster: int | None) -> None:
    if cluster is None:
        return
    auto = _automation(data)
    clusters = auto.setdefault(key, [])
    cluster_str = str(cluster)
    if cluster_str not in clusters:
        clusters.append(cluster_str)


def _cluster_processed(data: dict, key: str, cluster: int | None) -> bool:
    if cluster is None:
        return False
    return str(cluster) in _automation(data).get(key, [])


def _all_jobs_failed(data: dict) -> bool:
    checked = data.get("checked_outputs", {})
    checked_jobs = checked.get("checked_jobs")
    missing = checked.get("missing_outputs")
    return (
        isinstance(checked_jobs, int)
        and isinstance(missing, int)
        and checked_jobs > 0
        and missing >= checked_jobs
    )


def _automation_should_submit(data: dict, action: str) -> bool:
    auto = _automation(data)
    if auto.get("step") in {"blocked", "submitted_but_status_update_failed"}:
        return False
    if auto.get("unreconciled_submission"):
        return False
    last_action = auto.get("last_action")
    last_step = str(auto.get("step") or "")
    return not (last_action == action and last_step.startswith(action))


def advance_job_automation(
    job_dir: str,
    condor_status: dict | None,
    hadd_condor_status: dict | None,
) -> dict:
    data = load_status(job_dir, create=False)
    auto = _automation(data)
    if (
        data.get("all_done")
        or auto.get("step") in {"blocked", "submitted_but_status_update_failed"}
        or auto.get("unreconciled_submission")
    ):
        return data

    active_step = str(auto.get("step") or "")
    active_action = str(auto.get("last_action") or "")
    active_pid = _automation_active_pid(data)
    expected_fragment = str(auto.get("active_command") or "")
    if _automation_step_is_active(active_action, active_step):
        script_name = (
            "checkoutputs.py" if active_action == "checkoutputs" else "haddoutputs.py"
        )
        if active_pid is None:
            inferred_pid = _find_matching_automation_pid(job_dir, script_name)
            if inferred_pid is not None:
                auto["active_pid"] = inferred_pid
                auto["active_command"] = auto.get("active_command") or script_name
                save_status(job_dir, data)
                active_pid = inferred_pid
        if active_pid is not None and _process_alive(
            active_pid, expected_fragment=expected_fragment or script_name
        ):
            return data
        return record_automation_blocked(
            job_dir,
            f"Automation subprocess pid {active_pid if active_pid is not None else 'unknown'} is no longer running; the last {active_action or active_step} step likely died.",
        )

    if _all_jobs_failed(data):
        return record_automation_blocked(
            job_dir,
            "All checked jobs failed; skipping condor resubmission because this likely indicates a code/config bug.",
        )

    if condor_status:
        cluster = condor_status.get("cluster")
        cluster_state = condor_status.get("state")
        queue_empty = (
            cluster_state in {"done", "warn", "unsubmitted"}
            and _condor_queue_empty(cluster) is True
        )
        if cluster_state == "done" or queue_empty:
            if not _cluster_processed(data, "processed_job_clusters", cluster):
                _mark_cluster_processed(data, "processed_job_clusters", cluster)
                _set_automation_action(data, "checkoutputs", step="checking_outputs")
                save_status(job_dir, data)
                expect_submit = True
                return _run_automation_command(
                    job_dir,
                    [
                        sys.executable,
                        _script_path("checkoutputs.py"),
                        job_dir,
                        "--condor",
                    ],
                    "checkoutputs",
                    expect_condor_submit=expect_submit,
                )

    data = load_status(job_dir, create=False)
    auto = _automation(data)
    if auto.get("step") == "blocked" or data.get("all_done"):
        return data
    stale_checkoutputs = _resume_stale_checkoutputs(job_dir, data, condor_status)
    if stale_checkoutputs:
        return stale_checkoutputs
    checked = data.get("checked_outputs", {})
    if checked.get("status") != "All jobs done":
        return data

    if not _is_array_job(job_dir, data):
        return record_processing_complete(
            job_dir, "All non-array processing outputs done"
        )

    output_dir = _output_dir(job_dir, data)
    if not _has_array_outputs(output_dir):
        return record_automation_blocked(
            job_dir,
            "isArray is enabled, but no arrays_* root outputs were found; skipping hadd automation.",
        )

    hadd_checked = data.get("hadd_jobs_checked", {})
    if hadd_checked.get("status") == "All done":
        auto["step"] = "done"
        data["all_done"] = True
        save_status(job_dir, data)
        return data

    if data.get("hadd_jobs_submitted") is None:
        if not _automation_should_submit(data, "hadd_submit"):
            return data
        _set_automation_action(data, "hadd_submit", step="submitting_hadd")
        save_status(job_dir, data)
        return _run_automation_command(
            job_dir,
            [sys.executable, _script_path("haddoutputs.py"), job_dir, "--condor"],
            "hadd_submit",
            expect_condor_submit=True,
        )

    if hadd_condor_status:
        cluster = hadd_condor_status.get("cluster")
        cluster_state = hadd_condor_status.get("state")
        queue_empty = (
            cluster_state in {"done", "warn", "unsubmitted"}
            and _condor_queue_empty(cluster) is True
        )
        if (cluster_state == "done" or queue_empty) and not _cluster_processed(
            data, "processed_hadd_clusters", cluster
        ):
            _mark_cluster_processed(data, "processed_hadd_clusters", cluster)
            _set_automation_action(data, "hadd_check", step="checking_hadd")
            save_status(job_dir, data)
            return _run_automation_command(
                job_dir,
                [sys.executable, _script_path("haddoutputs.py"), job_dir, "--condor"],
                "hadd_check",
                expect_condor_submit=True,
            )

    return data


def _condor_status_text(status: dict | None) -> Text:
    if status is None:
        return Text(f"{ICONS['clock']} checking...", style="dim cyan")
    state = status.get("state")
    style = {
        "done": "bold green",
        "running": "bold cyan",
        "warn": "bold yellow",
        "unsubmitted": "dim",
        "error": "bold red",
    }.get(state, "cyan")
    icon = {
        "done": ICONS["ok"],
        "running": ICONS["clock"],
        "warn": ICONS["warn"],
        "unsubmitted": ICONS["info"],
        "error": ICONS["error"],
    }.get(state, ICONS["info"])
    done = status.get("done")
    total = status.get("total")
    job_kind = status.get("job_kind")
    if isinstance(done, int) and isinstance(total, int) and job_kind:
        text = Text(f"{icon} {done}/{total} done", style=style)
        if job_kind != "jobs":
            text.append(f"\n{job_kind}", style=style)
    else:
        text = Text(f"{icon} {status.get('label', 'checking...')}", style=style)
    detail = status.get("detail")
    if detail:
        text.append(f"\n{detail}", style="dim")
    return text


def _condor_progress(status: dict | None) -> Progress | None:
    if not RICH_AVAILABLE or status is None:
        return None
    done = status.get("done")
    total = status.get("total")
    if not isinstance(done, int) or not isinstance(total, int) or total <= 0:
        return None
    complete_style = "green" if done >= total else "cyan"
    progress = Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=14, complete_style=complete_style, finished_style="green"),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        expand=False,
    )
    progress.add_task("", total=total, completed=min(done, total))
    return progress


def _compact_condor_text(status: dict | None, label: str = "Condor") -> Text:
    if status is None:
        return Text(f"{ICONS['clock']} {label}: checking...", style="dim cyan")
    state = status.get("state")
    style = {
        "done": "bold green",
        "running": "bold cyan",
        "warn": "bold yellow",
        "unsubmitted": "dim",
        "error": "bold red",
    }.get(state, "cyan")
    icon = {
        "done": ICONS["ok"],
        "running": ICONS["clock"],
        "warn": ICONS["warn"],
        "unsubmitted": ICONS["info"],
        "error": ICONS["error"],
    }.get(state, ICONS["info"])
    done = status.get("done")
    total = status.get("total")
    job_kind = status.get("job_kind")
    if isinstance(done, int) and isinstance(total, int):
        suffix = ""
        if job_kind and job_kind not in ("jobs", "hadd jobs"):
            suffix = f" {job_kind.replace(' jobs', '')}"
        return Text(f"{icon} {label}: {done}/{total}{suffix}", style=style)
    return Text(f"{icon} {label}: {status.get('label', 'checking...')}", style=style)


def _processing_cell(data: dict, condor: dict | None):
    checked = data.get("checked_outputs", {})
    submitted = _fmt_value(data.get("submitted_jobs"), "-")
    checked_status = checked.get("status") or "-"
    checked_style = (
        "green"
        if checked_status == "All jobs done"
        else "yellow" if checked_status != "-" else "cyan"
    )
    auto = _automation(data)
    lines = [
        Text(f"{ICONS['submit']} Submitted: {submitted}", style="green"),
        _compact_condor_text(condor, "Condor"),
    ]
    progress = _condor_progress(condor)
    if progress:
        lines.append(progress)
    lines.append(
        Text(f"{ICONS['check']} Outputs: {checked_status}", style=checked_style)
    )
    resub = data.get("resubmitted_jobs", {})
    if resub.get("count"):
        lines.append(
            Text(f"{ICONS['resubmit']} Resubmitted: {resub['count']}", style="yellow")
        )
    if auto.get("blocked_reason"):
        lines.append(
            Text(f"{ICONS['warn']} {auto['blocked_reason']}", style="bold red")
        )
        output_summary = _automation_output_summary(data)
        if output_summary:
            lines.append(Text(output_summary, style="dim red"))
    elif auto.get("last_action"):
        lines.append(
            Text(
                f"{ICONS['info']} Auto: {auto.get('step') or auto['last_action']}",
                style="dim cyan",
            )
        )
    return Group(*lines)


def _hadd_cell(data: dict, hadd_condor: dict | None):
    hadd = data.get("hadd_jobs_checked", {})
    hadd_submitted = _fmt_value(data.get("hadd_jobs_submitted"), "-")
    hadd_checked = hadd.get("status") or "-"
    hadd_style = (
        "green"
        if data.get("all_done")
        else "yellow" if data.get("hadd_jobs_submitted") is not None else "cyan"
    )
    lines = [
        Text(f"{ICONS['hadd']} Submitted: {hadd_submitted}", style=hadd_style),
        _compact_condor_text(hadd_condor, "Condor"),
    ]
    progress = _condor_progress(hadd_condor)
    if progress:
        lines.append(progress)
    lines.append(Text(f"{ICONS['check']} Checked: {hadd_checked}", style=hadd_style))
    return Group(*lines)


def _compact_status(status: dict | None, unavailable: str = "-") -> str:
    if status is None:
        return "checking"
    done = status.get("done")
    total = status.get("total")
    if isinstance(done, int) and isinstance(total, int):
        label = f"{done}/{total}"
        job_kind = status.get("job_kind")
        if job_kind and job_kind not in ("jobs", "hadd jobs"):
            label += " resub"
        return label
    label = str(status.get("label") or unavailable)
    return "n/a" if label == "Unavailable" else label


def _short_status(value: str) -> str:
    return {
        "All jobs done": "done",
        "All done": "done",
        "Unavailable": "n/a",
        "unsubmitted": "unsubmitted",
    }.get(value, value)


def _compact_processing_value(data: dict, condor: dict | None) -> str:
    checked = data.get("checked_outputs", {})
    auto = _automation(data)
    pieces = [
        f"sub {_fmt_value(data.get('submitted_jobs'), '-')}",
        f"cond {_compact_status(condor)}",
        f"out {_short_status(checked.get('status') or '-')}",
    ]
    resub = data.get("resubmitted_jobs", {})
    if resub.get("count"):
        pieces.append(f"resub {resub['count']}")
    if auto.get("blocked_reason"):
        pieces.append("auto blocked")
    elif auto.get("last_action"):
        pieces.append(f"auto {auto.get('step') or auto['last_action']}")
    return " | ".join(pieces)


def _compact_hadd_value(data: dict, hadd_condor: dict | None) -> str:
    hadd = data.get("hadd_jobs_checked", {})
    return " | ".join(
        [
            f"sub {_fmt_value(data.get('hadd_jobs_submitted'), '-')}",
            f"cond {_compact_status(hadd_condor, unavailable='Unavailable')}",
            f"check {_short_status(hadd.get('status') or '-')}",
        ]
    )


def _build_detailed_table(
    rows: list[tuple[str, dict]], condor_known: dict, hadd_condor_known: dict
):
    table = Table(box=box.SIMPLE_HEAVY, expand=True, show_lines=True)
    table.add_column("Job", style="bold", ratio=1, min_width=10, max_width=20)
    table.add_column("State", justify="center", ratio=1, max_width=13)
    table.add_column("Processing", ratio=2, min_width=22, max_width=28, no_wrap=False)
    table.add_column("Hadd", ratio=2, min_width=20, no_wrap=False)
    table.add_column("Timeline", ratio=3, no_wrap=False)

    for job_dir, data in rows:
        condor = condor_known.get(job_dir)
        hadd_condor = hadd_condor_known.get(job_dir)
        overall_state = _job_overall_state(data)
        overall_label = f"{_state_icon(overall_state)} " + (
            "All done" if data.get("all_done") else "In progress"
        )
        table.add_row(
            Text(_job_title(data), style="bold white"),
            Text(overall_label, style=_state_style(overall_state)),
            _processing_cell(data, condor),
            _hadd_cell(data, hadd_condor),
            Text(_stage_summary(data), style="white"),
        )
    return table


def _build_compact_table(
    rows: list[tuple[str, dict]],
    condor_known: dict,
    hadd_condor_known: dict,
    title: str | None = None,
):
    table = Table(
        title=title,
        box=box.SIMPLE,
        expand=True,
        show_lines=False,
        header_style="bold cyan",
    )
    table.add_column("Job", style="bold", min_width=10, max_width=24)
    table.add_column("State", max_width=12)
    table.add_column("Processing", ratio=3, no_wrap=True, overflow="ellipsis")
    table.add_column("Hadd", ratio=2, no_wrap=True, overflow="ellipsis")

    for job_dir, data in rows:
        overall_state = _job_overall_state(data)
        state_label = "done" if data.get("all_done") else "active"
        table.add_row(
            Text(_job_title(data), style="bold white"),
            Text(
                f"{_state_icon(overall_state)} {state_label}",
                style=_state_style(overall_state),
            ),
            Text(
                _compact_processing_value(data, condor_known.get(job_dir)),
                style="white",
            ),
            Text(
                _compact_hadd_value(data, hadd_condor_known.get(job_dir)), style="white"
            ),
        )
    return table


def _stage_history_line(stages: Counter):
    if not stages:
        return None
    text = Text(
        f"{ICONS['info']} Stage history  ",
        style="bold cyan",
        no_wrap=True,
        overflow="ellipsis",
    )
    for idx, (stage, count) in enumerate(stages.most_common()):
        if idx:
            text.append("  •  ", style="dim")
        icon, color = _job_stage_label(stage)
        text.append(f"{icon} {stage}: {count}", style=color)
    return text


def _stage_label(stage: str, status: str) -> Text:
    label = Text(stage)
    if status == "done":
        label.stylize("bold green")
    elif status == "warn":
        label.stylize("bold yellow")
    elif status == "error":
        label.stylize("bold red")
    else:
        label.stylize("bold cyan")
    return label


def _progress_badge(
    label: str, value: str, state: str = "info", rich_markup: bool = True
) -> str:
    icon = _state_icon(state)
    if rich_markup and RICH_AVAILABLE:
        return f"[{_state_style(state)}]{icon} {label}:[/] {value}"
    return f"{icon} {label}: {value}"


def _automation_output_summary(data: dict) -> str:
    auto = _automation(data)
    pieces = []
    stdout = " ".join(str(auto.get("last_stdout") or "").split())
    stderr = " ".join(str(auto.get("last_stderr") or "").split())
    if stdout:
        pieces.append(f"stdout: {_trim_output(stdout, 180)}")
    if stderr:
        pieces.append(f"stderr: {_trim_output(stderr, 180)}")
    return " | ".join(pieces)


def build_dashboard(
    root: str = ".",
    condor_statuses: dict[str, dict] | None = None,
    hadd_condor_statuses: dict[str, dict] | None = None,
    terminal_height: int | None = None,
):
    job_dirs = _job_dirs(root)
    if not job_dirs:
        if RICH_AVAILABLE:
            return Panel.fit(
                "[bold yellow]No job directories found.[/]", title="Dashboard"
            )
        return "No job directories found."

    statuses = []
    for job_dir in job_dirs:
        data = load_status(job_dir)
        statuses.append((job_dir, data))

    condor_known = condor_statuses or {}
    hadd_condor_known = hadd_condor_statuses or {}
    statuses = sorted(
        statuses,
        key=lambda item: _job_progress_rank(
            item[1],
            condor_known.get(item[0]),
            hadd_condor_known.get(item[0]),
        ),
    )

    total = len(statuses)
    submitted = sum(1 for _, data in statuses if data.get("submitted_jobs") is not None)
    checked_done = sum(
        1
        for _, data in statuses
        if data.get("checked_outputs", {}).get("status") == "All jobs done"
    )
    resubmitted_pending = sum(
        1 for _, data in statuses if data.get("resubmitted_jobs", {}).get("count")
    )
    hadd_submitted = sum(
        1 for _, data in statuses if data.get("hadd_jobs_submitted") is not None
    )
    all_done = sum(1 for _, data in statuses if data.get("all_done"))
    automation_blocked = sum(
        1 for _, data in statuses if _automation(data).get("step") == "blocked"
    )
    condor_done = sum(
        1
        for job_dir, _ in statuses
        if condor_known.get(job_dir, {}).get("state") == "done"
    )
    condor_running = sum(
        1
        for job_dir, _ in statuses
        if condor_known.get(job_dir, {}).get("state") in ("running", "warn")
    )
    condor_unsubmitted = sum(
        1
        for job_dir, _ in statuses
        if condor_known.get(job_dir, {}).get("state") == "unsubmitted"
    )
    hadd_condor_done = sum(
        1
        for job_dir, _ in statuses
        if hadd_condor_known.get(job_dir, {}).get("state") == "done"
    )
    hadd_condor_running = sum(
        1
        for job_dir, _ in statuses
        if hadd_condor_known.get(job_dir, {}).get("state") in ("running", "warn")
    )
    stages = Counter()
    for _, data in statuses:
        for entry in data.get("stages_passed", []):
            if isinstance(entry, dict) and entry.get("stage"):
                stages[entry["stage"]] += 1

    if not RICH_AVAILABLE:
        lines = [f"{ICONS['info']} Condor Dashboard"]
        lines.append(
            f"{ICONS['submit']} Submitted: {submitted}/{total}   "
            f"{ICONS['check']} Checked done: {checked_done}   "
            f"{ICONS['resubmit']} Resubmitted: {resubmitted_pending}   "
            f"{ICONS['hadd']} Hadd submitted: {hadd_submitted}   "
            f"{ICONS['done']} All done: {all_done}   "
            f"{ICONS['warn']} Blocked: {automation_blocked}"
        )
        for job_dir, data in statuses:
            condor = condor_known.get(job_dir)
            hadd_condor = hadd_condor_known.get(job_dir)
            condor_value = condor.get("label") if condor else "checking..."
            hadd_condor_value = (
                hadd_condor.get("label") if hadd_condor else "checking..."
            )
            overall_state = _job_overall_state(data)
            checked = data.get("checked_outputs", {})
            hadd = data.get("hadd_jobs_checked", {})
            lines.append(f"{_state_icon(overall_state)} {_job_title(data)}")
            lines.append(
                f"  {_progress_badge('Submitted', _fmt_value(data.get('submitted_jobs'), '-'), 'done' if data.get('submitted_jobs') is not None else 'info', rich_markup=False)}"
            )
            lines.append(
                f"  {_progress_badge('Condor', condor_value, condor.get('state', 'info') if condor else 'info', rich_markup=False)}"
            )
            lines.append(
                f"  {_progress_badge('Hadd Condor', hadd_condor_value, hadd_condor.get('state', 'info') if hadd_condor else 'info', rich_markup=False)}"
            )
            lines.append(
                f"  {_progress_badge('Outputs', checked.get('status') or '-', 'done' if checked.get('status') == 'All jobs done' else 'warn' if checked.get('status') else 'info', rich_markup=False)}"
            )
            hadd_value = f"{_fmt_value(data.get('hadd_jobs_submitted'), '-')} | {hadd.get('status') or '-'}"
            hadd_state = (
                "done"
                if data.get("all_done")
                else "warn" if data.get("hadd_jobs_submitted") is not None else "info"
            )
            lines.append(
                f"  {_progress_badge('Hadd', hadd_value, hadd_state, rich_markup=False)}"
            )
            auto = _automation(data)
            if auto.get("blocked_reason"):
                lines.append(
                    f"  {_progress_badge('Automation', auto['blocked_reason'], 'error', rich_markup=False)}"
                )
                output_summary = _automation_output_summary(data)
                if output_summary:
                    lines.append(f"    {output_summary}")
            elif auto.get("last_action"):
                lines.append(
                    f"  {_progress_badge('Automation', auto.get('step') or auto['last_action'], 'info', rich_markup=False)}"
                )
            lines.append(f"  Timeline: {_stage_summary(data)}")
        if stages:
            pieces = []
            for stage, count in stages.most_common():
                icon, _ = _job_stage_label(stage)
                pieces.append(f"{icon} {stage}: {count}")
            lines.append(f"{ICONS['info']} Stage history: " + "  |  ".join(pieces))
        return "\n".join(lines)

    header_text = Text()
    header_text.append(f"{ICONS['info']} Condor Dashboard\n", style="bold white")
    header_text.append(
        f"{ICONS['submit']} Submitted: {submitted}/{total}   ", style="green"
    )
    header_text.append(
        f"{ICONS['check']} Checked done: {checked_done}   ", style="cyan"
    )
    header_text.append(
        f"{ICONS['resubmit']} Resubmitted: {resubmitted_pending}   ", style="yellow"
    )
    header_text.append(
        f"{ICONS['hadd']} Hadd submitted: {hadd_submitted}   ", style="magenta"
    )
    header_text.append(f"{ICONS['done']} All done: {all_done}   ", style="bold green")
    header_text.append(
        f"{ICONS['warn']} Blocked: {automation_blocked}\n", style="bold red"
    )
    header_text.append(
        f"{ICONS['clock']} Condor live: {condor_done} complete, {condor_running} running, {condor_unsubmitted} unsubmitted/logless",
        style="cyan",
    )
    header_text.append(
        f"   {ICONS['hadd']} Hadd condor: {hadd_condor_done} complete, {hadd_condor_running} running",
        style="magenta",
    )

    renderables = [
        Panel(
            header_text,
            title="Condor Dashboard",
            subtitle=f"Updated {_now_iso()}",
            border_style="cyan",
            box=box.ROUNDED,
        ),
    ]

    height = terminal_height or 32
    available_rows = max(8, height - 11)
    detailed_rows = 6
    compact_rows = 1
    active_statuses = [
        (job_dir, data) for job_dir, data in statuses if not data.get("all_done")
    ]
    done_statuses = [
        (job_dir, data) for job_dir, data in statuses if data.get("all_done")
    ]
    detail_estimate = len(statuses) * detailed_rows
    mixed_estimate = (
        len(done_statuses) * compact_rows + len(active_statuses) * detailed_rows + 4
    )
    view_mode = "compact"

    if detail_estimate <= available_rows:
        renderables.append(
            _build_detailed_table(statuses, condor_known, hadd_condor_known)
        )
        view_mode = "detailed"
    elif active_statuses and mixed_estimate <= available_rows:
        if done_statuses:
            renderables.append(
                _build_compact_table(
                    done_statuses,
                    condor_known,
                    hadd_condor_known,
                    title=f"{ICONS['done']} Done jobs ({len(done_statuses)})",
                )
            )
        renderables.append(
            _build_detailed_table(active_statuses, condor_known, hadd_condor_known)
        )
        view_mode = "mixed"
    else:
        renderables.append(
            _build_compact_table(
                statuses,
                condor_known,
                hadd_condor_known,
            )
        )

    stage_line = _stage_history_line(stages)
    if stage_line and view_mode != "compact":
        renderables.append(stage_line)
    return Group(*renderables)


def render_dashboard(
    root: str = ".",
    console: Console | None = None,
    condor_statuses: dict[str, dict] | None = None,
    hadd_condor_statuses: dict[str, dict] | None = None,
) -> None:
    console = console or Console() if RICH_AVAILABLE else None
    terminal_height = console.size.height if RICH_AVAILABLE and console else None
    dashboard = build_dashboard(
        root,
        condor_statuses=condor_statuses,
        hadd_condor_statuses=hadd_condor_statuses,
        terminal_height=terminal_height,
    )
    if RICH_AVAILABLE:
        console.print(dashboard)
    else:
        print(dashboard)


def watch_dashboard(root: str = ".", interval: float = 10.0, auto: bool = True) -> None:
    if not RICH_AVAILABLE:
        while True:
            job_dirs = _job_dirs(root)
            condor_statuses = {
                job_dir: read_condor_status(job_dir, load_status(job_dir, create=False))
                for job_dir in job_dirs
            }
            hadd_condor_statuses = {
                job_dir: read_hadd_condor_status(job_dir) for job_dir in job_dirs
            }
            if auto:
                for job_dir in job_dirs:
                    advance_job_automation(
                        job_dir,
                        condor_statuses.get(job_dir),
                        hadd_condor_statuses.get(job_dir),
                    )
            render_dashboard(
                root,
                condor_statuses=condor_statuses,
                hadd_condor_statuses=hadd_condor_statuses,
            )
            time.sleep(interval)

    console = Console()
    job_dirs = _job_dirs(root)
    condor_statuses: dict[str, dict] = {}
    hadd_condor_statuses: dict[str, dict] = {}
    pending: dict[tuple[str, str], concurrent.futures.Future] = {}
    max_workers = max(1, min(12, (len(job_dirs) or 1) * 3))

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        next_refresh = 0.0
        with Live(
            build_dashboard(
                root,
                condor_statuses,
                hadd_condor_statuses,
                terminal_height=console.size.height,
            ),
            console=console,
            refresh_per_second=1,
            transient=False,
        ) as live:
            try:
                while True:
                    now = time.monotonic()
                    job_dirs = _job_dirs(root)
                    if now >= next_refresh:
                        for job_dir in job_dirs:
                            data = load_status(job_dir, create=False)
                            job_key = ("job", job_dir)
                            hadd_key = ("hadd", job_dir)
                            if job_key not in pending:
                                pending[job_key] = executor.submit(
                                    read_condor_status, job_dir, data
                                )
                            if hadd_key not in pending:
                                pending[hadd_key] = executor.submit(
                                    read_hadd_condor_status, job_dir
                                )
                        next_refresh = now + interval

                    for key, future in list(pending.items()):
                        if not future.done():
                            continue
                        kind, job_dir = key
                        try:
                            if kind == "auto":
                                future.result()
                            elif kind == "hadd":
                                hadd_condor_statuses[job_dir] = future.result()
                            else:
                                condor_statuses[job_dir] = future.result()
                        except Exception as exc:
                            status = {
                                "state": "error",
                                "label": "status unavailable",
                                "detail": str(exc),
                                "cluster": None,
                                "done": None,
                                "total": None,
                            }
                            if kind == "hadd":
                                hadd_condor_statuses[job_dir] = status
                            elif kind == "job":
                                condor_statuses[job_dir] = status
                            else:
                                record_automation_warning(
                                    job_dir,
                                    f"Automation error: {exc}",
                                )
                        del pending[key]

                    if auto:
                        for job_dir in job_dirs:
                            auto_key = ("auto", job_dir)
                            if auto_key in pending:
                                continue
                            if job_dir not in condor_statuses:
                                continue
                            pending[auto_key] = executor.submit(
                                advance_job_automation,
                                job_dir,
                                condor_statuses.get(job_dir),
                                hadd_condor_statuses.get(job_dir),
                            )

                    live.update(
                        build_dashboard(
                            root,
                            condor_statuses,
                            hadd_condor_statuses,
                            terminal_height=console.size.height,
                        )
                    )
                    time.sleep(0.25)
            except KeyboardInterrupt:
                pass


def main() -> None:
    parser = argparse.ArgumentParser(description="Display condor job status dashboard")
    parser.add_argument(
        "--root",
        default=".",
        help="Directory to scan for job directories (default: current directory)",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=10.0,
        help="Seconds between condor log refreshes in live mode (default: %(default)s)",
    )
    parser.add_argument(
        "--once",
        action="store_true",
        help="Render once and exit instead of opening the live dashboard",
    )
    parser.add_argument(
        "--no-auto",
        action="store_true",
        help="Disable live automation; --once is always read-only.",
    )
    args = parser.parse_args()
    if not args.once and not args.no_auto:
        try:
            _validate_automation_environment()
        except RuntimeError as exc:
            print(f"[red][b]{exc}[/][/]", file=sys.stderr)
            raise SystemExit(1)
        resumed = resume_blocked_automation_jobs(args.root)
        if resumed:
            print(
                f"[green]Resumed {resumed} blocked automation job(s) after startup recheck.[/]"
            )
    if args.once:
        job_dirs = _job_dirs(args.root)
        condor_statuses = {
            job_dir: read_condor_status(job_dir, load_status(job_dir, create=False))
            for job_dir in job_dirs
        }
        hadd_condor_statuses = {
            job_dir: read_hadd_condor_status(job_dir) for job_dir in job_dirs
        }
        render_dashboard(
            args.root,
            condor_statuses=condor_statuses,
            hadd_condor_statuses=hadd_condor_statuses,
        )
    else:
        watch_dashboard(args.root, interval=args.interval, auto=not args.no_auto)


if __name__ == "__main__":
    main()

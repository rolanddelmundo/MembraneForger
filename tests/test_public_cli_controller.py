from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_cli(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    return subprocess.run(
        ["bash", str(ROOT / "run_pipeline.sh"), *args],
        cwd=cwd or ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=120,
    )


def test_help_lists_public_subcommands() -> None:
    result = run_cli("--help")
    assert result.returncode == 0
    for name in ("doctor", "init", "validate", "dry-run", "run", "resume", "stage", "status", "clean"):
        assert name in result.stdout


def test_init_creates_run_contract_with_path_spaces(tmp_path: Path) -> None:
    run_dir = tmp_path / "run with spaces"
    result = run_cli(
        "init",
        "--pdb",
        str(ROOT / "examples" / "minimal" / "inputs" / "minimal.pdb"),
        "--output-dir",
        str(run_dir),
    )
    assert result.returncode == 0, result.stderr
    assert (run_dir / "input" / "original.pdb").is_file()
    assert (run_dir / "config" / "config.yaml").is_file()
    status = json.loads((run_dir / "status.json").read_text(encoding="utf-8"))
    assert status["stages"]["stage1"]["state"] == "READY"
    provenance = json.loads((run_dir / "provenance" / "run.json").read_text(encoding="utf-8"))
    assert provenance["input"]["sha256"]


def test_doctor_json_reports_blocking_mstool_without_traceback() -> None:
    result = run_cli("doctor", "--json")
    assert result.returncode in {0, 1}
    rows = json.loads(result.stdout)
    assert any(row["name"] == "mstool" for row in rows)
    assert "Traceback" not in result.stderr


def test_clean_is_dry_run_by_default(tmp_path: Path) -> None:
    run_dir = ROOT / "runs" / "pytest_clean_dry_run"
    run_dir.mkdir(parents=True, exist_ok=True)
    result = run_cli("clean", "--run-dir", str(run_dir))
    assert result.returncode == 0, result.stderr
    assert run_dir.is_dir()
    assert "DRY RUN" in result.stdout
    run_dir.rmdir()


def test_clean_rejects_unsafe_confirm_targets(tmp_path: Path) -> None:
    targets = [
        ROOT,
        ROOT / "runs",
        ROOT.parent,
        ROOT / "membraneforger",
        tmp_path / "outside",
    ]
    (tmp_path / "outside").mkdir()
    for target in targets:
        result = run_cli("clean", "--run-dir", str(target), "--confirm")
        assert result.returncode == 2, target
        assert "refusing to clean" in result.stderr


def test_clean_confirm_removes_only_marked_run_directory() -> None:
    run_dir = ROOT / "runs" / "pytest_clean_valid"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / ".membraneforger-run.json").write_text(
        json.dumps({"schema_version": 1, "run_id": "pytest_clean_valid", "created_by": "MembraneForger"}) + "\n",
        encoding="utf-8",
    )
    (run_dir / "derived.txt").write_text("derived\n", encoding="utf-8")
    result = run_cli("clean", "--run-dir", str(run_dir), "--confirm")
    assert result.returncode == 0, result.stderr
    assert not run_dir.exists()


def test_clean_rejects_unmarked_and_symlink_run_targets(tmp_path: Path) -> None:
    unmarked = ROOT / "runs" / "pytest_clean_unmarked"
    unmarked.mkdir(parents=True, exist_ok=True)
    result = run_cli("clean", "--run-dir", str(unmarked), "--confirm")
    assert result.returncode == 2
    assert unmarked.is_dir()
    unmarked.rmdir()

    target = tmp_path / "target"
    target.mkdir()
    symlink = ROOT / "runs" / "pytest_clean_symlink"
    if symlink.exists() or symlink.is_symlink():
        symlink.unlink()
    symlink.symlink_to(target, target_is_directory=True)
    result = run_cli("clean", "--run-dir", str(symlink), "--confirm")
    assert result.returncode == 2
    assert target.is_dir()
    symlink.unlink()

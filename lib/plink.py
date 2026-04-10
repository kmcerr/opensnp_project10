"""
Project 10 — PLINK helper.

Provides a single run_plink() function used by all pipeline scripts
that call PLINK via subprocess.
"""

import subprocess
import sys
from pathlib import Path


def run_plink(cmd, label="", fatal=True):
    """
    Run a PLINK command via subprocess.

    - Prints informative lines from the PLINK .log file (variant/sample
      counts, genotyping rate, removals).
    - On non-zero return code, prints stderr and optionally exits.

    Parameters:
        cmd:    list of command-line arguments
        label:  human-readable label for console output
        fatal:  if True, exit on failure

    Returns:
        subprocess.CompletedProcess
    """
    print(f"\n{'─' * 60}")
    print(f"[{label}]")
    res = subprocess.run(cmd, capture_output=True, text=True)

    # Find the --out prefix to locate the .log file
    log_path = None
    for i, arg in enumerate(cmd):
        if arg == "--out" and i + 1 < len(cmd):
            log_path = Path(cmd[i + 1] + ".log")
            break

    if log_path and log_path.exists():
        log_keywords = [
            "variant", "sample", "people", "pass", "removed",
            "remaining", "loaded", "written", "marker", "snp",
            "total genotyping rate", "pruning complete",
        ]
        for line in log_path.read_text().strip().split("\n"):
            if any(kw in line.lower() for kw in log_keywords):
                print(f"  {line.strip()}")

    if res.returncode != 0:
        print(f"\n  FAILED (return code {res.returncode})")
        if res.stderr:
            print(f"  STDERR:\n{res.stderr[-2000:]}")
        if fatal:
            print("  Aborting pipeline.")
            sys.exit(1)
    else:
        print(f"  OK")

    return res

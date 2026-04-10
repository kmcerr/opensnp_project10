#!/usr/bin/env python3
"""
Run pipeline steps 01 → 07 in order.

Usage:
    python run_pipeline.py              # run all available steps
    python run_pipeline.py 7            # start from step 07
    python run_pipeline.py 3 6          # run steps 03 through 06

Environment variables:
    PROJECT10_LOG_LEVEL      Set log level (DEBUG, INFO, WARNING, ERROR)
    PROJECT10_DATA_DIR       Override data directory location
    PROJECT10_KEEP_INTERMEDIATE  Keep intermediate files (true/false)

Each step's main() is imported and called directly — no subprocesses.
"""

import importlib
import sys
import time
import logging

STEPS = [
    ("01", "01_create_mapping",      "Ancestry classification + corrections"),
    ("02", "02_ancestry_grouping",   "Map users to tier0/tier1"),
    ("03", "03_data_prep",           "Audit genotype files + filtering"),
    ("04", "04_build_verification",  "Verify genome builds + liftOver prep"),
    ("05", "05_plink_conversion",    "Convert to PLINK binary"),
    ("06", "06_merge",               "Build shared panel + merge"),
    ("07", "07_qc_ld_ibs",           "QC, LD pruning, pairwise IBS"),
    ("08", "08_visualizations",      "Generate presentation figures"),
    ("09", "09_network_analysis",    "Network construction + clustering analysis"),
]


def main():
    # Parse command line arguments
    start = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    end = int(sys.argv[2]) if len(sys.argv) > 2 else 9

    selected = [(n, mod, desc) for n, mod, desc in STEPS
                if start <= int(n) <= end]

    if not selected:
        print(f"No steps in range {start}–{end}")
        return 1

    # Setup logging and configuration
    import config
    from lib.validation import validate_config

    config.setup_logging()
    logger = logging.getLogger(__name__)

    logger.info("=" * 72)
    logger.info("OpenSNP Project 10 - Genetic Similarity Network Pipeline")
    logger.info("=" * 72)
    logger.info(f"Data directory: {config.DATA_DIR}")
    logger.info(f"Steps: {selected[0][0]} → {selected[-1][0]}")
    logger.info(f"Log level: {config.LOG_LEVEL}")

    # Validate configuration
    if not validate_config(config):
        logger.error("Configuration validation failed. Cannot proceed.")
        return 1

    logger.info("")

    # Execute pipeline steps
    for num, module_name, description in selected:
        logger.info("=" * 72)
        logger.info(f"Step {num}: {description}")
        logger.info("=" * 72)

        try:
            mod = importlib.import_module(module_name)
        except ModuleNotFoundError:
            logger.warning(f"⚠️  {module_name}.py not found — skipping\n")
            continue
        except Exception as e:
            logger.error(f"❌ Failed to import {module_name}: {e}")
            return 1

        t0 = time.time()
        try:
            mod.main()
            elapsed = time.time() - t0
            logger.info(f"✓ Step {num} completed in {elapsed:.1f}s\n")
        except Exception as e:
            elapsed = time.time() - t0
            logger.error(f"❌ Step {num} failed after {elapsed:.1f}s: {e}")
            logger.exception("Full traceback:")
            return 1

    logger.info("=" * 72)
    logger.info("✓ Pipeline complete")
    logger.info("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

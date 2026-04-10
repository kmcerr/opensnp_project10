"""
Validation utilities for data quality checks and configuration validation.
"""

import logging
import shutil
from pathlib import Path
from typing import List, Optional, Set
import pandas as pd

logger = logging.getLogger(__name__)


def validate_dataframe_columns(
    df: pd.DataFrame,
    required_cols: List[str],
    context: str = "DataFrame"
) -> None:
    """
    Validate that a DataFrame contains required columns.

    Parameters:
        df: DataFrame to validate
        required_cols: List of required column names
        context: Description for error message

    Raises:
        ValueError: If any required columns are missing
    """
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(
            f"{context} missing required columns: {missing}. "
            f"Available: {list(df.columns)}"
        )


def validate_ancestry_data(final: pd.DataFrame, expected_tier0: Set[str]) -> bool:
    """
    Run data quality checks on ancestry groupings DataFrame.

    Parameters:
        final: DataFrame with ancestry assignments
        expected_tier0: Set of valid tier0 values

    Returns:
        True if all checks pass, False otherwise
    """
    checks_passed = []
    checks_failed = []

    # Check 1: No duplicate user_ids
    if final['user_id'].duplicated().any():
        dup_count = final['user_id'].duplicated().sum()
        checks_failed.append(f"Found {dup_count} duplicate user_ids")
    else:
        checks_passed.append("✓ No duplicate user_ids")

    # Check 2: All tier0 values are valid
    invalid_tier0 = set(final['tier0'].dropna()) - expected_tier0
    if invalid_tier0:
        checks_failed.append(f"Invalid tier0 values: {invalid_tier0}")
    else:
        checks_passed.append("✓ All tier0 values valid")

    # Check 3: No missing values in required columns
    required = ['user_id', 'raw_ancestry', 'tier0', 'tier1']
    missing = final[required].isna().sum()
    if missing.any():
        checks_failed.append(
            f"Missing values in required columns:\n{missing[missing > 0]}"
        )
    else:
        checks_passed.append("✓ No missing values in required columns")

    # Print results
    logger.info("=" * 60)
    logger.info("DATA VALIDATION")
    logger.info("=" * 60)
    for check in checks_passed:
        logger.info(check)
    for check in checks_failed:
        logger.error(f"❌ {check}")

    return len(checks_failed) == 0


def validate_config(config_module) -> bool:
    """
    Validate configuration settings before running pipeline.

    Parameters:
        config_module: The config module to validate

    Returns:
        True if configuration is valid, False otherwise
    """
    errors = []
    warnings = []

    # Check PLINK is available
    if shutil.which(config_module.PLINK_EXEC) is None:
        errors.append(
            f"PLINK executable '{config_module.PLINK_EXEC}' not found on PATH. "
            f"Install PLINK 1.9 or set PLINK_EXEC environment variable."
        )

    # Validate thresholds
    if not 0 < config_module.SNP_PRESENCE_THRESHOLD <= 1.0:
        errors.append(
            f"SNP_PRESENCE_THRESHOLD must be (0,1], "
            f"got {config_module.SNP_PRESENCE_THRESHOLD}"
        )

    if not 0 < config_module.SAMPLE_MISS_THRESHOLD < 1.0:
        errors.append(
            f"SAMPLE_MISS_THRESHOLD must be (0,1), "
            f"got {config_module.SAMPLE_MISS_THRESHOLD}"
        )

    if not 0 < config_module.SNP_MISS_THRESHOLD < 1.0:
        errors.append(
            f"SNP_MISS_THRESHOLD must be (0,1), "
            f"got {config_module.SNP_MISS_THRESHOLD}"
        )

    if not 0 < config_module.MAF_THRESHOLD < 0.5:
        errors.append(
            f"MAF_THRESHOLD must be (0,0.5), got {config_module.MAF_THRESHOLD}"
        )

    # Check thread count
    if config_module.N_THREADS < 1:
        errors.append(f"N_THREADS must be >= 1, got {config_module.N_THREADS}")
    elif config_module.N_THREADS > 32:
        warnings.append(
            f"N_THREADS={config_module.N_THREADS} is very high. "
            f"Most systems perform best with 4-8 threads."
        )

    # Check required inputs exist (only for step 1+)
    if config_module.RAW_ANCESTRY_CSV.exists():
        logger.info(f"✓ Found input: {config_module.RAW_ANCESTRY_CSV}")
    else:
        warnings.append(
            f"Input file not found: {config_module.RAW_ANCESTRY_CSV}. "
            f"Required for step 01. See data/README.md for setup instructions."
        )

    # Validate LD pruning parameters
    if config_module.LD_R2 < 0 or config_module.LD_R2 > 1:
        errors.append(f"LD_R2 must be [0,1], got {config_module.LD_R2}")

    # Print results
    if errors:
        logger.error("❌ Configuration errors:")
        for e in errors:
            logger.error(f"  - {e}")

    if warnings:
        logger.warning("⚠️  Configuration warnings:")
        for w in warnings:
            logger.warning(f"  - {w}")

    if not errors and not warnings:
        logger.info("✓ Configuration validated successfully")
    elif not errors:
        logger.info("✓ Configuration valid (with warnings)")

    return len(errors) == 0


def validate_file_readable(filepath: Path, context: str = "File") -> bool:
    """
    Check if a file exists and is readable.

    Parameters:
        filepath: Path to file
        context: Description for error message

    Returns:
        True if file is readable, False otherwise
    """
    if not filepath.exists():
        logger.error(f"{context} not found: {filepath}")
        return False

    if not filepath.is_file():
        logger.error(f"{context} is not a regular file: {filepath}")
        return False

    try:
        with open(filepath, 'r') as f:
            f.read(1)
        return True
    except PermissionError:
        logger.error(f"{context} is not readable: {filepath}")
        return False
    except Exception as e:
        logger.error(f"{context} cannot be read: {filepath} ({e})")
        return False

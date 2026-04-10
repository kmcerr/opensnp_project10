"""
Project 10 — IBS computation and annotation helpers.

Consolidated from qc_ldpruning_ibs_v2, networkjoin_v2, step7_dedup_v2,
and step9_sensitivity.
"""

import pandas as pd
from pathlib import Path

from config import (
    GROUPINGS_CSV, STAGE4_INPUT_CSV,
    PIHAT_DUPLICATE, PIHAT_1ST_DEGREE, PIHAT_2ND_DEGREE, PIHAT_3RD_DEGREE,
)


# ═══════════════════════════════════════════════════════════════════
# IBS PROPORTIONS
# ═══════════════════════════════════════════════════════════════════

def compute_ibs_proportions(genome_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add IBS0_prop, IBS1_prop, IBS2_prop columns to a PLINK .genome
    DataFrame.

    PLINK 1.9 --genome full outputs raw COUNTS. This function
    converts them to proportions by dividing by total compared loci.

    Modifies genome_df in place and returns it.
    """
    for col in ("IBS0", "IBS1", "IBS2"):
        if col not in genome_df.columns:
            raise ValueError(
                f"Column '{col}' not found. Did you use --genome full?"
            )

    total = genome_df["IBS0"] + genome_df["IBS1"] + genome_df["IBS2"]

    # Guard against division by zero (no compared loci)
    n_zero = (total == 0).sum()
    if n_zero > 0:
        print(f"  WARNING: {n_zero} pair(s) with zero compared loci — "
              f"IBS proportions will be NaN for these.")
    total = total.replace(0, float("nan"))

    # Diagnostic: are values counts or already proportions?
    median_total = total.median()
    if median_total < 10:
        print(f"  NOTE: IBS values appear to be proportions "
              f"(median sum = {median_total:.4f}), not counts.")
    else:
        print(f"  IBS values are counts (median total loci: "
              f"{median_total:,.0f})")

    genome_df["IBS0_prop"] = genome_df["IBS0"] / total
    genome_df["IBS1_prop"] = genome_df["IBS1"] / total
    genome_df["IBS2_prop"] = genome_df["IBS2"] / total

    return genome_df


# ═══════════════════════════════════════════════════════════════════
# PAIR CATEGORIZATION
# ═══════════════════════════════════════════════════════════════════

def pair_category(pi_hat: float) -> str:
    """
    Classify a pairwise relationship degree from PI_HAT.

    Standard thresholds:
        >=0.9  → duplicate / MZ twin
        >=0.4  → first-degree (parent-child, full siblings)
        >=0.25 → second-degree (half-siblings, grandparent, avuncular)
        >=0.125 → third-degree (first cousins)
        else   → unrelated
    """
    if pi_hat >= PIHAT_DUPLICATE:
        return "duplicate_or_mz_twin"
    elif pi_hat >= PIHAT_1ST_DEGREE:
        return "first_degree"
    elif pi_hat >= PIHAT_2ND_DEGREE:
        return "second_degree"
    elif pi_hat >= PIHAT_3RD_DEGREE:
        return "third_degree"
    else:
        return "unrelated"


def add_pair_categories(df: pd.DataFrame) -> pd.DataFrame:
    """Add pair_category column to a DataFrame with PI_HAT."""
    df["pair_category"] = df["PI_HAT"].apply(pair_category)
    return df


# ═══════════════════════════════════════════════════════════════════
# METADATA LOADING
# ═══════════════════════════════════════════════════════════════════

def load_metadata(valid_ids: set = None,
                  groupings_path: Path = None,
                  stage4_path: Path = None) -> pd.DataFrame:
    """
    Load and merge ancestry (tier0/tier1/raw_ancestry) with genotype
    format metadata.

    Parameters:
        valid_ids:      if provided, restrict to these user_ids
        groupings_path: override path to groupings CSV
        stage4_path:    override path to stage4 manifest CSV

    Returns:
        DataFrame with columns: user_id, raw_ancestry, tier0, tier1,
        genotype_format
    """
    g_path = groupings_path or GROUPINGS_CSV
    s4_path = stage4_path or STAGE4_INPUT_CSV

    groupings = pd.read_csv(g_path)
    groupings["user_id"] = groupings["user_id"].astype(str)
    groupings_slim = groupings[
        ["user_id", "raw_ancestry", "tier0", "tier1"]
    ].copy()

    stage4 = pd.read_csv(s4_path)
    stage4["user_id"] = stage4["user_id"].astype(str)
    format_map = stage4[["user_id", "genotype_format"]].drop_duplicates()

    meta = groupings_slim.merge(format_map, on="user_id", how="left")

    if valid_ids is not None:
        valid_str = {str(x) for x in valid_ids}
        meta = meta[meta["user_id"].isin(valid_str)].copy()

    return meta


def load_fam_ids(fam_path: Path) -> set:
    """
    Read a PLINK .fam file and return the set of IID strings.
    """
    fam = pd.read_csv(
        fam_path, sep=r"\s+", header=None,
        names=["FID", "IID", "father", "mother", "sex", "pheno"],
    )
    fam["IID"] = fam["IID"].astype(str)
    return set(fam["IID"])


# ═══════════════════════════════════════════════════════════════════
# IBS ANNOTATION
# ═══════════════════════════════════════════════════════════════════

def annotate_ibs(ibs_df: pd.DataFrame,
                 meta_df: pd.DataFrame) -> pd.DataFrame:
    """
    Join IBS pairwise results with metadata for both individuals.

    Adds: tier0_1, tier1_1, raw_ancestry_1, genotype_format_1,
          tier0_2, tier1_2, raw_ancestry_2, genotype_format_2,
          same_tier0, same_tier1, same_format, pair_category.

    Parameters:
        ibs_df:  DataFrame with IID1, IID2, PI_HAT, IBS columns
        meta_df: DataFrame with user_id, tier0, tier1, raw_ancestry,
                 genotype_format

    Returns:
        Annotated DataFrame (same row count as ibs_df).
    """
    # Ensure string IDs in IBS dataframe
    for col in ("FID1", "IID1", "FID2", "IID2"):
        if col in ibs_df.columns:
            ibs_df[col] = ibs_df[col].astype(str)

    # Convert user_id to string in metadata before renaming
    meta_df = meta_df.copy()
    meta_df["user_id"] = meta_df["user_id"].astype(str)

    # Build rename dictionary for available columns only
    rename_map_1 = {"user_id": "IID1"}
    rename_map_2 = {"user_id": "IID2"}

    for col, suffix in [("tier0", "_1"), ("tier1", "_1"), ("raw_ancestry", "_1"),
                        ("genotype_format", "_1")]:
        if col in meta_df.columns:
            rename_map_1[col] = col + suffix

    for col, suffix in [("tier0", "_2"), ("tier1", "_2"), ("raw_ancestry", "_2"),
                        ("genotype_format", "_2")]:
        if col in meta_df.columns:
            rename_map_2[col] = col + suffix

    meta1 = meta_df.rename(columns=rename_map_1)
    meta2 = meta_df.rename(columns=rename_map_2)

    annot = ibs_df.merge(meta1, on="IID1", how="left") \
                   .merge(meta2, on="IID2", how="left")

    # Add comparison columns only if source columns exist
    if "tier0_1" in annot.columns and "tier0_2" in annot.columns:
        annot["same_tier0"] = annot["tier0_1"] == annot["tier0_2"]
    if "tier1_1" in annot.columns and "tier1_2" in annot.columns:
        annot["same_tier1"] = annot["tier1_1"] == annot["tier1_2"]
    if "genotype_format_1" in annot.columns and "genotype_format_2" in annot.columns:
        annot["same_format"] = annot["genotype_format_1"] == annot["genotype_format_2"]

    annot = add_pair_categories(annot)

    return annot


def load_genome_file(genome_path: Path) -> pd.DataFrame:
    """
    Load a PLINK .genome file, compute IBS proportions, and
    standardize ID columns to strings.
    """
    df = pd.read_csv(genome_path, sep=r"\s+")
    for col in ("FID1", "IID1", "FID2", "IID2"):
        df[col] = df[col].astype(str)
    df = compute_ibs_proportions(df)
    return df

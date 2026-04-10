import pandas as pd
from pathlib import Path

from config import (
    DATA_DIR, PIPELINE_MANIFEST_CSV, CONVERSION_READY_CSV,
    LIFTOVER_DIR, LIFTOVER_BED_DIR, LIFTOVER_OUT_DIR,
    LIFTOVER_UNMAPPED_DIR, LIFTOVER_EXEC, CHAIN_FILE,
    SENTINELS, AUTOSOMAL_CHROMS, ensure_dirs,
)
from lib.parsing import (
    split_flexible, parse_candidate_record, scan_file_for_positions,
    iter_genotype_records,
)

SENTINEL_RSIDS = set(SENTINELS.keys())


def extract_grch36_bed(filepath: Path, out_bed_path: Path):
    """Write autosomal rsID loci as 0-based BED for UCSC liftOver."""
    n_written = 0
    seen = set()

    with open(out_bed_path, "w") as out:
        for fields, _ in iter_genotype_records(filepath):
            rec = parse_candidate_record(fields)
            if rec is None:
                continue
            rsid, chrom, pos = rec
            if chrom not in AUTOSOMAL_CHROMS:
                continue
            key = (rsid, chrom, pos)
            if key in seen:
                continue
            seen.add(key)
            out.write(f"chr{chrom}\t{pos - 1}\t{pos}\t{rsid}\n")
            n_written += 1

    return n_written





def main():
    ensure_dirs()
    print(f"Pipeline manifest: {PIPELINE_MANIFEST_CSV}")
    print(f"liftOver dir:      {LIFTOVER_DIR}")

    pipeline = pd.read_csv(PIPELINE_MANIFEST_CSV)
    print(f"Pipeline rows: {len(pipeline)}")
    print("\nBy build_hint:")
    print(pipeline["build_hint"].value_counts(dropna=False).to_string())

    build_unknown = pipeline[pipeline["build_hint"].isna()].copy()
    grch36 = pipeline[pipeline["build_hint"] == "GRCh36"].copy()
    grch37 = pipeline[pipeline["build_hint"] == "GRCh37"].copy()

    print(f"\nGRCh37 known:  {len(grch37)}")
    print(f"GRCh36 known:  {len(grch36)}")
    print(f"Build unknown: {len(build_unknown)}")

    verification_rows = []

    for _, row in build_unknown.iterrows():
        user_id = row["user_id"]
        filepath = Path(row["filepath"])

        hits = scan_file_for_positions(filepath, SENTINEL_RSIDS)

        score37 = score36 = 0
        details = []

        for rsid, observed_list in hits.items():
            for chrom, pos in observed_list:
                expected_chr = SENTINELS[rsid]["chr"]
                g37 = SENTINELS[rsid]["GRCh37"]
                g36 = SENTINELS[rsid]["GRCh36"]

                if chrom == expected_chr and pos == g37:
                    score37 += 1
                    details.append(f"{rsid}: GRCh37")
                elif chrom == expected_chr and pos == g36:
                    score36 += 1
                    details.append(f"{rsid}: GRCh36")
                else:
                    details.append(f"{rsid}: chr{chrom}:{pos}=other")

        if score37 > 0 and score36 == 0:
            inferred_build = "GRCh37"
        elif score36 > 0 and score37 == 0:
            inferred_build = "GRCh36"
        elif score37 == 0 and score36 == 0:
            inferred_build = "Unresolved"
        else:
            inferred_build = "Conflicted"

        verification_rows.append({
            "user_id": user_id,
            "genotype_filename": row["genotype_filename"],
            "filepath": str(filepath),
            "score37": score37, "score36": score36,
            "inferred_build": inferred_build,
            "details": " | ".join(details),
        })

    verification_df = pd.DataFrame(verification_rows)
    print(f"Verified {len(verification_df)} build-unknown files")
    print(verification_df["inferred_build"].value_counts(dropna=False).to_string())

    build_verification_path = DATA_DIR / "build_verification_results.csv"
    verification_df.to_csv(build_verification_path, index=False)

    build_unknown_review = verification_df[
        verification_df["inferred_build"].isin(["Unresolved", "Conflicted"])
    ]
    build_unknown_review_path = DATA_DIR / "build_unknown_review.csv"
    build_unknown_review.to_csv(build_unknown_review_path, index=False)

    print(f"Saved build verification results to {build_verification_path}")
    print(f"Saved manual review list ({len(build_unknown_review)}) to "
          f"{build_unknown_review_path}")

    pipeline2 = pipeline.merge(
        verification_df[["user_id", "inferred_build"]],
        on="user_id", how="left",
    )

    pipeline2["resolved_build"] = pipeline2["build_hint"]
    mask_fill = pipeline2["resolved_build"].isna() & pipeline2["inferred_build"].notna()
    pipeline2.loc[mask_fill, "resolved_build"] = pipeline2.loc[mask_fill, "inferred_build"]

    print("Resolved build counts:")
    print(pipeline2["resolved_build"].value_counts(dropna=False).to_string())

    grch36_resolved = pipeline2[pipeline2["resolved_build"] == "GRCh36"].copy()
    liftover_rows = []

    for _, row in grch36_resolved.iterrows():
        user_id = row["user_id"]
        filepath = Path(row["filepath"])

        bed_in = LIFTOVER_BED_DIR / f"user_{user_id}.grch36.bed"
        bed_out = LIFTOVER_OUT_DIR / f"user_{user_id}.grch37.bed"
        unmapped = LIFTOVER_UNMAPPED_DIR / f"user_{user_id}.unmapped.bed"

        n_loci = extract_grch36_bed(filepath, bed_in)

        liftover_rows.append({
            "user_id": user_id,
            "genotype_filename": row["genotype_filename"],
            "filepath": str(filepath),
            "bed_in": str(bed_in), "bed_out": str(bed_out),
            "unmapped_out": str(unmapped), "n_loci": n_loci,
        })

    liftover_df = pd.DataFrame(liftover_rows)
    print(f"Prepared liftOver BEDs for {len(liftover_df)} files")

    liftover_manifest_path = DATA_DIR / "liftover_input_manifest.csv"
    liftover_df.to_csv(liftover_manifest_path, index=False)

    # Write shell commands
    liftover_cmds_path = DATA_DIR / "liftover_commands.sh"
    with open(liftover_cmds_path, "w") as f:
        f.write("#!/usr/bin/env bash\nset -euo pipefail\n\n")
        f.write(f'CHAIN="{CHAIN_FILE}"\nLIFTOVER="{LIFTOVER_EXEC}"\n\n')
        for _, row in liftover_df.iterrows():
            f.write(f'"$LIFTOVER" "{row["bed_in"]}" "$CHAIN" '
                    f'"{row["bed_out"]}" "{row["unmapped_out"]}"\n')

    print(f"Saved liftOver shell commands to {liftover_cmds_path}")

    conversion_ready = pipeline2[
        pipeline2["resolved_build"] == "GRCh37"
    ].copy()
    holdout_unresolved = pipeline2[
        pipeline2["resolved_build"].isin(["Unresolved", "Conflicted"])
    ].copy()
    holdout_grch36 = pipeline2[
        pipeline2["resolved_build"] == "GRCh36"
    ].copy()

    print(f"Conversion-ready now (GRCh37): {len(conversion_ready)}")
    print(f"Held out for liftOver (GRCh36): {len(holdout_grch36)}")
    print(f"Held out unresolved build:      {len(holdout_unresolved)}")

    conversion_ready.to_csv(CONVERSION_READY_CSV, index=False)
    holdout_grch36.to_csv(DATA_DIR / "holdout_grch36_manifest.csv", index=False)
    holdout_unresolved.to_csv(
        DATA_DIR / "holdout_unresolved_build_manifest.csv", index=False
    )

    print(f"\nSaved conversion-ready manifest to {CONVERSION_READY_CSV}")

    print(f"\nConversion-ready by genotype_format:")
    print(conversion_ready["genotype_format"].value_counts().to_string())

    print(f"\nConversion-ready by tier0:")
    print(conversion_ready["tier0"].value_counts().to_string())

    print("\n" + "=" * 72)
    print("SCRIPT 04 SUMMARY")
    print("=" * 72)
    print(f"\nInput pipeline manifest rows: {len(pipeline)}")
    print(f"Build-unknown files audited:  {len(build_unknown)}")
    print(f"GRCh36 files prepared:        {len(grch36_resolved)}")
    print(f"Conversion-ready now:         {len(conversion_ready)}")
    print(f"Held for liftOver:            {len(holdout_grch36)}")
    print(f"Held unresolved:              {len(holdout_unresolved)}")

    print("\nGenerated files:")
    for p in [
        build_verification_path, build_unknown_review_path,
        liftover_manifest_path, liftover_cmds_path,
        CONVERSION_READY_CSV,
        DATA_DIR / "holdout_grch36_manifest.csv",
        DATA_DIR / "holdout_unresolved_build_manifest.csv",
    ]:
        print(f"  - {p}")



if __name__ == "__main__":
    main()

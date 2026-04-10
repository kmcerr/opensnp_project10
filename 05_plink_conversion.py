import subprocess
import pandas as pd
from pathlib import Path

from config import (
    DATA_DIR, CONVERSION_READY_CSV, STAGE4_INPUT_CSV,
    PLINK_INDIVIDUAL_DIR, PLINK_EXEC, MIN_PARSED_SNPS,
    ensure_dirs,
)
from lib.parsing import parse_and_write_cleaned

def main():
    ensure_dirs()
    print(f"Conversion-ready manifest: {CONVERSION_READY_CSV}")
    print(f"PLINK output dir:          {PLINK_INDIVIDUAL_DIR}")

    manifest = pd.read_csv(CONVERSION_READY_CSV)
    print(f"Rows in conversion-ready manifest: {len(manifest)}")
    print("\nBy genotype format:")
    print(manifest["genotype_format"].value_counts().to_string())

    if "resolved_build" in manifest.columns:
        bad = manifest[manifest["resolved_build"] != "GRCh37"]
        if len(bad) > 0:
            print(f"WARNING: {len(bad)} non-GRCh37 rows present")

    conversion_logs = []
    failed_samples = []

    for i, row in manifest.iterrows():
        user_id = str(row["user_id"])
        filepath = Path(row["filepath"])
        fmt = row["format_detected"]
        out_prefix = PLINK_INDIVIDUAL_DIR / f"user_{user_id}"
        cleaned_path = PLINK_INDIVIDUAL_DIR / f"user_{user_id}_cleaned.txt"

        try:
            # Step 1: Parse and write cleaned 23andme-style file
            n_parsed, stats = parse_and_write_cleaned(filepath, fmt, cleaned_path)

            if n_parsed < MIN_PARSED_SNPS:
                failed_samples.append({
                    "user_id": user_id, "filepath": str(filepath),
                    "format_detected": fmt,
                    "reason": f"too_few_parsed_snps ({n_parsed})",
                    **stats,
                })
                if cleaned_path.exists():
                    cleaned_path.unlink()
                continue

            # Step 2: PLINK --23file import
            sex_raw = str(row.get("chrom_sex", "unknown")).upper()
            sex_flag = ["M"] if sex_raw == "XY" else (
                ["F"] if sex_raw == "XX" else []
            )

            cmd = [
                PLINK_EXEC,
                "--23file", str(cleaned_path), user_id, user_id,
            ] + sex_flag + [
                "--allow-no-sex",
                "--make-bed",
                "--out", str(out_prefix),
            ]

            res = subprocess.run(cmd, capture_output=True, text=True)

            if res.returncode != 0:
                failed_samples.append({
                    "user_id": user_id, "filepath": str(filepath),
                    "format_detected": fmt, "reason": "plink_23file_failed",
                    "stderr": res.stderr[:2000], **stats,
                })
                continue

            # Step 3: Clean up intermediate file
            if cleaned_path.exists():
                cleaned_path.unlink()

            conversion_logs.append({
                "user_id": user_id, "filepath": str(filepath),
                "genotype_filename": row["genotype_filename"],
                "genotype_format": row["genotype_format"],
                "format_detected": fmt,
                "tier0": row.get("tier0", ""),
                "tier1": row.get("tier1", ""),
                "parsed_snps": n_parsed,
                "bed_prefix": str(out_prefix),
                **stats,
            })

            if len(conversion_logs) % 25 == 0:
                print(f"  Converted {len(conversion_logs)} samples...")

        except Exception as e:
            failed_samples.append({
                "user_id": user_id, "filepath": str(filepath),
                "format_detected": fmt, "reason": f"exception: {e}",
            })

    print(f"\nSuccessful conversions: {len(conversion_logs)}")
    print(f"Failed conversions:     {len(failed_samples)}")

    conversion_df = pd.DataFrame(conversion_logs)
    failed_df = pd.DataFrame(failed_samples)

    if len(conversion_df) > 0:
        print("\nBy detected format:")
        print(conversion_df["format_detected"].value_counts().to_string())
        print(f"\nParsed SNPs: median={conversion_df['parsed_snps'].median():,.0f}, "
              f"min={conversion_df['parsed_snps'].min():,}, "
              f"max={conversion_df['parsed_snps'].max():,}")

    conversion_df.to_csv(DATA_DIR / "stage3_conversion_log.csv", index=False)
    failed_df.to_csv(DATA_DIR / "stage3_failed_samples.csv", index=False)

    stage4_cols = ["user_id", "genotype_filename", "genotype_format",
                   "format_detected", "tier0", "tier1", "parsed_snps",
                   "bed_prefix"]
    stage4_manifest = conversion_df[stage4_cols].copy()
    stage4_manifest.to_csv(STAGE4_INPUT_CSV, index=False)
    print(f"Saved Stage 4 input manifest ({len(stage4_manifest)} samples) to "
          f"{STAGE4_INPUT_CSV}")

    print("\n" + "=" * 60)
    print("SCRIPT 05 SUMMARY")
    print("=" * 60)
    print(f"Input manifest rows:    {len(manifest)}")
    print(f"Successful conversions: {len(conversion_df)}")
    print(f"Failed conversions:     {len(failed_df)}")
    if len(failed_df) > 0:
        print("\nFailed samples:")
        for _, r in failed_df.iterrows():
            print(f"  user_id={r['user_id']}: {r.get('reason', 'unknown')}")



if __name__ == "__main__":
    main()

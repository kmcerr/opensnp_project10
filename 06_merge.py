import subprocess
import sys
import pandas as pd
from pathlib import Path
from collections import Counter

from config import (
    DATA_DIR, STAGE4_INPUT_CSV, MERGE_DIR, FILTERED_DIR,
    MERGED_PREFIX, PLINK_EXEC, SNP_PRESENCE_THRESHOLD,
    MERGE_MAX_ATTEMPTS, ensure_dirs,
)
from lib.plink import run_plink

def write_snp_list(snps, path):
    with open(path, "w") as f:
        for rsid in sorted(snps):
            f.write(rsid + "\n")
    return path




def filter_all_samples(manifest_df, snp_list_path, output_dir, label="filt"):
    rows = []
    for _, row in manifest_df.iterrows():
        uid = row["user_id"]
        in_prefix = Path(str(row["bed_prefix"]))
        out_prefix = output_dir / f"user_{uid}_{label}"

        res = subprocess.run([
            PLINK_EXEC,
            "--bfile", str(in_prefix),
            "--allow-no-sex",
            "--extract", str(snp_list_path),
            "--make-bed",
            "--out", str(out_prefix),
        ], capture_output=True, text=True)

        rows.append({
            "user_id": uid,
            "out_prefix": str(out_prefix),
            "returncode": res.returncode,
        })

        if res.returncode != 0:
            print(f"    PLINK filter failed for user_{uid}: "
                  f"rc={res.returncode}")
            if res.stderr:
                print(f"      {res.stderr.strip()[:200]}")

    df = pd.DataFrame(rows)
    n_ok = (df["returncode"] == 0).sum()
    n_fail = (df["returncode"] != 0).sum()
    print(f"  Filtered OK: {n_ok}, Failed: {n_fail}")
    return df[df["returncode"] == 0].copy()




def attempt_merge(filtered_df, out_prefix, label="merge"):
    base = Path(filtered_df.iloc[0]["out_prefix"])
    merge_list = MERGE_DIR / f"{label}_merge_list.txt"

    with open(merge_list, "w") as f:
        for _, row in filtered_df.iloc[1:].iterrows():
            p = row["out_prefix"]
            # PLINK 1.9 --merge-list expects just the prefix (no extensions)
            f.write(f"{p}\n")

    run_plink([
        PLINK_EXEC,
        "--bfile", str(base),
        "--merge-list", str(merge_list),
        "--allow-no-sex",
        "--make-bed",
        "--out", str(out_prefix),
    ], label=label, fatal=False)

    missnp = Path(str(out_prefix) + "-merge.missnp")
    success = all(out_prefix.with_suffix(e).exists()
                  for e in [".bed", ".bim", ".fam"])
    return success, missnp if missnp.exists() else None





def main():
    ensure_dirs()
    print(f"Stage 4 input manifest: {STAGE4_INPUT_CSV}")
    print(f"Merge output dir:       {MERGE_DIR}")
    print(f"SNP presence threshold: {SNP_PRESENCE_THRESHOLD}")

    manifest = pd.read_csv(STAGE4_INPUT_CSV)
    print(f"Converted samples available: {len(manifest)}")
    print("\nBy genotype_format:")
    print(manifest["genotype_format"].value_counts().to_string())
    print("\nBy tier0:")
    print(manifest["tier0"].value_counts().to_string())

    # Verify BIM files exist
    missing_bims = []
    for _, row in manifest.iterrows():
        if not Path(str(row["bed_prefix"]) + ".bim").exists():
            missing_bims.append(row["user_id"])

    if missing_bims:
        print(f"WARNING: {len(missing_bims)} BIM files not found")
        manifest = manifest[~manifest["user_id"].isin(missing_bims)].copy()
        print(f"Proceeding with {len(manifest)} samples")

    snp_counts = Counter()
    n_samples = len(manifest)

    for idx, row in manifest.iterrows():
        bim_path = Path(str(row["bed_prefix"]) + ".bim")
        with open(bim_path) as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 2:
                    snp_counts[fields[1]] += 1
        if (idx + 1) % 50 == 0:
            print(f"  Scanned {idx + 1}/{n_samples} BIM files...")

    print(f"\nUnique SNPs: {len(snp_counts):,}")

    min_presence = int(round(n_samples * SNP_PRESENCE_THRESHOLD))
    shared_snps = {rsid for rsid, cnt in snp_counts.items()
                   if cnt >= min_presence}
    print(f"SNPs in ≥{min_presence}/{n_samples} samples "
          f"({SNP_PRESENCE_THRESHOLD:.0%}): {len(shared_snps):,}")

    merge_out = MERGED_PREFIX
    current_panel = shared_snps.copy()
    excluded_snps = set()

    for attempt in range(1, MERGE_MAX_ATTEMPTS + 1):
        print("\n" + "=" * 60)
        print(f"MERGE ATTEMPT {attempt}")
        print("=" * 60)

        panel_path = MERGE_DIR / f"snp_panel_v{attempt}.txt"
        write_snp_list(current_panel, panel_path)
        print(f"Panel size: {len(current_panel):,} SNPs")

        print(f"\nFiltering all samples to panel v{attempt}...")
        filtered = filter_all_samples(
            manifest, panel_path, FILTERED_DIR, label=f"v{attempt}"
        )

        if len(filtered) < 2:
            print("FATAL: fewer than 2 samples passed filtering")
            sys.exit(1)

        success, missnp_path = attempt_merge(
            filtered, merge_out, label=f"merge_v{attempt}"
        )

        if success:
            print(f"\nMerge succeeded on attempt {attempt}")
            break

        if missnp_path:
            n_new = sum(1 for _ in open(missnp_path))
            print(f"\n{n_new:,} SNPs in .missnp — excluding for next attempt")
            with open(missnp_path) as f:
                for line in f:
                    s = line.strip()
                    if s:
                        excluded_snps.add(s)
            current_panel = shared_snps - excluded_snps
            print(f"Cumulative excluded: {len(excluded_snps):,}")
        else:
            print("\nMerge failed without .missnp — cannot auto-recover")
            sys.exit(1)
    else:
        print(f"\nFATAL: Merge failed after {MERGE_MAX_ATTEMPTS} attempts")
        sys.exit(1)

    filtered.to_csv(MERGE_DIR / "final_filtering_log.csv", index=False)

    n_snps = sum(1 for _ in open(merge_out.with_suffix(".bim")))
    n_merged = sum(1 for _ in open(merge_out.with_suffix(".fam")))

    print("=" * 60)
    print("MERGE SUMMARY")
    print("=" * 60)
    print(f"  Samples:        {n_merged}")
    print(f"  SNPs:           {n_snps:,}")
    print(f"  Merge attempts: {attempt}")
    print(f"  SNPs excluded:  {len(excluded_snps):,}")
    print(f"  Output:         {merge_out}")

    manifest.to_csv(MERGE_DIR / "merged_sample_manifest.csv", index=False)
    print(f"\nTier0 distribution:")
    print(manifest["tier0"].value_counts().to_string())



if __name__ == "__main__":
    main()

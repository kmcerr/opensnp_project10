import sys
import pandas as pd
from pathlib import Path

from config import (
    MERGED_PREFIX, QC_DIR, PLINK_EXEC, N_THREADS,
    QC_MISSING_PREFIX, QC_STEP1_PREFIX, QC_STEP2_PREFIX,
    QC_STEP3_PREFIX, QC_STEP3B_PREFIX, QC_LD_PRUNE_PREFIX,
    QC_PRUNED_PREFIX, QC_IBS_PREFIX,
    SAMPLE_MISS_THRESHOLD, SNP_MISS_THRESHOLD, MAF_THRESHOLD,
    LD_WINDOW, LD_STEP, LD_R2,
    IBS2_RELATEDNESS_THRESHOLD, PIHAT_3RD_DEGREE,
    PIHAT_2ND_DEGREE, PIHAT_1ST_DEGREE,
    ensure_dirs,
)
from lib.plink import run_plink
from lib.ibs import compute_ibs_proportions

def main():
    ensure_dirs()
    print(f"Merged input prefix: {MERGED_PREFIX}")
    print(f"QC output dir:       {QC_DIR}")
    print(f"Threads for IBS:     {N_THREADS}")

    # Verify merged files exist
    for ext in (".bed", ".bim", ".fam"):
        path = Path(str(MERGED_PREFIX) + ext)
        if not path.exists():
            print(f"FATAL: {path.name} MISSING. Cannot proceed.")
            sys.exit(1)

    n_snps = sum(1 for _ in open(Path(str(MERGED_PREFIX) + ".bim")))
    n_samples = sum(1 for _ in open(Path(str(MERGED_PREFIX) + ".fam")))
    print(f"\nInitial merged dataset: {n_samples} samples, {n_snps:,} SNPs")

    run_plink([
        PLINK_EXEC, "--bfile", str(MERGED_PREFIX),
        "--autosome", "--allow-no-sex", "--missing",
        "--out", str(QC_MISSING_PREFIX),
    ], label="Missingness reports")

    imiss = pd.read_csv(Path(str(QC_MISSING_PREFIX) + ".imiss"), sep=r"\s+")
    lmiss = pd.read_csv(Path(str(QC_MISSING_PREFIX) + ".lmiss"), sep=r"\s+")

    print("Sample missingness summary (F_MISS):")
    print(imiss["F_MISS"].describe().round(4).to_string())

    print("\nSNP missingness summary (F_MISS):")
    print(lmiss["F_MISS"].describe().round(4).to_string())

    bad_samples = imiss[imiss["F_MISS"] > SAMPLE_MISS_THRESHOLD].copy()
    bad_samples_path = QC_DIR / "samples_fail_missingness.txt"

    if len(bad_samples) > 0:
        bad_samples[["FID", "IID"]].to_csv(
            bad_samples_path, sep="\t", header=False, index=False)
        print(f"Samples failing >{SAMPLE_MISS_THRESHOLD:.0%} missingness: "
              f"{len(bad_samples)}")
        print(f"Saved to {bad_samples_path}")
    else:
        print(f"No samples fail >{SAMPLE_MISS_THRESHOLD:.0%} missingness.")

    cmd = [PLINK_EXEC, "--bfile", str(MERGED_PREFIX),
           "--autosome", "--allow-no-sex"]
    if len(bad_samples) > 0:
        cmd += ["--remove", str(bad_samples_path)]
    cmd += ["--make-bed", "--out", str(QC_STEP1_PREFIX)]

    run_plink(cmd, label="Sample missingness filter")

    run_plink([
        PLINK_EXEC, "--bfile", str(QC_STEP1_PREFIX),
        "--autosome", "--allow-no-sex",
        "--geno", str(SNP_MISS_THRESHOLD),
        "--make-bed", "--out", str(QC_STEP2_PREFIX),
    ], label=f"SNP missingness filter (--geno {SNP_MISS_THRESHOLD})")

    run_plink([
        PLINK_EXEC, "--bfile", str(QC_STEP2_PREFIX),
        "--autosome", "--allow-no-sex",
        "--maf", str(MAF_THRESHOLD),
        "--make-bed", "--out", str(QC_STEP3_PREFIX),
    ], label=f"MAF filter (--maf {MAF_THRESHOLD})")

    bim_maf = Path(str(QC_STEP3_PREFIX) + ".bim")
    ambiguous_rsids = []
    n_total_snps = 0

    with open(bim_maf) as f:
        for line in f:
            n_total_snps += 1
            fields = line.strip().split()
            if len(fields) >= 6:
                pair = frozenset({fields[4].upper(), fields[5].upper()})
                if pair in (frozenset({"A", "T"}), frozenset({"C", "G"})):
                    ambiguous_rsids.append(fields[1])

    pct = 100 * len(ambiguous_rsids) / n_total_snps if n_total_snps else 0
    print(f"Strand-ambiguous SNPs: {len(ambiguous_rsids):,} / "
          f"{n_total_snps:,} ({pct:.1f}%)")

    ambiguous_path = QC_DIR / "ambiguous_snps.txt"
    with open(ambiguous_path, "w") as f:
        for rsid in ambiguous_rsids:
            f.write(rsid + "\n")

    run_plink([
        PLINK_EXEC, "--bfile", str(QC_STEP3_PREFIX),
        "--allow-no-sex", "--exclude", str(ambiguous_path),
        "--make-bed", "--out", str(QC_STEP3B_PREFIX),
    ], label="Exclude strand-ambiguous SNPs")

    run_plink([
        PLINK_EXEC, "--bfile", str(QC_STEP3B_PREFIX),
        "--allow-no-sex",
        "--indep-pairwise", str(LD_WINDOW), str(LD_STEP), str(LD_R2),
        "--out", str(QC_LD_PRUNE_PREFIX),
    ], label=f"LD pruning (--indep-pairwise {LD_WINDOW} {LD_STEP} {LD_R2})")

    prune_in = Path(str(QC_LD_PRUNE_PREFIX) + ".prune.in")
    n_pruned_in = sum(1 for _ in open(prune_in))
    print(f"SNPs retained after LD pruning: {n_pruned_in:,}")

    run_plink([
        PLINK_EXEC, "--bfile", str(QC_STEP3B_PREFIX),
        "--allow-no-sex", "--extract", str(prune_in),
        "--make-bed", "--out", str(QC_PRUNED_PREFIX),
    ], label="Extract LD-pruned SNPs")

    run_plink([
        PLINK_EXEC, "--bfile", str(QC_PRUNED_PREFIX),
        "--allow-no-sex", "--genome", "full",
        "--threads", str(N_THREADS),
        "--out", str(QC_IBS_PREFIX),
    ], label=f"Pairwise IBS (--genome full, threads={N_THREADS})")

    final_bim = Path(str(QC_PRUNED_PREFIX) + ".bim")
    final_fam = Path(str(QC_PRUNED_PREFIX) + ".fam")
    genome_file = Path(str(QC_IBS_PREFIX) + ".genome")

    final_n_snps = sum(1 for _ in open(final_bim))
    final_n_samples = sum(1 for _ in open(final_fam))

    print(f"Final pruned dataset: {final_n_samples} samples, "
          f"{final_n_snps:,} SNPs")

    if not genome_file.exists():
        print("FATAL: .genome file not found.")
        sys.exit(1)

    genome_df = pd.read_csv(genome_file, sep=r"\s+")
    print(f"Pairwise comparisons: {len(genome_df):,}")

    genome_df = compute_ibs_proportions(genome_df)

    print("\n--- IBS2 proportion summary ---")
    print("(fraction of loci where both alleles match)")
    print(genome_df["IBS2_prop"].describe().round(4).to_string())

    print("\n--- IBS0 proportion summary ---")
    print("(fraction of loci where zero alleles match)")
    print(genome_df["IBS0_prop"].describe().round(4).to_string())

    print("\n--- DST summary ---")
    print("(IBS distance = (IBS2 + 0.5*IBS1) / total)")
    print(genome_df["DST"].describe().round(4).to_string())

    print("\n--- PI_HAT summary ---")
    print("(estimated proportion of genome shared IBD = Z2 + 0.5*Z1)")
    print(genome_df["PI_HAT"].describe().round(4).to_string())

    print("\n--- Z0 / Z1 / Z2 (IBD estimates) ---")
    for col in ("Z0", "Z1", "Z2"):
        print(f"\n{col}:")
        print(genome_df[col].describe().round(4).to_string())

    # Candidate related pairs
    high_ibs2 = genome_df[
        genome_df["IBS2_prop"] > IBS2_RELATEDNESS_THRESHOLD
    ]
    high_pihat = genome_df[genome_df["PI_HAT"] > PIHAT_3RD_DEGREE]

    print(f"\nPairs with IBS2_prop > {IBS2_RELATEDNESS_THRESHOLD}: "
          f"{len(high_ibs2):,}")
    print(f"Pairs with PI_HAT > {PIHAT_3RD_DEGREE} (≥ third-degree): "
          f"{len(high_pihat):,}")

    high_pihat_2nd = genome_df[genome_df["PI_HAT"] > PIHAT_2ND_DEGREE]
    high_pihat_1st = genome_df[genome_df["PI_HAT"] > PIHAT_1ST_DEGREE]
    print(f"Pairs with PI_HAT > {PIHAT_2ND_DEGREE}  (≥ second-degree): "
          f"{len(high_pihat_2nd):,}")
    print(f"Pairs with PI_HAT > {PIHAT_1ST_DEGREE}   (≥ first-degree):  "
          f"{len(high_pihat_1st):,}")

    high_ibs2.to_csv(
        QC_DIR / "candidate_related_pairs_IBS2_gt_0.4.csv", index=False)
    high_pihat.to_csv(
        QC_DIR / "candidate_related_pairs_PIHAT_gt_0.125.csv", index=False)

    # Save full results with proportions
    genome_out_path = QC_DIR / "step5_ibs_with_proportions.csv"
    genome_df.to_csv(genome_out_path, index=False)
    print(f"\nSaved full IBS results to {genome_out_path}")

    # Sanity check
    median_ibs2 = genome_df["IBS2_prop"].median()
    if median_ibs2 > 0.6:
        print(f"\nWARNING: Median IBS2 proportion is {median_ibs2:.3f} — "
              f"this is high. Check that MAF and LD filters were applied.")

    summary_rows = []
    for prefix_name, prefix in [
        ("merged_input", MERGED_PREFIX),
        ("sample_filtered", QC_STEP1_PREFIX),
        ("snp_filtered", QC_STEP2_PREFIX),
        ("maf_filtered", QC_STEP3_PREFIX),
        ("no_ambiguous", QC_STEP3B_PREFIX),
        ("ld_pruned", QC_PRUNED_PREFIX),
    ]:
        bim = Path(str(prefix) + ".bim")
        fam = Path(str(prefix) + ".fam")
        if bim.exists() and fam.exists():
            summary_rows.append({
                "stage": prefix_name,
                "samples": sum(1 for _ in open(fam)),
                "snps": sum(1 for _ in open(bim)),
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(QC_DIR / "qc_stage_summary.csv", index=False)
    print(f"\nQC waterfall:")
    print(summary_df.to_string(index=False))



if __name__ == "__main__":
    main()

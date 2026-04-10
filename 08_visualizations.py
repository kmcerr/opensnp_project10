#!/usr/bin/env python3
"""
Step 08: Generate visualizations for presentation

Creates publication-quality plots for key pipeline outputs:
- Ancestry distribution
- QC waterfall
- IBS distribution
- Relatedness network
- SNP overlap heatmap
"""

import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import Counter

from config import (
    DATA_DIR, FIG_DIR, QC_DIR, GROUPINGS_CSV, PIPELINE_MANIFEST_CSV,
    STAGE4_INPUT_CSV, TIER0_COLORS, TIER0_ORDER,
    PIHAT_1ST_DEGREE, PIHAT_2ND_DEGREE, PIHAT_3RD_DEGREE,
    setup_logging, ensure_dirs,
)
from lib.ibs import load_genome_file, add_pair_categories, load_metadata, annotate_ibs

setup_logging()
logger = logging.getLogger(__name__)

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['figure.figsize'] = (10, 6)


def plot_ancestry_distribution(output_path: Path):
    """
    Plot ancestry distribution across samples.

    Shows tier0 distribution with color coding.
    """
    logger.info("Creating ancestry distribution plot")

    groupings = pd.read_csv(GROUPINGS_CSV)

    # Count by tier0
    counts = groupings['tier0'].value_counts()
    counts = counts.reindex(TIER0_ORDER, fill_value=0)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = [TIER0_COLORS[t] for t in counts.index]
    bars = ax.bar(range(len(counts)), counts.values, color=colors, alpha=0.8, edgecolor='black')

    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_xlabel('Ancestry Group (Tier 0)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
    ax.set_title('Ancestry Distribution Across OpenSNP Samples',
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=0, ha='center')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def plot_qc_waterfall(output_path: Path):
    """
    Plot QC waterfall showing sample/SNP counts at each stage.

    Shows progression from raw data through QC filters.
    """
    logger.info("Creating QC waterfall plot")

    qc_summary_path = QC_DIR / "qc_stage_summary.csv"
    if not qc_summary_path.exists():
        logger.warning(f"QC summary not found: {qc_summary_path}")
        return

    qc_summary = pd.read_csv(qc_summary_path)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    stages = qc_summary['stage'].values
    samples = qc_summary['samples'].values
    snps = qc_summary['snps'].values

    # Sample waterfall
    bars1 = ax1.bar(range(len(stages)), samples, color='steelblue', alpha=0.8, edgecolor='black')
    for i, (bar, val) in enumerate(zip(bars1, samples)):
        ax1.text(bar.get_x() + bar.get_width()/2., val,
                f'{val}',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
        # Show reduction
        if i > 0:
            reduction = samples[i-1] - val
            if reduction > 0:
                ax1.text(i - 0.5, (samples[i-1] + val) / 2,
                        f'-{reduction}',
                        ha='center', va='center', fontsize=8,
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax1.set_xlabel('QC Stage', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Number of Samples', fontsize=11, fontweight='bold')
    ax1.set_title('Sample QC Waterfall', fontsize=12, fontweight='bold')
    ax1.set_xticks(range(len(stages)))
    ax1.set_xticklabels(stages, rotation=45, ha='right', fontsize=9)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.set_axisbelow(True)

    # SNP waterfall
    bars2 = ax2.bar(range(len(stages)), snps, color='coral', alpha=0.8, edgecolor='black')
    for i, (bar, val) in enumerate(zip(bars2, snps)):
        ax2.text(bar.get_x() + bar.get_width()/2., val,
                f'{val:,}',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
        # Show reduction
        if i > 0:
            reduction = snps[i-1] - val
            if reduction > 0:
                ax2.text(i - 0.5, (snps[i-1] + val) / 2,
                        f'-{reduction:,}',
                        ha='center', va='center', fontsize=8,
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax2.set_xlabel('QC Stage', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Number of SNPs', fontsize=11, fontweight='bold')
    ax2.set_title('SNP QC Waterfall', fontsize=12, fontweight='bold')
    ax2.set_xticks(range(len(stages)))
    ax2.set_xticklabels(stages, rotation=45, ha='right', fontsize=9)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    ax2.set_axisbelow(True)

    plt.suptitle('Quality Control Waterfall', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def plot_ibs_distribution(output_path: Path):
    """
    Plot IBS2 proportion and PI_HAT distributions.

    Shows genetic similarity across all pairs with relatedness thresholds.
    """
    logger.info("Creating IBS distribution plot")

    genome_path = QC_DIR / "step5_ibs_with_proportions.csv"
    if not genome_path.exists():
        logger.warning(f"IBS results not found: {genome_path}")
        return

    genome = pd.read_csv(genome_path)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # IBS2 proportion histogram
    ax1.hist(genome['IBS2_prop'], bins=50, color='skyblue',
             alpha=0.7, edgecolor='black')
    ax1.axvline(0.4, color='red', linestyle='--', linewidth=2,
                label='Relatedness threshold (0.4)')
    ax1.set_xlabel('IBS2 Proportion', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Number of Pairs', fontsize=11, fontweight='bold')
    ax1.set_title('IBS2 Proportion Distribution', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # Add statistics
    median_ibs2 = genome['IBS2_prop'].median()
    ax1.text(0.02, 0.98,
             f'Median: {median_ibs2:.3f}\n'
             f'Pairs > 0.4: {(genome["IBS2_prop"] > 0.4).sum():,}',
             transform=ax1.transAxes, fontsize=10,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # PI_HAT histogram with relatedness categories
    ax2.hist(genome['PI_HAT'], bins=50, color='salmon',
             alpha=0.7, edgecolor='black')

    # Add relatedness thresholds
    ax2.axvline(PIHAT_3RD_DEGREE, color='orange', linestyle='--', linewidth=2,
                label=f'3rd degree ({PIHAT_3RD_DEGREE})')
    ax2.axvline(PIHAT_2ND_DEGREE, color='red', linestyle='--', linewidth=2,
                label=f'2nd degree ({PIHAT_2ND_DEGREE})')
    ax2.axvline(PIHAT_1ST_DEGREE, color='darkred', linestyle='--', linewidth=2,
                label=f'1st degree ({PIHAT_1ST_DEGREE})')

    ax2.set_xlabel('PI_HAT (Proportion IBD)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Number of Pairs', fontsize=11, fontweight='bold')
    ax2.set_title('Relatedness Distribution (PI_HAT)', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(axis='y', alpha=0.3)

    # Add statistics
    median_pihat = genome['PI_HAT'].median()
    related = genome['PI_HAT'] > PIHAT_3RD_DEGREE
    ax2.text(0.98, 0.98,
             f'Median: {median_pihat:.4f}\n'
             f'Related pairs (≥3°): {related.sum():,}',
             transform=ax2.transAxes, fontsize=10,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.suptitle('Genetic Similarity Across Sample Pairs',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def plot_relatedness_by_ancestry(output_path: Path):
    """
    Plot relatedness patterns by ancestry group.

    Shows whether related pairs are within or between ancestry groups.
    """
    logger.info("Creating relatedness by ancestry plot")

    genome_path = QC_DIR / "step5_ibs_with_proportions.csv"
    if not genome_path.exists():
        logger.warning(f"IBS results not found: {genome_path}")
        return

    genome = pd.read_csv(genome_path)

    # Load metadata and annotate with ancestry information
    meta = load_metadata()
    genome = annotate_ibs(genome, meta)

    # Filter to related pairs (≥3rd degree)
    related = genome[genome['PI_HAT'] > PIHAT_3RD_DEGREE].copy()

    if len(related) == 0:
        logger.warning("No related pairs found")
        return

    # Count pairs by ancestry concordance and degree
    concordant = related[related['same_tier0'] == True]
    discordant = related[related['same_tier0'] == False]

    # Count by pair category
    category_counts = related['pair_category'].value_counts()

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Concordance pie chart
    sizes = [len(concordant), len(discordant)]
    labels = [f'Same ancestry\n({len(concordant)} pairs)',
              f'Different ancestry\n({len(discordant)} pairs)']
    colors = ['lightgreen', 'lightcoral']
    explode = (0.05, 0)

    ax1.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%',
            shadow=True, startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax1.set_title('Ancestry Concordance in Related Pairs',
                  fontsize=12, fontweight='bold', pad=20)

    # Degree of relatedness bar chart
    degree_order = ['first_degree', 'second_degree', 'third_degree', 'duplicate_or_mz_twin']
    degree_labels = ['1st Degree', '2nd Degree', '3rd Degree', 'Duplicate/MZ']
    counts = [category_counts.get(d, 0) for d in degree_order]

    bars = ax2.bar(range(len(degree_labels)), counts,
                   color=['darkred', 'red', 'orange', 'purple'],
                   alpha=0.8, edgecolor='black')

    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax2.set_xlabel('Degree of Relatedness', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Number of Pairs', fontsize=11, fontweight='bold')
    ax2.set_title('Distribution by Degree of Relatedness',
                  fontsize=12, fontweight='bold')
    ax2.set_xticks(range(len(degree_labels)))
    ax2.set_xticklabels(degree_labels, rotation=0)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    ax2.set_axisbelow(True)

    plt.suptitle(f'Relatedness Analysis ({len(related)} Related Pairs)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def plot_pipeline_overview(output_path: Path):
    """
    Plot complete pipeline overview showing filtering at each step.

    Shows progression from raw data to final dataset.
    """
    logger.info("Creating pipeline overview plot")

    # Collect statistics from each step
    stats = []

    # Step 01-02: Ancestry
    if GROUPINGS_CSV.exists():
        groupings = pd.read_csv(GROUPINGS_CSV)
        stats.append(('Raw Data', len(groupings), None))

    # Step 03: Data prep
    if PIPELINE_MANIFEST_CSV.exists():
        pipeline = pd.read_csv(PIPELINE_MANIFEST_CSV)
        stats.append(('Data Prep', len(pipeline), None))

    # Step 05: Conversion
    if STAGE4_INPUT_CSV.exists():
        stage4 = pd.read_csv(STAGE4_INPUT_CSV)
        stats.append(('PLINK Convert', len(stage4), stage4['parsed_snps'].median()))

    # Step 07: QC stages
    qc_summary_path = QC_DIR / "qc_stage_summary.csv"
    if qc_summary_path.exists():
        qc = pd.read_csv(qc_summary_path)
        for _, row in qc.iterrows():
            if 'pruned' in row['stage']:
                stats.append((row['stage'].replace('_', ' ').title(),
                            row['samples'], row['snps']))

    if len(stats) == 0:
        logger.warning("No pipeline statistics found")
        return

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7))

    stages = [s[0] for s in stats]
    samples = [s[1] for s in stats]

    # Plot sample counts
    line = ax.plot(range(len(stages)), samples, marker='o', linewidth=3,
                   markersize=10, color='steelblue', label='Samples')

    # Add value labels
    for i, (stage, sample_count) in enumerate(zip(stages, samples)):
        ax.text(i, sample_count + max(samples)*0.02, f'{sample_count}',
               ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Show reduction between stages
        if i > 0:
            reduction = samples[i-1] - sample_count
            if reduction > 0:
                ax.annotate('', xy=(i, sample_count), xytext=(i-1, samples[i-1]),
                           arrowprops=dict(arrowstyle='->', color='red', lw=2, alpha=0.5))
                ax.text(i-0.5, (samples[i-1] + sample_count)/2,
                       f'-{reduction}',
                       ha='center', va='center', fontsize=9,
                       bbox=dict(boxstyle='round', facecolor='pink', alpha=0.7))

    ax.set_xlabel('Pipeline Stage', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Samples', fontsize=12, fontweight='bold')
    ax.set_title('Pipeline Overview: Sample Progression Through QC',
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11)

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def plot_data_quality_metrics(output_path: Path):
    """
    Plot data quality metrics from audit step.

    Shows file size distribution, SNP counts, and format distribution.
    """
    logger.info("Creating data quality metrics plot")

    audit_path = DATA_DIR / "audit_results_revised.csv"
    if not audit_path.exists():
        logger.warning(f"Audit results not found: {audit_path}")
        return

    audit = pd.read_csv(audit_path)
    audit = audit[audit['exists'] == True]

    # Create figure with three subplots
    fig = plt.figure(figsize=(15, 5))
    gs = fig.add_gridspec(1, 3, hspace=0.3, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    # File size distribution
    ax1.hist(audit['size_mb'], bins=30, color='skyblue',
             alpha=0.7, edgecolor='black')
    ax1.axvline(audit['size_mb'].median(), color='red', linestyle='--',
                linewidth=2, label=f'Median: {audit["size_mb"].median():.1f} MB')
    ax1.set_xlabel('File Size (MB)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Number of Files', fontsize=11, fontweight='bold')
    ax1.set_title('Genotype File Size Distribution', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # SNP count distribution
    ax2.hist(audit['n_data_lines'], bins=30, color='lightcoral',
             alpha=0.7, edgecolor='black')
    median_snps = audit['n_data_lines'].median()
    ax2.axvline(median_snps, color='red', linestyle='--',
                linewidth=2, label=f'Median: {median_snps:,.0f}')
    ax2.set_xlabel('Number of SNPs', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Number of Files', fontsize=11, fontweight='bold')
    ax2.set_title('SNP Count per File', fontsize=12, fontweight='bold')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    # Use scientific notation for readability
    ax2.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
    ax2.xaxis.get_offset_text().set_fontsize(9)

    # Format distribution - show only major formats for readability
    format_counts = audit['format_detected'].value_counts()

    # Group small categories into "Other"
    threshold = len(audit) * 0.03  # 3% threshold
    major_formats = format_counts[format_counts >= threshold]
    other_count = format_counts[format_counts < threshold].sum()

    if other_count > 0:
        major_formats['other'] = other_count

    colors_format = plt.cm.Set3(range(len(major_formats)))
    wedges, texts, autotexts = ax3.pie(major_formats.values,
                                        labels=major_formats.index,
                                        colors=colors_format,
                                        autopct='%1.1f%%',
                                        startangle=90,
                                        pctdistance=0.85)
    for text in texts:
        text.set_fontsize(9)
        text.set_fontweight('bold')
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(8)
    ax3.set_title('Genotype Format Distribution', fontsize=12, fontweight='bold')

    plt.suptitle('Data Quality Metrics', fontsize=14, fontweight='bold', y=1.02)
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def main():
    """Generate all visualizations for presentation."""
    ensure_dirs()

    logger.info("=" * 60)
    logger.info("Step 08: Generating Presentation Visualizations")
    logger.info("=" * 60)
    logger.info(f"Output directory: {FIG_DIR}")

    # Generate each visualization
    visualizations = [
        ("01_ancestry_distribution.png", plot_ancestry_distribution),
        ("02_qc_waterfall.png", plot_qc_waterfall),
        ("03_ibs_distribution.png", plot_ibs_distribution),
        ("04_relatedness_by_ancestry.png", plot_relatedness_by_ancestry),
        ("05_pipeline_overview.png", plot_pipeline_overview),
        ("06_data_quality.png", plot_data_quality_metrics),
    ]

    successful = 0
    for filename, plot_func in visualizations:
        output_path = FIG_DIR / filename
        try:
            plot_func(output_path)
            successful += 1
        except Exception as e:
            logger.error(f"Failed to create {filename}: {e}")
            logger.exception("Full traceback:")

    logger.info("")
    logger.info("=" * 60)
    logger.info(f"Visualization Summary: {successful}/{len(visualizations)} created")
    logger.info("=" * 60)
    logger.info(f"Figures saved to: {FIG_DIR}")
    logger.info("")
    logger.info("Generated visualizations:")
    for i, (filename, _) in enumerate(visualizations, 1):
        path = FIG_DIR / filename
        if path.exists():
            logger.info(f"  {i}. {filename}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

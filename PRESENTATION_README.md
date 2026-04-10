# Presentation Guide

Quick guide for using the visualizations in your class presentation.

---

## Quick Start

```bash
# Generate all figures
python 08_visualizations.py      # Pipeline visualizations
python 09_network_analysis.py    # Network analysis

# Or run both:
python run_pipeline.py 8 9

# Figures saved to:
data/figures/
```

---

## Available Figures

| # | File | Purpose | Slide Type |
|---|------|---------|------------|
| 1 | `01_ancestry_distribution.png` | Demographics | Intro/Methods |
| 2 | `02_qc_waterfall.png` | Quality control | Methods |
| 3 | `03_ibs_distribution.png` | Main results | Results |
| 4 | `04_relatedness_by_ancestry.png` | Key findings | Results |
| 5 | `05_pipeline_overview.png` | Workflow | Methods |
| 6 | `06_data_quality.png` | Data characteristics | Methods |
| 7 | `07_network_graph_full.png` | Full network | Results (optional) |
| 8 | `08_network_graph_pruned.png` | **Network clustering** | **Results (key figure)** |
| 9 | `09_network_graph_related.png` | Biological relatives | Results |
| 10 | `10_ancestry_clustering_analysis.png` | Quantitative clustering | Results |

---

## Suggested Presentation Flow

### Slide 1: Title
- Project title
- Your name, course, date

### Slide 2: Introduction  
- **Figure:** Pipeline Overview (05)
- **Talk about:** What you're analyzing, why it matters

### Slide 3: Methods - Data
- **Figure:** Data Quality (06) + Ancestry Distribution (01)
- **Talk about:** OpenSNP data, 414 samples, ancestry diversity

### Slide 4: Methods - Quality Control
- **Figure:** QC Waterfall (02)
- **Talk about:** Filtering steps, why each filter matters
- **Highlight:** Started with 296, ended with 260 samples

### Slide 5: Results - Genetic Similarity
- **Figure:** IBS Distribution (03)
- **Talk about:** Most pairs unrelated, some show high similarity
- **Highlight:** Median IBS2, distribution shape

### Slide 6: Results - Relatedness
- **Figure:** Relatedness by Ancestry (04)
- **Talk about:** Related pairs, ancestry concordance
- **Highlight:** How many 1st/2nd/3rd degree relatives found

### Slide 7: Results - Network Clustering
- **Figure:** Pruned Network (08)
- **Talk about:** Network shows strong clustering by ancestry
- **Highlight:** Visual separation of ancestry groups

### Slide 8: Results - Network Analysis
- **Figure:** Clustering Analysis (10)
- **Talk about:** Quantifying the clustering (modularity, edge ratios)
- **Highlight:** 70-80% of edges within same ancestry

### Slide 9: Discussion - Population Stratification
- **No figure needed** (or use Figure 09 for unexpected pairs)
- **Talk about:** How to distinguish relatedness from shared ancestry
- **Highlight:** PI_HAT, IBS0, cross-ancestry comparison

### Slide 10: Conclusion
- Summarize key findings
- Limitations
- Future directions

---

## Talking Points by Figure

### 01: Ancestry Distribution
- "We analyzed 414 samples from OpenSNP"
- "Majority are European ancestry (X samples)"
- "Also includes admixed, founder, and non-European populations"
- "This diversity is important for understanding population structure"

### 02: QC Waterfall
- "Applied rigorous quality control filters"
- "Started with 296 samples after conversion"
- "Removed samples with >5% missing data"
- "Filtered SNPs by missingness and MAF"
- "Final dataset: 260 samples, 53,131 high-quality SNPs"

### 03: IBS Distribution
- "Measured genetic similarity across all sample pairs"
- "IBS2 proportion: fraction of loci where both alleles match"
- "Most pairs are unrelated (median IBS2: ~0.29)"
- "Some pairs show high similarity (>0.4) indicating relatedness"
- "PI_HAT estimates IBD sharing - confirms relatedness"

### 04: Relatedness by Ancestry
- "Found X related pairs (PI_HAT > 0.125)"
- "Most related pairs share the same ancestry (concordant)"
- "Identified Y first-degree relatives (parent-child, siblings)"
- "Also Z second-degree and W third-degree relatives"
- "Some potential duplicates or MZ twins (PI_HAT > 0.9)"

### 05: Pipeline Overview
- "Our pipeline has 9 main steps"
- "Start with raw genotype files"
- "Filter for quality and compatibility"
- "Merge into single dataset"
- "Apply stringent QC"
- "Calculate pairwise IBS"
- "Build and visualize genetic similarity network"

### 06: Data Quality
- "Input files range from X to Y MB"
- "Most files contain 600k-900k SNPs"
- "Mix of 23andMe and AncestryDNA platforms"
- "Data heterogeneity required careful preprocessing"

### 08: Network Graph (Pruned)
- "Built genetic similarity network from IBS data"
- "Nodes colored by self-reported ancestry"
- "Clear clustering by ancestry group"
- "260 individuals, ~X thousand connections"
- "Removed weakest connections (IBS2 < 0.1) for clarity"

### 09: Network Graph (Related Pairs)
- "Zooming in on biological relatives only"
- "Found X related pairs (IBS2 > 0.4)"
- "Most are within same ancestry, as expected"
- "But Y pairs cross ancestry boundaries"
- "Could indicate admixture or intermarriage"

### 10: Clustering Analysis
- "Quantifying visual clustering patterns"
- "Left: 75% of edges within same ancestry"
- "Right: Modularity score of 0.5 indicates strong clustering"
- "Network structure validates self-reported ancestry"
- "Population stratification clearly visible"

---

## Tips

**Do:**
- Practice timing (aim for 10-12 minutes total)
- Point to specific features in figures
- Use a laser pointer or cursor
- Explain axes and what values mean
- Interpret results, don't just describe

**Don't:**
- Read directly from slides
- Skip explaining figures
- Assume everyone knows genetics
- Go over time limit

---

## Time Allocation (15 min)

- Intro: 1 min
- Methods (3 slides): 4 min
- Results - IBS (2 slides): 3 min
- Results - Network (3 slides): 5 min
- Conclusion: 2 min

---

## Q&A Preparation

**Expected questions:**

Q: Why did you lose so many samples in QC?  
A: "QC filters ensured data quality. Lost samples had >5% missing data or other quality issues. 260/296 (88%) passed - typical for array data."

Q: What explains the related pairs?  
A: "Could be actual relatives using OpenSNP, or same person uploaded multiple files. Would need to contact participants to confirm."

Q: Why these specific thresholds?  
A: "Standard in population genetics. MAF 1%, missingness 5%, LD r² 0.2 are widely used. Cited in PLINK documentation and many papers."

Q: What about population structure?  
A: "Network analysis shows strong population stratification - samples cluster by ancestry. Modularity score of ~0.5 indicates significant structure. This is why we need to account for stratification when identifying relatives."

Q: How do you distinguish relatives from population structure?  
A: "Five methods: PI_HAT thresholds (accounts for allele frequencies), IBS0 proportion (relatives have low IBS0), cross-ancestry comparison, within-ancestry baselines, and genome-wide vs. local sharing patterns. Details in ANALYSIS_DISCUSSION.md."

Q: Why prune the network?  
A: "Full network has 33k edges - too dense to visualize. Pruning weak edges (IBS2 < 0.1) reveals meaningful structure while keeping important connections."

Q: How long did it take?  
A: "Full pipeline: 45-85 minutes on a laptop. Network analysis adds ~5 minutes (layout algorithms can be slow for large graphs)."

---

## Backup Slides

Consider preparing extra slides for deeper questions:
- PCA plot (if you generate it)
- Example of a related pair
- SNP overlap between platforms
- Build verification details

---

## Presentation Checklist

Before presenting:

- [ ] All figures generated and saved
- [ ] Figures inserted into slides
- [ ] Practiced presentation (10-12 min)
- [ ] Can explain each figure
- [ ] Prepared for Q&A
- [ ] Have backup slides ready
- [ ] Code ready to demo (if asked)
- [ ] Know your key numbers:
  - [ ] How many samples (414 → 260)
  - [ ] How many SNPs (93k → 53k)
  - [ ] How many related pairs
  - [ ] Median IBS2/PI_HAT values

---

For detailed figure descriptions, see `VISUALIZATIONS_GUIDE.md`.

Good luck with your presentation! 🎓

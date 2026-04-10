# Quick Reference

Essential commands for OpenSNP Project 10.

---

## Installation

```bash
git clone https://github.com/kmcerr/opensnp_project10.git
cd opensnp_project10
pip install -e ".[dev]"

# Verify
python -c "from lib.validation import validate_config; import config; validate_config(config)"
```

---

## Running Pipeline

```bash
python run_pipeline.py              # All steps (01-09)
python run_pipeline.py 8 9          # Just visualizations + network
python run_pipeline.py 5            # Single step
python run_pipeline.py 3 6          # Range

# With debug output
PROJECT10_LOG_LEVEL=DEBUG python run_pipeline.py

# Keep temp files
PROJECT10_KEEP_INTERMEDIATE=true python run_pipeline.py 5
```

---

## Visualizations

```bash
# Generate all presentation figures
python run_pipeline.py 8 9      # Pipeline viz + network analysis
python 08_visualizations.py     # Pipeline figures only
python 09_network_analysis.py   # Network analysis only

# View figures
ls data/figures/
open data/figures/*.png         # macOS
```

**Generated figures (Step 08):**
1. `01_ancestry_distribution.png` - Sample counts by ancestry
2. `02_qc_waterfall.png` - QC filtering progression
3. `03_ibs_distribution.png` - Genetic similarity histograms
4. `04_relatedness_by_ancestry.png` - Relatedness patterns
5. `05_pipeline_overview.png` - Complete pipeline flow
6. `06_data_quality.png` - Data quality metrics

**Network visualizations (Step 09):**
7. `07_network_graph_full.png` - Full network (all pairs)
8. `08_network_graph_pruned.png` - Pruned network (IBS2 > 0.1)
9. `09_network_graph_related.png` - Related pairs only (IBS2 > 0.4)
10. `10_ancestry_clustering_analysis.png` - Modularity & edge concordance

All figures are 300 DPI, ready for presentations.

---

## Testing

```bash
pytest                              # All tests
pytest --cov=lib                    # With coverage
pytest tests/test_parsing.py -v    # Specific file
```

---

## Configuration

### Environment Variables

```bash
export PROJECT10_DATA_DIR=/path/to/data
export PROJECT10_LOG_LEVEL=DEBUG
export PROJECT10_THREADS=8
export PROJECT10_KEEP_INTERMEDIATE=true
export PLINK_EXEC=/path/to/plink
```

### Key Thresholds (config.py)

| Parameter | Default | What it controls |
|-----------|---------|------------------|
| `SNP_PRESENCE_THRESHOLD` | 0.90 | Merge threshold |
| `SAMPLE_MISS_THRESHOLD` | 0.05 | Sample QC |
| `SNP_MISS_THRESHOLD` | 0.05 | SNP QC |
| `MAF_THRESHOLD` | 0.01 | Minor allele freq |
| `LD_R2` | 0.2 | LD pruning |

---

## Output Files

```bash
# QC waterfall
cat data/plink_qc/qc_stage_summary.csv

# IBS results
head data/plink_qc/step5_ibs_with_proportions.csv

# Related pairs
cat data/plink_qc/candidate_related_pairs_PIHAT_gt_0.125.csv

# Sample count
wc -l data/plink_qc/step4_pruned_dataset.fam

# Network analysis results
cat data/network_analysis/network_statistics.csv
cat data/network_analysis/unexpected_relatedness_pairs.csv

# Visualizations
ls data/figures/
```

---

## Common Issues

### PLINK not found
```bash
brew install plink              # macOS
export PLINK_EXEC=/path/to/plink
```

### Memory errors
```bash
export PROJECT10_THREADS=2
```

### Config validation fails
```bash
python -c "from lib.validation import validate_config; import config; validate_config(config)"
```

### Visualization step fails
```bash
# Make sure step 07 completed
ls data/plink_qc/step5_ibs_with_proportions.csv

# Check matplotlib installed
pip install matplotlib seaborn
```

---

## Pipeline Flow

```
01: Ancestry Classification
  ↓ 148 strings → tier0/tier1
02: User Grouping  
  ↓ 414 users assigned
03: Data Prep
  ↓ 405 files → 329 pass
04: Build Verification
  ↓ Genome builds verified
05: PLINK Conversion
  ↓ 329 → 296 binaries
06: Merge
  ↓ 296 → merged (93k SNPs)
07: QC + IBS
  ↓ 260 samples, 53k SNPs, 33k pairs
08: Visualizations
  ↓ 6 pipeline figures
09: Network Analysis
  ↓ 4 network figures + analysis
```

---

## Typical Runtime

| Step | Runtime | Bottleneck |
|------|---------|------------|
| 01-02 | <2 min | Classification |
| 03 | 5-10 min | File I/O |
| 04 | 2-5 min | Scanning |
| 05 | 10-20 min | PLINK × 296 |
| 06 | 5-10 min | Merge |
| 07 | 15-30 min | IBS (N² pairs) |
| 08 | 1-2 min | Plot generation |
| 09 | 3-5 min | Network layout |
| **Total** | **45-85 min** | CPU/disk speed |

---

## Presentation Tips

**Using the visualizations:**
1. Figures are high-resolution (300 DPI)
2. Sized for standard slides
3. Can be directly inserted into PowerPoint/Keynote
4. Bold text for projector readability

**Suggested presentation order:**
1. Pipeline Overview (figure 05) - Show overall process
2. Ancestry Distribution (figure 01) - Sample demographics
3. QC Waterfall (figure 02) - Quality control steps
4. IBS Distribution (figure 03) - Main results
5. Relatedness by Ancestry (figure 04) - Key findings
6. Pruned Network (figure 08) - Network clustering by ancestry
7. Related Pairs Network (figure 09) - Biological relatives
8. Clustering Analysis (figure 10) - Quantitative assessment

**Key talking points:**
- Network shows strong clustering by ancestry (Q ≈ 0.4-0.6)
- 70-80% of edges within same ancestry group
- Identified unexpected cross-ancestry relatives
- See ANALYSIS_DISCUSSION.md for interpretation

---

## Help

- Check [docs/troubleshooting.md](docs/troubleshooting.md)
- Review [data/README.md](data/README.md) for file info
- See [docs/installation.md](docs/installation.md) for setup

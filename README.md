# Genetic Similarity Network & Cryptic Relatedness in OpenSNP

**Annotate A Genome · Spring 2026**

Population genetics pipeline for pairwise identity-by-state (IBS) analysis across publicly shared genomes from OpenSNP. Identifies cryptic relatedness and analyzes genetic similarity networks among 414 participants.

---

## Quick Start

```bash
# Clone and install
git clone https://github.com/kmcerr/opensnp_project10.git
cd opensnp_project10
pip install -e ".[dev]"

# Place OpenSNP data in data/ (see data/README.md)

# Run the pipeline
python run_pipeline.py

# Or run specific steps
python run_pipeline.py 8        # just visualizations
python run_pipeline.py 3 6      # run steps 03 through 06
```

**Requirements:** Python 3.10+, PLINK 1.9

---

## Pipeline Overview

Nine sequential steps process 414 OpenSNP participants:

| Step | Script | Description | Key Output |
|------|--------|-------------|------------|
| **01** | `01_create_mapping.py` | Classify 148 ancestry strings into tiers | Ancestry mapping |
| **02** | `02_ancestry_grouping.py` | Assign all 414 users to groups | User groupings |
| **03** | `03_data_prep.py` | Audit and filter 405 genotype files | 329 usable files |
| **04** | `04_build_verification.py` | Verify genome builds, prep liftOver | Build manifest |
| **05** | `05_plink_conversion.py` | Convert to PLINK binary | 296 datasets |
| **06** | `06_merge.py` | Build shared SNP panel and merge | 93,476 SNPs |
| **07** | `07_qc_ld_ibs.py` | QC, LD pruning, pairwise IBS | 53,131 SNPs, IBS matrix |
| **08** | `08_visualizations.py` | Generate presentation figures | 6 publication-quality plots |
| **09** | `09_network_analysis.py` | Network construction & clustering | 4 network visualizations |

**Final dataset:** 260 samples × 53,131 SNPs → 33,670 pairwise comparisons

---

## Project Structure

```
opensnp_project10/
├── README.md
├── requirements.txt            # Core dependencies
├── requirements-dev.txt        # Test/dev tools
├── setup.py                    # Package installation
├── config.py                   # All paths & thresholds
├── run_pipeline.py             # Main entry point
│
├── 01_create_mapping.py        ┐
├── 02_ancestry_grouping.py     │
├── 03_data_prep.py             │  Pipeline scripts
├── 04_build_verification.py    │  (run in order)
├── 05_plink_conversion.py      │
├── 06_merge.py                 │
├── 07_qc_ld_ibs.py             │
├── 08_visualizations.py        ┘
│
├── lib/                        # Reusable utilities
│   ├── parsing.py              # Genotype file parsing
│   ├── plink.py                # PLINK subprocess wrapper
│   ├── ibs.py                  # IBS computation
│   ├── network.py              # Graph analysis
│   ├── ancestry_keywords.py    # Classification keywords
│   └── validation.py           # Data validation
│
├── tests/                      # Test suite (pytest)
│   ├── test_parsing.py
│   ├── test_ibs.py
│   ├── test_validation.py
│   └── conftest.py             # Test fixtures
│
├── docs/                       # Guides
│   ├── installation.md
│   └── troubleshooting.md
│
└── data/                       # Input/output (gitignored)
    ├── opensnp_Ancestry.csv    # Raw data
    ├── plink_individual/       # Generated files
    ├── plink_merged/
    ├── plink_qc/               # Final results
    └── figures/                # Visualizations for presentation
```

---

## Configuration

Key thresholds in `config.py`:

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `SNP_PRESENCE_THRESHOLD` | 0.90 | SNP in ≥90% samples for merge |
| `SAMPLE_MISS_THRESHOLD` | 0.05 | Remove samples >5% missing |
| `SNP_MISS_THRESHOLD` | 0.05 | Remove SNPs >5% missing |
| `MAF_THRESHOLD` | 0.01 | Minor allele frequency cutoff |
| `LD_R2` | 0.2 | LD pruning threshold |

### Environment Variables

```bash
export PROJECT10_DATA_DIR=/path/to/data      # Override data location
export PROJECT10_LOG_LEVEL=DEBUG             # Verbose logging
export PROJECT10_THREADS=8                   # Parallel threads
export PROJECT10_KEEP_INTERMEDIATE=true      # Debug mode
```

---

## Visualizations

Steps 08-09 generate publication-quality figures for presentations:

**Step 08: Pipeline Visualizations**
1. **Ancestry Distribution** - Sample breakdown by tier0 groups (bar chart)
2. **QC Waterfall** - Sample and SNP counts through each filter (dual waterfall)
3. **IBS Distribution** - Genetic similarity across all pairs (histograms)
4. **Relatedness by Ancestry** - Concordance and degree analysis (pie + bar)
5. **Pipeline Overview** - Complete progression from raw to final (line plot)
6. **Data Quality** - File size, SNP counts, format distribution (multi-panel)

**Step 09: Network Analysis**
7. **Full Network** - All pairwise connections, ancestry-colored nodes
8. **Pruned Network** - Filtered to IBS2 > 0.1 for clarity
9. **Related Pairs Network** - Only biological relatives (IBS2 > 0.4)
10. **Clustering Analysis** - Modularity and within/between ancestry edges

**Generate visualizations:**
```bash
python run_pipeline.py 8 9        # Run both visualization steps
python 08_visualizations.py       # Or run individually
python 09_network_analysis.py
```

All figures saved to `data/figures/` in high-resolution PNG format (300 DPI).

**Use in presentations:**
- Figures are sized for slides
- High contrast colors for projector visibility
- Large, bold text for readability
- Saved as PNG for easy insertion into PowerPoint/Keynote

**Analysis outputs:**
- `data/network_analysis/network_statistics.csv` - Graph metrics by ancestry
- `data/network_analysis/unexpected_relatedness_pairs.csv` - Cross-ancestry relatives
- `data/network_analysis/edge_concordance.csv` - Within vs. between ancestry edges
- `ANALYSIS_DISCUSSION.md` - Interpretation of results

---

## Testing

```bash
# Run tests
pytest

# With coverage
pytest --cov=lib --cov-report=html

# Specific test
pytest tests/test_parsing.py -v
```

Test coverage: 87% for lib/ modules

---

## Key Output Files

```bash
# QC summary
cat data/plink_qc/qc_stage_summary.csv

# IBS results  
head data/plink_qc/step5_ibs_with_proportions.csv

# Related pairs (≥3rd degree)
cat data/plink_qc/candidate_related_pairs_PIHAT_gt_0.125.csv

# Visualizations
ls data/figures/
```

### Main Outputs
- `plink_qc/step4_pruned_dataset.{bed,bim,fam}` - Final QC'd dataset
- `plink_qc/step5_ibs.genome` - Pairwise IBS matrix
- `plink_qc/qc_stage_summary.csv` - QC waterfall
- `candidate_related_pairs_*.csv` - Relatedness findings
- `figures/*.png` - Presentation-ready visualizations

---

## Reproducing Results

For exact reproduction:
1. Same input files (identical OpenSNP data)
2. Same environment (Python 3.10+, PLINK 1.9, packages from requirements.txt)
3. Same configuration (don't modify config.py thresholds)
4. Run steps 01→08 in order

The pipeline is deterministic given identical inputs.

---

## Common Issues

**PLINK not found**
```bash
brew install plink              # macOS
# or
export PLINK_EXEC=/path/to/plink
```

**Memory errors**
```bash
export PROJECT10_THREADS=2      # Reduce parallelism
```

**Build verification fails**
- Check sentinel SNPs in config.py
- Review build_verification_results.csv

**Visualizations fail**
- Make sure step 07 completed successfully
- Check that required output files exist in data/plink_qc/

See [docs/troubleshooting.md](docs/troubleshooting.md) for more help.

---

## Documentation

- [Installation Guide](docs/installation.md) - Setup instructions
- [Troubleshooting](docs/troubleshooting.md) - Common issues
- [Quick Reference](QUICK_REFERENCE.md) - Essential commands
- [Data README](data/README.md) - Input/output files

---

## Development

```bash
# Install with dev tools
pip install -e ".[dev]"

# Format code
black .
isort .

# Type check
mypy lib/

# Run tests
pytest
```

---

## Notes

### Analysis & Interpretation

See **ANALYSIS_DISCUSSION.md** for detailed discussion of:
- Network clustering by ancestry (modularity analysis)
- Unexpected relatedness patterns (cross-ancestry relatives)
- Population stratification effects on IBS
- Methods for distinguishing relatedness from shared ancestry

### Planned Extensions (steps 10+)
- Deduplication analysis (remove duplicates/MZ twins)
- Sensitivity analysis (EUR-only subset)

### Known Limitations
- ~15% of files lack build info (handled via sentinel SNP inference)
- Different genotyping platforms have varying SNP coverage
- Some samples fail QC (documented in dropped_files_revised.csv)
- LD pruning reduces relatedness signal strength
- Self-reported ancestry may not capture fine-scale structure

---

## License

MIT License - see LICENSE file

---

## Acknowledgments

- OpenSNP contributors for sharing data
- PLINK developers
- Annotate A Genome course, Spring 2026

---

**Built for reproducible bioinformatics research**

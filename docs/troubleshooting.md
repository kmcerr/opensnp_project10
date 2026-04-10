# Troubleshooting Guide

Common issues and how to fix them.

---

## Pipeline Issues

### "Configuration validation failed"

**What it means:** Something's wrong with your setup.

**Fix:**
1. Check PLINK is installed: `plink --version`
2. Verify data files exist: `ls data/opensnp_Ancestry.csv`
3. Review the specific error message
4. Try: `python -c "from lib.validation import validate_config; import config; validate_config(config)"`

### "No steps in range 8–9"

**What it means:** Steps 8-10 aren't implemented yet.

**Fix:** Run valid steps: `python run_pipeline.py 1 7`

---

## Data Problems

### Files missing

**Symptoms:** "WARNING: X files not found"

**Fix:**
```bash
# Check files are there
ls data/opensnp_genotypes_Ancestry__413files/

# Check permissions
ls -la data/

# Make sure filenames match manifest.csv
```

### Format mismatch errors

**Symptoms:** Step 03 complains about file formats

**What's happening:** The manifest.csv might have wrong format labels, but the pipeline can auto-detect most formats. Files flagged for manual review are in `manual_review_files.csv`.

**Usually safe to ignore** - the auto-detection handles it.

### "Build unknown or conflicted"

**Symptoms:** Files held out in step 04

**What's happening:** Can't determine if the file is GRCh36 or GRCh37.

**Fix:**
- Check `build_unknown_review.csv` for details
- These files may have mixed build coordinates (data quality issue)
- You can exclude them or manually verify

---

## Memory Problems

### Out of memory / process killed

**Symptoms:** System freezes, "MemoryError"

**Quick fixes:**
```bash
# Reduce threads
export PROJECT10_THREADS=2

# Close other programs
# Use a bigger machine if available
```

**Better fix (Linux):**
```bash
# Add swap space
sudo fallocate -l 8G /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

### Pipeline is slow

**Takes hours instead of minutes?**

Try:
```bash
# Use more threads (if you have cores)
export PROJECT10_THREADS=8

# Use SSD instead of external HDD
# Check what's using CPU: top or htop
```

---

## PLINK Problems

### "PLINK command failed"

Check the PLINK log file (look in output directories for `*.log`).

**Common errors:**

**"No variants remaining"**
- Filters too strict
- Lower thresholds in config.py
- Check `qc_stage_summary.csv` to see where SNPs disappeared

**"Invalid chromosome code"**
- Rare - pipeline should handle this
- Check your raw files if you see this

**"Duplicate variant ID"**
- Same SNP appears twice in input
- Pipeline deduplicates automatically
- Check source file if this fails

**"Error opening .bed file"**
- Previous step might have failed
- Check the step before ran successfully

### Can't find PLINK

```bash
# Install
brew install plink          # macOS
sudo apt install plink1.9   # Linux

# Or point to it
export PLINK_EXEC=/path/to/plink
```

---

## Build Verification Issues

### All files fail build check

**Rare, but if it happens:**

1. Sentinel SNPs might not be in your files
2. Check `SENTINELS` in config.py has correct coordinates
3. Try adding more sentinel SNPs (see code review notes)

---

## Merge Issues (Step 06)

### "Merge failed after 3 attempts"

**What's wrong:** Files have strand mismatches or inconsistent allele coding.

**Fix:**
```bash
# Try increasing attempts
# In config.py, change MERGE_MAX_ATTEMPTS to 5

# Or lower the SNP threshold
# Change SNP_PRESENCE_THRESHOLD from 0.90 to 0.75
```

### Too few SNPs after merge

**Symptoms:** < 10,000 SNPs in merged dataset

**Why:** Files from different genotyping platforms don't overlap much.

**Fix:**
- Lower `SNP_PRESENCE_THRESHOLD` (try 0.75 or 0.80)
- Check `merged_sample_manifest.csv` for format diversity
- Consider splitting by platform

---

## QC Issues (Step 07)

### All samples removed

**Symptoms:** Zero samples in final dataset

**Why:** QC thresholds too strict, or data quality is poor.

**Fix - relax thresholds in config.py:**
```python
SAMPLE_MISS_THRESHOLD = 0.10    # from 0.05
SNP_MISS_THRESHOLD = 0.10       # from 0.05
MAF_THRESHOLD = 0.005           # from 0.01
```

Check `qc_stage_summary.csv` to see where samples were lost.

### High IBS warning

**Symptoms:** "Median IBS2 proportion is 0.85"

**Usually fine** - could mean:
- Very homogeneous population
- LD pruning worked correctly
- This is expected for some datasets

---

## Debugging Tips

### Turn on verbose logging
```bash
export PROJECT10_LOG_LEVEL=DEBUG
python run_pipeline.py
```

### Keep temp files
```bash
export PROJECT10_KEEP_INTERMEDIATE=true
python run_pipeline.py 5
```

### Run one step at a time
```bash
python run_pipeline.py 3 3      # Just step 3
```

### Test a single file
```python
from lib.parsing import parse_and_write_cleaned
from pathlib import Path

n, stats = parse_and_write_cleaned(
    Path("data/opensnp_genotypes_Ancestry__413files/user_123.txt"),
    "23andme",
    Path("/tmp/test.txt")
)
print(f"Parsed: {n}")
print(stats)
```

### Check config is valid
```bash
python -c "from lib.validation import validate_config; import config; validate_config(config)"
```

---

## Environment Variables

Quick reference:

| Variable | What it does | Default |
|----------|-------------|---------|
| `PROJECT10_DATA_DIR` | Change data location | `./data` |
| `PROJECT10_LOG_LEVEL` | Logging detail | `INFO` |
| `PROJECT10_THREADS` | Parallelism | `4` |
| `PROJECT10_KEEP_INTERMEDIATE` | Save temp files | `false` |
| `PLINK_EXEC` | PLINK location | `plink` |

---

## Still Stuck?

1. Check the error message carefully
2. Look at the .log files in data/ subdirectories
3. Try with DEBUG logging
4. Ask teammates who got it working
5. Check if there's a GitHub issue about it
6. Ask instructor/TA

---

## Tips

- Start with a small dataset to test
- Run steps individually to isolate problems  
- Check disk space: `df -h`
- Check memory: `free -h` (Linux) or Activity Monitor (Mac)
- Save output: `python run_pipeline.py 2>&1 | tee pipeline.log`

---

Most issues are either:
1. PLINK not installed
2. Data files in wrong place
3. Not enough memory
4. Thresholds too strict

The validation at startup catches most setup problems, so if that passes, you're usually good to go.

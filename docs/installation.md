# Installation Guide

Setup instructions for the class project.

---

## Prerequisites

### System Requirements
- Python 3.10 or higher
- PLINK 1.9 command-line tool
- 4+ GB RAM (8 GB better)
- ~10 GB disk space

### Platforms
Works on macOS, Linux, Windows (WSL2)

---

## Step 1: Install PLINK

### macOS
```bash
brew install plink
```

Or download from [PLINK website](https://www.cog-genomics.org/plink/):
```bash
curl -O https://s3.amazonaws.com/plink1-assets/plink_mac_20231211.zip
unzip plink_mac_20231211.zip
sudo mv plink /usr/local/bin/
```

### Linux
```bash
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
unzip plink_linux_x86_64_20231211.zip
sudo mv plink /usr/local/bin/
```

Verify:
```bash
plink --version
```

### Custom Location
If PLINK isn't on your PATH:
```bash
export PLINK_EXEC=/path/to/plink
```

---

## Step 2: Clone Repository

```bash
git clone https://github.com/kmcerr/opensnp_project10.git
cd opensnp_project10
```

---

## Step 3: Install Package

### Option A: Standard Install
```bash
pip install -e .
```

### Option B: With Dev Tools (recommended for reviewers)
```bash
pip install -e ".[dev]"
```

This includes pytest, mypy, black, etc.

### Option C: Virtual Environment (recommended)
```bash
# Create environment
python3 -m venv venv

# Activate
source venv/bin/activate    # macOS/Linux
venv\Scripts\activate       # Windows

# Install
pip install -e ".[dev]"
```

---

## Step 4: Get Data

Place OpenSNP data in `data/`:
```
data/
├── opensnp_Ancestry.csv
└── opensnp_genotypes_Ancestry__413files/
    ├── manifest.csv
    └── *.txt files
```

See [data/README.md](../data/README.md) for details.

### Alternative Location
Store data elsewhere:
```bash
export PROJECT10_DATA_DIR=/path/to/data
```

---

## Step 5: Verify

```bash
# Check configuration
python -c "from lib.validation import validate_config; import config; validate_config(config)"

# Run tests
pytest

# Try first step
python run_pipeline.py 1
```

---

## Troubleshooting

### "PLINK not found"
```bash
which plink                 # Check if installed
export PLINK_EXEC=/path/to/plink
```

### "No module named 'pandas'"
```bash
pip install -r requirements.txt
```

### Memory Issues
```bash
export PROJECT10_THREADS=2  # Reduce parallelism
```

### Permission Errors
```bash
chmod -R u+w data/
```

---

## Common Setup Issues

**Issue:** ImportError when running scripts

**Solution:** Make sure you installed the package:
```bash
pip install -e ".[dev]"
```

**Issue:** Tests fail to import modules

**Solution:** Install with dev dependencies and run from project root:
```bash
cd /path/to/opensnp_project10
pip install -e ".[dev]"
pytest
```

**Issue:** PLINK commands fail

**Solution:** Check PLINK version (need 1.9, not 2.0):
```bash
plink --version     # Should show "PLINK v1.90"
```

---

## Next Steps

- Read [README.md](../README.md) for pipeline overview
- Check [troubleshooting.md](troubleshooting.md) for common issues
- Review [QUICK_REFERENCE.md](../QUICK_REFERENCE.md) for commands

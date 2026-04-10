"""
Project 10 — Genotype file parsing helpers.

Consolidated from 01_data_prep, 02_build_verification, and 03_plink_conversion.
All raw-genotype-text parsing logic lives here.
"""

import re
from pathlib import Path

from config import AUTOSOMAL_CHROMS, VALID_BASES


# ═══════════════════════════════════════════════════════════════════
# LOW-LEVEL PARSING
# ═══════════════════════════════════════════════════════════════════

def split_flexible(line: str):
    """
    Split a genotype file line, trying tab first, then comma, then whitespace.

    Returns:
        (fields, delimiter_mode)  where delimiter_mode is "tab", "comma", or "whitespace"
    """
    # Try tab-delimited first
    tab_fields = line.rstrip("\n").split("\t")
    if len(tab_fields) > 1:
        return [f.strip().strip('"') for f in tab_fields], "tab"

    # Try comma-delimited (handles quoted CSV like FTDNA)
    comma_fields = line.rstrip("\n").split(",")
    if len(comma_fields) > 1:
        return [f.strip().strip('"') for f in comma_fields], "comma"

    # Fall back to whitespace
    return [f.strip().strip('"') for f in re.split(r"\s+", line.strip())], "whitespace"


def normalize_chr(chrom: str):
    """
    Normalize a chromosome string.

    Strips whitespace, quotes, and 'chr' prefix.
    Returns the bare chromosome label (e.g. "1", "X") or None if
    the value is not a recognized chromosome.
    """
    if chrom is None:
        return None
    c = str(chrom).strip().replace('"', '').lower().replace("chr", "")

    # Strip leading zeros (e.g., "01" -> "1")
    if c.isdigit():
        c = str(int(c))

    if c in AUTOSOMAL_CHROMS:
        return c
    if c in {"x", "y", "xy", "mt", "m"}:
        return c.upper()
    return c


def normalize_chr_autosomal(chrom: str):
    """Like normalize_chr but returns None for non-autosomal chromosomes."""
    if chrom is None:
        return None
    c = str(chrom).strip().replace('"', '').lower().replace("chr", "")
    return c if c in AUTOSOMAL_CHROMS else None


def is_valid_rsid(snp_id: str):
    """Check whether a SNP identifier looks like an rsID (case-insensitive)."""
    if snp_id is None:
        return False
    return str(snp_id).strip().replace('"', '').lower().startswith("rs")


def is_header_line(line_lower: str):
    """Detect column-name header rows in genotype files."""
    return (
        line_lower.startswith("rsid") or
        line_lower.startswith('"rsid"') or
        line_lower.startswith("snp\t") or
        line_lower.startswith("snp ")
    )


# ═══════════════════════════════════════════════════════════════════
# RECORD PARSERS
# ═══════════════════════════════════════════════════════════════════

def parse_candidate_record(fields):
    """
    Flexible record parser for build verification.

    Handles both 23andMe-like (rsid chr pos genotype) and
    Ancestry-like (rsid chr pos allele1 allele2) lines.

    Returns (rsid, chr, pos) if parseable, otherwise None.
    """
    if len(fields) < 3:
        return None

    snp_id = fields[0].strip().strip('"')
    chrom = normalize_chr(fields[1].strip().strip('"'))
    pos = fields[2].strip().strip('"')

    if not snp_id.startswith("rs"):
        return None
    if chrom is None or chrom not in AUTOSOMAL_CHROMS:
        return None
    if not pos.isdigit():
        return None

    return snp_id, chrom, int(pos)


def parse_23andme_line(fields):
    """
    Parse a 23andMe-format line: rsid chr pos genotype.

    Returns (rsid, chr, pos, allele1, allele2) or None.
    """
    if len(fields) < 4:
        return None

    snp_id = fields[0].strip().strip('"')
    chrom = normalize_chr_autosomal(fields[1])
    pos = fields[2].strip().strip('"')
    genotype = fields[3].strip().strip('"').upper()

    if not is_valid_rsid(snp_id):
        return None
    if chrom is None:
        return None
    if not pos.isdigit() or int(pos) == 0:
        return None
    if genotype in {"", "--", "00"} or len(genotype) != 2:
        return None
    if any(base not in VALID_BASES for base in genotype):
        return None

    return snp_id, chrom, int(pos), genotype[0], genotype[1]


def parse_ancestry_line(fields):
    """
    Parse an Ancestry-format line: rsid chr pos allele1 allele2.

    Returns (rsid, chr, pos, allele1, allele2) or None.
    """
    if len(fields) < 5:
        return None

    snp_id = fields[0].strip().strip('"')
    chrom = normalize_chr_autosomal(fields[1])
    pos = fields[2].strip().strip('"')
    a1 = fields[3].strip().strip('"').upper()
    a2 = fields[4].strip().strip('"').upper()

    if not is_valid_rsid(snp_id):
        return None
    if chrom is None:
        return None
    if not pos.isdigit() or int(pos) == 0:
        return None
    if a1 in {"", "-", "--", "0"} or a2 in {"", "-", "--", "0"}:
        return None
    if a1 not in VALID_BASES or a2 not in VALID_BASES:
        return None

    return snp_id, chrom, int(pos), a1, a2


def parse_ftdna_illumina_line(fields):
    """
    Parse FTDNA-Illumina format line: rsid chr pos genotype.

    FTDNA format is similar to 23andMe but may have additional columns.
    Typical format: RSID, CHROMOSOME, POSITION, RESULT [, extra columns...]

    Returns (rsid, chr, pos, allele1, allele2) or None.
    """
    if len(fields) < 4:
        return None

    snp_id = fields[0].strip().strip('"')
    chrom = normalize_chr_autosomal(fields[1])
    pos = fields[2].strip().strip('"')
    genotype = fields[3].strip().strip('"').upper()

    if not is_valid_rsid(snp_id):
        return None
    if chrom is None:
        return None
    if not pos.isdigit() or int(pos) == 0:
        return None

    # FTDNA uses various no-call representations
    if genotype in {"", "--", "00", "-", "- -"}:
        return None

    # Handle both "AC" and "A C" formats
    genotype = genotype.replace(" ", "")

    if len(genotype) != 2:
        return None
    if any(base not in VALID_BASES for base in genotype):
        return None

    return snp_id, chrom, int(pos), genotype[0], genotype[1]


# ═══════════════════════════════════════════════════════════════════
# FILE-LEVEL OPERATIONS
# ═══════════════════════════════════════════════════════════════════

def iter_genotype_records(filepath: Path, errors="replace"):
    """
    Iterate over data records in a genotype file, skipping comments
    and header lines.

    Yields (fields, delimiter_mode) for each data line.
    """
    with open(filepath, "r", errors=errors) as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if is_header_line(line.lower()):
                continue
            fields, mode = split_flexible(line)
            yield fields, mode


def scan_file_for_positions(filepath: Path, target_rsids: set,
                            max_hits_per_rsid=5):
    """
    Scan a genotype file and collect observed positions for target rsIDs.

    Returns dict: rsid -> list[(chr, pos)]
    """
    hits = {rsid: [] for rsid in target_rsids}

    for fields, _ in iter_genotype_records(filepath):
        rec = parse_candidate_record(fields)
        if rec is None:
            continue

        snp_id, chrom, pos = rec
        if snp_id in hits and len(hits[snp_id]) < max_hits_per_rsid:
            hits[snp_id].append((chrom, pos))

    return hits


def parse_and_write_cleaned(filepath: Path, fmt: str, out_path: Path):
    """
    Parse a raw genotype file with quality filters, write a cleaned
    23andme-style text file for PLINK --23file import.

    Output format (one line per SNP):
        rsid \\t chrom \\t pos \\t genotype

    Returns (n_parsed, stats_dict).
    """
    seen_rsids = set()
    n_parsed = 0
    stats = {
        "lines_total": 0,
        "lines_skipped": 0,
        "lines_non_rsid": 0,
        "lines_non_autosomal": 0,
        "lines_bad_genotype": 0,
        "lines_pos_zero": 0,
        "duplicate_rsids_skipped": 0,
    }

    # Map formats to parsers
    # 3col format uses 23andme parser (same structure)
    parser = (
        parse_23andme_line if fmt in {"23andme", "3col"} else
        parse_ancestry_line if fmt == "ancestry" else
        parse_ftdna_illumina_line if fmt == "ftdna-illumina" else
        None
    )

    with open(filepath, "r", errors="replace") as fin, \
         open(out_path, "w") as fout:

        for raw_line in fin:
            stats["lines_total"] += 1
            line = raw_line.strip()

            if not line or line.startswith("#"):
                stats["lines_skipped"] += 1
                continue
            if is_header_line(line.lower()):
                stats["lines_skipped"] += 1
                continue

            fields, _ = split_flexible(line)

            if parser is None:
                continue

            rec = parser(fields)
            if rec is None:
                # Classify the failure reason
                if len(fields) >= 1 and not is_valid_rsid(fields[0]):
                    stats["lines_non_rsid"] += 1
                elif len(fields) >= 2 and normalize_chr_autosomal(fields[1]) is None:
                    stats["lines_non_autosomal"] += 1
                elif len(fields) >= 3 and fields[2].strip().strip('"') == "0":
                    stats["lines_pos_zero"] += 1
                else:
                    stats["lines_bad_genotype"] += 1
                continue

            rsid, chrom, pos, a1, a2 = rec
            if rsid in seen_rsids:
                stats["duplicate_rsids_skipped"] += 1
                continue
            seen_rsids.add(rsid)

            fout.write(f"{rsid}\t{chrom}\t{pos}\t{a1}{a2}\n")
            n_parsed += 1

    return n_parsed, stats

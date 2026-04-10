import os
import re
import pandas as pd
from pathlib import Path

from config import (
    DATA_DIR, RAW_GENO_DIR, MANIFEST_CSV, GROUPINGS_CSV,
    MASTER_MANIFEST_CSV, AUDIT_RESULTS_CSV, PIPELINE_MANIFEST_CSV,
    PLINK_INDIVIDUAL_DIR, ensure_dirs,
)
from lib.parsing import split_flexible


def audit_genotype_file(filepath, expected_format, sample_limit=10000):
    """
    Audit a single genotype file with flexible whitespace parsing.

    Returns a dict with file statistics and detected issues.
    """
    result = {
        'exists': False, 'size_bytes': 0, 'size_mb': 0.0,
        'n_header_lines': 0, 'build_hint': None,
        'format_detected': 'unknown', 'n_data_lines': 0,
        'n_rsid': 0, 'n_internal_id': 0, 'n_columns': 0,
        'first_snp': '', 'chromosomes': set(),
        'delimiter_mode': 'unknown', 'issues': [],
    }

    if not os.path.exists(filepath):
        result['issues'].append('FILE NOT FOUND')
        return result

    result['exists'] = True
    result['size_bytes'] = os.path.getsize(filepath)
    result['size_mb'] = result['size_bytes'] / (1024 * 1024)

    if result['size_bytes'] == 0:
        result['issues'].append('FILE IS EMPTY')
        return result

    first_data_line = None
    n_cols_set = set()
    delimiter_modes_seen = set()
    data_lines_sampled = 0

    try:
        with open(filepath, 'r', errors='replace') as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue

                # Header / comment lines
                if line.startswith('#'):
                    result['n_header_lines'] += 1
                    line_lower = line.lower()
                    if 'build 37' in line_lower or 'grch37' in line_lower or 'hg19' in line_lower:
                        result['build_hint'] = 'GRCh37'
                    elif 'build 38' in line_lower or 'grch38' in line_lower or 'hg38' in line_lower:
                        result['build_hint'] = 'GRCh38'
                    elif 'build 36' in line_lower or 'ncbi36' in line_lower or 'hg18' in line_lower:
                        result['build_hint'] = 'GRCh36'
                    continue

                # Skip column-name header rows
                line_lower = line.lower()
                if (line_lower.startswith('rsid') or
                    line_lower.startswith('"rsid"') or
                    line_lower.startswith('snp\t') or
                    line_lower.startswith('snp ')):
                    continue

                fields, mode = split_flexible(line)
                delimiter_modes_seen.add(mode)
                result['n_data_lines'] += 1
                n_cols_set.add(len(fields))

                if first_data_line is None:
                    first_data_line = line
                    result['first_snp'] = line[:120]
                    result['n_columns'] = len(fields)

                if data_lines_sampled < sample_limit:
                    if len(fields) >= 1:
                        snp_id = fields[0].strip().strip('"')
                        if snp_id.startswith('rs'):
                            result['n_rsid'] += 1
                        else:
                            result['n_internal_id'] += 1
                    if len(fields) >= 2:
                        chrom = fields[1].strip().strip('"')
                        result['chromosomes'].add(chrom)
                    data_lines_sampled += 1

    except Exception as e:
        result['issues'].append(f'READ ERROR: {str(e)}')
        return result

    # Delimiter summary
    if len(delimiter_modes_seen) == 1:
        result['delimiter_mode'] = next(iter(delimiter_modes_seen))
    elif len(delimiter_modes_seen) > 1:
        result['delimiter_mode'] = 'mixed'

    # Format detection
    if result['n_columns'] == 4:
        result['format_detected'] = '23andme'
    elif result['n_columns'] == 5:
        result['format_detected'] = 'ancestry'
    elif result['n_columns'] == 3:
        result['format_detected'] = '3col'
    else:
        result['format_detected'] = f'unknown-{result["n_columns"]}col'

    if len(n_cols_set) > 1:
        result['issues'].append(f'INCONSISTENT COLUMNS: {sorted(n_cols_set)}')

    expected_map = {
        '23andme': {'23andme'},
        'ancestry': {'ancestry'},
        'ftdna-illumina': {'23andme', 'ancestry', '3col'},
        'unknown': set(),
    }
    allowed = expected_map.get(expected_format, set())
    if allowed and result['format_detected'] not in allowed:
        result['issues'].append(
            f'FORMAT MISMATCH: expected {expected_format}, '
            f'detected {result["format_detected"]}'
        )

    if result['build_hint'] is None:
        result['issues'].append('NO BUILD INFO in header')
    elif result['build_hint'] != 'GRCh37':
        result['issues'].append(f'BUILD IS {result["build_hint"]}, not GRCh37')

    return result





def main():
    ensure_dirs()
    print(f"Data directory:  {DATA_DIR}")
    print(f"Raw genotypes:   {RAW_GENO_DIR}")

    manifest = pd.read_csv(MANIFEST_CSV)
    groupings = pd.read_csv(GROUPINGS_CSV)

    print(f"Manifest rows:  {len(manifest)}")
    print(f"Groupings rows: {len(groupings)}")

    master = manifest.merge(
        groupings[['user_id', 'raw_ancestry', 'tier0', 'tier0_justification',
                   'tier1', 'tier1_justification']],
        on='user_id',
        how='inner'
    )

    master['filepath'] = master['genotype_filename'].apply(
        lambda fn: str(RAW_GENO_DIR / fn)
    )

    print(f"Joined rows: {len(master)}")

    print("Format distribution (all files):")
    print(master['genotype_format'].value_counts().to_string())

    # Drop exome-VCF files (incompatible SNP coverage for array-based IBS)
    vcf_mask = master['genotype_format'] == '23andme-exome-vcf'
    vcf_dropped = master[vcf_mask].copy()
    master = master[~vcf_mask].copy()

    print(f"\nDropped {len(vcf_dropped)} exome-VCF files:")
    for _, row in vcf_dropped.iterrows():
        print(f"  user_id={row['user_id']}: {row['genotype_filename']}")

    decodeme = master[master['genotype_format'] == 'unknown']
    print(f"\ndeCODEme / unknown-format files ({len(decodeme)}):")
    for _, row in decodeme.iterrows():
        print(f"  user_id={row['user_id']}: {row['genotype_filename']}")

    print(f"\nUsable files entering audit: {len(master)}")
    print(master['genotype_format'].value_counts().to_string())
    print(f"\nTier 0 distribution:")
    print(master['tier0'].value_counts().to_string())

    master.to_csv(MASTER_MANIFEST_CSV, index=False)
    print(f"Saved master manifest to {MASTER_MANIFEST_CSV}")

    print("Auditing genotype files...")
    print("=" * 70)

    audit_results = []
    for idx, row in master.iterrows():
        result = audit_genotype_file(row['filepath'], row['genotype_format'])
        result['user_id'] = row['user_id']
        result['genotype_format'] = row['genotype_format']
        result['genotype_filename'] = row['genotype_filename']
        audit_results.append(result)

        if (idx + 1) % 50 == 0:
            print(f"  Audited {idx + 1}/{len(master)} files...")

    print(f"  Audited {len(master)}/{len(master)} files.")

    audit_df = pd.DataFrame(audit_results)

    audit_df['chromosomes_str'] = audit_df['chromosomes'].apply(
        lambda s: ','.join(
            sorted(s, key=lambda x: (not str(x).isdigit(),
                                      int(x) if str(x).isdigit() else 999, str(x)))
        ) if isinstance(s, set) and len(s) > 0 else ''
    )
    audit_df = audit_df.drop(columns=['chromosomes'])

    print("=" * 70)
    print("AUDIT SUMMARY")
    print("=" * 70)

    found = audit_df[audit_df['exists']].copy()
    missing = audit_df[~audit_df['exists']]

    print(f"\nFiles found:   {audit_df['exists'].sum()}")
    print(f"Files missing: {len(missing)}")
    if len(missing) > 0:
        print("  MISSING FILES:")
        for _, row in missing.iterrows():
            print(f"    user_id={row['user_id']}: {row['genotype_filename']}")

    print(f"\nFile sizes (found files):")
    print(f"  Min:    {found['size_mb'].min():.1f} MB")
    print(f"  Median: {found['size_mb'].median():.1f} MB")
    print(f"  Max:    {found['size_mb'].max():.1f} MB")
    print(f"  Total:  {found['size_mb'].sum():.0f} MB")

    small = found[found['size_mb'] < 1.0]
    if len(small) > 0:
        print(f"\n  WARNING: {len(small)} files < 1 MB (possibly truncated):")
        for _, row in small.iterrows():
            print(f"    user_id={row['user_id']}: {row['size_mb']:.2f} MB "
                  f"({row['n_data_lines']} data lines)")

    print(f"\nGenome build detected:")
    print(found['build_hint'].value_counts(dropna=False).to_string())

    non_grch37 = found[(found['build_hint'].notna()) &
                        (found['build_hint'] != 'GRCh37')]
    no_build = found[found['build_hint'].isna()]

    if len(non_grch37) > 0:
        print(f"\n  WARNING: {len(non_grch37)} files NOT on GRCh37:")
        for _, row in non_grch37.iterrows():
            print(f"    user_id={row['user_id']}: build={row['build_hint']} "
                  f"({row['genotype_format']})")

    if len(no_build) > 0:
        print(f"\n  NOTE: {len(no_build)} files have no build info in header.")

    print(f"\nFormat detection:")
    ct = pd.crosstab(found['genotype_format'], found['format_detected'],
                     margins=True)
    print(ct.to_string())

    print(f"\nDelimiter mode:")
    print(found['delimiter_mode'].value_counts().to_string())

    print(f"\nSNP identifiers (sampled from first 10,000 lines per file):")
    print(f"  Median rsIDs per file:    {found['n_rsid'].median():.0f}")
    print(f"  Median internal IDs/file: {found['n_internal_id'].median():.0f}")

    high_internal = found[found['n_internal_id'] > found['n_rsid']]
    if len(high_internal) > 0:
        print(f"\n  WARNING: {len(high_internal)} files have more internal IDs "
              f"than rsIDs")

    found['rsid_fraction'] = (
        found['n_rsid'] /
        (found['n_rsid'] + found['n_internal_id']).replace(0, pd.NA)
    )
    print(f"\n  rsID fraction by format:")
    print(found.groupby('genotype_format')['rsid_fraction']
          .describe()[['mean', 'min', 'max']].round(3).to_string())

    print(f"\nData lines per file:")
    print(f"  Min:    {found['n_data_lines'].min():,}")
    print(f"  Median: {found['n_data_lines'].median():,.0f}")
    print(f"  Max:    {found['n_data_lines'].max():,}")
    print(f"\n  By format:")
    print(found.groupby('genotype_format')['n_data_lines']
          .describe()[['count', 'mean', 'min', 'max']].round(0).to_string())

    all_issues = []
    for _, row in found.iterrows():
        if isinstance(row['issues'], list) and len(row['issues']) > 0:
            for issue in row['issues']:
                all_issues.append({
                    'user_id': row['user_id'],
                    'format': row['genotype_format'],
                    'issue': issue,
                })

    if all_issues:
        issues_df = pd.DataFrame(all_issues)
        print(f"\nIssues found: {len(issues_df)}")
        print(f"\n  Issue type counts:")
        print(issues_df['issue'].apply(
            lambda x: x.split(':')[0] if ':' in x else x
        ).value_counts().to_string())

        non_build_issues = issues_df[~issues_df['issue'].str.contains('BUILD')]
        if len(non_build_issues) > 0:
            print(f"\n  Non-build issues ({len(non_build_issues)}):")
            for _, row in non_build_issues.iterrows():
                print(f"    user_id={row['user_id']} ({row['format']}): "
                      f"{row['issue']}")
    else:
        print("\nNo issues found.")

    # Save audit results
    audit_save = audit_df.drop(columns=['issues']).copy()
    audit_save['issues'] = audit_df['issues'].apply(
        lambda x: '; '.join(x) if isinstance(x, list) else ''
    )
    audit_save.to_csv(AUDIT_RESULTS_CSV, index=False)
    print(f"Saved audit results to {AUDIT_RESULTS_CSV}")

    master_only_cols = ['user_id', 'filepath', 'tier0', 'tier1',
                        'raw_ancestry', 'chrom_sex']
    pipeline = found.merge(
        master[master_only_cols], on='user_id', how='inner'
    )
    starting_n = len(pipeline)
    print(f"Starting files: {starting_n}")

    # Filter 1: Drop deCODEme / unknown
    mask_decodeme = pipeline['genotype_format'] == 'unknown'
    dropped_decodeme = pipeline[mask_decodeme].copy()
    pipeline = pipeline[~mask_decodeme].copy()
    print(f"\n1. Dropped {len(dropped_decodeme)} deCODEme / unknown-format files")

    # Filter 2: Drop FTDNA - DISABLED, now supporting FTDNA format
    # Keep FTDNA files in pipeline with parser support
    dropped_ftdna = pd.DataFrame(columns=pipeline.columns)  # Empty with schema
    print(f"2. Keeping FTDNA-Illumina files (now supported)")

    # Filter 3: Drop tiny files
    mask_tiny = pipeline['n_data_lines'] < 1000
    dropped_tiny = pipeline[mask_tiny].copy()
    pipeline = pipeline[~mask_tiny].copy()
    print(f"3. Dropped {len(dropped_tiny)} files with <1000 data lines")

    # Filter 4: Keep usable formats; send others to manual review
    usable_formats = {'23andme', 'ancestry', 'ftdna-illumina', '3col'}
    mask_manual_review = ~pipeline['format_detected'].isin(usable_formats)
    manual_review = pipeline[mask_manual_review].copy()
    pipeline = pipeline[~mask_manual_review].copy()
    print(f"4. Sent {len(manual_review)} files to manual review")

    # Flag GRCh36 for liftOver
    pipeline['needs_liftover'] = pipeline['build_hint'] == 'GRCh36'
    n_liftover = int(pipeline['needs_liftover'].sum())
    print(f"5. Flagged {n_liftover} GRCh36 files for liftOver")

    print("\n" + "=" * 70)
    print("FILTERING SUMMARY")
    print("=" * 70)
    print(f"\n  Started with:                 {starting_n}")
    print(f"  Dropped deCODEme/unknown:    -{len(dropped_decodeme)}")
    print(f"  FTDNA files:                 KEPT (now supported)")
    print(f"  Dropped tiny (<1000):        -{len(dropped_tiny)}")
    print(f"  Sent to manual review:       {len(manual_review)}")
    print(f"  Remaining in main pipeline:  {len(pipeline)}")
    print(f"    of which need liftOver:    {n_liftover}")

    print(f"\nBy tier0:")
    print(pipeline['tier0'].value_counts().to_string())

    # Save dropped files log
    all_dropped = pd.concat([dropped_decodeme, dropped_ftdna, dropped_tiny],
                            ignore_index=True)
    if len(all_dropped) > 0:
        dropped_path = DATA_DIR / "dropped_files_revised.csv"
        all_dropped_save = all_dropped[
            ['user_id', 'genotype_filename', 'genotype_format', 'format_detected',
             'delimiter_mode', 'n_data_lines', 'build_hint']
        ].copy()
        all_dropped_save['drop_reason'] = ''
        all_dropped_save.loc[
            all_dropped_save['user_id'].isin(dropped_decodeme['user_id']),
            'drop_reason'] = 'decodeme_or_unknown'
        all_dropped_save.loc[
            all_dropped_save['user_id'].isin(dropped_ftdna['user_id']),
            'drop_reason'] = 'ftdna_main_pipeline_excluded'
        all_dropped_save.loc[
            all_dropped_save['user_id'].isin(dropped_tiny['user_id']),
            'drop_reason'] = 'too_few_lines'
        all_dropped_save.to_csv(dropped_path, index=False)
        print(f"\nSaved dropped file log to {dropped_path}")

    # Save manual review list
    if len(manual_review) > 0:
        manual_review_path = DATA_DIR / "manual_review_files.csv"
        manual_review_cols = [
            'user_id', 'genotype_filename', 'filepath', 'genotype_format',
            'format_detected', 'delimiter_mode', 'n_data_lines',
            'n_rsid', 'n_internal_id', 'build_hint', 'first_snp', 'tier0', 'tier1'
        ]
        manual_review[[c for c in manual_review_cols
                        if c in manual_review.columns]].to_csv(
            manual_review_path, index=False
        )
        print(f"Saved manual-review file list to {manual_review_path}")

    # Save pipeline manifest
    pipeline_cols = [
        'user_id', 'genotype_filename', 'filepath',
        'genotype_format', 'format_detected', 'delimiter_mode',
        'build_hint', 'needs_liftover', 'n_data_lines',
        'n_rsid', 'n_internal_id', 'chrom_sex', 'raw_ancestry', 'tier0', 'tier1'
    ]
    pipeline_cols = [c for c in pipeline_cols if c in pipeline.columns]
    pipeline[pipeline_cols].to_csv(PIPELINE_MANIFEST_CSV, index=False)
    print(f"Saved pipeline manifest ({len(pipeline)} files) to "
          f"{PIPELINE_MANIFEST_CSV}")



if __name__ == "__main__":
    main()

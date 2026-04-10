import pandas as pd

from config import RAW_ANCESTRY_CSV, MAPPING_CSV, GROUPINGS_CSV, EXPECTED_TIER0



def main():
    raw = pd.read_csv(RAW_ANCESTRY_CSV)
    raw = raw.rename(columns={"value": "raw_ancestry"})
    raw["raw_ancestry"] = raw["raw_ancestry"].str.strip()

    print(f"Raw ancestry rows: {len(raw)}")
    print(f"Unique user_ids:   {raw['user_id'].nunique()}")
    print(f"Unique raw values: {raw['raw_ancestry'].nunique()}")

    mapping = pd.read_csv(MAPPING_CSV)
    mapping["raw_ancestry"] = mapping["raw_ancestry"].str.strip()

    print(f"\nMapping table rows: {len(mapping)}")
    print(f"Mapping columns:    {list(mapping.columns)}")

    # No duplicate raw_ancestry entries (would cause row explosion)
    dup_mappings = mapping[mapping.duplicated(subset=["raw_ancestry"], keep=False)]
    if len(dup_mappings) > 0:
        print(f"WARNING: {len(dup_mappings)} duplicate raw_ancestry entries:")
        print(dup_mappings[["raw_ancestry", "tier0", "tier1"]].to_string())
    else:
        print("Mapping table: no duplicate raw_ancestry entries")

    # Tier0 vocabulary check
    actual_tier0 = set(mapping["tier0"].dropna().unique())
    unexpected_tier0 = actual_tier0 - EXPECTED_TIER0
    if unexpected_tier0:
        print(f"WARNING: unexpected tier0 values: {unexpected_tier0}")
    else:
        print(f"Tier0 vocabulary check passed: {sorted(actual_tier0)}")

    print(f"\nTier1 values in mapping ({mapping['tier1'].nunique()}):")
    print(mapping["tier1"].value_counts().to_string())

    raw_values = set(raw["raw_ancestry"].unique())
    mapped_values = set(mapping["raw_ancestry"].unique())

    not_in_mapping = raw_values - mapped_values
    not_in_data = mapped_values - raw_values

    if not_in_mapping:
        print(f"WARNING: {len(not_in_mapping)} raw values have NO mapping:")
        for v in sorted(not_in_mapping):
            n = raw[raw["raw_ancestry"] == v].shape[0]
            print(f"  \"{v}\" ({n} rows)")
    else:
        print("All raw ancestry values have a mapping entry")

    if not_in_data:
        print(f"\nNOTE: {len(not_in_data)} mapping entries unused (harmless)")

    lower_map = {}
    for val in raw_values:
        key = val.lower()
        lower_map.setdefault(key, []).append(val)

    case_variants = {k: v for k, v in lower_map.items() if len(v) > 1}
    if case_variants:
        print(f"Case variants detected ({len(case_variants)} groups):")
        for key, variants in sorted(case_variants.items()):
            print(f"  {variants}")
    else:
        print("No case variants detected")

    final = raw.merge(mapping, on="raw_ancestry", how="left")

    unmapped = final[final["tier0"].isna()]
    print(f"Total rows:    {len(final)}")
    print(f"Unmapped rows: {len(unmapped)}")

    if len(unmapped) > 0:
        print("\nWARNING — unmapped raw_ancestry values:")
        for val in unmapped["raw_ancestry"].unique():
            n = unmapped[unmapped["raw_ancestry"] == val].shape[0]
            print(f"  \"{val}\" ({n} rows)")
        print("\nAdd these to the mapping table and re-run.")

    if len(final) != len(raw):
        print(f"\nWARNING: row count changed ({len(raw)} → {len(final)})")
    else:
        print(f"Row count preserved: {len(final)}")

    dup_users = final[final.duplicated(subset=["user_id"], keep=False)]
    if len(dup_users) > 0:
        print(f"\nWARNING: {len(dup_users)} duplicate user_ids in output")
    else:
        print("No duplicate user_ids")

    print("\n" + "=" * 60)
    print("GROUPING SUMMARY")
    print("=" * 60)

    print(f"\nTier 0 counts:")
    print(final["tier0"].value_counts(dropna=False).to_string())

    print(f"\nTier 0 × Tier 1:")
    cross = (final.groupby(["tier0", "tier1"])
             .size()
             .reset_index(name="n")
             .sort_values(["tier0", "n"], ascending=[True, False]))
    print(cross.to_string(index=False))

    final.to_csv(GROUPINGS_CSV, index=False)
    print(f"\nSaved: {GROUPINGS_CSV}")
    print(f"Columns: {list(final.columns)}")


if __name__ == "__main__":
    main()

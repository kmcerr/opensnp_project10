"""
Project 10 — Network construction and analysis helpers.

Consolidated from networkjoin_v2 and step7_dedup_v2.
"""

import pandas as pd
import networkx as nx


def build_graph(node_df: pd.DataFrame,
                edge_df: pd.DataFrame,
                weight_col: str,
                label: str = "") -> nx.Graph:
    """
    Build a networkx graph from node and edge DataFrames.

    Parameters:
        node_df:     DataFrame with node_id, tier0, tier1, raw_ancestry,
                     genotype_format
        edge_df:     DataFrame with IID1, IID2, and weight_col
        weight_col:  column to use as edge weight ("IBS2_prop" or "PI_HAT")
        label:       human-readable label for printing

    Returns:
        nx.Graph with node attributes and weighted edges.
    """
    G = nx.Graph()

    # Add all nodes with attributes
    for _, row in node_df.iterrows():
        G.add_node(
            str(row["node_id"]),
            tier0=row.get("tier0", ""),
            tier1=row.get("tier1", ""),
            raw_ancestry=row.get("raw_ancestry", ""),
            genotype_format=row.get("genotype_format", ""),
        )

    # Add edges
    for _, row in edge_df.iterrows():
        G.add_edge(
            str(row["IID1"]),
            str(row["IID2"]),
            weight=row[weight_col],
        )

    # Print summary
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    components = list(nx.connected_components(G))
    n_components = len(components)
    largest_cc = max((len(c) for c in components), default=0)
    n_isolates = sum(1 for _ in nx.isolates(G))

    print(f"\nGraph ({label}, weight={weight_col}):")
    print(f"  Nodes:                {n_nodes}")
    print(f"  Edges:                {n_edges}")
    print(f"  Connected components: {n_components}")
    print(f"  Largest component:    {largest_cc}")
    print(f"  Isolates:             {n_isolates}")

    # Tier0 breakdown
    isolate_ids = set(nx.isolates(G))
    connected_ids = set(G.nodes()) - isolate_ids

    tier0_connected = [G.nodes[n].get("tier0", "") for n in connected_ids]
    tier0_isolate = [G.nodes[n].get("tier0", "") for n in isolate_ids]

    if tier0_connected:
        print(f"\n  Connected nodes by tier0:")
        for val, cnt in pd.Series(tier0_connected).value_counts().items():
            print(f"    {val}: {cnt}")

    if tier0_isolate:
        print(f"\n  Isolates by tier0:")
        for val, cnt in pd.Series(tier0_isolate).value_counts().items():
            print(f"    {val}: {cnt}")

    return G


def concordance_summary(edge_df: pd.DataFrame, label: str):
    """
    Summarize ancestry concordance for an edge set.

    Prints the fraction of edges connecting same-ancestry vs
    different-ancestry pairs, plus a cross-tabulation of the top
    tier0 pair types.
    """
    if len(edge_df) == 0:
        print(f"\n{label}: no edges")
        return

    n = len(edge_df)
    n_same = int(edge_df["same_tier0"].sum())
    n_diff = int((~edge_df["same_tier0"]).sum())
    n_na = int(edge_df["same_tier0"].isna().sum())

    print(f"\n{label}:")
    print(f"  Total edges:              {n:,}")
    print(f"  Same tier0 (concordant):  {n_same:,} ({100 * n_same / n:.1f}%)")
    print(f"  Diff tier0 (discordant):  {n_diff:,} ({100 * n_diff / n:.1f}%)")
    if n_na > 0:
        print(f"  Missing tier0 (NA):       {n_na:,}")

    # Cross-tabulation
    if "tier0_1" in edge_df.columns and "tier0_2" in edge_df.columns:
        pair_labels = edge_df.apply(
            lambda r: " × ".join(
                sorted([str(r["tier0_1"]), str(r["tier0_2"])])
            ),
            axis=1,
        )
        print(f"\n  Top tier0 pair types:")
        for pair, cnt in pair_labels.value_counts().head(10).items():
            print(f"    {pair}: {cnt:,}")

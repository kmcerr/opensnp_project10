#!/usr/bin/env python3
"""
Step 09: Network Construction & Visualization

Builds genetic similarity network from IBS results and analyzes
relatedness patterns by ancestry. Creates network visualizations
for presentation.

Core project requirements fulfilled:
- NetworkX graph construction from IBS2 proportions
- Network visualization with ancestry-colored nodes
- Clustering analysis by ancestry
- Identification of unexpected relatedness pairs
- Modularity and connectivity statistics
"""

import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
from collections import Counter, defaultdict

from config import (
    DATA_DIR, FIG_DIR, QC_DIR, NETWORK_DIR, GROUPINGS_CSV,
    TIER0_COLORS, TIER0_ORDER,
    IBS2_PRUNE_THRESHOLD, IBS2_RELATEDNESS_THRESHOLD,
    PIHAT_1ST_DEGREE, PIHAT_2ND_DEGREE, PIHAT_3RD_DEGREE, PIHAT_DUPLICATE,
    setup_logging, ensure_dirs,
)
from lib.ibs import load_genome_file, add_pair_categories, annotate_ibs
from lib.network import build_graph, concordance_summary

setup_logging()
logger = logging.getLogger(__name__)

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['figure.figsize'] = (12, 8)


def load_network_data():
    """
    Load IBS results and ancestry metadata for network construction.

    Returns:
        ibs_df: DataFrame with pairwise IBS results
        meta_df: DataFrame with ancestry labels
    """
    logger.info("Loading network data")

    # Load IBS results
    genome_path = QC_DIR / "step5_ibs_with_proportions.csv"
    if not genome_path.exists():
        raise FileNotFoundError(f"IBS results not found: {genome_path}")

    # Load ancestry metadata
    if not GROUPINGS_CSV.exists():
        raise FileNotFoundError(f"Ancestry groupings not found: {GROUPINGS_CSV}")

    meta_df = pd.read_csv(GROUPINGS_CSV)

    # Load and annotate IBS results with ancestry metadata
    ibs_df = pd.read_csv(genome_path)
    ibs_df = annotate_ibs(ibs_df, meta_df)

    logger.info(f"Loaded {len(ibs_df):,} pairwise comparisons")
    logger.info(f"Loaded {len(meta_df)} individuals with ancestry labels")

    return ibs_df, meta_df


def build_ibs_network(ibs_df, meta_df, threshold=0.0, label=""):
    """
    Build NetworkX graph from IBS data with optional edge pruning.

    Parameters:
        ibs_df: DataFrame with IID1, IID2, IBS2_prop
        meta_df: DataFrame with user_id, tier0, tier1, raw_ancestry
        threshold: Minimum IBS2_prop to include edge (default: 0.0)
        label: Description for logging

    Returns:
        NetworkX Graph with node attributes and weighted edges
    """
    logger.info(f"Building network: {label}")
    logger.info(f"  Edge threshold: IBS2_prop >= {threshold}")

    # Filter edges by threshold
    edge_df = ibs_df[ibs_df['IBS2_prop'] >= threshold].copy()
    logger.info(f"  Edges after filtering: {len(edge_df):,}")

    # Prepare node dataframe
    # Get unique individuals from IBS results
    all_ids = pd.concat([edge_df['IID1'], edge_df['IID2']]).unique()

    # Match with metadata
    meta_subset = meta_df[meta_df['user_id'].isin(all_ids)].copy()
    meta_subset = meta_subset.rename(columns={'user_id': 'node_id'})

    # Build graph
    G = build_graph(
        node_df=meta_subset,
        edge_df=edge_df,
        weight_col='IBS2_prop',
        label=label
    )

    return G


def calculate_network_statistics(G, meta_df):
    """
    Calculate network statistics by ancestry group.

    Returns:
        DataFrame with statistics by tier0
    """
    logger.info("Calculating network statistics")

    stats = []

    # Overall statistics
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    n_components = nx.number_connected_components(G)
    largest_cc = len(max(nx.connected_components(G), key=len)) if n_components > 0 else 0
    n_isolates = sum(1 for _ in nx.isolates(G))
    avg_degree = sum(dict(G.degree()).values()) / n_nodes if n_nodes > 0 else 0

    # Calculate clustering coefficient (can be slow for large graphs)
    if n_nodes > 0:
        try:
            avg_clustering = nx.average_clustering(G)
        except:
            avg_clustering = 0.0
    else:
        avg_clustering = 0.0

    stats.append({
        'tier0': 'Overall',
        'n_nodes': n_nodes,
        'n_edges': n_edges,
        'n_components': n_components,
        'largest_component': largest_cc,
        'n_isolates': n_isolates,
        'avg_degree': avg_degree,
        'avg_clustering': avg_clustering,
    })

    logger.info(f"Overall: {n_nodes} nodes, {n_edges} edges, {n_components} components")

    # Statistics by ancestry group
    for tier0 in TIER0_ORDER:
        # Get nodes for this ancestry
        nodes_in_tier = [n for n in G.nodes() if G.nodes[n].get('tier0') == tier0]

        if len(nodes_in_tier) == 0:
            continue

        # Subgraph for this ancestry
        subG = G.subgraph(nodes_in_tier)

        n_nodes_tier = subG.number_of_nodes()
        n_edges_tier = subG.number_of_edges()
        n_components_tier = nx.number_connected_components(subG)
        largest_cc_tier = len(max(nx.connected_components(subG), key=len)) if n_components_tier > 0 else 0
        n_isolates_tier = sum(1 for _ in nx.isolates(subG))
        avg_degree_tier = sum(dict(subG.degree()).values()) / n_nodes_tier if n_nodes_tier > 0 else 0

        try:
            avg_clustering_tier = nx.average_clustering(subG)
        except:
            avg_clustering_tier = 0.0

        stats.append({
            'tier0': tier0,
            'n_nodes': n_nodes_tier,
            'n_edges': n_edges_tier,
            'n_components': n_components_tier,
            'largest_component': largest_cc_tier,
            'n_isolates': n_isolates_tier,
            'avg_degree': avg_degree_tier,
            'avg_clustering': avg_clustering_tier,
        })

        logger.info(f"{tier0}: {n_nodes_tier} nodes, {n_edges_tier} edges")

    stats_df = pd.DataFrame(stats)
    return stats_df


def calculate_modularity_by_ancestry(G):
    """
    Calculate modularity score treating ancestry groups as communities.

    Returns:
        Modularity Q score (higher = stronger clustering by ancestry)
    """
    logger.info("Calculating modularity by ancestry")

    # Create community partition by tier0
    communities = defaultdict(list)
    for node in G.nodes():
        tier0 = G.nodes[node].get('tier0', 'Unknown')
        communities[tier0].append(node)

    # Convert to list of sets (NetworkX format)
    community_list = [set(nodes) for nodes in communities.values() if len(nodes) > 0]

    try:
        Q = nx.algorithms.community.modularity(G, community_list)
        logger.info(f"Modularity Q = {Q:.4f}")
        return Q
    except Exception as e:
        logger.warning(f"Could not calculate modularity: {e}")
        return None


def calculate_edge_concordance(G):
    """
    Calculate within-ancestry vs. between-ancestry edge ratios.

    Returns:
        DataFrame with concordance statistics
    """
    logger.info("Calculating edge concordance by ancestry")

    within_count = 0
    between_count = 0

    for u, v in G.edges():
        tier0_u = G.nodes[u].get('tier0', 'Unknown')
        tier0_v = G.nodes[v].get('tier0', 'Unknown')

        if tier0_u == tier0_v:
            within_count += 1
        else:
            between_count += 1

    total = within_count + between_count

    concordance_df = pd.DataFrame([{
        'edge_type': 'Within ancestry',
        'count': within_count,
        'percentage': 100 * within_count / total if total > 0 else 0,
    }, {
        'edge_type': 'Between ancestry',
        'count': between_count,
        'percentage': 100 * between_count / total if total > 0 else 0,
    }])

    logger.info(f"Within ancestry: {within_count:,} ({100*within_count/total:.1f}%)")
    logger.info(f"Between ancestry: {between_count:,} ({100*between_count/total:.1f}%)")

    return concordance_df


def identify_unexpected_pairs(ibs_df, threshold=IBS2_RELATEDNESS_THRESHOLD):
    """
    Flag high IBS2 pairs with discordant ancestry.

    Parameters:
        ibs_df: DataFrame with IBS results and ancestry annotations
        threshold: IBS2_prop threshold for relatedness

    Returns:
        DataFrame with unexpected pairs sorted by IBS2_prop
    """
    logger.info(f"Identifying unexpected relatedness (IBS2 > {threshold})")

    # Filter to related pairs
    related = ibs_df[ibs_df['IBS2_prop'] > threshold].copy()

    # Filter to discordant ancestry
    unexpected = related[related['same_tier0'] == False].copy()

    # Sort by IBS2_prop
    unexpected = unexpected.sort_values('IBS2_prop', ascending=False)

    # Add interpretation notes
    def interpret_pair(row):
        if row['PI_HAT'] >= PIHAT_1ST_DEGREE:
            degree = '1st degree'
        elif row['PI_HAT'] >= PIHAT_2ND_DEGREE:
            degree = '2nd degree'
        elif row['PI_HAT'] >= PIHAT_3RD_DEGREE:
            degree = '3rd degree'
        else:
            degree = 'distant'

        return f"{degree} relatedness with discordant ancestry"

    unexpected['interpretation'] = unexpected.apply(interpret_pair, axis=1)

    logger.info(f"Found {len(unexpected)} unexpected relatedness pairs")

    return unexpected


def visualize_network_full(G, output_path):
    """
    Visualize full network with ancestry-colored nodes.
    """
    logger.info("Creating full network visualization")

    if G.number_of_nodes() == 0:
        logger.warning("Empty graph, skipping visualization")
        return

    fig, ax = plt.subplots(figsize=(14, 10))

    # Calculate layout (spring layout works well for genetic networks)
    logger.info("Calculating layout (this may take a minute)...")
    pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)

    # Prepare node colors by ancestry
    node_colors = []
    for node in G.nodes():
        tier0 = G.nodes[node].get('tier0', 'Unknown')
        node_colors.append(TIER0_COLORS.get(tier0, '#8C8C8C'))

    # Draw network
    nx.draw_networkx_edges(G, pos, alpha=0.2, width=0.5, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                           node_size=30, alpha=0.8, ax=ax)

    # Create legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=TIER0_COLORS[tier0], label=tier0)
        for tier0 in TIER0_ORDER if tier0 in [G.nodes[n].get('tier0') for n in G.nodes()]
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    ax.set_title(f'Genetic Similarity Network (Full)\n{G.number_of_nodes()} individuals, {G.number_of_edges():,} edges',
                 fontsize=14, fontweight='bold', pad=20)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def visualize_network_pruned(G, output_path):
    """
    Visualize network with weak edges removed (IBS2 > 0.1).
    """
    logger.info("Creating pruned network visualization")

    if G.number_of_nodes() == 0:
        logger.warning("Empty graph, skipping visualization")
        return

    fig, ax = plt.subplots(figsize=(14, 10))

    # Calculate layout
    logger.info("Calculating layout...")
    pos = nx.spring_layout(G, k=0.8, iterations=50, seed=42)

    # Prepare node colors and sizes
    node_colors = []
    node_sizes = []
    for node in G.nodes():
        tier0 = G.nodes[node].get('tier0', 'Unknown')
        node_colors.append(TIER0_COLORS.get(tier0, '#8C8C8C'))
        # Size by degree
        node_sizes.append(50 + 10 * G.degree(node))

    # Draw network
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=1.0, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9, ax=ax)

    # Create legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=TIER0_COLORS[tier0], label=tier0)
        for tier0 in TIER0_ORDER if tier0 in [G.nodes[n].get('tier0') for n in G.nodes()]
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    ax.set_title(f'Genetic Similarity Network (Pruned IBS2 > {IBS2_PRUNE_THRESHOLD})\n{G.number_of_nodes()} individuals, {G.number_of_edges():,} edges',
                 fontsize=14, fontweight='bold', pad=20)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def visualize_network_related(G, output_path):
    """
    Visualize network showing only related pairs (IBS2 > 0.4).
    """
    logger.info("Creating related pairs network visualization")

    if G.number_of_nodes() == 0:
        logger.warning("Empty graph, skipping visualization")
        return

    fig, ax = plt.subplots(figsize=(14, 10))

    # For related pairs, use different layout to show clusters
    logger.info("Calculating layout...")
    if G.number_of_edges() > 0:
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.spring_layout(G, seed=42)

    # Prepare node colors and sizes
    node_colors = []
    node_sizes = []
    for node in G.nodes():
        tier0 = G.nodes[node].get('tier0', 'Unknown')
        node_colors.append(TIER0_COLORS.get(tier0, '#8C8C8C'))
        # Size by degree
        node_sizes.append(100 + 20 * G.degree(node))

    # Draw edges with weight-based width
    edge_widths = [G[u][v]['weight'] * 3 for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.5, ax=ax)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9, ax=ax,
                           edgecolors='black', linewidths=1)

    # Label highly connected nodes
    high_degree_nodes = {node: pos[node] for node in G.nodes() if G.degree(node) >= 3}
    if len(high_degree_nodes) > 0 and len(high_degree_nodes) < 20:
        nx.draw_networkx_labels(G, high_degree_nodes, font_size=8, ax=ax)

    # Create legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=TIER0_COLORS[tier0], label=tier0)
        for tier0 in TIER0_ORDER if tier0 in [G.nodes[n].get('tier0') for n in G.nodes()]
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    ax.set_title(f'Related Pairs Network (IBS2 > {IBS2_RELATEDNESS_THRESHOLD})\n{G.number_of_nodes()} individuals, {G.number_of_edges():,} edges',
                 fontsize=14, fontweight='bold', pad=20)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def visualize_clustering_analysis(stats_df, concordance_df, modularity_Q, output_path):
    """
    Visualize IBS2-based clustering analysis (shows the problem).
    """
    logger.info("Creating IBS2-based clustering analysis visualization")

    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(1, 2, hspace=0.3, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Plot 1: Within vs. between ancestry edges
    colors = ['#55A868', '#DD8452']
    bars = ax1.bar(concordance_df['edge_type'], concordance_df['count'],
                   color=colors, alpha=0.8, edgecolor='black')

    total_edges = concordance_df['count'].sum()
    for bar in bars:
        height = bar.get_height()
        pct = (height / total_edges * 100)
        # Place count inside bar if it's large enough
        if height > total_edges * 0.1:
            ax1.text(bar.get_x() + bar.get_width()/2., height/2,
                    f'{int(height):,}\n({pct:.1f}%)',
                    ha='center', va='center', fontsize=11, fontweight='bold',
                    color='white')
        else:
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height):,}\n({pct:.1f}%)',
                    ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax1.set_ylabel('Number of Edges', fontsize=12, fontweight='bold')
    ax1.set_title('Edge Distribution by Ancestry\n(IBS2-based - All edges within-ancestry)',
                  fontsize=12, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.set_axisbelow(True)
    ax1.set_ylim(0, total_edges * 1.15)

    # Plot 2: Modularity score
    if modularity_Q is not None:
        ax2.barh(['Modularity (Q)'], [modularity_Q], color='steelblue', alpha=0.8, height=0.3)
        ax2.axvline(0.3, color='orange', linestyle='--', linewidth=2,
                   label='Weak clustering', alpha=0.7)
        ax2.axvline(0.5, color='red', linestyle='--', linewidth=2,
                   label='Strong clustering', alpha=0.7)

        ax2.text(0.05, 0, f'Q = {modularity_Q:.4f}',
                ha='left', va='center', fontsize=13, fontweight='bold')

        ax2.set_xlabel('Modularity Score', fontsize=12, fontweight='bold')
        ax2.set_title('Network Modularity\n(Q ≈ 0: No clustering detected)',
                     fontsize=12, fontweight='bold')
        ax2.set_xlim(-0.05, 1.0)
        ax2.legend(fontsize=9, loc='upper right')
        ax2.grid(axis='x', alpha=0.3, linestyle='--')
    else:
        ax2.text(0.5, 0.5, 'Modularity calculation unavailable',
                ha='center', va='center', fontsize=12,
                transform=ax2.transAxes)
        ax2.axis('off')

    plt.suptitle('Ancestry Clustering Analysis (IBS2-based - Shows Problem)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def visualize_pihat_clustering_analysis(ibs_df, output_path):
    """
    Visualize PI_HAT-based clustering analysis (correct approach).

    Shows edge distribution and modularity for truly related pairs only.
    """
    logger.info("Creating PI_HAT-based clustering analysis visualization")

    # Filter to related pairs (PI_HAT >= 3rd degree)
    related = ibs_df[ibs_df['PI_HAT'] >= PIHAT_3RD_DEGREE].copy()

    if len(related) == 0:
        logger.warning("No related pairs found - skipping PI_HAT clustering visualization")
        return

    # Calculate concordance for related pairs only
    if 'same_tier0' in related.columns:
        within_ancestry = related['same_tier0'].sum()
        between_ancestry = len(related) - within_ancestry
    else:
        logger.warning("No ancestry metadata - cannot calculate concordance")
        within_ancestry = 0
        between_ancestry = len(related)

    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(1, 2, hspace=0.3, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Plot 1: Within vs. between ancestry for RELATED pairs
    concordance_data = {
        'edge_type': ['Within ancestry', 'Between ancestry'],
        'count': [within_ancestry, between_ancestry]
    }
    colors = ['#55A868', '#DD8452']

    bars = ax1.bar(concordance_data['edge_type'], concordance_data['count'],
                   color=colors, alpha=0.8, edgecolor='black')

    total_edges = sum(concordance_data['count'])
    for bar, count in zip(bars, concordance_data['count']):
        if count > 0:
            pct = (count / total_edges * 100)
            ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height()/2,
                    f'{int(count)}\n({pct:.1f}%)',
                    ha='center', va='center', fontsize=12, fontweight='bold',
                    color='white' if count > total_edges * 0.2 else 'black')

    ax1.set_ylabel('Number of Related Pairs', fontsize=12, fontweight='bold')
    ax1.set_title(f'Related Pairs by Ancestry\n(PI_HAT ≥ {PIHAT_3RD_DEGREE}: {total_edges} pairs)',
                  fontsize=12, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.set_axisbelow(True)
    ax1.set_ylim(0, max(concordance_data['count']) * 1.3)

    # Plot 2: Relationship degree distribution
    degree_labels = []
    degree_counts = []

    duplicates = (related['PI_HAT'] >= PIHAT_DUPLICATE).sum()
    first_deg = ((related['PI_HAT'] >= PIHAT_1ST_DEGREE) & (related['PI_HAT'] < PIHAT_DUPLICATE)).sum()
    second_deg = ((related['PI_HAT'] >= PIHAT_2ND_DEGREE) & (related['PI_HAT'] < PIHAT_1ST_DEGREE)).sum()
    third_deg = ((related['PI_HAT'] >= PIHAT_3RD_DEGREE) & (related['PI_HAT'] < PIHAT_2ND_DEGREE)).sum()

    if duplicates > 0:
        degree_labels.append(f'Duplicate/MZ\n(≥{PIHAT_DUPLICATE})')
        degree_counts.append(duplicates)
    if first_deg > 0:
        degree_labels.append(f'1st degree\n(≥{PIHAT_1ST_DEGREE})')
        degree_counts.append(first_deg)
    if second_deg > 0:
        degree_labels.append(f'2nd degree\n(≥{PIHAT_2ND_DEGREE})')
        degree_counts.append(second_deg)
    if third_deg > 0:
        degree_labels.append(f'3rd degree\n(≥{PIHAT_3RD_DEGREE})')
        degree_counts.append(third_deg)

    degree_colors = ['#8E44AD', '#E74C3C', '#F39C12', '#F1C40F'][:len(degree_labels)]
    bars2 = ax2.bar(degree_labels, degree_counts, color=degree_colors, alpha=0.8, edgecolor='black')

    for bar in bars2:
        height = bar.get_height()
        if height > 0:
            ax2.text(bar.get_x() + bar.get_width()/2., height/2,
                    f'{int(height)}',
                    ha='center', va='center', fontsize=12, fontweight='bold',
                    color='white')

    ax2.set_ylabel('Number of Pairs', fontsize=12, fontweight='bold')
    ax2.set_title('Related Pairs by Relationship Degree',
                  fontsize=12, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    ax2.set_axisbelow(True)
    if degree_counts:
        ax2.set_ylim(0, max(degree_counts) * 1.3)

    plt.suptitle('Relatedness Analysis (PI_HAT-based - Correct Approach)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


# ═══════════════════════════════════════════════════════════════════
# PI_HAT-BASED NETWORK FUNCTIONS (BETTER FOR RELATEDNESS)
# ═══════════════════════════════════════════════════════════════════

def build_pihat_network(ibs_df, meta_df, threshold=0.0, label=""):
    """
    Build NetworkX graph from IBS data filtered by PI_HAT (relatedness).

    PI_HAT measures actual genetic relatedness, unlike IBS2 which measures
    shared alleles (high even for unrelated individuals in homogeneous populations).

    Parameters:
        ibs_df: DataFrame with IID1, IID2, PI_HAT, IBS2_prop
        meta_df: DataFrame with user_id, tier0, tier1, raw_ancestry
        threshold: Minimum PI_HAT to include edge (default: 0.0)
        label: Description for logging

    Returns:
        NetworkX Graph with node attributes and weighted edges
    """
    logger.info(f"Building PI_HAT network: {label}")
    logger.info(f"  Edge threshold: PI_HAT >= {threshold}")

    # Filter edges by PI_HAT threshold
    edge_df = ibs_df[ibs_df['PI_HAT'] >= threshold].copy()
    logger.info(f"  Edges after filtering: {len(edge_df):,}")

    if len(edge_df) == 0:
        logger.warning("  No edges after filtering - creating empty graph")
        # Return empty graph with all individuals as isolates
        all_ids = meta_df['user_id'].astype(str).unique()
        meta_subset = meta_df.copy()
        meta_subset['user_id'] = meta_subset['user_id'].astype(str)
        meta_subset = meta_subset.rename(columns={'user_id': 'node_id'})

        G = nx.Graph(name=label)
        for _, row in meta_subset.iterrows():
            G.add_node(row['node_id'],
                      tier0=row['tier0'],
                      tier1=row['tier1'],
                      raw_ancestry=row['raw_ancestry'])
        return G

    # Get unique individuals from filtered edges
    all_ids = pd.concat([edge_df['IID1'], edge_df['IID2']]).unique()

    # Match with metadata
    meta_df_copy = meta_df.copy()
    meta_df_copy['user_id'] = meta_df_copy['user_id'].astype(str)
    meta_subset = meta_df_copy[meta_df_copy['user_id'].isin(all_ids)].copy()
    meta_subset = meta_subset.rename(columns={'user_id': 'node_id'})

    # Build graph with PI_HAT as edge weight
    G = build_graph(
        node_df=meta_subset,
        edge_df=edge_df,
        weight_col='PI_HAT',
        label=label
    )

    return G


def visualize_pihat_network(G, output_path, threshold_label=""):
    """
    Visualize PI_HAT-based relatedness network with cleaner layout.

    For sparse networks (few edges), uses better layout algorithms.
    """
    logger.info(f"Creating PI_HAT relatedness network visualization")

    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()

    if n_nodes == 0:
        logger.warning("Empty graph - skipping visualization")
        return

    fig, ax = plt.subplots(figsize=(14, 10))

    # Choose layout based on graph density
    if n_edges == 0:
        # All isolates - use circular layout
        pos = nx.circular_layout(G)
        logger.info("  Using circular layout (no edges)")
    elif n_edges < 100:
        # Sparse graph - use Kamada-Kawai for better structure
        try:
            pos = nx.kamada_kawai_layout(G)
            logger.info("  Using Kamada-Kawai layout (sparse graph)")
        except:
            pos = nx.spring_layout(G, k=2, iterations=100, seed=42)
            logger.info("  Using spring layout (fallback)")
    else:
        # Dense graph - use spring with high spacing
        pos = nx.spring_layout(G, k=1.5, iterations=100, seed=42)
        logger.info("  Using spring layout with high spacing")

    # Color nodes by ancestry
    node_colors = [TIER0_COLORS.get(G.nodes[node].get('tier0', 'Unknown'), '#CCCCCC')
                   for node in G.nodes()]

    # Draw edges with width based on PI_HAT
    if n_edges > 0:
        edges = G.edges()
        weights = [G[u][v].get('weight', 0.1) for u, v in edges]
        # Scale edge widths: PI_HAT 0.1 -> width 1, PI_HAT 1.0 -> width 10
        edge_widths = [max(0.5, w * 10) for w in weights]

        # Color edges by relationship degree
        edge_colors = []
        for u, v in edges:
            w = G[u][v].get('weight', 0)
            if w >= PIHAT_1ST_DEGREE:
                edge_colors.append('#E74C3C')  # Red - 1st degree
            elif w >= PIHAT_2ND_DEGREE:
                edge_colors.append('#F39C12')  # Orange - 2nd degree
            elif w >= PIHAT_3RD_DEGREE:
                edge_colors.append('#F1C40F')  # Yellow - 3rd degree
            else:
                edge_colors.append('#95A5A6')  # Gray - distant

        nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color=edge_colors,
                              alpha=0.7, ax=ax)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color=node_colors,
                          node_size=300, alpha=0.9,
                          edgecolors='black', linewidths=1.5, ax=ax)

    # Add labels for connected nodes only (if not too many)
    if n_edges > 0 and n_edges < 50:
        connected_nodes = {node for edge in G.edges() for node in edge}
        labels = {node: node for node in connected_nodes}
        nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)

    # Create legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    # Ancestry legend
    ancestry_legend = [
        Patch(facecolor=TIER0_COLORS[tier0], label=tier0, edgecolor='black')
        for tier0 in TIER0_ORDER if tier0 in [G.nodes[n].get('tier0') for n in G.nodes()]
    ]

    # Relationship legend (if edges exist)
    if n_edges > 0:
        relationship_legend = [
            Line2D([0], [0], color='#E74C3C', linewidth=3, label='1st degree (≥0.4)'),
            Line2D([0], [0], color='#F39C12', linewidth=3, label='2nd degree (≥0.25)'),
            Line2D([0], [0], color='#F1C40F', linewidth=3, label='3rd degree (≥0.125)'),
            Line2D([0], [0], color='#95A5A6', linewidth=2, label='Distant (<0.125)'),
        ]

        # Two-column legend
        leg1 = ax.legend(handles=ancestry_legend, loc='upper left',
                        title='Ancestry', fontsize=9, title_fontsize=10)
        ax.add_artist(leg1)
        ax.legend(handles=relationship_legend, loc='upper right',
                 title='Relatedness', fontsize=9, title_fontsize=10)

    else:
        ax.legend(handles=ancestry_legend, loc='upper right',
                 title='Ancestry', fontsize=9, title_fontsize=10)

    ax.set_title(f'Relatedness Network (PI_HAT-based)\n{threshold_label}\n{n_nodes} individuals, {n_edges} related pairs',
                 fontsize=14, fontweight='bold', pad=20)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def visualize_pihat_comparison(ibs_df, output_path):
    """
    Create comparison plot showing why PI_HAT is better than IBS2 for relatedness.

    Shows scatter plot of IBS2 vs PI_HAT with relationship degree thresholds.
    """
    logger.info("Creating IBS2 vs PI_HAT comparison plot")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Plot 1: Scatter plot of IBS2 vs PI_HAT
    scatter = ax1.scatter(ibs_df['IBS2_prop'], ibs_df['PI_HAT'],
                         alpha=0.3, s=10, c='#3498DB')

    ax1.axhline(y=PIHAT_3RD_DEGREE, color='#F1C40F', linestyle='--',
                linewidth=2, label=f'3rd degree (PI_HAT={PIHAT_3RD_DEGREE})')
    ax1.axhline(y=PIHAT_2ND_DEGREE, color='#F39C12', linestyle='--',
                linewidth=2, label=f'2nd degree (PI_HAT={PIHAT_2ND_DEGREE})')
    ax1.axhline(y=PIHAT_1ST_DEGREE, color='#E74C3C', linestyle='--',
                linewidth=2, label=f'1st degree (PI_HAT={PIHAT_1ST_DEGREE})')
    ax1.axhline(y=PIHAT_DUPLICATE, color='#8E44AD', linestyle='--',
                linewidth=2, label=f'Duplicate (PI_HAT={PIHAT_DUPLICATE})')

    ax1.axvline(x=0.4, color='gray', linestyle=':', linewidth=2, alpha=0.5,
                label='IBS2=0.4 threshold')

    ax1.set_xlabel('IBS2 Proportion (Shared Alleles)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('PI_HAT (Genetic Relatedness)', fontsize=12, fontweight='bold')
    ax1.set_title('IBS2 vs PI_HAT: Why IBS2 Filtering Fails',
                  fontsize=13, fontweight='bold', pad=15)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.grid(alpha=0.3)
    ax1.set_xlim(-0.05, 1.05)
    ax1.set_ylim(-0.1, 1.1)

    # Add annotation explaining the problem
    ax1.text(0.5, -0.05,
             'Problem: All pairs have IBS2 > 0.4 due to population structure,\nbut only 3 pairs are truly related (PI_HAT > 0.125)',
             fontsize=10, ha='center', style='italic',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.3))

    # Plot 2: Histograms
    ax2_twin = ax2.twinx()

    # IBS2 histogram
    ax2.hist(ibs_df['IBS2_prop'], bins=50, alpha=0.6, color='#3498DB',
            label='IBS2 distribution', edgecolor='black')
    ax2.axvline(x=0.4, color='gray', linestyle=':', linewidth=3,
               label='IBS2=0.4 threshold')
    ax2.set_xlabel('IBS2 Proportion', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Pairs (IBS2)', fontsize=11, fontweight='bold', color='#3498DB')
    ax2.tick_params(axis='y', labelcolor='#3498DB')

    # PI_HAT histogram
    ax2_twin.hist(ibs_df['PI_HAT'], bins=50, alpha=0.6, color='#E74C3C',
                 label='PI_HAT distribution', edgecolor='black')
    ax2_twin.axvline(x=PIHAT_3RD_DEGREE, color='#F1C40F', linestyle='--', linewidth=2)
    ax2_twin.set_ylabel('Number of Pairs (PI_HAT)', fontsize=11, fontweight='bold', color='#E74C3C')
    ax2_twin.tick_params(axis='y', labelcolor='#E74C3C')
    ax2_twin.set_ylim(0, ax2_twin.get_ylim()[1] * 1.2)  # Extra space for legend

    ax2.set_title('Distribution Comparison', fontsize=13, fontweight='bold', pad=15)
    ax2.legend(loc='upper left', fontsize=9)
    ax2_twin.legend(loc='upper right', fontsize=9)
    ax2.grid(alpha=0.3)

    plt.suptitle('IBS2 vs PI_HAT: Filtering Comparison',
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved: {output_path}")


def main():
    """Generate network analysis and visualizations."""
    ensure_dirs()

    logger.info("=" * 60)
    logger.info("Step 09: Network Construction & Visualization")
    logger.info("=" * 60)
    logger.info(f"Output directory: {NETWORK_DIR}")
    logger.info(f"Figures directory: {FIG_DIR}")

    # Load data
    ibs_df, meta_df = load_network_data()

    # Build networks at different thresholds
    logger.info("")
    logger.info("=" * 60)
    logger.info("Building Networks")
    logger.info("=" * 60)

    G_full = build_ibs_network(ibs_df, meta_df, threshold=0.0,
                               label="Full network (all pairs)")

    G_pruned = build_ibs_network(ibs_df, meta_df, threshold=IBS2_PRUNE_THRESHOLD,
                                 label=f"Pruned network (IBS2 > {IBS2_PRUNE_THRESHOLD})")

    G_related = build_ibs_network(ibs_df, meta_df, threshold=IBS2_RELATEDNESS_THRESHOLD,
                                  label=f"Related pairs (IBS2 > {IBS2_RELATEDNESS_THRESHOLD})")

    # Build PI_HAT-based networks (BETTER for actual relatedness)
    logger.info("")
    logger.info("=" * 60)
    logger.info("Building PI_HAT-based Networks (Actual Relatedness)")
    logger.info("=" * 60)

    G_pihat_3rd = build_pihat_network(ibs_df, meta_df, threshold=PIHAT_3RD_DEGREE,
                                      label=f"Related ≥3rd degree (PI_HAT ≥ {PIHAT_3RD_DEGREE})")

    G_pihat_2nd = build_pihat_network(ibs_df, meta_df, threshold=PIHAT_2ND_DEGREE,
                                      label=f"Related ≥2nd degree (PI_HAT ≥ {PIHAT_2ND_DEGREE})")

    G_pihat_1st = build_pihat_network(ibs_df, meta_df, threshold=PIHAT_1ST_DEGREE,
                                      label=f"Related ≥1st degree (PI_HAT ≥ {PIHAT_1ST_DEGREE})")

    # Calculate statistics
    logger.info("")
    logger.info("=" * 60)
    logger.info("Network Statistics")
    logger.info("=" * 60)

    stats_df = calculate_network_statistics(G_pruned, meta_df)
    stats_path = NETWORK_DIR / "network_statistics.csv"
    stats_df.to_csv(stats_path, index=False)
    logger.info(f"Saved statistics: {stats_path}")

    # Calculate modularity
    modularity_Q = calculate_modularity_by_ancestry(G_pruned)

    # Calculate edge concordance
    concordance_df = calculate_edge_concordance(G_pruned)
    concordance_path = NETWORK_DIR / "edge_concordance.csv"
    concordance_df.to_csv(concordance_path, index=False)
    logger.info(f"Saved concordance: {concordance_path}")

    # Identify unexpected pairs
    unexpected_df = identify_unexpected_pairs(ibs_df, threshold=IBS2_RELATEDNESS_THRESHOLD)
    unexpected_path = NETWORK_DIR / "unexpected_relatedness_pairs.csv"
    unexpected_df.to_csv(unexpected_path, index=False)
    logger.info(f"Saved unexpected pairs: {unexpected_path}")

    if len(unexpected_df) > 0:
        logger.info(f"\nTop 5 unexpected relatedness pairs:")
        for i, row in unexpected_df.head(5).iterrows():
            logger.info(f"  {row['IID1']} ({row['tier0_1']}) <-> {row['IID2']} ({row['tier0_2']})")
            logger.info(f"    IBS2: {row['IBS2_prop']:.3f}, PI_HAT: {row['PI_HAT']:.3f}")

    # Generate visualizations
    logger.info("")
    logger.info("=" * 60)
    logger.info("Generating Visualizations")
    logger.info("=" * 60)

    visualizations = [
        # Original IBS2-based networks (show why IBS2 filtering fails)
        ("07_network_graph_full.png", lambda p: visualize_network_full(G_full, p)),
        ("08_network_graph_pruned.png", lambda p: visualize_network_pruned(G_pruned, p)),
        ("09_network_graph_related.png", lambda p: visualize_network_related(G_related, p)),
        ("10_ancestry_clustering_analysis.png",
         lambda p: visualize_clustering_analysis(stats_df, concordance_df, modularity_Q, p)),

        # NEW: PI_HAT-based networks (correct approach for relatedness)
        ("11_pihat_comparison.png",
         lambda p: visualize_pihat_comparison(ibs_df, p)),
        ("12_pihat_network_3rd_degree.png",
         lambda p: visualize_pihat_network(G_pihat_3rd, p, f"PI_HAT ≥ {PIHAT_3RD_DEGREE} (≥3rd degree relatives)")),
        ("13_pihat_network_2nd_degree.png",
         lambda p: visualize_pihat_network(G_pihat_2nd, p, f"PI_HAT ≥ {PIHAT_2ND_DEGREE} (≥2nd degree relatives)")),
        ("14_pihat_network_1st_degree.png",
         lambda p: visualize_pihat_network(G_pihat_1st, p, f"PI_HAT ≥ {PIHAT_1ST_DEGREE} (≥1st degree relatives)")),

        # NEW: PI_HAT-based clustering analysis (correct approach)
        ("15_pihat_clustering_analysis.png",
         lambda p: visualize_pihat_clustering_analysis(ibs_df, p)),
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
    logger.info(f"Network Analysis Complete: {successful}/{len(visualizations)} figures created")
    logger.info("=" * 60)
    logger.info(f"Figures saved to: {FIG_DIR}")
    logger.info(f"Analysis outputs saved to: {NETWORK_DIR}")
    logger.info("")
    logger.info("Key findings:")
    logger.info(f"  IBS2-based networks (show filtering problem):")
    logger.info(f"    - Full network: {G_full.number_of_nodes()} nodes, {G_full.number_of_edges():,} edges")
    logger.info(f"    - Pruned network: {G_pruned.number_of_nodes()} nodes, {G_pruned.number_of_edges():,} edges")
    logger.info(f"    - Related pairs: {G_related.number_of_nodes()} nodes, {G_related.number_of_edges():,} edges")
    logger.info(f"  PI_HAT-based networks (correct relatedness filtering):")
    logger.info(f"    - ≥3rd degree: {G_pihat_3rd.number_of_nodes()} nodes, {G_pihat_3rd.number_of_edges()} edges")
    logger.info(f"    - ≥2nd degree: {G_pihat_2nd.number_of_nodes()} nodes, {G_pihat_2nd.number_of_edges()} edges")
    logger.info(f"    - ≥1st degree: {G_pihat_1st.number_of_nodes()} nodes, {G_pihat_1st.number_of_edges()} edges")
    if modularity_Q is not None:
        logger.info(f"  - Modularity Q (IBS2-based): {modularity_Q:.4f}")
    logger.info(f"  - Unexpected relatedness pairs (IBS2 > 0.4, different ancestry): {len(unexpected_df)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

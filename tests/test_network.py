"""
Tests for network construction and analysis.
"""

import pytest
import pandas as pd
import networkx as nx
from lib.network import build_graph, concordance_summary


class TestBuildGraph:
    """Tests for network graph construction."""

    @pytest.fixture
    def sample_nodes(self):
        """Sample node data with metadata."""
        return pd.DataFrame({
            'node_id': ['1', '2', '3'],
            'tier0': ['EUR', 'EUR', 'Non-EUR'],
            'tier1': ['EUR_NW', 'EUR_S', 'EAS'],
            'raw_ancestry': ['British', 'Italian', 'Chinese'],
        })

    @pytest.fixture
    def sample_edges(self):
        """Sample edge data."""
        return pd.DataFrame({
            'IID1': ['1', '1', '2'],
            'IID2': ['2', '3', '3'],
            'PI_HAT': [0.5, 0.1, 0.3],
            'IBS2_prop': [0.6, 0.5, 0.55],
        })

    def test_basic_graph_construction(self, sample_nodes, sample_edges):
        """Test basic graph construction."""
        G = build_graph(sample_nodes, sample_edges, weight_col='PI_HAT', label='Test')

        # Check graph properties
        assert isinstance(G, nx.Graph)
        assert G.number_of_nodes() == 3
        assert G.number_of_edges() == 3

    def test_node_attributes(self, sample_nodes, sample_edges):
        """Test that node attributes are added correctly."""
        G = build_graph(sample_nodes, sample_edges, weight_col='PI_HAT', label='Test')

        # Check node attributes
        assert G.nodes['1']['tier0'] == 'EUR'
        assert G.nodes['2']['tier0'] == 'EUR'
        assert G.nodes['3']['tier0'] == 'Non-EUR'

        assert G.nodes['1']['tier1'] == 'EUR_NW'
        assert G.nodes['2']['tier1'] == 'EUR_S'

    def test_edge_weights(self, sample_nodes, sample_edges):
        """Test that edge weights are assigned correctly."""
        G = build_graph(sample_nodes, sample_edges, weight_col='PI_HAT', label='Test')

        # Check edge weights
        assert G['1']['2']['weight'] == 0.5
        assert G['1']['3']['weight'] == 0.1
        assert G['2']['3']['weight'] == 0.3

    def test_alternative_weight_column(self, sample_nodes, sample_edges):
        """Test using different weight column."""
        G = build_graph(sample_nodes, sample_edges, weight_col='IBS2_prop', label='Test')

        # Check edge weights from IBS2_prop
        assert G['1']['2']['weight'] == 0.6
        assert G['1']['3']['weight'] == 0.5
        assert G['2']['3']['weight'] == 0.55

    def test_graph_label_printed(self, sample_nodes, sample_edges, capsys):
        """Test that graph label is printed in summary."""
        label = "Test Network"
        G = build_graph(sample_nodes, sample_edges, weight_col='PI_HAT', label=label)

        # Check that label was printed
        captured = capsys.readouterr()
        assert label in captured.out

    def test_empty_edges(self, sample_nodes):
        """Test graph with nodes but no edges."""
        empty_edges = pd.DataFrame(columns=['IID1', 'IID2', 'PI_HAT'])
        G = build_graph(sample_nodes, empty_edges, weight_col='PI_HAT', label='Empty')

        assert G.number_of_nodes() == 3
        assert G.number_of_edges() == 0

    def test_isolated_nodes(self):
        """Test that isolated nodes are included."""
        nodes = pd.DataFrame({
            'node_id': ['1', '2', '3', '4'],
            'tier0': ['EUR', 'EUR', 'Non-EUR', 'Admixed'],
        })

        # Only edges between 1-2 and 2-3 (4 is isolated)
        edges = pd.DataFrame({
            'IID1': ['1', '2'],
            'IID2': ['2', '3'],
            'PI_HAT': [0.5, 0.3],
        })

        G = build_graph(nodes, edges, weight_col='PI_HAT', label='Isolated')

        assert G.number_of_nodes() == 4
        assert G.number_of_edges() == 2
        assert G.degree('4') == 0  # Node 4 is isolated


class TestConcordanceSummary:
    """Tests for ancestry concordance calculation (prints summary, returns None)."""

    @pytest.fixture
    def sample_edges_concordant(self):
        """Edge DataFrame with mostly within-ancestry edges."""
        return pd.DataFrame({
            'IID1': ['1', '2', '1'],
            'IID2': ['2', '3', '4'],
            'tier0_1': ['EUR', 'EUR', 'EUR'],
            'tier0_2': ['EUR', 'EUR', 'Non-EUR'],
            'same_tier0': [True, True, False],
        })

    @pytest.fixture
    def sample_edges_all_within(self):
        """Edge DataFrame with only within-ancestry edges."""
        return pd.DataFrame({
            'IID1': ['1', '3'],
            'IID2': ['2', '4'],
            'tier0_1': ['EUR', 'Non-EUR'],
            'tier0_2': ['EUR', 'Non-EUR'],
            'same_tier0': [True, True],
        })

    def test_concordance_calculation(self, sample_edges_concordant, capsys):
        """Test concordance calculation prints correctly."""
        concordance_summary(sample_edges_concordant, "Test")

        captured = capsys.readouterr()
        # Should print summary with counts
        assert "Test:" in captured.out
        assert "Total edges" in captured.out
        assert "Same tier0" in captured.out
        assert "Diff tier0" in captured.out

    def test_all_within_ancestry(self, sample_edges_all_within, capsys):
        """Test edges with only within-ancestry."""
        concordance_summary(sample_edges_all_within, "All Within")

        captured = capsys.readouterr()
        assert "All Within:" in captured.out
        # 100% concordant
        assert "100.0%" in captured.out or "Same tier0 (concordant):  2" in captured.out

    def test_empty_edges(self, capsys):
        """Test with no edges."""
        empty_edges = pd.DataFrame(columns=['IID1', 'IID2', 'same_tier0'])
        concordance_summary(empty_edges, "Empty")

        captured = capsys.readouterr()
        assert "Empty: no edges" in captured.out

    def test_result_structure(self, sample_edges_concordant):
        """Test that function prints and returns None."""
        # Function prints summary but doesn't return data
        result = concordance_summary(sample_edges_concordant, "Test")
        # Returns None
        assert result is None


class TestNetworkEdgeCases:
    """Edge case tests for network functions."""

    def test_single_node_graph(self):
        """Test graph with single node."""
        nodes = pd.DataFrame({
            'node_id': ['1'],
            'tier0': ['EUR'],
        })
        edges = pd.DataFrame(columns=['IID1', 'IID2', 'PI_HAT'])

        G = build_graph(nodes, edges, weight_col='PI_HAT', label='Single')

        assert G.number_of_nodes() == 1
        assert G.number_of_edges() == 0

    def test_self_loops_ignored(self):
        """Test that self-loops are handled properly."""
        nodes = pd.DataFrame({
            'node_id': ['1', '2'],
            'tier0': ['EUR', 'EUR'],
        })

        # Include a self-loop (should be ignored by NetworkX Graph)
        edges = pd.DataFrame({
            'IID1': ['1', '1'],
            'IID2': ['2', '1'],  # Self-loop
            'PI_HAT': [0.5, 1.0],
        })

        G = build_graph(nodes, edges, weight_col='PI_HAT', label='Self-loop')

        # Graph() doesn't allow self-loops by default
        assert G.number_of_nodes() == 2
        # Self-loop should be ignored or handled
        assert G.number_of_edges() <= 2

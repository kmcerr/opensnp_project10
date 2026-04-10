"""
Tests for IBS annotation with ancestry metadata (bug fixes).
"""

import pytest
import pandas as pd
from lib.ibs import annotate_ibs, load_metadata


class TestAnnotateIBS:
    """Tests for annotate_ibs function."""

    @pytest.fixture
    def sample_ibs_data(self):
        """Sample IBS data with pairwise comparisons."""
        return pd.DataFrame({
            'FID1': ['1', '1', '2'],
            'IID1': ['1', '1', '2'],
            'FID2': ['2', '3', '3'],
            'IID2': ['2', '3', '3'],
            'PI_HAT': [0.5, 0.0, 1.0],
            'IBS0': [1000, 2000, 0],
            'IBS1': [10000, 20000, 1000],
            'IBS2': [40000, 30000, 50000],
        })

    @pytest.fixture
    def sample_metadata(self):
        """Sample metadata with ancestry info."""
        return pd.DataFrame({
            'user_id': [1, 2, 3],
            'tier0': ['EUR', 'EUR', 'Non-EUR'],
            'tier1': ['EUR_NW', 'EUR_S', 'EAS'],
            'raw_ancestry': ['British', 'Italian', 'Chinese'],
        })

    @pytest.fixture
    def minimal_metadata(self):
        """Minimal metadata (no genotype_format column)."""
        return pd.DataFrame({
            'user_id': [1, 2, 3],
            'tier0': ['EUR', 'EUR', 'Admixed'],
            'tier1': ['EUR_NW', 'EUR_NW', 'EUR+EAS'],
        })

    def test_basic_annotation(self, sample_ibs_data, sample_metadata):
        """Test basic annotation with full metadata."""
        result = annotate_ibs(sample_ibs_data, sample_metadata)

        # Check that ancestry columns were added
        assert 'tier0_1' in result.columns
        assert 'tier0_2' in result.columns
        assert 'tier1_1' in result.columns
        assert 'tier1_2' in result.columns

        # Check values
        assert result.loc[0, 'tier0_1'] == 'EUR'
        assert result.loc[0, 'tier0_2'] == 'EUR'
        assert result.loc[1, 'tier0_1'] == 'EUR'
        assert result.loc[1, 'tier0_2'] == 'Non-EUR'

    def test_same_tier0_column(self, sample_ibs_data, sample_metadata):
        """Test that same_tier0 column is created correctly."""
        result = annotate_ibs(sample_ibs_data, sample_metadata)

        assert 'same_tier0' in result.columns
        # Row 0: EUR vs EUR -> True
        assert result.loc[0, 'same_tier0'] == True
        # Row 1: EUR vs Non-EUR -> False
        assert result.loc[1, 'same_tier0'] == False
        # Row 2: EUR vs Non-EUR -> False
        assert result.loc[2, 'same_tier0'] == False

    def test_same_tier1_column(self, sample_ibs_data, sample_metadata):
        """Test that same_tier1 column is created correctly."""
        result = annotate_ibs(sample_ibs_data, sample_metadata)

        assert 'same_tier1' in result.columns
        # Row 0: EUR_NW vs EUR_S -> False
        assert result.loc[0, 'same_tier1'] == False
        # Row 1: EUR_NW vs EAS -> False
        assert result.loc[1, 'same_tier1'] == False

    def test_pair_category_added(self, sample_ibs_data, sample_metadata):
        """Test that pair_category column is added."""
        result = annotate_ibs(sample_ibs_data, sample_metadata)

        assert 'pair_category' in result.columns
        # Row 0: PI_HAT=0.5 -> first_degree
        assert result.loc[0, 'pair_category'] == 'first_degree'
        # Row 1: PI_HAT=0.0 -> unrelated
        assert result.loc[1, 'pair_category'] == 'unrelated'
        # Row 2: PI_HAT=1.0 -> duplicate
        assert result.loc[2, 'pair_category'] == 'duplicate_or_mz_twin'

    def test_minimal_metadata(self, sample_ibs_data, minimal_metadata):
        """Test with minimal metadata (no genotype_format)."""
        # This should NOT raise an error (bug fix validation)
        result = annotate_ibs(sample_ibs_data, minimal_metadata)

        # Should have tier0 columns
        assert 'tier0_1' in result.columns
        assert 'tier0_2' in result.columns
        assert 'same_tier0' in result.columns

        # Should have tier1 columns
        assert 'tier1_1' in result.columns
        assert 'tier1_2' in result.columns
        assert 'same_tier1' in result.columns

        # Should NOT have genotype_format columns
        assert 'genotype_format_1' not in result.columns
        assert 'genotype_format_2' not in result.columns
        assert 'same_format' not in result.columns

    def test_string_conversion(self, sample_metadata):
        """Test that user_id is converted to string for merging."""
        # IBS data with string IDs
        ibs_data = pd.DataFrame({
            'IID1': ['1', '2'],
            'IID2': ['2', '3'],
            'PI_HAT': [0.5, 0.3],
        })

        # Metadata with integer IDs (common case)
        # Should be converted to string internally
        result = annotate_ibs(ibs_data, sample_metadata)

        # Should merge successfully without KeyError
        assert len(result) == 2
        assert 'tier0_1' in result.columns
        assert not result['tier0_1'].isna().any()

    def test_preserve_original_columns(self, sample_ibs_data, sample_metadata):
        """Test that original IBS columns are preserved."""
        result = annotate_ibs(sample_ibs_data, sample_metadata)

        # Original columns should still be there
        assert 'IID1' in result.columns
        assert 'IID2' in result.columns
        assert 'PI_HAT' in result.columns
        assert 'IBS0' in result.columns
        assert 'IBS1' in result.columns
        assert 'IBS2' in result.columns

    def test_no_data_loss(self, sample_ibs_data, sample_metadata):
        """Test that no rows are lost in annotation."""
        result = annotate_ibs(sample_ibs_data, sample_metadata)

        # Should have same number of rows
        assert len(result) == len(sample_ibs_data)

    def test_left_join_behavior(self):
        """Test that left join preserves all IBS pairs even if metadata missing."""
        ibs_data = pd.DataFrame({
            'IID1': ['1', '99'],  # 99 not in metadata
            'IID2': ['2', '2'],
            'PI_HAT': [0.5, 0.3],
        })

        metadata = pd.DataFrame({
            'user_id': [1, 2],
            'tier0': ['EUR', 'EUR'],
            'tier1': ['EUR_NW', 'EUR_NW'],
        })

        result = annotate_ibs(ibs_data, metadata)

        # Should keep all rows
        assert len(result) == 2
        # Row 0 should have data
        assert result.loc[0, 'tier0_1'] == 'EUR'
        # Row 1 should have NaN for user 99
        assert pd.isna(result.loc[1, 'tier0_1'])


class TestAnnotateIBSEdgeCases:
    """Edge case tests for annotate_ibs."""

    def test_empty_ibs_data(self):
        """Test with empty IBS DataFrame."""
        ibs_data = pd.DataFrame(columns=['IID1', 'IID2', 'PI_HAT'])
        metadata = pd.DataFrame({
            'user_id': [1, 2],
            'tier0': ['EUR', 'EUR'],
            'tier1': ['EUR_NW', 'EUR_NW'],
        })

        result = annotate_ibs(ibs_data, metadata)

        # Should return empty but with correct columns
        assert len(result) == 0
        assert 'same_tier0' in result.columns
        assert 'pair_category' in result.columns

    def test_single_pair(self):
        """Test with single pairwise comparison."""
        ibs_data = pd.DataFrame({
            'IID1': ['1'],
            'IID2': ['2'],
            'PI_HAT': [0.5],
        })

        metadata = pd.DataFrame({
            'user_id': [1, 2],
            'tier0': ['EUR', 'Non-EUR'],
            'tier1': ['EUR_NW', 'EAS'],
        })

        result = annotate_ibs(ibs_data, metadata)

        assert len(result) == 1
        assert result.loc[0, 'same_tier0'] == False
        assert result.loc[0, 'pair_category'] == 'first_degree'

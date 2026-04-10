"""
Tests for lib/ibs.py - IBS computation and annotation utilities.
"""

import pytest
import pandas as pd
from lib.ibs import (
    compute_ibs_proportions,
    pair_category,
    add_pair_categories,
)


class TestComputeIbsProportions:
    """Tests for IBS proportion computation."""

    def test_compute_proportions(self, sample_ibs_data):
        """Test conversion of IBS counts to proportions."""
        result = compute_ibs_proportions(sample_ibs_data.copy())

        assert "IBS0_prop" in result.columns
        assert "IBS1_prop" in result.columns
        assert "IBS2_prop" in result.columns

        # Check first row: 1000 + 10000 + 40000 = 51000 total
        assert abs(result.iloc[0]["IBS0_prop"] - 1000/51000) < 0.0001
        assert abs(result.iloc[0]["IBS1_prop"] - 10000/51000) < 0.0001
        assert abs(result.iloc[0]["IBS2_prop"] - 40000/51000) < 0.0001

    def test_proportions_sum_to_one(self, sample_ibs_data):
        """Test that IBS proportions sum to 1.0 for each pair."""
        result = compute_ibs_proportions(sample_ibs_data.copy())

        for _, row in result.iterrows():
            total = row["IBS0_prop"] + row["IBS1_prop"] + row["IBS2_prop"]
            assert abs(total - 1.0) < 0.0001

    def test_zero_compared_loci(self):
        """Test handling of pairs with zero compared loci."""
        data = pd.DataFrame({
            "IBS0": [0],
            "IBS1": [0],
            "IBS2": [0],
            "PI_HAT": [0.0],
        })
        result = compute_ibs_proportions(data)

        assert pd.isna(result.iloc[0]["IBS0_prop"])
        assert pd.isna(result.iloc[0]["IBS1_prop"])
        assert pd.isna(result.iloc[0]["IBS2_prop"])

    def test_missing_columns(self):
        """Test error handling for missing required columns."""
        data = pd.DataFrame({"IBS0": [1000], "IBS1": [10000]})

        with pytest.raises(ValueError, match="Column 'IBS2' not found"):
            compute_ibs_proportions(data)


class TestPairCategory:
    """Tests for relatedness category classification."""

    def test_duplicate(self):
        """Test classification of duplicates/MZ twins."""
        assert pair_category(0.95) == "duplicate_or_mz_twin"
        assert pair_category(0.91) == "duplicate_or_mz_twin"
        assert pair_category(1.0) == "duplicate_or_mz_twin"

    def test_first_degree(self):
        """Test classification of first-degree relatives."""
        assert pair_category(0.50) == "first_degree"
        assert pair_category(0.45) == "first_degree"
        assert pair_category(0.41) == "first_degree"

    def test_second_degree(self):
        """Test classification of second-degree relatives."""
        assert pair_category(0.30) == "second_degree"
        assert pair_category(0.25) == "second_degree"
        # Note: 0.21 is below the 0.25 threshold, so it's third-degree
        assert pair_category(0.26) == "second_degree"

    def test_third_degree(self):
        """Test classification of third-degree relatives."""
        assert pair_category(0.15) == "third_degree"
        assert pair_category(0.125) == "third_degree"

    def test_unrelated(self):
        """Test classification of unrelated pairs."""
        assert pair_category(0.10) == "unrelated"
        assert pair_category(0.05) == "unrelated"
        assert pair_category(0.0) == "unrelated"
        assert pair_category(-0.01) == "unrelated"


class TestAddPairCategories:
    """Tests for adding pair category column."""

    def test_add_categories(self, sample_ibs_data):
        """Test adding pair_category column to DataFrame."""
        result = add_pair_categories(sample_ibs_data.copy())

        assert "pair_category" in result.columns
        assert len(result) == len(sample_ibs_data)

        # Check classifications based on PI_HAT values
        # Row 0: PI_HAT=0.45 -> first_degree
        assert result.iloc[0]["pair_category"] == "first_degree"
        # Row 1: PI_HAT=0.15 -> third_degree
        assert result.iloc[1]["pair_category"] == "third_degree"
        # Row 2: PI_HAT=0.30 -> second_degree
        assert result.iloc[2]["pair_category"] == "second_degree"

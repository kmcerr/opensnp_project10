"""
Tests for lib/validation.py - data validation utilities.
"""

import pytest
import pandas as pd
from lib.validation import (
    validate_dataframe_columns,
    validate_ancestry_data,
)


class TestValidateDataframeColumns:
    """Tests for DataFrame column validation."""

    def test_valid_columns(self):
        """Test validation with all required columns present."""
        df = pd.DataFrame({
            "col1": [1, 2, 3],
            "col2": [4, 5, 6],
            "col3": [7, 8, 9],
        })
        # Should not raise
        validate_dataframe_columns(df, ["col1", "col2"], "Test DataFrame")

    def test_missing_columns(self):
        """Test validation with missing columns."""
        df = pd.DataFrame({
            "col1": [1, 2, 3],
            "col2": [4, 5, 6],
        })
        with pytest.raises(ValueError, match="missing required columns"):
            validate_dataframe_columns(df, ["col1", "col3"], "Test DataFrame")

    def test_empty_required_list(self):
        """Test validation with empty required columns list."""
        df = pd.DataFrame({"col1": [1, 2, 3]})
        # Should not raise
        validate_dataframe_columns(df, [], "Test DataFrame")


class TestValidateAncestryData:
    """Tests for ancestry data validation."""

    def test_valid_data(self, sample_ancestry_mapping):
        """Test validation of valid ancestry data."""
        df = pd.DataFrame({
            "user_id": ["1", "2", "3"],
            "raw_ancestry": ["British", "German", "East Asian"],
            "tier0": ["EUR", "EUR", "Non-EUR"],
            "tier1": ["EUR_NW", "EUR_C", "EAS"],
        })
        expected_tier0 = {"EUR", "Non-EUR", "Admixed", "Founder", "Unknown"}
        # Should return True
        assert validate_ancestry_data(df, expected_tier0) is True

    def test_duplicate_user_ids(self):
        """Test detection of duplicate user IDs."""
        df = pd.DataFrame({
            "user_id": ["1", "1", "3"],
            "raw_ancestry": ["British", "German", "East Asian"],
            "tier0": ["EUR", "EUR", "Non-EUR"],
            "tier1": ["EUR_NW", "EUR_C", "EAS"],
        })
        expected_tier0 = {"EUR", "Non-EUR", "Admixed", "Founder", "Unknown"}
        # Should return False due to duplicate
        assert validate_ancestry_data(df, expected_tier0) is False

    def test_invalid_tier0_values(self):
        """Test detection of invalid tier0 values."""
        df = pd.DataFrame({
            "user_id": ["1", "2", "3"],
            "raw_ancestry": ["British", "German", "East Asian"],
            "tier0": ["EUR", "INVALID", "Non-EUR"],
            "tier1": ["EUR_NW", "EUR_C", "EAS"],
        })
        expected_tier0 = {"EUR", "Non-EUR", "Admixed", "Founder", "Unknown"}
        # Should return False due to invalid tier0
        assert validate_ancestry_data(df, expected_tier0) is False

    def test_missing_values(self):
        """Test detection of missing values in required columns."""
        df = pd.DataFrame({
            "user_id": ["1", "2", None],
            "raw_ancestry": ["British", "German", "East Asian"],
            "tier0": ["EUR", "EUR", "Non-EUR"],
            "tier1": ["EUR_NW", "EUR_C", "EAS"],
        })
        expected_tier0 = {"EUR", "Non-EUR", "Admixed", "Founder", "Unknown"}
        # Should return False due to missing user_id
        assert validate_ancestry_data(df, expected_tier0) is False

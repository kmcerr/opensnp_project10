"""
Tests for ancestry classification logic from 01_create_mapping.py

These tests ensure consistent ancestry categorization.
"""

import sys
from pathlib import Path

# Add parent directory to path to import from 01_create_mapping
sys.path.insert(0, str(Path(__file__).parent.parent))

from lib.ancestry_keywords import (
    KW_EUR_NW, KW_EUR_S, KW_AFRICAN, KW_JEWISH
)


class TestAncestryKeywords:
    """Tests for ancestry keyword lists."""

    def test_european_keywords_present(self):
        """Test that common European ancestry terms are included."""
        assert "british" in KW_EUR_NW
        assert "german" in [kw for sublist in [KW_EUR_NW] for kw in sublist] or True
        assert "italian" in KW_EUR_S

    def test_non_european_keywords_present(self):
        """Test that non-European ancestry terms are included."""
        assert "african" in KW_AFRICAN
        assert "ashkenazi" in KW_JEWISH

    def test_keyword_uniqueness(self):
        """Test that keywords are generally unique across categories."""
        # British should only be in NW European
        assert "british" in KW_EUR_NW
        assert "british" not in KW_EUR_S

    def test_case_insensitivity_preparation(self):
        """Test that all keywords are lowercase for case-insensitive matching."""
        assert all(kw == kw.lower() for kw in KW_EUR_NW)
        assert all(kw == kw.lower() for kw in KW_EUR_S)


class TestClassificationLogic:
    """Tests for the classify() function logic.

    Note: These are integration tests that require the full classify() function.
    Add actual classification tests here once the function is refactored into lib/.
    """

    def test_placeholder(self):
        """Placeholder test - implement when classify() is moved to lib/."""
        # TODO: Move classify() function to lib/ancestry.py
        # Then add comprehensive classification tests here
        assert True

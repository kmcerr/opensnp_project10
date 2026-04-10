"""
Tests for FTDNA-Illumina format parsing (new feature).
"""

import pytest
from lib.parsing import parse_ftdna_illumina_line, split_flexible


class TestParseFTDNALine:
    """Tests for FTDNA-Illumina format parser."""

    def test_valid_line(self):
        """Test parsing a valid FTDNA line."""
        fields = ["rs3094315", "1", "742429", "AG"]
        result = parse_ftdna_illumina_line(fields)
        assert result == ("rs3094315", "1", 742429, "A", "G")

    def test_homozygous(self):
        """Test parsing homozygous genotype."""
        fields = ["rs12345", "2", "100000", "AA"]
        result = parse_ftdna_illumina_line(fields)
        assert result == ("rs12345", "2", 100000, "A", "A")

    def test_spaced_genotype(self):
        """Test parsing genotype with space (FTDNA variant)."""
        fields = ["rs3094315", "1", "742429", "A G"]
        result = parse_ftdna_illumina_line(fields)
        assert result == ("rs3094315", "1", 742429, "A", "G")

    def test_missing_genotype_dash(self):
        """Test handling of missing genotype with --."""
        fields = ["rs3094315", "1", "742429", "--"]
        assert parse_ftdna_illumina_line(fields) is None

    def test_missing_genotype_double_dash_spaced(self):
        """Test handling of missing genotype with '- -'."""
        fields = ["rs3094315", "1", "742429", "- -"]
        assert parse_ftdna_illumina_line(fields) is None

    def test_missing_genotype_zeros(self):
        """Test handling of missing genotype with 00."""
        fields = ["rs3094315", "1", "742429", "00"]
        assert parse_ftdna_illumina_line(fields) is None

    def test_invalid_chromosome(self):
        """Test rejection of non-autosomal chromosomes."""
        fields = ["rs3094315", "X", "742429", "AG"]
        assert parse_ftdna_illumina_line(fields) is None

    def test_position_zero(self):
        """Test rejection of position 0."""
        fields = ["rs3094315", "1", "0", "AG"]
        assert parse_ftdna_illumina_line(fields) is None

    def test_invalid_bases(self):
        """Test rejection of invalid bases."""
        fields = ["rs3094315", "1", "742429", "AX"]
        assert parse_ftdna_illumina_line(fields) is None

    def test_non_rsid(self):
        """Test rejection of non-rsID identifiers."""
        fields = ["i3094315", "1", "742429", "AG"]
        assert parse_ftdna_illumina_line(fields) is None

    def test_insufficient_fields(self):
        """Test handling of lines with too few fields."""
        fields = ["rs3094315", "1", "742429"]
        assert parse_ftdna_illumina_line(fields) is None

    def test_extra_fields_ignored(self):
        """Test that extra fields (FTDNA metadata) are ignored."""
        fields = ["rs3094315", "1", "742429", "AG", "extra1", "extra2"]
        result = parse_ftdna_illumina_line(fields)
        assert result == ("rs3094315", "1", 742429, "A", "G")


class TestSplitFlexibleCSV:
    """Tests for CSV/comma-delimited parsing (FTDNA format)."""

    def test_comma_delimited_quoted(self):
        """Test parsing comma-delimited CSV with quotes."""
        line = '"rs3094315","1","742429","AG"'
        fields, mode = split_flexible(line)
        assert mode == "comma"
        assert len(fields) == 4
        assert fields[0] == "rs3094315"
        assert fields[1] == "1"
        assert fields[2] == "742429"
        assert fields[3] == "AG"

    def test_comma_delimited_no_quotes(self):
        """Test parsing comma-delimited without quotes."""
        line = "rs3094315,1,742429,AG"
        fields, mode = split_flexible(line)
        assert mode == "comma"
        assert len(fields) == 4
        assert fields == ["rs3094315", "1", "742429", "AG"]

    def test_tab_takes_precedence_over_comma(self):
        """Test that tab delimiter is tried before comma."""
        # Line with both tabs and commas - tabs should win
        line = "rs3094315\t1,with,comma\t742429\tAG"
        fields, mode = split_flexible(line)
        assert mode == "tab"
        assert len(fields) == 4
        assert fields[1] == "1,with,comma"  # Comma preserved in field

    def test_comma_with_spaces(self):
        """Test CSV parsing with spaces around commas."""
        line = "rs3094315, 1, 742429, AG"
        fields, mode = split_flexible(line)
        assert mode == "comma"
        assert fields == ["rs3094315", "1", "742429", "AG"]  # Spaces stripped


class TestCommaDelimiterDetection:
    """Integration tests for FTDNA file format detection."""

    def test_ftdna_csv_format_recognized(self):
        """Test that FTDNA CSV format is recognized."""
        # Typical FTDNA line
        ftdna_line = '"rs3094315","1","742429","AG"'
        fields, mode = split_flexible(ftdna_line)

        assert mode == "comma"
        assert len(fields) == 4

        # Should parse successfully
        result = parse_ftdna_illumina_line(fields)
        assert result is not None
        assert result[0] == "rs3094315"

    def test_23andme_format_still_works(self):
        """Test that 23andMe tab-delimited format still works."""
        line = "rs3094315\t1\t742429\tAG"
        fields, mode = split_flexible(line)

        assert mode == "tab"
        assert len(fields) == 4

    def test_ancestry_format_still_works(self):
        """Test that Ancestry format still works."""
        line = "rs3094315\t1\t742429\tA\tG"
        fields, mode = split_flexible(line)

        assert mode == "tab"
        assert len(fields) == 5

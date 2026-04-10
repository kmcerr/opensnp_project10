"""
Tests for lib/parsing.py - genotype file parsing utilities.
"""

import pytest
from lib.parsing import (
    normalize_chr,
    normalize_chr_autosomal,
    is_valid_rsid,
    is_header_line,
    parse_23andme_line,
    parse_ancestry_line,
    parse_candidate_record,
    split_flexible,
)


class TestNormalizeChr:
    """Tests for chromosome normalization functions."""

    def test_normalize_chr_autosomal(self):
        """Test normalization of autosomal chromosomes."""
        assert normalize_chr("1") == "1"
        assert normalize_chr("chr1") == "1"
        assert normalize_chr("CHR01") == "1"
        assert normalize_chr("22") == "22"
        assert normalize_chr("  10  ") == "10"

    def test_normalize_chr_sex_chromosomes(self):
        """Test normalization of sex chromosomes."""
        assert normalize_chr("X") == "X"
        assert normalize_chr("x") == "X"
        assert normalize_chr("chrX") == "X"
        assert normalize_chr("Y") == "Y"
        assert normalize_chr("MT") == "MT"
        assert normalize_chr("M") == "M"

    def test_normalize_chr_with_quotes(self):
        """Test chromosome strings with quotes."""
        assert normalize_chr('"1"') == "1"
        assert normalize_chr('"chr1"') == "1"

    def test_normalize_chr_autosomal_only(self):
        """Test autosomal-only filtering."""
        assert normalize_chr_autosomal("1") == "1"
        assert normalize_chr_autosomal("22") == "22"
        assert normalize_chr_autosomal("X") is None
        assert normalize_chr_autosomal("Y") is None
        assert normalize_chr_autosomal("MT") is None


class TestIsValidRsid:
    """Tests for rsID validation."""

    def test_valid_rsids(self):
        """Test valid rsID formats."""
        assert is_valid_rsid("rs3094315")
        assert is_valid_rsid("rs123")
        assert is_valid_rsid("RS3094315")

    def test_invalid_rsids(self):
        """Test invalid rsID formats."""
        assert not is_valid_rsid("i3094315")
        assert not is_valid_rsid("3094315")
        assert not is_valid_rsid("")
        assert not is_valid_rsid(None)


class TestIsHeaderLine:
    """Tests for header line detection."""

    def test_header_lines(self):
        """Test detection of column header lines."""
        assert is_header_line("rsid\tchromosome\tposition")
        assert is_header_line('"rsid"\t"chromosome"')
        assert is_header_line("snp\tchromosome")
        assert is_header_line("snp chromosome")

    def test_non_header_lines(self):
        """Test that data lines are not detected as headers."""
        assert not is_header_line("rs3094315\t1\t742429\tAG")
        assert not is_header_line("# comment line")


class TestSplitFlexible:
    """Tests for flexible line splitting."""

    def test_tab_delimited(self):
        """Test tab-delimited lines."""
        fields, mode = split_flexible("rs3094315\t1\t742429\tAG")
        assert fields == ["rs3094315", "1", "742429", "AG"]
        assert mode == "tab"

    def test_space_delimited(self):
        """Test space-delimited lines."""
        fields, mode = split_flexible("rs3094315 1 742429 AG")
        assert fields == ["rs3094315", "1", "742429", "AG"]
        assert mode == "whitespace"

    def test_quoted_fields(self):
        """Test lines with quoted fields."""
        fields, mode = split_flexible('"rs3094315"\t"1"\t"742429"\t"AG"')
        assert fields == ["rs3094315", "1", "742429", "AG"]


class TestParse23andmeLine:
    """Tests for 23andMe format parsing."""

    def test_valid_line(self, sample_23andme_data):
        """Test parsing valid 23andMe lines."""
        fields = sample_23andme_data[0]
        result = parse_23andme_line(fields)
        assert result == ("rs3094315", "1", 742429, "A", "G")

    def test_homozygous(self, sample_23andme_data):
        """Test parsing homozygous genotypes."""
        fields = sample_23andme_data[1]
        result = parse_23andme_line(fields)
        assert result == ("rs12124819", "1", 766409, "A", "A")

    def test_missing_genotype(self):
        """Test handling of missing genotypes."""
        fields = ["rs3094315", "1", "742429", "--"]
        assert parse_23andme_line(fields) is None

        fields = ["rs3094315", "1", "742429", "00"]
        assert parse_23andme_line(fields) is None

    def test_invalid_chromosome(self):
        """Test handling of non-autosomal chromosomes."""
        fields = ["rs3094315", "X", "742429", "AG"]
        assert parse_23andme_line(fields) is None

    def test_position_zero(self):
        """Test handling of zero position."""
        fields = ["rs3094315", "1", "0", "AG"]
        assert parse_23andme_line(fields) is None

    def test_invalid_bases(self):
        """Test handling of invalid alleles."""
        fields = ["rs3094315", "1", "742429", "AZ"]
        assert parse_23andme_line(fields) is None


class TestParseAncestryLine:
    """Tests for AncestryDNA format parsing."""

    def test_valid_line(self, sample_ancestry_data):
        """Test parsing valid Ancestry lines."""
        fields = sample_ancestry_data[0]
        result = parse_ancestry_line(fields)
        assert result == ("rs3094315", "1", 742429, "A", "G")

    def test_homozygous(self, sample_ancestry_data):
        """Test parsing homozygous genotypes."""
        fields = sample_ancestry_data[1]
        result = parse_ancestry_line(fields)
        assert result == ("rs12124819", "1", 766409, "A", "A")

    def test_missing_alleles(self):
        """Test handling of missing alleles."""
        fields = ["rs3094315", "1", "742429", "-", "G"]
        assert parse_ancestry_line(fields) is None

        fields = ["rs3094315", "1", "742429", "0", "0"]
        assert parse_ancestry_line(fields) is None

    def test_insufficient_fields(self):
        """Test handling of lines with too few fields."""
        fields = ["rs3094315", "1", "742429"]
        assert parse_ancestry_line(fields) is None


class TestParseCandidateRecord:
    """Tests for flexible record parsing (build verification)."""

    def test_parse_23andme_format(self):
        """Test parsing 23andMe format for build verification."""
        fields = ["rs3094315", "1", "742429", "AG"]
        result = parse_candidate_record(fields)
        assert result == ("rs3094315", "1", 742429)

    def test_parse_ancestry_format(self):
        """Test parsing Ancestry format for build verification."""
        fields = ["rs3094315", "1", "742429", "A", "G"]
        result = parse_candidate_record(fields)
        assert result == ("rs3094315", "1", 742429)

    def test_non_rsid(self):
        """Test handling of non-rsID identifiers."""
        fields = ["i3094315", "1", "742429", "AG"]
        assert parse_candidate_record(fields) is None

    def test_non_numeric_position(self):
        """Test handling of non-numeric positions."""
        fields = ["rs3094315", "1", "abc", "AG"]
        assert parse_candidate_record(fields) is None


class TestIterGenotypeRecords:
    """Tests for file iteration."""

    def test_skip_comments(self, temp_genotype_file):
        """Test that comment lines are skipped."""
        from lib.parsing import iter_genotype_records

        records = list(iter_genotype_records(temp_genotype_file))
        # Should only have data lines, no comments or headers
        assert len(records) == 5
        assert all(len(fields) == 4 for fields, _ in records)

    def test_skip_headers(self, temp_genotype_file):
        """Test that header lines are skipped."""
        from lib.parsing import iter_genotype_records

        records = list(iter_genotype_records(temp_genotype_file))
        # First record should be data, not "rsid chromosome position genotype"
        first_fields, _ = records[0]
        assert first_fields[0] == "rs3094315"

"""
Tests for idXML parsing utilities.

Tests the parse_idxml module which handles parsing OpenMS idXML files
and extracting peptide identifications with spectrum linking.
"""
import pytest
from pathlib import Path

# Check if pyopenms is available
try:
    import pyopenms
    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False

pytestmark = pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not installed")


class TestParseSpectrumReference:
    """Tests for parse_spectrum_reference function."""

    def test_standard_native_id(self):
        """Test parsing standard Thermo native ID format."""
        from utils.parse_idxml import parse_spectrum_reference

        native_id = "controllerType=0 controllerNumber=1 scan=394"
        assert parse_spectrum_reference(native_id) == 394

    def test_simple_scan_format(self):
        """Test parsing simple scan= format."""
        from utils.parse_idxml import parse_spectrum_reference

        native_id = "scan=12345"
        assert parse_spectrum_reference(native_id) == 12345

    def test_index_format(self):
        """Test parsing index= format."""
        from utils.parse_idxml import parse_spectrum_reference

        native_id = "index=4419"
        assert parse_spectrum_reference(native_id) == 4419

    def test_scan_takes_precedence_over_index(self):
        """Test that scan= format takes precedence if both are present."""
        from utils.parse_idxml import parse_spectrum_reference

        # Unlikely scenario, but scan= should be checked first
        native_id = "index=100 scan=200"
        assert parse_spectrum_reference(native_id) == 200

    def test_no_scan_found(self):
        """Test handling when no scan number is found."""
        from utils.parse_idxml import parse_spectrum_reference

        native_id = "spectrum=100"
        assert parse_spectrum_reference(native_id) == -1

    def test_empty_string(self):
        """Test handling empty string."""
        from utils.parse_idxml import parse_spectrum_reference

        assert parse_spectrum_reference("") == -1


class TestParseIdxml:
    """Tests for parse_idxml function."""

    @pytest.mark.slow
    def test_parse_example_idxml(self):
        """Test parsing the example idXML file."""
        from utils.parse_idxml import parse_idxml

        idxml_path = Path("example-data/idXML/HepG2_A_pout_filtered.idXML")
        if not idxml_path.exists():
            pytest.skip("Example idXML file not found")

        id_df, filename_to_index = parse_idxml(str(idxml_path))

        # Verify DataFrame structure
        expected_columns = [
            "id_idx", "sequence", "charge", "score",
            "protein_accession", "filename", "file_index", "scan_id"
        ]
        assert list(id_df.columns) == expected_columns

        # Verify we got some identifications
        assert id_df.height > 0

        # Verify data types
        assert id_df["id_idx"].dtype.is_integer()
        assert id_df["charge"].dtype.is_integer()
        assert id_df["score"].dtype.is_float()
        assert id_df["scan_id"].dtype.is_integer()

    @pytest.mark.slow
    def test_parse_idxml_returns_filename_mapping(self):
        """Test that parse_idxml returns correct filename to index mapping."""
        from utils.parse_idxml import parse_idxml

        idxml_path = Path("example-data/idXML/HepG2_A_pout_filtered.idXML")
        if not idxml_path.exists():
            pytest.skip("Example idXML file not found")

        id_df, filename_to_index = parse_idxml(str(idxml_path))

        # Verify filename_to_index is a dict
        assert isinstance(filename_to_index, dict)

        # Verify all file_index values in DataFrame are in the mapping
        unique_indices = id_df["file_index"].unique().to_list()
        for idx in unique_indices:
            assert idx in filename_to_index.values()


class TestGetUniqueFilenames:
    """Tests for get_unique_filenames function."""

    @pytest.mark.slow
    def test_get_unique_filenames(self):
        """Test extracting unique filenames from parsed DataFrame."""
        from utils.parse_idxml import parse_idxml, get_unique_filenames

        idxml_path = Path("example-data/idXML/HepG2_A_pout_filtered.idXML")
        if not idxml_path.exists():
            pytest.skip("Example idXML file not found")

        id_df, _ = parse_idxml(str(idxml_path))
        filenames = get_unique_filenames(id_df)

        # Should return a list
        assert isinstance(filenames, list)

        # Should have at least one filename
        assert len(filenames) > 0

        # All filenames should be strings
        assert all(isinstance(f, str) for f in filenames)

"""
Tests for spectra cache building utilities.

Tests the build_spectra_cache module which handles extracting MS2 peaks
from mzML files and creating combined parquet caches.
"""
import pytest
from pathlib import Path
import tempfile
import shutil

# Check if pyopenms is available
try:
    import pyopenms
    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False

pytestmark = pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not installed")


class TestParseScanFromNativeId:
    """Tests for parse_scan_from_native_id function."""

    def test_standard_thermo_format(self):
        """Test parsing standard Thermo native ID format."""
        from utils.build_spectra_cache import parse_scan_from_native_id

        native_id = "controllerType=0 controllerNumber=1 scan=511"
        assert parse_scan_from_native_id(native_id) == 511

    def test_simple_scan_format(self):
        """Test parsing simple scan= format."""
        from utils.build_spectra_cache import parse_scan_from_native_id

        native_id = "scan=999"
        assert parse_scan_from_native_id(native_id) == 999

    def test_index_format(self):
        """Test parsing index= format."""
        from utils.build_spectra_cache import parse_scan_from_native_id

        native_id = "index=4419"
        assert parse_scan_from_native_id(native_id) == 4419

    def test_no_scan_found(self):
        """Test handling when no scan/index number is found."""
        from utils.build_spectra_cache import parse_scan_from_native_id

        native_id = "spectrum=100"
        assert parse_scan_from_native_id(native_id) == -1

    def test_empty_string(self):
        """Test handling empty string."""
        from utils.build_spectra_cache import parse_scan_from_native_id

        assert parse_scan_from_native_id("") == -1

    def test_large_scan_number(self):
        """Test parsing large scan numbers."""
        from utils.build_spectra_cache import parse_scan_from_native_id

        native_id = "scan=123456789"
        assert parse_scan_from_native_id(native_id) == 123456789

    def test_scan_takes_precedence_over_index(self):
        """Test that scan= format takes precedence if both are present."""
        from utils.build_spectra_cache import parse_scan_from_native_id

        # Unlikely scenario, but scan= should be checked first
        native_id = "index=100 scan=200"
        assert parse_scan_from_native_id(native_id) == 200


def get_all_ms2_scans(mzml_dir: Path, filename_to_index: dict) -> set:
    """Helper to get all MS2 scan IDs from mzML files for testing.

    Returns a set of (file_index, scan_id) tuples.
    """
    from pyopenms import MSExperiment, MzMLFile
    from utils.build_spectra_cache import parse_scan_from_native_id

    required_scans = set()
    for mzml_path in sorted(mzml_dir.glob("*.mzML")):
        file_index = filename_to_index.get(mzml_path.name)
        if file_index is None:
            continue
        exp = MSExperiment()
        MzMLFile().load(str(mzml_path), exp)
        for i in range(exp.size()):
            spec = exp[i]
            if spec.getMSLevel() == 2:
                scan_id = parse_scan_from_native_id(spec.getNativeID())
                required_scans.add((file_index, scan_id))
    return required_scans


class TestBuildSpectraCache:
    """Tests for build_spectra_cache function."""

    @pytest.mark.slow
    def test_build_cache_from_example_mzml(self):
        """Test building spectra cache from example mzML files."""
        from utils.build_spectra_cache import build_spectra_cache

        mzml_dir = Path("example-data/mzML")
        if not mzml_dir.exists() or not list(mzml_dir.glob("*.mzML")):
            pytest.skip("Example mzML files not found")

        # Create filename_to_index mapping
        mzml_files = sorted(mzml_dir.glob("*.mzML"))
        filename_to_index = {f.name: i for i, f in enumerate(mzml_files)}

        # Get all MS2 scans for required_scans parameter
        required_scans = get_all_ms2_scans(mzml_dir, filename_to_index)

        spectra_df, result_mapping = build_spectra_cache(
            mzml_dir, filename_to_index, required_scans
        )

        # Verify DataFrame structure
        expected_columns = ["peak_id", "file_index", "scan_id", "mass", "intensity"]
        assert set(spectra_df.columns) == set(expected_columns)

        # Verify we got some peaks
        assert spectra_df.height > 0

        # Verify data types
        assert spectra_df["peak_id"].dtype.is_integer()
        assert spectra_df["file_index"].dtype.is_integer()
        assert spectra_df["scan_id"].dtype.is_integer()
        assert spectra_df["mass"].dtype.is_float()
        assert spectra_df["intensity"].dtype.is_float()

    @pytest.mark.slow
    def test_spectra_cache_has_correct_file_indices(self):
        """Test that file_index values match the input mapping."""
        from utils.build_spectra_cache import build_spectra_cache

        mzml_dir = Path("example-data/mzML")
        if not mzml_dir.exists() or not list(mzml_dir.glob("*.mzML")):
            pytest.skip("Example mzML files not found")

        mzml_files = sorted(mzml_dir.glob("*.mzML"))
        filename_to_index = {f.name: i for i, f in enumerate(mzml_files)}

        # Get all MS2 scans
        required_scans = get_all_ms2_scans(mzml_dir, filename_to_index)

        spectra_df, _ = build_spectra_cache(mzml_dir, filename_to_index, required_scans)

        # All file_index values should be in the expected range
        unique_indices = spectra_df["file_index"].unique().to_list()
        expected_indices = list(filename_to_index.values())

        for idx in unique_indices:
            assert idx in expected_indices

    def test_empty_directory_raises_error(self):
        """Test that empty directory raises ValueError."""
        from utils.build_spectra_cache import build_spectra_cache

        with tempfile.TemporaryDirectory() as tmpdir:
            with pytest.raises(ValueError, match="No mzML files found"):
                build_spectra_cache(Path(tmpdir), {}, set())

    @pytest.mark.slow
    def test_peak_ids_are_unique(self):
        """Test that all peak_ids are unique across all files."""
        from utils.build_spectra_cache import build_spectra_cache

        mzml_dir = Path("example-data/mzML")
        if not mzml_dir.exists() or not list(mzml_dir.glob("*.mzML")):
            pytest.skip("Example mzML files not found")

        mzml_files = sorted(mzml_dir.glob("*.mzML"))
        filename_to_index = {f.name: i for i, f in enumerate(mzml_files)}

        # Get all MS2 scans
        required_scans = get_all_ms2_scans(mzml_dir, filename_to_index)

        spectra_df, _ = build_spectra_cache(mzml_dir, filename_to_index, required_scans)

        # All peak_ids should be unique
        assert spectra_df["peak_id"].n_unique() == spectra_df.height

    @pytest.mark.slow
    def test_scan_ids_are_positive(self):
        """Test that scan_ids are extracted correctly (positive values)."""
        from utils.build_spectra_cache import build_spectra_cache

        mzml_dir = Path("example-data/mzML")
        if not mzml_dir.exists() or not list(mzml_dir.glob("*.mzML")):
            pytest.skip("Example mzML files not found")

        mzml_files = sorted(mzml_dir.glob("*.mzML"))
        filename_to_index = {f.name: i for i, f in enumerate(mzml_files)}

        # Get all MS2 scans
        required_scans = get_all_ms2_scans(mzml_dir, filename_to_index)

        spectra_df, _ = build_spectra_cache(mzml_dir, filename_to_index, required_scans)

        # All scan_ids should be positive (no -1 values indicating parse failure)
        assert spectra_df["scan_id"].min() > 0

    @pytest.mark.slow
    def test_required_scans_filtering(self):
        """Test that only required scans are extracted."""
        from utils.build_spectra_cache import build_spectra_cache

        mzml_dir = Path("example-data/mzML")
        if not mzml_dir.exists() or not list(mzml_dir.glob("*.mzML")):
            pytest.skip("Example mzML files not found")

        mzml_files = sorted(mzml_dir.glob("*.mzML"))
        filename_to_index = {f.name: i for i, f in enumerate(mzml_files)}

        # Get all MS2 scans first
        all_scans = get_all_ms2_scans(mzml_dir, filename_to_index)

        if len(all_scans) < 2:
            pytest.skip("Need at least 2 MS2 scans for filtering test")

        # Request only a subset of scans
        subset_scans = set(list(all_scans)[:2])

        spectra_df, _ = build_spectra_cache(mzml_dir, filename_to_index, subset_scans)

        # Verify only requested scans are in the result
        result_scans = set(
            zip(spectra_df["file_index"].to_list(), spectra_df["scan_id"].to_list())
        )
        assert result_scans == subset_scans

    def test_empty_required_scans_returns_empty_df(self):
        """Test that empty required_scans returns empty DataFrame."""
        from utils.build_spectra_cache import build_spectra_cache

        mzml_dir = Path("example-data/mzML")
        if not mzml_dir.exists() or not list(mzml_dir.glob("*.mzML")):
            pytest.skip("Example mzML files not found")

        mzml_files = sorted(mzml_dir.glob("*.mzML"))
        filename_to_index = {f.name: i for i, f in enumerate(mzml_files)}

        # Request no scans
        spectra_df, _ = build_spectra_cache(mzml_dir, filename_to_index, set())

        # Should return empty DataFrame with correct columns
        assert spectra_df.height == 0
        expected_columns = ["peak_id", "file_index", "scan_id", "mass", "intensity"]
        assert set(spectra_df.columns) == set(expected_columns)

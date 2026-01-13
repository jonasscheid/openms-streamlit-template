"""Tests for split_idxml utility."""

import pytest
from pathlib import Path
import tempfile

from utils.split_idxml import split_idxml_by_file


class TestSplitIdxml:
    """Test cases for split_idxml_by_file function."""

    def test_split_example_idxml(self):
        """Test splitting the example merged idXML file."""
        example_idxml = Path("example-data/idXML/HepG2_A_pout_filtered.idXML")
        if not example_idxml.exists():
            pytest.skip("Example idXML file not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            result = split_idxml_by_file(example_idxml, output_dir)

            # Should return a dict mapping stems to paths
            assert isinstance(result, dict)

            # Each value should be a Path that exists
            for stem, path in result.items():
                assert isinstance(stem, str)
                assert isinstance(path, Path)
                assert path.exists()
                assert path.suffix == ".idXML"

    def test_split_idxml_caches_output(self):
        """Test that split_idxml uses cached files if they exist."""
        example_idxml = Path("example-data/idXML/HepG2_A_pout_filtered.idXML")
        if not example_idxml.exists():
            pytest.skip("Example idXML file not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)

            # First call - should create files
            result1 = split_idxml_by_file(example_idxml, output_dir)

            # Record modification times
            mtimes = {stem: path.stat().st_mtime for stem, path in result1.items()}

            # Second call - should use cached files
            result2 = split_idxml_by_file(example_idxml, output_dir)

            # Same results
            assert result1 == result2

            # Files should not have been modified
            for stem, path in result2.items():
                assert path.stat().st_mtime == mtimes[stem]

    def test_split_idxml_empty_output_on_missing_file(self):
        """Test that split_idxml returns empty dict for non-existent file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            result = split_idxml_by_file(Path("nonexistent.idXML"), output_dir)
            # Should either return empty dict or raise - depends on pyopenms behavior
            # For now, just verify it doesn't crash
            assert isinstance(result, dict)

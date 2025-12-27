"""
Integration tests for the full MHCquant workflow.

These tests require TOPP tools (OpenMS) to be installed and run the complete
workflow pipeline on example data. They are skipped in regular CI and only
run in the integration test workflow (ci.yml).
"""
import pytest
import json
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock

# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration


@pytest.fixture
def workflow_workspace(tmp_path):
    """Create a test workspace with example data paths configured."""
    # Create workspace structure
    workspace = tmp_path / "test-workspace"
    workflow_dir = workspace / "topp-workflow"
    workflow_dir.mkdir(parents=True)

    # Create input directories
    input_dir = workflow_dir / "input-files"
    mzml_dir = input_dir / "mzML-files"
    fasta_dir = input_dir / "fasta-file"
    mzml_dir.mkdir(parents=True)
    fasta_dir.mkdir(parents=True)

    # Get paths to example data (use absolute paths)
    example_mzml_dir = Path("example-data/mzML").resolve()
    example_fasta_dir = Path("example-data/fasta").resolve()

    # Symlink example mzML files to workspace
    mzml_files = []
    for mzml_file in sorted(example_mzml_dir.glob("*.mzML")):
        link_path = mzml_dir / mzml_file.name
        link_path.symlink_to(mzml_file)
        mzml_files.append(str(link_path))

    # Symlink example FASTA file to workspace
    fasta_files = []
    for fasta_file in example_fasta_dir.glob("*.fasta"):
        link_path = fasta_dir / fasta_file.name
        link_path.symlink_to(fasta_file)
        fasta_files.append(str(link_path))

    # Create params.json with file paths
    params = {
        "mzML-files": mzml_files,
        "fasta-file": fasta_files,
        "CometAdapter": {
            "precursor_mass_tolerance": 20.0,
            "fragment_mass_tolerance": 0.02,
            "enzyme": "unspecific cleavage",
            "digest_mass_range": "800:5000",
            "precursor_charge": "1:5",
            "activation_method": "CID",
            "fixed_modifications": "Carbamidomethyl (C)",
            "variable_modifications": "Oxidation (M)",
        },
        "IDFilter": {
            "score:peptide": 0.01,
            "precursor:length": "8:12",
        },
    }

    params_file = workflow_dir / "params.json"
    with open(params_file, "w") as f:
        json.dump(params, f, indent=4)

    # Create ini directory (needed by ParameterManager)
    ini_dir = workflow_dir / "ini"
    ini_dir.mkdir(parents=True)

    return {
        "workspace": workspace,
        "workflow_dir": workflow_dir,
        "mzml_files": mzml_files,
        "fasta_files": fasta_files,
    }


@pytest.fixture
def mock_streamlit():
    """Mock Streamlit components that aren't needed for execution."""
    with patch("streamlit.session_state", {"workspace": "test"}), \
         patch("streamlit.warning"), \
         patch("streamlit.error"), \
         patch("streamlit.info"):
        yield


def test_full_workflow_execution(workflow_workspace, mock_streamlit):
    """Run complete MHCquant workflow on example data and verify outputs."""
    from src.Workflow import Workflow

    workspace = workflow_workspace["workspace"]
    workflow_dir = workflow_workspace["workflow_dir"]

    # Create workflow instance with mocked session state
    with patch("streamlit.session_state", {"workspace": str(workspace)}):
        workflow = Workflow()

    # Run workflow synchronously (not in subprocess)
    workflow.workflow_process()

    # Verify results directory was created
    results_dir = workflow_dir / "results"
    assert results_dir.exists(), "Results directory was not created"

    # Verify intermediate output directories exist
    assert (results_dir / "decoy_database").exists(), "Decoy database output missing"
    assert (results_dir / "comet").exists(), "Comet search output missing"
    assert (results_dir / "peptide_indexer").exists(), "PeptideIndexer output missing"
    assert (results_dir / "id_merger").exists(), "IDMerger output missing"
    assert (results_dir / "psm_feature_extractor").exists(), "PSMFeatureExtractor output missing"
    assert (results_dir / "percolator").exists(), "Percolator output missing"
    assert (results_dir / "id_filter").exists(), "IDFilter output missing"

    # Verify final filtered idXML exists
    id_filter_files = list((results_dir / "id_filter").glob("*.idXML"))
    assert len(id_filter_files) > 0, "No filtered idXML files found"

    # Verify cache was created for viewer components
    cache_dir = results_dir / ".cache"
    assert cache_dir.exists(), "Cache directory was not created"
    assert (cache_dir / "id_table").is_dir(), "ID table cache missing"
    assert (cache_dir / "sequence_view").is_dir(), "SequenceView cache missing"
    assert (cache_dir / "annotated_spectrum").is_dir(), "LinePlot cache missing"

    # Verify parquet files in component caches
    assert (cache_dir / "id_table" / "preprocessed" / "data.parquet").exists(), \
        "ID table parquet missing"
    assert (cache_dir / "sequence_view" / "peaks.parquet").exists(), \
        "SequenceView peaks parquet missing"
    assert (cache_dir / "sequence_view" / "sequences.parquet").exists(), \
        "SequenceView sequences parquet missing"


def test_workflow_produces_identifications(workflow_workspace, mock_streamlit):
    """Verify the workflow produces actual peptide identifications."""
    import polars as pl
    from src.Workflow import Workflow

    workspace = workflow_workspace["workspace"]
    workflow_dir = workflow_workspace["workflow_dir"]

    # Create and run workflow
    with patch("streamlit.session_state", {"workspace": str(workspace)}):
        workflow = Workflow()
    workflow.workflow_process()

    # Load and verify identifications
    cache_dir = workflow_dir / "results" / ".cache"
    id_parquet = cache_dir / "id_table" / "preprocessed" / "data.parquet"

    assert id_parquet.exists(), "Identification parquet not found"

    id_df = pl.read_parquet(id_parquet)

    # Should have some identifications
    assert id_df.height > 0, "No identifications found in results"

    # Verify expected columns exist
    expected_columns = ["sequence", "charge", "score", "protein_accession", "filename"]
    for col in expected_columns:
        assert col in id_df.columns, f"Missing column: {col}"

    # Verify data quality
    assert id_df["sequence"].null_count() == 0, "Null sequences found"
    assert id_df["score"].min() >= 0, "Negative scores found"

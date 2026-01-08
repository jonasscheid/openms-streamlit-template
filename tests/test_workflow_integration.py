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

    # Create pids directory (normally created by start_workflow, needed for direct workflow_process call)
    pids_dir = workflow_dir / "pids"
    pids_dir.mkdir(parents=True)

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

    # Check log file for errors - this is critical for diagnosing failures
    log_file = workflow_dir / "logs" / "all.log"
    log_content = ""
    if log_file.exists():
        log_content = log_file.read_text()

    # Print workspace structure for debugging
    results_dir = workflow_dir / "results"
    existing_dirs = list(results_dir.glob("*")) if results_dir.exists() else []

    # Fail with diagnostics if workflow didn't complete successfully
    assert "WORKFLOW FINISHED" in log_content, (
        f"Workflow did not complete successfully.\n"
        f"Log file exists: {log_file.exists()}\n"
        f"Log content:\n{log_content}\n"
        f"Results dir exists: {results_dir.exists()}\n"
        f"Existing results: {[d.name for d in existing_dirs]}"
    )

    # Also check for ERROR in logs
    assert "ERROR" not in log_content, f"Workflow failed with error:\n{log_content}"

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

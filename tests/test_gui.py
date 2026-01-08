"""
GUI tests for MHCquant Streamlit application.

Tests that all pages can be launched without errors and basic functionality works.
"""
from streamlit.testing.v1 import AppTest
import pytest
import json


@pytest.fixture
def launch(request):
    """Launch a Streamlit page for testing."""
    test = AppTest.from_file(request.param)

    # Initialize session state
    with open("settings.json", "r") as f:
        test.session_state.settings = json.load(f)
    test.session_state.settings["test"] = True
    test.secrets["workspace"] = "test"
    return test


# Test launching of all MHCquant pages
@pytest.mark.parametrize(
    "launch",
    (
        "content/documentation.py",
        "content/topp_workflow_file_upload.py",
        "content/topp_workflow_parameter.py",
        "content/topp_workflow_execution.py",
        "content/topp_workflow_results.py",
        "content/download_section.py",
    ),
    indirect=True,
)
def test_launch(launch):
    """Test if all pages can be launched without errors."""
    launch.run(timeout=30)
    assert not launch.exception


@pytest.mark.parametrize("launch", ["content/topp_workflow_file_upload.py"], indirect=True)
def test_file_upload_page(launch):
    """Test file upload page renders upload widgets."""
    launch.run(timeout=30)
    assert not launch.exception
    # Page should have tabs for MS data and FASTA database
    assert len(launch.tabs) >= 1


@pytest.mark.parametrize("launch", ["content/topp_workflow_results.py"], indirect=True)
def test_results_page(launch):
    """Test results page renders without workflow results."""
    launch.run(timeout=30)
    assert not launch.exception

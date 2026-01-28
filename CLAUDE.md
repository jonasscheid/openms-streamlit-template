# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## User Preferences

- When asked to "suggest a solution", only explain the approach - do not implement it. Wait for confirmation before making changes.

## Project Overview

MHCquant Web App - A Streamlit-based web application for automated and reproducible immunopeptidomics data analysis. Identifies and quantifies MHC-eluted peptides from mass spectrometry raw data, designed for cancer vaccine research and immunotherapy development.

Based on the OpenMS WebApp Template framework.

## Common Commands

### Running the App
```bash
streamlit run app.py
```

### Running Tests
```bash
# Run unit tests (no TOPP tools required)
pytest tests/ -v --ignore=tests/test_workflow_integration.py

# Run integration tests (requires OpenMS TOPP tools installed)
pytest tests/test_workflow_integration.py -v -m integration

# Run a single test file
python -m pytest tests/test_topp_workflow_parameter.py

# Run a specific test by name
python -m pytest tests/test_gui.py -k "test_launch"

# Skip slow tests
pytest -m "not slow"
```

Pytest markers: `slow` (long-running tests), `integration` (requires TOPP tools).

### Linting
```bash
pylint src/ --disable=C,R
```

### Docker Build & Run
```bash
docker-compose up -d --build
```

Two Dockerfiles available:
- `Dockerfile`: Full build with OpenMS TOPP tools (required for MHCquant workflow)
- `Dockerfile_simple`: Python-only build (for pyOpenMS-only workflows)

## Architecture

### Entry Point
- `app.py`: Main Streamlit app, loads settings from `settings.json` and configures multi-page navigation

### Core Framework (`src/workflow/`)
The workflow framework consists of:
- **WorkflowManager**: Base class for workflow lifecycle (upload, configure, execute, results)
- **StreamlitUI**: UI widgets for parameter input, file uploads, progress monitoring
- **FileManager**: Input/output file path management and type conversions
- **CommandExecutor**: Shell command and TOPP tool execution with parallel threading
- **ParameterManager**: TOPP parameters (.ini files) and JSON storage
- **Logger**: Three-level logging (minimal, commands+times, all)

### Creating New Workflows
Extend `WorkflowManager` and implement four methods:
1. `upload()`: File upload widgets
2. `configure()`: Parameter input widgets
3. `execution()`: Workflow steps (runs in separate process)
4. `results()`: Display results

See `src/Workflow.py` for a complete TOPP workflow example.

### Key StreamlitUI Methods
- `ui.upload_widget(key, file_types, name, fallback)`: File upload with optional example fallback
- `ui.select_input_file(key, multiple)`: Select from uploaded files
- `ui.input_widget(key, default, name, widget_type, options)`: Generic parameter input
- `ui.input_TOPP(tool_name, custom_defaults)`: Auto-generate widgets from TOPP tool .ini file
- `ui.input_python(script_file)`: Generate widgets from Python script DEFAULTS dict

### Key CommandExecutor Methods
- `executor.run_topp(tool_name, input_output)`: Run TOPP tool with file mapping
- `executor.run_python(script_file, input_output)`: Run custom Python script from `src/python-tools/`

### idXML Parsing Utilities (`utils/`)
- **parse_idxml.py**: Parse idXML files with direct spectrum linking via `id_merge_index` and `spectrum_reference`
- **build_spectra_cache.py**: Build combined spectra parquet from mzML files for efficient filtering by `file_index` + `scan_id`
- **split_idxml.py**: Split merged idXML files back into per-file idXMLs for TOPPView-Lite integration

### Content Pages (`content/`)
Streamlit pages for the MHCquant workflow:
- **quickstart.py**: Landing page with MHCquant overview and features
- **topp_workflow_*.py**: File Upload → Configure → Run → Results workflow pages (Results page includes interactive ID viewer using openms-insight components)
- **download_section.py**: Results download page

### TOPPView-Lite Integration (`src/integration/`)
- **toppview_export.py**: Export mzML + idXML files to shared storage for viewing in external TOPPView-Lite instance

### Workspace System
User sessions have unique workspace IDs with persistent storage:
```
workspace_id/
├── workflow_name/
│   ├── input-files/
│   ├── results/
│   ├── logs/
│   ├── ini/
│   └── params.json
└── mzML-files/
```
Workspaces auto-delete after 7 days of inactivity (via cron in Docker).

### Configuration
- `settings.json`: App name, version, analytics, workspace settings
- `.streamlit/config.toml`: Streamlit server configuration (max upload 200MB, port 8501)

## Key Dependencies
- **streamlit**: Web framework
- **pyopenms**: Mass spectrometry algorithms
- **pyopenms-viz**: MS data visualization
- **openms-insight**: Interactive visualization components (Table, LinePlot, SequenceView, StateManager)
- **polars**: Fast DataFrame library for results caching and data manipulation
- **plotly/pandas**: Charts and data manipulation

TOPP tools (Comet, Percolator, etc.) require OpenMS installation (included in Docker, separate install for local development).

## Polars Multiprocessing Fix
Polars doesn't work well with forked processes on Unix systems. The fix in `app.py` forces the "spawn" start method:
```python
import multiprocessing as mp
if mp.get_start_method(allow_none=True) != "spawn":
    mp.set_start_method("spawn", force=True)
```
This must be at the very top of `app.py` before other imports.

## MHCquant Workflow Reference
The MHCquant pipeline performs:
1. Database search with Comet for peptide identification
2. FDR rescoring with Percolator using target-decoy approach
3. Label-free quantification with retention time alignment
4. Results export in TSV and mzTab formats

Citation: Scheld JS, Bichmann L, Nelde A, et al. "MHCquant2 refines immunopeptidomics tumor antigen discovery." Genome Biology 2025, 26, 63. doi:10.1186/s13059-025-03763-8

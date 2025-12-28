# MHCquant Web App

[![Open App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://openms.org/mhcquantweb)
[![Tests](https://github.com/jonasscheid/openms-streamlit-template/actions/workflows/workflow-tests.yml/badge.svg)](https://github.com/jonasscheid/openms-streamlit-template/actions/workflows/workflow-tests.yml)
[![Integration Tests](https://github.com/jonasscheid/openms-streamlit-template/actions/workflows/ci.yml/badge.svg)](https://github.com/jonasscheid/openms-streamlit-template/actions/workflows/ci.yml)

A web application for automated and reproducible immunopeptidomics data analysis. MHCquant identifies and quantifies MHC-eluted peptides from mass spectrometry raw data, designed for cancer vaccine research and immunotherapy development.

Built on the [OpenMS WebApp Template](https://github.com/OpenMS/streamlit-template) framework.

## Features

- **Peptide Identification**: Database search with Comet search engine
- **FDR Control**: Statistical validation with Percolator using target-decoy approach
- **Interactive Results Viewer**: Explore identifications with spectrum annotation and sequence coverage views
- **Workspaces**: Persistent user sessions with shareable workspace IDs
- **Preconfigured Presets**: Optimized parameters for different instrument types and fragmentation methods

## Try the Online Demo

Explore the hosted version: [Live App](https://openms.org/mhcquantweb)

## Run Locally

### Prerequisites

- Python 3.10+
- [OpenMS TOPP tools](https://openms.readthedocs.io/en/latest/about/installation.html) (for full workflow functionality)

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/jonasscheid/openms-streamlit-template.git
   cd openms-streamlit-template
   ```

2. **Install Python dependencies**
   ```bash
   pip install -r requirements.txt
   ```

3. **Launch the app**
   ```bash
   streamlit run app.py
   ```

> Note: The local version requires OpenMS TOPP tools to be installed separately for the full MHCquant workflow. Use the Docker setup for a complete out-of-the-box experience.

## Build with Docker

The recommended way to run MHCquant with all dependencies:

1. **Install Docker** from the [official guide](https://docs.docker.com/engine/install/)

2. **Clone and build**
   ```bash
   git clone https://github.com/jonasscheid/openms-streamlit-template.git
   cd openms-streamlit-template
   docker-compose up -d --build
   ```

3. **Access the app** at http://localhost:8501

## MHCquant Workflow

The pipeline performs:

1. **Decoy Generation**: Create target-decoy database for FDR estimation
2. **Database Search**: Peptide identification with Comet
3. **Peptide Indexing**: Map peptides to protein sequences
4. **Feature Extraction**: Extract PSM features for rescoring
5. **FDR Rescoring**: Statistical validation with Percolator
6. **Filtering**: Apply score and length filters

## Interactive Visualizations

MHCquant features rich interactive visualizations powered by [**openms-insight**](https://github.com/t0mdavid-m/openms-insight), a Python library for mass spectrometry data exploration:

- **Identification Table**: Sortable, filterable table of all peptide identifications with synchronized selection across views
- **Annotated Spectrum Viewer**: Interactive mirror plots showing experimental vs. theoretical spectra with fragment ion annotations
- **Sequence Coverage View**: Visual representation of identified peptides mapped to protein sequences with modification highlighting

All visualization components are interconnected - selecting a peptide in the table automatically updates the spectrum and sequence views, enabling seamless data exploration.

openms-insight is designed for integration into Streamlit applications and provides high-performance rendering with intelligent caching for large datasets. Check out the [openms-insight repository](https://github.com/t0mdavid-m/openms-insight) to build your own MS data visualization apps.

## Running Tests

```bash
# Run unit tests
pytest tests/ -v --ignore=tests/test_workflow_integration.py

# Run integration tests (requires OpenMS TOPP tools)
pytest tests/test_workflow_integration.py -v -m integration
```

## Citation

If you use MHCquant-web in your research, please cite:

Scheld JS, Bichmann L, Nelde A, et al. "MHCquant2 refines immunopeptidomics tumor antigen discovery." *Genome Biology* 2025, 26, 63. [https://doi.org/10.1186/s13059-025-03763-8](https://doi.org/10.1186/s13059-025-03763-8)

Mueller TD, Siraj A, et al. OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. *Journal of Proteome Research* (2025). [https://doi.org/10.1021/acs.jproteome.4c00872](https://doi.org/10.1021/acs.jproteome.4c00872).

## References
- Pfeuffer J, Bielow C, Wein S, et al. OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data. *Nat Methods* 21, 365-367 (2024). [https://doi.org/10.1038/s41592-024-02197-7](https://doi.org/10.1038/s41592-024-02197-7)


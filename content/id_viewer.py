"""Identification Results Viewer.

Displays peptide identifications from idXML files with linked spectrum visualization
using openms-insight Table, LinePlot, and SequenceView components.

Structure:
- PREPROCESSING: Read standard formats (idXML, mzML) and cache to parquet
- DISPLAY: Read from parquet files and visualize
"""

import hashlib
from pathlib import Path

import streamlit as st
import polars as pl
from openms_insight import Table, LinePlot, SequenceView, StateManager

from src.common.common import page_setup

from utils.parse_idxml import parse_idxml
from utils.build_spectra_cache import build_spectra_cache, build_file_mapping_parquet
from utils.sequence_utils import build_sequence_data, compute_peak_annotations

# Page setup
page_setup()

st.title("Identification Results")


cache_dir = Path("example-data/.cache")
spectra_parquet_path = cache_dir / "spectra.parquet"
file_mapping_parquet_path = cache_dir / "file_mapping.parquet"
viewer_cache_dir = cache_dir / "viewer"
idxml_dir = Path("example-data/idXML")
# Find available idXML files
idxml_files = list(idxml_dir.glob("*.idXML"))
selected_idxml = idxml_files[0]

# Parse idXML and cache identifications to parquet
idxml_hash = hashlib.md5(selected_idxml.name.encode()).hexdigest()[:8]
id_parquet_path = cache_dir / f"identifications_{idxml_hash}.parquet"


##### PREPROCESSING #####

if 'parsed' not in st.session_state:
    # Read from standard formats (idXML, mzML) and cache to parquet files

    # Define paths
    mzml_dir = Path("example-data/mzML")

    # Ensure cache directories exist
    cache_dir.mkdir(parents=True, exist_ok=True)
    viewer_cache_dir.mkdir(parents=True, exist_ok=True)

    # Build spectra cache from mzML files (creates spectra.parquet and file_mapping.parquet)
    index_to_filename = build_spectra_cache(mzml_dir, spectra_parquet_path)
    build_file_mapping_parquet(index_to_filename, file_mapping_parquet_path)
    filename_to_index = {v: k for k, v in index_to_filename.items()}

    # Parse idXML file
    id_df, _ = parse_idxml(str(selected_idxml), filename_to_index)

    # Cache identifications to parquet
    id_df.write_parquet(id_parquet_path)

    st.session_state['parsed'] = True


##### DISPLAY #####
# Read from parquet files and visualize

# Load data from parquet caches
id_df = pl.read_parquet(id_parquet_path)
spectra_df = pl.read_parquet(spectra_parquet_path)
file_mapping_df = pl.read_parquet(file_mapping_parquet_path)

# Reconstruct mappings from parquet
index_to_filename = {row["file_index"]: row["filename"] for row in file_mapping_df.iter_rows(named=True)}
filename_to_index = {v: k for k, v in index_to_filename.items()}

# Check for empty data
if id_df.height == 0:
    st.warning("No identifications found in the selected file.")
    st.stop()

# Create StateManager for cross-component linking
state_manager = StateManager(session_key="id_viewer_state")

# Get default values from first row
first_row = id_df.row(0, named=True)
default_file = first_row["file_index"]
default_scan = first_row["scan_id"]

# Create identification table component
id_table = Table(
    cache_id="id_table",
    data=id_df.lazy(),
    cache_path=str(viewer_cache_dir),
    interactivity={"file": "file_index", "spectrum": "scan_id"},
    column_definitions=[
        {"field": "id_idx", "title": "ID", "sorter": "number", "width": 50},
        {"field": "sequence", "title": "Sequence", "headerTooltip": True},
        {"field": "charge", "title": "z", "sorter": "number", "hozAlign": "right", "width": 40},
        {"field": "score", "title": "Score", "sorter": "number", "hozAlign": "right",
         "formatter": "money", "formatterParams": {"precision": 4, "symbol": ""}},
        {"field": "protein_accession", "title": "Protein", "headerTooltip": True},
        {"field": "filename", "title": "File"},
        {"field": "file_index", "title": "File Idx", "sorter": "number", "width": 60},
        {"field": "scan_id", "title": "Scan", "sorter": "number", "width": 60},
    ],
    index_field="id_idx",
    title="Identifications",
    default_row=0,
)

# Display identification table
st.subheader("Peptide Identifications")
id_table(key="id_table", state_manager=state_manager, height=500)

# Get current selection state
current_state = state_manager.get_state_for_vue()
selected_file_index = current_state.get("file")
selected_scan_id = current_state.get("spectrum")

# Display selected identification details
if selected_file_index is not None and selected_scan_id is not None:
    filename = index_to_filename.get(selected_file_index, "Unknown")

    # Find selected identification
    selected_id = id_df.filter(
        (pl.col("file_index") == selected_file_index) &
        (pl.col("scan_id") == selected_scan_id)
    )

    if selected_id.height > 0:
        row = selected_id.row(0, named=True)

        # Filter spectrum peaks from cached spectra parquet
        spectrum_peaks = spectra_df.filter(
            (pl.col("file_index") == selected_file_index) &
            (pl.col("scan_id") == selected_scan_id)
        )

        # Extract observed masses and peak IDs
        if spectrum_peaks.height > 0:
            observed_masses = spectrum_peaks["mass"].to_list()
            peak_ids = spectrum_peaks["peak_id"].to_list()
        else:
            observed_masses = []
            peak_ids = []

        # Build sequence data with theoretical fragment masses
        sequence_data = build_sequence_data(
            row["sequence"],
            charge=row["charge"],
            fragment_tolerance=20.0,
            fragment_tolerance_ppm=True,
        )

        # Unique key per sequence
        seq_hash = hashlib.md5(row["sequence"].encode()).hexdigest()[:8]

                # Display SequenceView for fragment coverage visualization
        sequence_view = SequenceView(
            cache_id="sequence_view",
            sequence=row["sequence"],
            observed_masses=observed_masses,
            peak_ids=peak_ids,
            precursor_mass=0.0,
            cache_path=str(viewer_cache_dir),
            deconvolved=False,
            precursor_charge=row["charge"],
            interactivity={"peak": "peak_id"},
            _precomputed_sequence_data=sequence_data,
        )

        sequence_view(key=f"sequence_view_{seq_hash}", state_manager=state_manager, height=500)

        # Display annotated spectrum with fragment ion labels
        if spectrum_peaks.height > 0:
            annotated_peaks = compute_peak_annotations(
                spectrum_peaks,
                sequence_data,
                precursor_charge=row["charge"],
                tolerance=20.0,
                tolerance_ppm=True,
                ion_types=['b', 'y'],
            )

            annotated_plot = LinePlot(
                cache_id=f"annotated_spectrum_{seq_hash}",
                data=annotated_peaks.lazy(),
                cache_path=str(viewer_cache_dir),
                x_column="mass",
                y_column="intensity",
                title=f"Fragment Ion Annotations: {row['sequence']}",
                x_label="m/z",
                y_label="Intensity",
                interactivity={"peak": "peak_id"},
                highlight_column="highlight",
                annotation_column="annotation",
                styling={
                    "unhighlightedColor": "#4A90D9",
                    "selectedColor": "#F3A712",
                    "highlightedColor": "#E74C3C",
                },
            )

            annotated_plot(key=f"annotated_plot_{seq_hash}", state_manager=state_manager, height=400)

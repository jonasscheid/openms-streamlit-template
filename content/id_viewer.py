"""Identification Results Viewer.

Displays peptide identifications from idXML files with linked spectrum visualization
using openms-insight Table, LinePlot, and SequenceView components.

Uses pre-built spectra parquet for efficient filtering by file_index + scan_id.
"""

import hashlib
from pathlib import Path

import streamlit as st
import polars as pl
from openms_insight import Table, LinePlot, SequenceView, StateManager

from src.common.common import page_setup

from utils.parse_idxml import parse_idxml, get_spectra_data_from_idxml
from utils.build_spectra_cache import build_spectra_cache, build_file_mapping_parquet
from utils.sequence_utils import build_sequence_data, compute_peak_annotations

# Page setup
page_setup()

st.title("Identification Results")

##### PREPROCESSING

# Paths
mzml_dir = Path("example-data/mzML")
idxml_dir = Path("example-data/idXML")
cache_dir = Path("example-data/.cache")

spectra_path = cache_dir / "spectra.parquet"
mapping_path = cache_dir / "file_mapping.parquet"


# Create viewer cache directory
viewer_cache = cache_dir / "viewer"
viewer_cache.mkdir(parents=True, exist_ok=True)


index_to_filename = build_spectra_cache(mzml_dir, spectra_path)
build_file_mapping_parquet(index_to_filename, mapping_path)
filename_to_index = {v: k for k, v in index_to_filename.items()}

# Find available idXML files
idxml_files = list(idxml_dir.glob("*.idXML"))

if not idxml_files:
    st.warning("No idXML files found.")
    st.stop()

# File selection
selected_idxml = st.selectbox(
    "Select identification file",
    idxml_files,
    format_func=lambda p: p.name,
    key="id_viewer_file",
)

id_df, spectra_data = parse_idxml(str(selected_idxml), filename_to_index)

# Get default values from first row for initial display
first_row = id_df.row(0, named=True)
default_file = first_row["file_index"]
default_scan = first_row["scan_id"]


id_table = Table(
    cache_id="id_table",
    data=id_df.lazy(),
    cache_path=str(viewer_cache),
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

# Always render the LinePlot - filtering handles showing correct data
spectrum_plot = LinePlot(
    cache_id="spectrum_plot",
    data_path=str(spectra_path),
    cache_path=str(viewer_cache),
    filters={"file": "file_index", "spectrum": "scan_id"},
    filter_defaults={"file": default_file, "spectrum": default_scan},
    x_column="mass",
    y_column="intensity",
    title="Mass Spectrum",
    x_label="m/z",
    y_label="Intensity",
    interactivity={"peak": "peak_id"},
    styling={
        "unhighlightedColor": "#4A90D9",
        "selectedColor": "#F3A712",
    },
)

if id_df.height == 0:
    st.warning("No identifications found in the selected file.")
    st.stop()

# Display summary
col1, col2, col3 = st.columns(3)
col1.metric("Identifications", id_df.height)
col2.metric("Source Files", len(spectra_data))
col3.metric("Linked", id_df.filter(pl.col("file_index") >= 0).height)

# Show file mapping
with st.expander("File Index Mapping"):
    st.markdown("**idXML spectra_data order:**")
    for i, fname in enumerate(spectra_data):
        mapped_idx = filename_to_index.get(fname, -1)
        st.write(f"  id_merge_index={i} → `{fname}` → file_index={mapped_idx}")

    st.markdown("**Spectra cache file order:**")
    for idx, fname in sorted(index_to_filename.items()):
        st.write(f"  file_index={idx} → `{fname}`")

# Create StateManager for cross-component linking
state_manager = StateManager(session_key="id_viewer_state")

st.subheader("Peptide Identifications")

id_table(key="id_table", state_manager=state_manager, height=500)

# Show current selection info below
current_state = state_manager.get_state_for_vue()
selected_file_index = current_state.get("file")
selected_scan_id = current_state.get("spectrum")

if selected_file_index is not None and selected_scan_id is not None:
    filename = index_to_filename.get(selected_file_index, "Unknown")

    # Find the selected identification for additional info
    selected_id = id_df.filter(
        (pl.col("file_index") == selected_file_index) &
        (pl.col("scan_id") == selected_scan_id)
    )

    if selected_id.height > 0:
        row = selected_id.row(0, named=True)
        st.info(f"**Selected:** {row['sequence']} (z={row['charge']}, score={row['score']:.4f}) | File: `{filename}` | Scan: {selected_scan_id}")

        # SequenceView for fragment ion visualization
        st.markdown("---")
        st.subheader(f"Peptide: {row['sequence']}")

        # Load spectra data to get observed masses
        spectra_df = pl.read_parquet(spectra_path)
        spectrum_peaks = spectra_df.filter(
            (pl.col("file_index") == selected_file_index) &
            (pl.col("scan_id") == selected_scan_id)
        )

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

        # Unique key per sequence (used for cache_id and component keys)
        seq_hash = hashlib.md5(row["sequence"].encode()).hexdigest()[:8]

        # Compute peak annotations for fragment ions
        if spectrum_peaks.height > 0:
            annotated_peaks = compute_peak_annotations(
                spectrum_peaks,
                sequence_data,
                precursor_charge=row["charge"],
                tolerance=20.0,
                tolerance_ppm=True,
                ion_types=['b', 'y'],
            )

            # Show annotated spectrum plot
            st.subheader("Annotated Spectrum")

            annotated_plot = LinePlot(
                cache_id=f"annotated_spectrum_{seq_hash}",
                data=annotated_peaks.lazy(),
                cache_path=str(viewer_cache),
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

        # Create SequenceView
        sequence_view = SequenceView(
            cache_id="sequence_view",
            sequence=row["sequence"],
            observed_masses=observed_masses,
            peak_ids=peak_ids,
            precursor_mass=0.0,  # Not needed for display
            cache_path=str(viewer_cache),
            deconvolved=False,  # Using m/z data
            precursor_charge=row["charge"],
            interactivity={"peak": "peak_id"},
            _precomputed_sequence_data=sequence_data,
        )

        sequence_view(key=f"sequence_view_{seq_hash}", state_manager=state_manager, height=500)
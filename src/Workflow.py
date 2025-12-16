import json
import streamlit as st
from src.workflow.WorkflowManager import WorkflowManager
from src.workflow.ParameterManager import ParameterManager

# for result section:
import hashlib
from pathlib import Path
import polars as pl

from utils.parse_idxml import parse_idxml
from utils.build_spectra_cache import build_spectra_cache
from utils.sequence_utils import build_sequence_data, compute_peak_annotations

from openms_insight import Table, LinePlot, SequenceView, StateManager


class Workflow(WorkflowManager):
    # Setup pages for upload, parameter, execution and results.
    # For layout use any streamlit components such as tabs (as shown in example), columns, or even expanders.
    def __init__(self) -> None:
        # Initialize the parent class with the workflow name.
        super().__init__("TOPP Workflow", st.session_state["workspace"])

    def upload(self) -> None:
        t = st.tabs(["MS data", "FASTA database"])
        with t[0]:
            # Use the upload method from StreamlitUI to handle mzML file uploads.
            self.ui.upload_widget(
                key="mzML-files",
                name="MS data",
                file_types="mzML",
                fallback=[str(f) for f in Path("example-data", "mzML").glob("*.mzML")],
            )
        with t[1]:
            self.ui.upload_widget(
                key="fasta-file",
                name="FASTA database",
                file_types="fasta",
                fallback=[str(f) for f in Path("example-data", "fasta").glob("*.fasta")],
            )

    @st.fragment
    def configure(self) -> None:
        # Allow users to select mzML files for the analysis.
        self.ui.select_input_file("mzML-files", multiple=True)
        self.ui.select_input_file("fasta-file", multiple=False)

        # Create tabs for different analysis steps.
        t = st.tabs(
            ["**Search Parameters**", "**Filter Parameters**"]
        )
        with t[0]:
            # Load Presets
            presets_path = Path("src", "assets", "comet_presets.json")
            with open(presets_path, "r") as f:
                presets = json.load(f)
            
            st.markdown("##### Load Presets")
            cols = st.columns(len(presets))
            for i, (name, params) in enumerate(presets.items()):
                if cols[i].button(name, use_container_width=True):
                    # Update params
                    if "CometAdapter" not in self.params:
                        self.params["CometAdapter"] = {}
                    
                    tool_prefix = f"{self.parameter_manager.topp_param_prefix}CometAdapter:1:"
                    
                    for k, v in params.items():
                        # Update self.params
                        self.params["CometAdapter"][k] = v
                        # Update st.session_state
                        full_key = f"{tool_prefix}{k}"
                        st.session_state[full_key] = v
                        
                    self.parameter_manager.save_parameters()
                    st.rerun()

            # Parameters for CometAdapter
            self.ui.input_TOPP(
                "CometAdapter",
                exclude_parameters=["binary_modifications","missed_cleavages"]
            )
        with t[1]:
            # Parameters for IDFilter TOPP tool.
            st.markdown("#### ID Filter Settings")
            st.caption("Configure FDR filter threshold (score:peptide) and Peptide length threshold (precursor:length).")
            self.ui.input_TOPP(
                "IDFilter",
                include_parameters=["score:peptide", "precursor:length"],
                custom_defaults={
                    "score:peptide": 0.01,
                    "precursor:length": "8:12"
                },
                display_tool_name=False
            )
            

    def execution(self) -> None:
        # Reload parameters to ensure we have the latest from the UI
        self.params = self.parameter_manager.get_parameters_from_json()

        # Any parameter checks
        #if "mzML-files" not in self.params or not self.params["mzML-files"]:
        #    self.logger.log("ERROR: No mzML files selected.")
        #    return
        #if "fasta-file" not in self.params or not self.params["fasta-file"]:
        #    self.logger.log("ERROR: No FASTA file selected.")
        #    return

        # 1. Input Data
        in_mzML = self.file_manager.get_files(self.params["mzML-files"])
        in_fasta = self.file_manager.get_files(self.params["fasta-file"])
        
        self.logger.log(f"Number of input mzML files: {len(in_mzML)}")

        # 2. Decoy Generation
        # Use first FASTA file for decoy generation
        self.logger.log("Generating decoys...")
        out_decoy_db = self.file_manager.get_files(
            in_fasta[0:1], 
            set_file_type="fasta", 
            set_results_dir="decoy_database", 
        )
        self.executor.run_topp(
            "DecoyDatabase",
            input_output={"in": in_fasta[0:1], "out": out_decoy_db},
            custom_params={
                "decoy_string": "DECOY_", 
                "decoy_string_position": "prefix",
                "enzyme": "no cleavage"
            }
        )

        # 3. Search Engine (Comet)
        # Input: Link mzML files with the (single) Target-Decoy DB
        self.logger.log("Running CometAdapter...")
        out_comet = self.file_manager.get_files(
            in_mzML, set_file_type="idXML", set_results_dir="comet"
        )
        self.executor.run_topp(
            "CometAdapter",
            input_output={"in": in_mzML, "out": out_comet, "database": out_decoy_db},
            #custom_params={
            #    "missed_cleavages": "0",
            #}
        )

        # 4. PeptideIndexer
        self.logger.log("Running PeptideIndexer...")
        out_indexer = self.file_manager.get_files(
            out_comet, set_file_type="idXML", set_results_dir="peptide_indexer"
        )
        self.executor.run_topp(
            "PeptideIndexer",
            input_output={"in": out_comet, "out": out_indexer, "fasta": out_decoy_db},
            custom_params={
                "decoy_string": "DECOY",
                "enzyme:specificity": "none"
            }
        )

        # 5. IDMerger
        self.logger.log("Merging idXML files...")
        out_merged = self.file_manager.get_files(
            "merged.idXML", 
            set_results_dir="id_merger", 
        )
        # Pass list of files as nested list to indicate merging (all inputs to one command)
        self.executor.run_topp(
            "IDMerger",
            input_output={"in": [out_indexer], "out": out_merged},
            custom_params={
                "annotate_file_origin": "true",
                "merge_proteins_add_PSMs": ""
            }
        )

        # 5.5 MS2Rescore
        self.logger.log("Running ms2rescore...")
        in_idxml = out_merged[0]
        mzml_dir = str(Path(in_mzML[0]).parent)
        
        # Prepare content for ms2rescore
        out_ms2rescore = self.file_manager.get_files(
            "merged_ms2rescore.idXML", 
            set_results_dir="ms2rescore", 
        )
        out_ms2rescore_path = out_ms2rescore[0]
        
        # Determine output stem (remove .idXML extension for ms2rescore argument)
        out_stem = str(Path(out_ms2rescore_path).with_suffix(""))
        
        # Calculate tolerance (default 0.02 -> 0.04)
        frag_tol = float(self.params.get("CometAdapter", {}).get("fragment_mass_tolerance", 0.02))
        ms2_tol = 2 * frag_tol

        # Use the wrapper script in src/python-tools
        wrapper_script = Path("src", "python-tools", "ms2rescore_wrapper.py")
        
        cmd = [
            "python", str(wrapper_script),
            "--psm_file", str(in_idxml),
            "--spectrum_path", mzml_dir,
            "--output_path", out_stem + ".idXML", # Wrapper writes to this file
            "--processes", "1",
            "--ms2_tolerance", str(ms2_tol),
            "--ms2pip_model", "Immuno-HCD",
            "--feature_generators", "deeplc,ms2pip",
            "--rescoring_engine", "percolator"
        ]
        
        self.executor.run_command(cmd)

        # Parse feature names
        # The wrapper doesn't explicitly write feature names to a separate file, 
        # but ms2rescore library might. 
        # If not, we might need to assume a name or extract from output idXML?
        # mhcquant expects '*_feature_names.tsv'. 
        # Let's hope ms2rescore writes it.
        feature_file = Path(out_stem + "_feature_names.tsv")
        extra_features = []
        if feature_file.exists():
            with open(feature_file, "r") as f:
                # mhcquant approach: feature names one per line or TSV?
                # Check lines
                lines = f.readlines()
                for line in lines:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        # feature_generator (0), feature_name (1)
                        if "psm_file" not in parts[0]: 
                             extra_features.append(parts[1])
        else:
            self.logger.log(f"Warning: Feature file {feature_file} not found. Proceeding without extra features.")

        # 6. PSMFeatureExtractor
        self.logger.log("Running PSMFeatureExtractor...")
        out_psm = self.file_manager.get_files(
            out_ms2rescore, 
            set_file_type="idXML", 
            set_results_dir="psm_feature_extractor", 
        )
        self.executor.run_topp(
            "PSMFeatureExtractor",
            input_output={"in": out_ms2rescore, "out": out_psm},
            custom_params={
                "extra": extra_features
            }
        )

        # 7. PercolatorAdapter (Rescoring)
        self.logger.log("Running PercolatorAdapter...")
        out_percolator = self.file_manager.get_files(
            out_psm, 
            set_file_type="idXML", 
            set_results_dir="percolator", 
        )
        self.executor.run_topp(
            "PercolatorAdapter",
            input_output={"in": out_psm, "out": out_percolator},
            custom_params={
                "subset_max_train": "0",
                "seed": 4711,
                "trainFDR": 0.05,
                "testFDR": 0.05,
                "enzyme": "no_enzyme",
                "post_processing_tdc": ""
            }
        )

        # 8. IDFilter
        self.logger.log("Running IDFilter...")
        out_filtered = self.file_manager.get_files(
            out_percolator, 
            set_file_type="idXML", 
            set_results_dir="id_filter", 
        )
        self.executor.run_topp(
            "IDFilter",
            input_output={"in": out_percolator, "out": out_filtered},
            custom_params={
                "delete_unreferenced_peptide_hits": "",
                "remove_decoys": ""
            }
        )

        # Postprocessing
        self.logger.log("Postprocessing...")
        # Create cache directories
        cache_dir = self.file_manager.workflow_dir / '.cache'
        cache_dir.mkdir(parents=True, exist_ok=True)
        viewer_cache_dir = cache_dir / "viewer"
        viewer_cache_dir.mkdir(parents=True, exist_ok=True)
        self.logger.log(self.params['mzML-files'])
        self.logger.log(Path(in_mzML[0]).parent)
        # Build spectra cache from mzML files
        index_to_filename = build_spectra_cache(
            Path(in_mzML[0]).parent, cache_dir / "spectra.parquet"
        )
        filename_to_index = {v: k for k, v in index_to_filename.items()}

        # Parse idXML file
        id_df, _ = parse_idxml(out_filtered[0], filename_to_index)

        # Cache identifications to parquet
        id_df.write_parquet(cache_dir / 'identifications.parquet')



    @st.fragment
    def results(self) -> None:
        cache_dir = self.file_manager.workflow_dir / '.cache'
        viewer_cache_dir = cache_dir / "viewer"

        # Load data from parquet caches
        id_df = pl.read_parquet(cache_dir / 'identifications.parquet')
        spectra_df = pl.read_parquet(cache_dir / "spectra.parquet")

        # Check for empty data
        if id_df.height == 0:
            st.warning("No identifications found in the selected file.")
            st.stop()

        # Create StateManager for cross-component linking
        state_manager = StateManager(session_key="id_viewer_state")

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


import json
import streamlit as st
from src.workflow.WorkflowManager import WorkflowManager

# for result section:
from pathlib import Path
import polars as pl

from utils.parse_idxml import parse_idxml
from utils.build_spectra_cache import build_spectra_cache
from utils.split_idxml import split_idxml_by_file

from src.integration import render_toppview_button, is_toppview_available

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
                exclude_parameters=["binary_modifications","missed_cleavages"],
                custom_defaults={"fixed_modifications": []}
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

        # # 5.5 MS2Rescore
        # self.logger.log("Running ms2rescore...")
        # in_idxml = out_merged[0]
        # mzml_dir = str(Path(in_mzML[0]).parent)
        
        # # Prepare content for ms2rescore
        # out_ms2rescore = self.file_manager.get_files(
        #     "merged_ms2rescore.idXML", 
        #     set_results_dir="ms2rescore", 
        # )
        # out_ms2rescore_path = out_ms2rescore[0]
        
        # # Determine output stem (remove .idXML extension for ms2rescore argument)
        # out_stem = str(Path(out_ms2rescore_path).with_suffix(""))
        
        # # Calculate tolerance (default 0.02 -> 0.04)
        # frag_tol = float(self.params.get("CometAdapter", {}).get("fragment_mass_tolerance", 0.02))
        # ms2_tol = 2 * frag_tol

        # # Use the wrapper script in src/python-tools
        # wrapper_script = Path("src", "python-tools", "ms2rescore_wrapper.py")
        
        # cmd = [
        #     "python", str(wrapper_script),
        #     "--psm_file", str(in_idxml),
        #     "--spectrum_path", mzml_dir,
        #     "--output_path", out_stem + ".idXML", # Wrapper writes to this file
        #     "--processes", "1",
        #     "--ms2_tolerance", str(ms2_tol),
        #     "--ms2pip_model", "Immuno-HCD",
        #     "--feature_generators", "deeplc,ms2pip",
        #     "--rescoring_engine", "percolator"
        # ]
        
        # self.executor.run_command(cmd)

        # Parse feature names
        # The wrapper doesn't explicitly write feature names to a separate file, 
        # but ms2rescore library might. 
        # If not, we might need to assume a name or extract from output idXML?
        # mhcquant expects '*_feature_names.tsv'. 
        # Let's hope ms2rescore writes it.
        # feature_file = Path(out_stem + "_feature_names.tsv")
        extra_features = []
        # if feature_file.exists():
        #     with open(feature_file, "r") as f:
        #         # mhcquant approach: feature names one per line or TSV?
        #         # Check lines
        #         lines = f.readlines()
        #         for line in lines:
        #             parts = line.strip().split("\t")
        #             if len(parts) >= 2:
        #                 # feature_generator (0), feature_name (1)
        #                 if "psm_file" not in parts[0]: 
        #                      extra_features.append(parts[1])
        # else:
        #     self.logger.log(f"Warning: Feature file {feature_file} not found. Proceeding without extra features.")

        # 6. PSMFeatureExtractor
        self.logger.log("Running PSMFeatureExtractor...")
        out_psm = self.file_manager.get_files(
            # out_ms2rescore, 
            out_merged, 
            set_file_type="idXML", 
            set_results_dir="psm_feature_extractor", 
        )
        self.executor.run_topp(
            "PSMFeatureExtractor",
            # input_output={"in": out_ms2rescore, "out": out_psm},
            input_output={"in": out_merged, "out": out_psm},
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

        # Create cache directory
        cache_dir = self.file_manager.workflow_dir / 'results' / '.cache'
        cache_dir.mkdir(parents=True, exist_ok=True)

        # Parse idXML file
        id_df, filename_to_index = parse_idxml(out_filtered[0])

        # Extract required scans from identifications
        required_scans = set(
            zip(id_df["file_index"].to_list(), id_df["scan_id"].to_list())
        )

        # Build spectra cache from mzML files (only for required scans)
        spectra_df, filename_to_index = build_spectra_cache(
            Path(in_mzML[0]).parent, filename_to_index, required_scans
        )

        # Create identification table component
        Table(
            cache_id="id_table",
            data=id_df.lazy(),
            cache_path=str(cache_dir),
            interactivity={"file": "file_index", "spectrum": "scan_id", "identification": "id_idx"},
            column_definitions=[
                {"field": "sequence", "title": "Sequence", "headerTooltip": True},
                {"field": "charge", "title": "Charge", "sorter": "number", "hozAlign": "right"},
                {"field": "score", "title": "Score", "sorter": "number", "hozAlign": "right",
                "formatter": "money", "formatterParams": {"precision": 4, "symbol": ""}},
                {"field": "protein_accession", "title": "Protein", "headerTooltip": True},
                {"field": "filename", "title": "File"},
            ],
            initial_sort=[{'column': 'score', 'dir': 'desc'}],
            index_field="id_idx",
            title="Identifications",
            default_row=0,
        )

        # Create SequenceView
        # - Uses sequence_data from id_df (filtered by identification selection)
        # - Uses peaks_data from spectra_df (filtered by file + spectrum selection)
        # - Fragment matching happens in Vue, annotations returned to Python
        comet_params = self.parameter_manager.get_topp_parameters("CometAdapter")
        frag_tol = comet_params.get("fragment_mass_tolerance")
        absolute_error = comet_params.get("fragment_error_units") == 'Da'

        sequence_view = SequenceView(
            cache_id="sequence_view",
            sequence_data=id_df.lazy().select(["id_idx", "sequence", "charge"]).rename({
                "id_idx": "sequence_id",
                "charge": "precursor_charge",
            }),
            peaks_data=spectra_df.lazy(),
            filters={
                "identification": "sequence_id",  # Filter sequence by selected identification
                "file": "file_index",              # Filter peaks by selected file
                "spectrum": "scan_id",             # Filter peaks by selected spectrum
            },
            interactivity={"peak": "peak_id"},
            deconvolved=False,
            annotation_config={
                "ion_types": ["b", "y"],
                "neutral_losses": True,
                "proton_loss_addition": True,
                "tolerance": frag_tol,
                "tolerance_ppm": not absolute_error,
            },
            cache_path=str(cache_dir),
        )

        # Create linked LinePlot using factory method
        # Automatically gets annotations from SequenceView
        LinePlot.from_sequence_view(
            sequence_view,
            cache_id="annotated_spectrum",
            cache_path=str(cache_dir),
            title="Annotated Spectrum",
            styling={
                "unhighlightedColor": "#CCCCCC",
                "highlightColor": "#E74C3C",
                "selectedColor": "#F3A712",
            },
        )


    @st.fragment
    def results(self) -> None:
        cache_dir = self.file_manager.workflow_dir / 'results' / '.cache'

        # Check if workflow has been run
        if not (cache_dir / 'id_table').is_dir():
            st.warning("Please run a workflow to display results.")
            st.stop()

        # TOPPView-Lite integration button
        if is_toppview_available():
            # Get mzML paths from params
            mzml_paths = [Path(p) for p in self.params.get("mzML-files", [])]

            # Get merged idXML and split it
            id_filter_dir = self.file_manager.workflow_dir / 'results' / 'id_filter'
            merged_idxmls = list(id_filter_dir.glob("*.idXML")) if id_filter_dir.exists() else []
            merged_idxml = merged_idxmls[0] if merged_idxmls else None

            if merged_idxml and mzml_paths:
                # Split idXML by source file (cached)
                split_cache_dir = cache_dir / 'split_idxml'
                split_mapping = split_idxml_by_file(merged_idxml, split_cache_dir)

                # Match idXML files to mzML files by stem
                idxml_paths = [
                    split_mapping[p.stem] for p in mzml_paths
                    if p.stem in split_mapping
                ]

                render_toppview_button(
                    mzml_paths=mzml_paths,
                    idxml_paths=idxml_paths,
                    app_name="MHCquant",
                )

        # Create StateManager for cross-component linking
        state_manager = StateManager(session_key="id_viewer_state")

        # Load components from cache
        id_table = Table(cache_id="id_table", cache_path=str(cache_dir))
        sequence_view = SequenceView(cache_id="sequence_view", cache_path=str(cache_dir))
        annotated_plot = LinePlot(cache_id="annotated_spectrum", cache_path=str(cache_dir))

        # Display identification table
        st.subheader("Peptide Identifications")
        id_table(key="id_table", state_manager=state_manager, height=400)

        sv_result = sequence_view(key="sequence_view", state_manager=state_manager, height=800)

        if sv_result.annotations is not None and sv_result.annotations.height > 0:
            st.caption(f"Matched {sv_result.annotations.height} fragments")

        annotated_plot(
            key="annotated_spectrum",
            state_manager=state_manager,
            height=450,
            sequence_view_key="sequence_view",
        )


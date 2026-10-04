import sys

import pyopenms as poms
import streamlit as st
from src.workflow.WorkflowManager import WorkflowManager

# for result section:
from pathlib import Path
import polars as pl

from utils.parse_idxml import parse_idxml
from utils.build_spectra_cache import build_spectra_cache
from utils.feature_names import read_extra_features
# from utils.split_idxml import split_idxml_by_file

# from src.integration import render_toppview_button, is_toppview_available

from openms_insight import Table, LinePlot, SequenceView, StateManager


# Identification settings of nf-core/mhcquant 3.3.0 (nextflow.config, conf/modules.config)
MHCQUANT_VERSION = "3.3.0"

COMET_DEFAULTS = {
    "instrument": "high_res",
    "enzyme": "unspecific cleavage",
    "activation_method": "ALL",
    "digest_mass_range": "800:2500",
    "precursor_charge": "2:3",
    "precursor_mass_tolerance": 5.0,
    "precursor_error_units": "ppm",
    "fragment_mass_tolerance": 0.01,
    "fragment_bin_offset": 0.0,
    "max_variable_mods_in_peptide": 3,
    "fixed_modifications": [],
    "variable_modifications": ["Oxidation (M)"],
    "num_hits": 1,
    "missed_cleavages": 0,
    "spectrum_batch_size": 0,
    "use_A_ions": "false",
    "use_C_ions": "false",
    "use_X_ions": "false",
    "use_Z_ions": "false",
    "use_NL_ions": "false",
    "remove_precursor_peak": "no",
}

# Visible in the nf-core/mhcquant parameter schema; hidden there -> advanced here; the rest is not shown
COMET_SHOWN = [
    "instrument", "activation_method", "digest_mass_range", "precursor_charge",
    "precursor_mass_tolerance", "precursor_error_units", "fragment_mass_tolerance",
    "fragment_bin_offset", "max_variable_mods_in_peptide", "fixed_modifications",
    "variable_modifications", "num_hits", "use_X_ions", "use_Z_ions", "use_A_ions",
    "use_C_ions", "use_NL_ions", "remove_precursor_peak",
]
COMET_ADVANCED = ["enzyme", "spectrum_batch_size"]

IDFILTER_DEFAULTS = {"score:peptide": 0.01, "precursor:length": "8:12"}
IDFILTER_SHOWN = ["score:peptide", "precursor:length"]

# Valueless TOPP switches; run_topp emits them bare when set to True
TOPP_FLAGS = {
    "IDMerger": ["merge_proteins_add_PSMs"],
    "PercolatorAdapter": ["post_processing_tdc", "peptide_level_fdrs"],
    "IDFilter": ["remove_decoys", "delete_unreferenced_peptide_hits"],
    "TextExporter": ["id:peptides_only"],
}

SUMMARY_COLUMNS = ",".join([
    "peptidoform", "sequence", "score", "score_type", "psm", "rt", "mz", "charge",
    "accessions", "aa_before", "aa_after", "start", "end",
    "COMET:deltaCn", "COMET:deltaLCn", "COMET:lnExpect", "COMET:xcorr",
    "rt_diff_best", "observed_retention_time_best", "predicted_retention_time_best",
    "spec_pearson", "std_abs_diff", "ccs_predicted_im2deep", "ccs_error_im2deep", "ion_mobility",
])


class Workflow(WorkflowManager):
    # Setup pages for upload, parameter, execution and results.
    # For layout use any streamlit components such as tabs (as shown in example), columns, or even expanders.
    def __init__(self) -> None:
        # Initialize the parent class with the workflow name.
        super().__init__("TOPP Workflow", st.session_state["workspace"])

    def upload(self) -> None:
        st.caption(
            "MHCquant takes centroided or profile mzML files. Convert Thermo .raw files with "
            "ThermoRawFileParser and Bruker .d folders with tdf2mzml first."
        )
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
        st.info(
            "The selected mzML files are processed as replicates of one sample. "
            "Use a separate workspace for each sample."
        )
        preset = self.params.get("mhcquant-preset")
        st.caption(
            f"Last applied preset: **{preset}**" if preset
            else f"No preset applied: nf-core/mhcquant {MHCQUANT_VERSION} defaults."
        )

        # Create tabs for different analysis steps.
        t = st.tabs(
            ["**Search Parameters**", "**Rescoring**", "**Filter Parameters**", "**Preprocessing**"]
        )
        with t[0]:
            self.ui.input_widget(
                "skip-decoy-generation", False, "Skip decoy generation", widget_type="checkbox",
                help="Use when the FASTA file already contains decoys with the prefix DECOY_.",
            )
            self._input_topp_curated(
                "CometAdapter", COMET_SHOWN, COMET_ADVANCED, custom_defaults=COMET_DEFAULTS,
            )
        with t[1]:
            self.ui.input_python("ms2rescore_wrapper")
            if st.session_state.get("advanced"):
                self.ui.input_widget(
                    "subset-max-train", 0, "Percolator subset_max_train", widget_type="number",
                    min_value=0, help="Maximum number of PSMs used for Percolator training (0: all).",
                )
        with t[2]:
            # Parameters for IDFilter TOPP tool.
            st.markdown("#### ID Filter Settings")
            st.caption("Configure FDR filter threshold (score:peptide) and Peptide length threshold (precursor:length).")
            self.ui.input_widget(
                "fdr-level", "peptide_level_fdrs", "FDR level", widget_type="selectbox",
                options=["peptide_level_fdrs", "psm_level_fdrs"],
                help="Level at which Percolator computes the false discovery rate.",
            )
            self._input_topp_curated(
                "IDFilter", IDFILTER_SHOWN, custom_defaults=IDFILTER_DEFAULTS,
                flag_parameters=TOPP_FLAGS["IDFilter"], display_tool_name=False,
            )
        with t[3]:
            self.ui.input_widget(
                "run-centroidisation", False, "Centroid spectra (PeakPickerHiRes)", widget_type="checkbox",
                help="Comet needs centroided spectra; enable for profile-mode mzML files.",
            )
            self.ui.input_widget(
                "pick-ms-levels", 2, "MS levels to centroid", widget_type="number", min_value=1, max_value=2,
            )
            self.ui.input_widget(
                "filter-mzml", False, "Clean up mzML files (FileFilter)", widget_type="checkbox",
                help="Removes MS2 spectra without precursor charge (FileFilter -peak_options:rm_pc_charge 0).",
            )

    def _input_topp_curated(self, tool: str, shown: list, advanced: list = (), **kwargs) -> None:
        """input_TOPP limited to `shown` (+ `advanced` behind the toggle); all other parameters stay hidden."""
        if self.parameter_manager.create_ini(tool):
            ini_path = Path(self.parameter_manager.ini_dir, f"{tool}.ini")
            param = poms.Param()
            poms.ParamXMLFile().load(str(ini_path), param)
            prefix = f"{tool}:1:"
            keys = [k.decode().split(prefix, 1)[1] for k in param.keys() if prefix in k.decode()]
            # OpenMS tags some schema-visible parameters as advanced; the tag only affects display
            untagged = False
            for key in shown:
                full_key = f"{prefix}{key}".encode()
                tags = param.getTags(full_key) if param.exists(full_key) else []
                if b"advanced" in tags:
                    param.clearTags(full_key)
                    for tag in tags:
                        if tag != b"advanced":
                            param.addTag(full_key, tag)
                    untagged = True
            if untagged:
                poms.ParamXMLFile().store(str(ini_path), param)
            kwargs["exclude_parameters"] = [k for k in keys if k not in set(shown) | set(advanced)]
        self.ui.input_TOPP(tool, include_parameters=list(shown), **kwargs)

    def _sync_tool_settings(self) -> None:
        """Write pipeline defaults and flag definitions to params.json; the configure page may not have rendered."""
        params = self.parameter_manager.get_parameters_from_json()
        params.setdefault("_defaults", {}).update({"CometAdapter": COMET_DEFAULTS, "IDFilter": IDFILTER_DEFAULTS})
        params.setdefault("_flag_params", {}).update(TOPP_FLAGS)
        self.parameter_manager.write_parameters(params)

    def _link_spectra(self, mzml_files: list) -> Path:
        """Directory with one link per selected mzML, for tools that take a spectrum folder."""
        spectra_dir = Path(self.workflow_dir, "results", "spectra")
        spectra_dir.mkdir(parents=True, exist_ok=True)
        for mzml in mzml_files:
            link = spectra_dir / Path(mzml).name
            if not link.exists():
                link.symlink_to(Path(mzml).resolve())
        return spectra_dir

    def execution(self) -> bool:
        self.params = self.parameter_manager.get_parameters_from_json()
        if not self.params.get("mzML-files"):
            self.logger.log("ERROR: No mzML files selected.")
            return False
        if not self.params.get("fasta-file"):
            self.logger.log("ERROR: No FASTA file selected.")
            return False

        self._sync_tool_settings()
        comet_params = self.parameter_manager.get_merged_params("CometAdapter")

        # 1. Input Data
        in_mzML = self.file_manager.get_files(self.params["mzML-files"])
        in_fasta = self.file_manager.get_files(self.params["fasta-file"])
        self.logger.log(f"Number of input mzML files: {len(in_mzML)}")

        # 2. Optional preprocessing
        mzml = in_mzML
        if self.params.get("run-centroidisation", False):
            self.logger.log("Centroiding spectra...")
            out_picked = self.file_manager.get_files(mzml, set_results_dir="centroided")
            if not self.executor.run_topp(
                "PeakPickerHiRes",
                input_output={"in": mzml, "out": out_picked},
                custom_params={"algorithm:ms_levels": int(self.params.get("pick-ms-levels", 2))},
            ):
                return False
            mzml = out_picked
        if self.params.get("filter-mzml", False):
            self.logger.log("Filtering mzML files...")
            out_filtered_mzml = self.file_manager.get_files(mzml, set_results_dir="filtered_mzml")
            if not self.executor.run_topp(
                "FileFilter",
                input_output={"in": mzml, "out": out_filtered_mzml},
                custom_params={"peak_options:rm_pc_charge": 0},
            ):
                return False
            mzml = out_filtered_mzml
        spectra_dir = self._link_spectra(mzml)

        # 3. Decoy Generation
        if self.params.get("skip-decoy-generation", False):
            database = in_fasta[0:1]
        else:
            self.logger.log("Generating decoys...")
            database = self.file_manager.get_files(
                in_fasta[0:1],
                set_file_type="fasta",
                set_results_dir="decoy_database",
            )
            if not self.executor.run_topp(
                "DecoyDatabase",
                input_output={"in": in_fasta[0:1], "out": database},
                custom_params={
                    "decoy_string": "DECOY_",
                    "decoy_string_position": "prefix",
                    "enzyme": "no cleavage"
                }
            ):
                return False

        # 4. Search Engine (Comet)
        self.logger.log("Running CometAdapter...")
        out_comet = self.file_manager.get_files(
            mzml, set_file_type="idXML", set_results_dir="comet"
        )
        comet_io = {"in": mzml, "out": out_comet, "database": database}
        # run_topp drops empty lists, which would leave Comet on its own defaults; a bare option means empty
        for key in ("fixed_modifications", "variable_modifications"):
            if not comet_params.get(key):
                comet_io[key] = [[]]
        if not self.executor.run_topp("CometAdapter", input_output=comet_io):
            return False

        # 5. PeptideIndexer
        self.logger.log("Running PeptideIndexer...")
        out_indexer = self.file_manager.get_files(
            out_comet, set_file_type="idXML", set_results_dir="peptide_indexer"
        )
        if not self.executor.run_topp(
            "PeptideIndexer",
            input_output={"in": out_comet, "out": out_indexer, "fasta": database},
            custom_params={
                "decoy_string": "DECOY",
                "enzyme:specificity": "none"
            }
        ):
            return False

        # 6. IDMerger
        self.logger.log("Merging idXML files...")
        out_merged = self.file_manager.get_files(
            "merged.idXML",
            set_results_dir="id_merger",
        )
        # Pass list of files as nested list to indicate merging (all inputs to one command)
        if not self.executor.run_topp(
            "IDMerger",
            input_output={"in": [out_indexer], "out": out_merged},
            custom_params={
                "annotate_file_origin": "true",
                "merge_proteins_add_PSMs": True
            }
        ):
            return False

        # 7. MS²Rescore feature generation
        self.logger.log("Running MS²Rescore...")
        out_ms2rescore = self.file_manager.get_files(
            "merged_ms2rescore.idXML",
            set_results_dir="ms2rescore",
        )
        if not self.executor.run_python(
            "ms2rescore_wrapper",
            input_output={
                "in": out_merged[0],
                "spectrum_path": str(spectra_dir),
                "out": out_ms2rescore[0],
                "ms2_tolerance": 2 * float(comet_params.get("fragment_mass_tolerance", 0.01)),
                "processes": self.executor._get_max_threads(),
            }
        ):
            return False
        extra_features = read_extra_features(Path(out_ms2rescore[0]).with_suffix(".feature_names.tsv"))
        self.logger.log(f"Loaded {len(extra_features)} extra features from MS²Rescore")

        # 8. PSMFeatureExtractor
        self.logger.log("Running PSMFeatureExtractor...")
        out_psm = self.file_manager.get_files(
            out_ms2rescore,
            set_file_type="idXML",
            set_results_dir="psm_feature_extractor",
        )
        if not self.executor.run_topp(
            "PSMFeatureExtractor",
            input_output={"in": out_ms2rescore, "out": out_psm},
            custom_params={
                "extra": extra_features
            }
        ):
            return False

        # 9. PercolatorAdapter (Rescoring)
        self.logger.log("Running PercolatorAdapter...")
        out_percolator = self.file_manager.get_files(
            out_psm,
            set_file_type="idXML",
            set_results_dir="percolator",
        )
        if not self.executor.run_topp(
            "PercolatorAdapter",
            input_output={"in": out_psm, "out": out_percolator},
            custom_params={
                "seed": 4711,
                "trainFDR": 0.05,
                "testFDR": 0.05,
                "enzyme": "no_enzyme",
                "subset_max_train": int(self.params.get("subset-max-train", 0)),
                "post_processing_tdc": True,
                "peptide_level_fdrs": self.params.get("fdr-level", "peptide_level_fdrs") == "peptide_level_fdrs",
                "weights": str(Path(out_percolator[0]).parent / "percolator_feature_weights.tsv"),
            }
        ):
            return False

        # 10. IDFilter
        self.logger.log("Running IDFilter...")
        out_filtered = self.file_manager.get_files(
            out_percolator,
            set_file_type="idXML",
            set_results_dir="id_filter",
        )
        if not self.executor.run_topp(
            "IDFilter",
            input_output={"in": out_percolator, "out": out_filtered},
            custom_params={
                "remove_decoys": True,
                "delete_unreferenced_peptide_hits": True,
            }
        ):
            return False

        # 11. mhcquant TSV export
        self.logger.log("Exporting identifications...")
        out_text = self.file_manager.get_files(
            out_filtered, set_file_type="tsv", set_results_dir="text_exporter"
        )
        if not self.executor.run_topp(
            "TextExporter",
            input_output={"in": out_filtered, "out": out_text},
            custom_params={"id:peptides_only": True, "id:add_hit_metavalues": 0, "id:add_metavalues": 0},
        ):
            return False
        tsv_dir = Path(self.workflow_dir, "results", "mhcquant_tsv")
        tsv_dir.mkdir(parents=True, exist_ok=True)
        if not self.executor.run_command([
            sys.executable, str(Path("src", "mhcquant", "summarize_results.py")),
            "--input", out_text[0],
            "--out_prefix", str(tsv_dir / Path(self.workflow_dir).parent.name),
            "--columns", SUMMARY_COLUMNS,
        ]):
            return False

        # Postprocessing
        self.logger.log("Postprocessing...")

        # Parse idXML file
        id_df, filename_to_index = parse_idxml(out_filtered[0])
        if id_df.height == 0:
            self.logger.log(
                "WARNING: No peptides passed FDR filtering; an empty TSV was written. "
                "Consider raising the FDR threshold."
            )
            return True

        # Create cache directory
        cache_dir = self.file_manager.workflow_dir / 'results' / '.cache'
        cache_dir.mkdir(parents=True, exist_ok=True)

        # Extract required scans from identifications
        required_scans = set(
            zip(id_df["file_index"].to_list(), id_df["scan_id"].to_list())
        )

        # Build spectra cache from mzML files (only for required scans)
        spectra_df, filename_to_index = build_spectra_cache(
            spectra_dir, filename_to_index, required_scans
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
            initial_sort=[{'column': 'score', 'dir': 'asc'}],
            index_field="id_idx",
            title="Identifications",
            default_row=0,
        )

        # Create SequenceView
        # - Uses sequence_data from id_df (filtered by identification selection)
        # - Uses peaks_data from spectra_df (filtered by file + spectrum selection)
        # - Fragment matching happens in Vue, annotations returned to Python
        frag_tol = float(comet_params.get("fragment_mass_tolerance", 0.01))
        absolute_error = comet_params.get("fragment_error_units", "Da") == 'Da'

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
        return True


    @st.fragment
    def results(self) -> None:
        results_dir = self.file_manager.workflow_dir / 'results'
        cache_dir = results_dir / '.cache'
        tsv_files = sorted((results_dir / 'mhcquant_tsv').glob('*.tsv'))

        if tsv_files:
            st.download_button(
                "⬇️ Download identifications (mhcquant TSV)",
                data=tsv_files[0].read_bytes(),
                file_name=tsv_files[0].name,
                mime="text/tab-separated-values",
            )

        # Check if workflow has been run
        if not (cache_dir / 'id_table').is_dir():
            if tsv_files:
                st.warning("No peptides passed FDR filtering. Consider raising the FDR threshold.")
                st.stop()
            st.warning("Please run a workflow to display results.")
            st.stop()

        # TOPPView-Lite integration button
        # if is_toppview_available():
        #     # Get mzML paths from params
        #     mzml_paths = [Path(p) for p in self.params.get("mzML-files", [])]

        #     # Get merged idXML and split it
        #     id_filter_dir = self.file_manager.workflow_dir / 'results' / 'id_filter'
        #     merged_idxmls = list(id_filter_dir.glob("*.idXML")) if id_filter_dir.exists() else []
        #     merged_idxml = merged_idxmls[0] if merged_idxmls else None

        #     if merged_idxml and mzml_paths:
        #         # Split idXML by source file (cached)
        #         split_cache_dir = cache_dir / 'split_idxml'
        #         split_mapping = split_idxml_by_file(merged_idxml, split_cache_dir)

        #         # Match idXML files to mzML files by stem
        #         idxml_paths = [
        #             split_mapping[p.stem] for p in mzml_paths
        #             if p.stem in split_mapping
        #         ]

        #         render_toppview_button(
        #             mzml_paths=mzml_paths,
        #             idxml_paths=idxml_paths,
        #             app_name="MHCquant",
        #         )

        # Create StateManager for cross-component linking
        state_manager = StateManager(session_key="id_viewer_state")

        # Load components from cache
        id_table = Table(cache_id="id_table", cache_path=str(cache_dir))
        sequence_view = SequenceView(cache_id="sequence_view", cache_path=str(cache_dir))
        annotated_plot = LinePlot(cache_id="annotated_spectrum", cache_path=str(cache_dir))

        # Display identification table
        st.subheader("Peptide Identifications")
        st.info("Scores are q-values (FDR). Lower scores indicate more confident identifications.")
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

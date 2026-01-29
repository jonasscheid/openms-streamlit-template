"""
MS2Rescore wrapper script for the MHCquant workflow.

Generates PSM features using DeepLC and MS2PIP for improved rescoring with Percolator.
Follows the nf-core/mhcquant approach using psm_utils for idXML I/O.
"""
import importlib.resources
import json
import sys
from pathlib import Path

from psm_utils.io.idxml import IdXMLReader, IdXMLWriter
from psm_utils.psm_list import PSMList
from ms2rescore import package_data, rescore

############################
# default parameter values #
############################

DEFAULTS = [
    {"key": "in", "value": "", "help": "Input idXML file.", "hide": True},
    {"key": "spectrum_path", "value": "", "help": "Path to mzML directory.", "hide": True},
    {"key": "out", "value": "", "help": "Output idXML file.", "hide": True},
    {"key": "processes", "value": 8, "help": "Number of parallel processes.", "hide": True},
    {"key": "ms2_tolerance", "value": 0.04, "help": "MS2 fragment tolerance in Da.", "hide": True},
]


def get_params():
    """Load parameters from JSON file passed as command line argument."""
    if len(sys.argv) > 1:
        with open(sys.argv[1], "r") as f:
            return json.load(f)
    else:
        return {}


def filter_incomplete_psms(psm_list, peptide_ids):
    """
    Filter out PSMs that don't have complete feature sets.

    MS2Rescore may fail to generate features for some PSMs (e.g., if
    spectrum not found). This removes PSMs with incomplete rescoring_features
    and removes corresponding PeptideIdentification entries.

    IMPORTANT: Modifies peptide_ids in-place to preserve pyopenms type
    (PeptideIdentificationList), which is required by IdXMLFile.store().

    Returns:
        tuple: (filtered_psm_list, peptide_ids) - peptide_ids is modified in-place
    """
    if not psm_list:
        return psm_list, peptide_ids

    # Get maximum feature count as the expected standard
    max_features = max(
        (len(psm.rescoring_features) for psm in psm_list if psm.rescoring_features),
        default=0
    )

    if max_features == 0:
        return psm_list, peptide_ids

    # Filter PSMs with incomplete feature sets
    filtered = [
        psm for psm in psm_list
        if psm.rescoring_features and len(psm.rescoring_features) == max_features
    ]

    removed_count = len(psm_list) - len(filtered)
    if removed_count > 0:
        print(f"Filtered out {removed_count} PSMs with incomplete features")

        # Get spectrum_ids that are kept
        kept_spectrum_ids = set(psm.spectrum_id for psm in filtered)

        # Filter peptide_ids in-place using clear() and push_back()
        # to preserve the PeptideIdentificationList type
        kept_pep_ids = [
            pep_id for pep_id in peptide_ids
            if pep_id.getMetaValue("spectrum_reference") in kept_spectrum_ids
        ]
        removed_count = len(peptide_ids) - len(kept_pep_ids)

        peptide_ids.clear()
        for pep_id in kept_pep_ids:
            peptide_ids.push_back(pep_id)

        if removed_count > 0:
            print(f"Removed {removed_count} PeptideIdentification entries")

    # Return as PSMList to preserve type for IdXMLWriter
    return PSMList(psm_list=filtered), peptide_ids


def main():
    params = get_params()

    in_idxml = params["in"]
    spectrum_path = params["spectrum_path"]
    out_idxml = params["out"]
    processes = int(params.get("processes", 8))
    ms2_tolerance = float(params.get("ms2_tolerance", 0.04))

    print(f"Input idXML: {in_idxml}")
    print(f"Spectrum path: {spectrum_path}")
    print(f"Output idXML: {out_idxml}")
    print(f"Processes: {processes}")
    print(f"MS2 tolerance: {ms2_tolerance}")

    # Read PSMs from idXML
    print("Reading idXML file...")
    reader = IdXMLReader(in_idxml, spectrum_path=spectrum_path)
    psm_list = reader.read_file()
    print(f"Read {len(psm_list)} PSMs")

    # Load default config from ms2rescore package (ensures all required keys)
    config = json.load(importlib.resources.open_text(package_data, "config_default.json"))

    # Override with our values (matching nf-core/mhcquant)
    config["ms2rescore"]["psm_file"] = in_idxml
    config["ms2rescore"]["spectrum_path"] = spectrum_path
    config["ms2rescore"]["output_path"] = str(Path(out_idxml).with_suffix(""))
    config["ms2rescore"]["processes"] = processes
    config["ms2rescore"]["rename_to_usi"] = False
    config["ms2rescore"]["fasta_file"] = None
    config["ms2rescore"]["id_decoy_pattern"] = "^DECOY_"
    config["ms2rescore"]["lower_score_is_better"] = True
    config["ms2rescore"]["log_level"] = "info"

    # Configure feature generators (reset to only use deeplc + ms2pip)
    config["ms2rescore"]["feature_generators"] = {
        "deeplc": {
            "deeplc_retrain": False,
            "calibration_set_size": 0.15,
        },
        "ms2pip": {
            "model": "Immuno-HCD",
            "ms2_tolerance": ms2_tolerance,
        },
    }

    # Configure rescoring engine (generate features for external Percolator)
    config["ms2rescore"]["rescoring_engine"] = {"percolator": {}}

    print("Running ms2rescore feature generation...")
    print(f"Feature generators: deeplc, ms2pip (model: Immuno-HCD)")

    # Run ms2rescore - this generates features and writes feature_names.tsv
    rescore(config, psm_list)

    # Filter PSMs with incomplete features (also filters peptide_ids to match)
    psm_list, peptide_ids = filter_incomplete_psms(psm_list, reader.peptide_ids)
    print(f"PSMs after filtering: {len(psm_list)}")

    # Write output idXML with rescoring features
    print(f"Writing output idXML: {out_idxml}")
    writer = IdXMLWriter(out_idxml, protein_ids=reader.protein_ids, peptide_ids=peptide_ids)
    writer.write_file(psm_list)

    print("MS2Rescore complete")


if __name__ == "__main__":
    main()

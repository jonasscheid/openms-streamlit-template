"""
MS²Rescore feature generation for the MHCquant workflow.

Port of nf-core/mhcquant 3.3.0 bin/ms2rescore_cli.py: adds DeepLC/MS²PIP features to a
merged idXML for a separate Percolator run, without rescoring inside MS²Rescore.
"""
import importlib.resources
import json
import logging
import sys
from pathlib import Path

import pyopenms as oms

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

MS2PIP_MODELS = [
    "Immuno-HCD", "HCD", "HCD2021", "HCD2019", "HCDch2", "CID", "CIDch2", "timsTOF",
    "timsTOF2024", "timsTOF2023", "TTOF5600", "TMT", "CID-TMT", "iTRAQ", "iTRAQphospho",
]

DEFAULTS = [
    {"key": "in", "value": "", "hide": True},
    {"key": "spectrum_path", "value": "", "hide": True},
    {"key": "out", "value": "", "hide": True},
    {"key": "processes", "value": 1, "hide": True},
    {"key": "ms2_tolerance", "value": 0.02, "hide": True},
    {
        "key": "feature_generators",
        "value": ["deeplc", "ms2pip"],
        "name": "Feature generators",
        "help": "MS²Rescore feature generators added to the Percolator features.",
        "widget_type": "multiselect",
        "options": ["basic", "deeplc", "ms2pip", "im2deep"],
    },
    {
        "key": "ms2pip_model",
        "value": "Immuno-HCD",
        "name": "MS²PIP model",
        "help": "Fragment intensity model; choose the one matching instrument and fragmentation.",
        "widget_type": "selectbox",
        "options": MS2PIP_MODELS,
    },
    {
        "key": "calibration_set_size",
        "value": 0.15,
        "name": "DeepLC calibration set size",
        "help": "Fraction of PSMs used to calibrate DeepLC retention time predictions.",
        "widget_type": "number",
        "min": 0.01,
        "max": 1.0,
        "step_size": 0.01,
        "advanced": True,
    },
]


def build_config(params: dict) -> dict:
    """MS²Rescore config as built by the pipeline for the Percolator engine (features only)."""
    from ms2rescore import package_data

    config = json.load(importlib.resources.open_text(package_data, "config_default.json"))
    ms2rescore_config = config["ms2rescore"]

    generators = params["feature_generators"]
    if isinstance(generators, str):
        generators = [g for g in generators.split(",") if g]
    ms2rescore_config["feature_generators"] = {}
    if "basic" in generators:
        ms2rescore_config["feature_generators"]["basic"] = {}
    if "ms2pip" in generators:
        ms2rescore_config["feature_generators"]["ms2pip"] = {
            "model": params["ms2pip_model"],
            "ms2_tolerance": float(params["ms2_tolerance"]),
            "model_dir": None,
        }
    if "deeplc" in generators:
        ms2rescore_config["feature_generators"]["deeplc"] = {
            "deeplc_retrain": False,
            "calibration_set_size": float(params["calibration_set_size"]),
        }
    if "im2deep" in generators:
        ms2rescore_config["feature_generators"]["im2deep"] = {}

    # Empty engine: MS²Rescore only adds features, Percolator runs as a separate TOPP step
    ms2rescore_config["rescoring_engine"] = {}
    ms2rescore_config.update(
        {
            "psm_file": params["in"],
            "spectrum_path": params["spectrum_path"],
            "output_path": str(Path(params["out"]).with_suffix("")),
            "log_level": "info",
            "processes": int(params["processes"]),
            "fasta_file": None,
            "id_decoy_pattern": "^DECOY_",
            "lower_score_is_better": True,
        }
    )
    return config


def filter_out_artifact_psms(psm_list, peptide_ids):
    """Drop PeptideHits whose peptidoform could not be processed by all feature generators."""
    num_mandatory_features = max(len(psm.rescoring_features) for psm in psm_list)
    complete = [psm for psm in psm_list if len(psm.rescoring_features) == num_mandatory_features]

    all_peptides = {next(iter(psm.provenance_data.items()))[1] for psm in psm_list}
    complete_peptides = {next(iter(psm.provenance_data.items()))[1] for psm in complete}
    not_supported = all_peptides - complete_peptides
    if not not_supported:
        return peptide_ids

    kept = oms.PeptideIdentificationList()
    for peptide_id in peptide_ids:
        hits = [hit for hit in peptide_id.getHits() if hit.getSequence().toString() not in not_supported]
        if hits:
            peptide_id.setHits(hits)
            kept.push_back(peptide_id)
    logging.info(f"Removed {len(not_supported)} peptides not supported by all feature generators: {not_supported}")
    return kept


def write_idxml(path, protein_ids, peptide_ids, psm_list) -> None:
    """Store the identifications with rescoring features (psm-utils 1.4 IdXMLWriter update semantics)."""
    # psm-utils' IdXMLWriter updates copies when pyOpenMS 3.5 returns them from PeptideIdentificationList
    from psm_utils.io.idxml import RESCORING_FEATURE_LIST

    runs = [Path(r.decode() if isinstance(r, bytes) else r).stem for r in protein_ids[0].getMetaValue("spectra_data")]
    psms_by_spectrum = psm_list.get_psm_dict()[None]
    updated = oms.PeptideIdentificationList()
    for peptide_id in peptide_ids:
        run = runs[peptide_id.getMetaValue("id_merge_index")] if len(runs) > 1 else runs[0]
        psms = psms_by_spectrum[run][peptide_id.getMetaValue("spectrum_reference")]
        psm_by_sequence = {psm.provenance_data[str(psm.peptidoform)]: psm for psm in psms}
        hits = []
        for hit in peptide_id.getHits():
            psm = psm_by_sequence[hit.getSequence().toString()]
            if psm.score is not None:
                hit.setScore(psm.score)
            if psm.rank is not None:
                hit.setRank(psm.rank - 1)
            if psm.qvalue is not None:
                hit.setMetaValue("q-value", psm.qvalue)
            if psm.pep is not None:
                hit.setMetaValue("PEP", psm.pep)
            for feature, value in psm.rescoring_features.items():
                if feature not in RESCORING_FEATURE_LIST:
                    hit.setMetaValue(feature, float(value))
            hits.append(hit)
        peptide_id.setHits(hits)
        updated.push_back(peptide_id)
    oms.IdXMLFile().store(str(path), protein_ids, updated)


def main():
    # Heavy imports stay here: the app imports this module to read DEFAULTS
    from ms2rescore import rescore
    from psm_utils.io.idxml import IdXMLReader

    with open(sys.argv[1], encoding="utf-8") as f:
        params = json.load(f)

    config = build_config(params)
    logging.info(f"MS²Rescore config: {config}")

    reader = IdXMLReader(params["in"])
    psm_list = reader.read_file()
    rescore(config, psm_list)

    peptide_ids = filter_out_artifact_psms(psm_list, reader.peptide_ids)
    write_idxml(params["out"], reader.protein_ids, peptide_ids, psm_list)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
# Adapted from mhcquant/bin/ms2rescore_cli.py
# Written by Jonas Scheid (mhcquant)

import sys
import click
import importlib.resources
import json
import logging
from typing import List

import pandas as pd

from ms2rescore import rescore, package_data
from psm_utils.io.idxml import IdXMLReader, IdXMLWriter
from psm_utils import PSMList
import pyopenms as oms

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def parse_cli_arguments_to_config(**kwargs):
    """Update default MS²Rescore config with CLI arguments"""
    try:
        config = json.load(importlib.resources.open_text(package_data, "config_default.json"))
    except Exception as e:
        logging.warning(f"Could not load config_default.json from package_data: {e}. Trying empty dict.")
        config = {"ms2rescore": {}}

    for key, value in kwargs.items():
        if key in ["ms2pip_model", "ms2_tolerance", "test_fdr", "calibration_set_size"]:
            continue

        elif key == "feature_generators":
            if not value:
                feature_generators = []
            else:
                feature_generators = value.split(",")
            
            config["ms2rescore"]["feature_generators"] = {}
            if "basic" in feature_generators:
                config["ms2rescore"]["feature_generators"]["basic"] = {}
            if "ms2pip" in feature_generators:
                config["ms2rescore"]["feature_generators"]["ms2pip"] = {
                    "model": kwargs["ms2pip_model"],
                    "ms2_tolerance": kwargs["ms2_tolerance"],
                    "model_dir": kwargs["ms2pip_model_dir"],
                }
            if "deeplc" in feature_generators:
                config["ms2rescore"]["feature_generators"]["deeplc"] = {
                    "deeplc_retrain": False,
                    "calibration_set_size": kwargs["calibration_set_size"],
                }

        elif key == "rescoring_engine":
            config["ms2rescore"]["rescoring_engine"] = {}
            if value == "mokapot":
                config["ms2rescore"]["rescoring_engine"]["mokapot"] = {
                    "write_weights": True,
                    "write_txt": False,
                    "write_flashlfq": False,
                    "max_workers": kwargs["processes"],
                    "test_fdr" : kwargs["test_fdr"]
                }
            if value == "percolator":
                logging.info(
                    "Percolator rescoring engine has been specified. Use the idXML containing rescoring features and run Percolator in a separate step."
                )
                continue

        else:
            config["ms2rescore"][key] = value

    return config


def rescore_idxml(input_file, output_file, config) -> None:
    """Rescore PSMs in an idXML file and keep other information unchanged."""
    reader = IdXMLReader(input_file)
    psm_list = reader.read_file()

    rescore(config, psm_list)

    peptide_ids_filtered = filter_out_artifact_psms(psm_list, reader.peptide_ids)

    writer = IdXMLWriter(output_file, reader.protein_ids, peptide_ids_filtered)
    writer.write_file(psm_list)


def filter_out_artifact_psms(
    psm_list: PSMList, peptide_ids: List[oms.PeptideIdentification]
) -> List[oms.PeptideIdentification]:
    """Filter out PeptideHits that could not be processed by all feature generators"""
    if not psm_list:
        return peptide_ids
        
    num_mandatory_features = max([len(psm.rescoring_features) for psm in psm_list])
    new_psm_list = PSMList(psm_list=[psm for psm in psm_list if len(psm.rescoring_features) == num_mandatory_features])

    psm_list_peptides = set([next(iter(psm.provenance_data.items()))[1] for psm in psm_list])
    new_psm_list_peptides = set([next(iter(psm.provenance_data.items()))[1] for psm in new_psm_list])
    not_supported_peptides = psm_list_peptides - new_psm_list_peptides

    if len(not_supported_peptides) == 0:
        return peptide_ids

    new_peptide_ids = []
    for peptide_id in peptide_ids:
        new_hits = []
        for hit in peptide_id.getHits():
            if hit.getSequence().toString() in not_supported_peptides:
                continue
            new_hits.append(hit)
        if len(new_hits) == 0:
            continue
        peptide_id.setHits(new_hits)
        new_peptide_ids.append(peptide_id)
    logging.info(
        f"Removed {len(psm_list_peptides) - len(new_psm_list_peptides)} PSMs. Peptides not supported: {not_supported_peptides}"
    )
    return new_peptide_ids


@click.command()
@click.option("-p", "--psm_file", help="Path to PSM file", required=True)
@click.option("-s", "--spectrum_path", help="Path to spectrum file or directory", required=True)
@click.option("-o", "--output_path", help="Path and stem for output file names")
@click.option("-l", "--log_level", help="Logging level", default="info")
@click.option("-n", "--processes", help="Number of parallel processes", type=int, default=16)
@click.option("-f", "--fasta_file", help="Path to FASTA file")
@click.option("-fg", "--feature_generators", help="Comma-separated list of feature generators", default="")
@click.option("-pipm", "--ms2pip_model", help="MS2PIP model", default="Immuno-HCD")
@click.option("-pipmdir", "--ms2pip_model_dir", help="Path to MS2PIP models", default=None)
@click.option("-ms2tol", "--ms2_tolerance", help="Fragment mass tolerance [Da]", type=float, default=0.02)
@click.option("-cs", "--calibration_set_size", help="Percentage of number of calibration set for DeepLC", default=0.15)
@click.option("-re", "--rescoring_engine", help="Either mokapot or percolator", default="mokapot")
@click.option("--test_fdr", help="Test FDR for Mokapot", type=float, default=0.05)
@click.option("-d", "--id_decoy_pattern", help="Regex decoy pattern", default="^DECOY_")
@click.option("-lsb", "--lower_score_is_better", help="Interpretation of primary search engine score", default=True)
def main(**kwargs):
    config = parse_cli_arguments_to_config(**kwargs)
    logging.info("MS2Rescore config:")
    logging.info(config)
    rescore_idxml(kwargs["psm_file"], kwargs["output_path"], config)


if __name__ == "__main__":
    sys.exit(main())

"""
Tests for src/python-tools/ms2rescore_wrapper.py.

Builds a small merged idXML (two runs that share spectrum references) and checks that
PSMs without a full feature set are removed without breaking the idXML write.
"""
import importlib.util
from pathlib import Path

import pytest

oms = pytest.importorskip("pyopenms")
pytest.importorskip("psm_utils")
pytest.importorskip("ms2rescore")

from psm_utils.io.idxml import IdXMLReader  # noqa: E402

WRAPPER = Path(__file__).resolve().parents[1] / "src" / "python-tools" / "ms2rescore_wrapper.py"


def load_wrapper():
    spec = importlib.util.spec_from_file_location("ms2rescore_wrapper", WRAPPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_merged_idxml(path: Path, psms: list[tuple[int, str, str, float]]) -> None:
    """psms: (run index, spectrum reference, sequence, e-value)."""
    protein_id = oms.ProteinIdentification()
    protein_id.setIdentifier("run")
    protein_id.setSearchEngine("Comet")
    protein_id.setScoreType("expect")
    protein_id.setMetaValue("spectra_data", [b"run_a.mzML", b"run_b.mzML"])
    protein_hit = oms.ProteinHit()
    protein_hit.setAccession("P1")
    protein_hit.setMetaValue("target_decoy", "target")
    protein_id.setHits([protein_hit])

    peptide_ids = oms.PeptideIdentificationList()
    for run, spectrum_ref, sequence, evalue in psms:
        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(sequence))
        hit.setCharge(2)
        hit.setScore(evalue)
        hit.setRank(0)
        hit.setMetaValue("target_decoy", "target")
        hit.setMetaValue("MS:1002252", 2.5)
        evidence = oms.PeptideEvidence()
        evidence.setProteinAccession("P1")
        hit.setPeptideEvidences([evidence])
        peptide_id = oms.PeptideIdentification()
        peptide_id.setIdentifier("run")
        peptide_id.setScoreType("expect")
        peptide_id.setHigherScoreBetter(False)
        peptide_id.setRT(100.0 + run)
        peptide_id.setMZ(500.0)
        peptide_id.setMetaValue("spectrum_reference", spectrum_ref)
        peptide_id.setMetaValue("id_merge_index", run)
        peptide_id.setHits([hit])
        peptide_ids.push_back(peptide_id)
    oms.IdXMLFile().store(str(path), [protein_id], peptide_ids)


def test_incomplete_psm_is_removed_across_runs_sharing_scan_numbers(tmp_path):
    wrapper = load_wrapper()
    in_idxml = tmp_path / "merged.idXML"
    write_merged_idxml(
        in_idxml,
        [
            (0, "scan=1", "SIINFEKL", 0.01),
            (1, "scan=1", "SIINFEKLL", 0.02),
            (0, "scan=2", "AAAWYLWEV", 0.03),
            (1, "scan=2", "YLLPAIVHI", 0.04),
        ],
    )
    reader = IdXMLReader(in_idxml)
    psm_list = reader.read_file()
    for psm in psm_list:
        incomplete = psm.peptidoform.sequence == "SIINFEKL"
        psm.rescoring_features = {"f1": 1.0} if incomplete else {"f1": 1.0, "f2": 2.0}

    peptide_ids = wrapper.filter_out_artifact_psms(psm_list, reader.peptide_ids)
    out_idxml = tmp_path / "out.idXML"
    wrapper.write_idxml(out_idxml, reader.protein_ids, peptide_ids, psm_list)

    protein_ids, written = [], oms.PeptideIdentificationList()
    oms.IdXMLFile().load(str(out_idxml), protein_ids, written)
    kept = sorted(
        (written.at(i).getMetaValue("id_merge_index"), written.at(i).getHits()[0].getSequence().toString())
        for i in range(written.size())
    )
    assert kept == [(0, "AAAWYLWEV"), (1, "SIINFEKLL"), (1, "YLLPAIVHI")]
    # pyOpenMS 3.5 lists return copies; features must still reach the stored hits
    assert all(written.at(i).getHits()[0].getMetaValue("f2") == 2.0 for i in range(written.size()))


def test_feature_names_skip_header_and_search_engine_features(tmp_path):
    from utils.feature_names import read_extra_features

    names = tmp_path / "x.feature_names.tsv"
    names.write_text(
        "feature_generator\tfeature_name\n"
        "psm_file\tCOMET:lnExpect\n"
        "ms2pip\tspec_pearson\n"
        "deeplc\trt_diff\n"
    )
    assert read_extra_features(names) == ["spec_pearson", "rt_diff"]

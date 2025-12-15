"""Parse idXML files with direct spectrum linking.

Uses id_merge_index and spectrum_reference for linking identifications
to source mzML files and scan numbers - no RT/mz tolerance matching needed.
"""

import re
from pathlib import Path
from typing import Tuple, List, Dict, Optional

import polars as pl
from pyopenms import IdXMLFile


def parse_spectrum_reference(spectrum_ref: str) -> int:
    """Extract scan number from spectrum_reference string.

    Args:
        spectrum_ref: Native ID string, e.g. "controllerType=0 controllerNumber=1 scan=394"

    Returns:
        Scan number as integer, or -1 if not found
    """
    match = re.search(r'scan=(\d+)', spectrum_ref)
    if match:
        return int(match.group(1))
    return -1


def get_spectra_data_from_idxml(idxml_path: str) -> List[str]:
    """Extract spectra_data list from idXML file.

    Args:
        idxml_path: Path to the idXML file

    Returns:
        List of source mzML filenames from spectra_data UserParam
    """
    protein_ids = []
    peptide_ids = []
    IdXMLFile().load(str(idxml_path), protein_ids, peptide_ids)

    spectra_data = []
    if protein_ids:
        try:
            sd = protein_ids[0].getMetaValue("spectra_data")
            if sd:
                # Handle both list and string formats
                if isinstance(sd, (list, tuple)):
                    # List of strings or bytes
                    spectra_data = [
                        s.decode("utf-8") if isinstance(s, bytes) else str(s)
                        for s in sd
                    ]
                elif isinstance(sd, str):
                    # Parse "[file1.mzML,file2.mzML]" format
                    sd = sd.strip("[]")
                    spectra_data = [s.strip() for s in sd.split(",")]
                elif isinstance(sd, bytes):
                    sd = sd.decode("utf-8").strip("[]")
                    spectra_data = [s.strip() for s in sd.split(",")]
        except Exception as e:
            print(f"Warning: Could not extract spectra_data: {e}")

    return spectra_data


def parse_idxml(
    idxml_path: str,
    filename_to_index: Optional[Dict[str, int]] = None,
) -> Tuple[pl.DataFrame, List[str]]:
    """Parse idXML to DataFrame with direct spectrum linking.

    Extracts identification data and links to source files via:
    - spectra_data: list of source mzML filenames in ProteinIdentification
    - id_merge_index: index into spectra_data per PeptideHit
    - spectrum_reference: native ID containing scan number

    Args:
        idxml_path: Path to the idXML file
        filename_to_index: Optional dict mapping filename to file_index.
                          If provided, file_index column will use these indices.
                          If None, file_index will equal id_merge_index.

    Returns:
        Tuple of:
        - DataFrame with columns: id_idx, sequence, charge, score,
          protein_accession, filename, file_index, scan_id
        - spectra_data list (source mzML filenames)
    """
    protein_ids = []
    peptide_ids = []
    IdXMLFile().load(str(idxml_path), protein_ids, peptide_ids)

    # Get spectra_data from first ProteinIdentification
    spectra_data = []
    if protein_ids:
        try:
            sd = protein_ids[0].getMetaValue("spectra_data")
            if sd:
                # Handle both list and string formats
                if isinstance(sd, (list, tuple)):
                    # List of strings or bytes
                    spectra_data = [
                        s.decode("utf-8") if isinstance(s, bytes) else str(s)
                        for s in sd
                    ]
                elif isinstance(sd, str):
                    # Parse "[file1.mzML,file2.mzML]" format
                    sd = sd.strip("[]")
                    spectra_data = [s.strip() for s in sd.split(",")]
                elif isinstance(sd, bytes):
                    sd = sd.decode("utf-8").strip("[]")
                    spectra_data = [s.strip() for s in sd.split(",")]
        except Exception as e:
            print(f"Warning: Could not extract spectra_data: {e}")

    print(f"Found spectra_data: {spectra_data}")

    # Extract identification data
    id_data = []
    for idx, pep_id in enumerate(peptide_ids):
        hits = pep_id.getHits()
        if not hits:
            continue

        best_hit = hits[0]

        # Get id_merge_index from PeptideIdentification (not PeptideHit!)
        id_merge_index = -1
        try:
            id_merge_index = int(pep_id.getMetaValue("id_merge_index"))
        except Exception:
            pass

        # Determine filename from id_merge_index
        filename = ""
        if 0 <= id_merge_index < len(spectra_data):
            filename = spectra_data[id_merge_index]

        # Map filename to file_index
        file_index = -1
        if filename_to_index is not None and filename:
            file_index = filename_to_index.get(filename, -1)
        else:
            # Fallback: use id_merge_index as file_index
            file_index = id_merge_index

        # Parse scan number from spectrum_reference
        # Use getSpectrumReference() as the primary API method
        spectrum_ref = ""
        try:
            spectrum_ref = pep_id.getSpectrumReference()
            if isinstance(spectrum_ref, bytes):
                spectrum_ref = spectrum_ref.decode("utf-8")
        except Exception:
            pass

        # Fallback to meta value if getSpectrumReference() didn't work
        if not spectrum_ref:
            try:
                spectrum_ref = pep_id.getMetaValue("spectrum_reference")
                if isinstance(spectrum_ref, bytes):
                    spectrum_ref = spectrum_ref.decode("utf-8")
            except Exception:
                pass

        scan_id = parse_spectrum_reference(str(spectrum_ref) if spectrum_ref else "")

        # Extract sequence and other data
        sequence = best_hit.getSequence().toString()
        charge = best_hit.getCharge()
        score = best_hit.getScore()

        # Get protein accessions
        protein_accessions = []
        for evidence in best_hit.getPeptideEvidences():
            acc = evidence.getProteinAccession()
            if isinstance(acc, bytes):
                acc = acc.decode("utf-8")
            protein_accessions.append(acc)

        id_data.append({
            "id_idx": idx,
            "sequence": sequence,
            "charge": charge,
            "score": score,
            "protein_accession": ";".join(protein_accessions) if protein_accessions else "",
            "filename": filename,
            "file_index": file_index,
            "scan_id": scan_id,
        })

    # Create DataFrame
    if id_data:
        df = pl.DataFrame(id_data).cast({
            "id_idx": pl.Int64,
            "sequence": pl.Utf8,
            "charge": pl.Int32,
            "score": pl.Float64,
            "protein_accession": pl.Utf8,
            "filename": pl.Utf8,
            "file_index": pl.Int32,
            "scan_id": pl.Int32,
        })
    else:
        df = pl.DataFrame({
            "id_idx": pl.Series([], dtype=pl.Int64),
            "sequence": pl.Series([], dtype=pl.Utf8),
            "charge": pl.Series([], dtype=pl.Int32),
            "score": pl.Series([], dtype=pl.Float64),
            "protein_accession": pl.Series([], dtype=pl.Utf8),
            "filename": pl.Series([], dtype=pl.Utf8),
            "file_index": pl.Series([], dtype=pl.Int32),
            "scan_id": pl.Series([], dtype=pl.Int32),
        })

    return df, spectra_data


def get_unique_filenames(id_df: pl.DataFrame) -> List[str]:
    """Get list of unique source filenames from parsed ID DataFrame."""
    return id_df["filename"].unique().to_list()

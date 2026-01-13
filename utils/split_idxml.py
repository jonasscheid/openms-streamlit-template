"""Split merged idXML files into separate files per source mzML.

Used for TOPPView-Lite integration which expects separate idXML files
rather than the merged output from IDMerger.
"""

from pathlib import Path
from typing import Dict, List
from collections import defaultdict

from pyopenms import IdXMLFile, ProteinIdentification, PeptideIdentification, PeptideIdentificationList


def split_idxml_by_file(
    merged_idxml_path: Path,
    output_dir: Path,
) -> Dict[str, Path]:
    """Split a merged idXML into separate files per source mzML.

    Takes a merged idXML file (from IDMerger with annotate_file_origin=true)
    and splits it into individual idXML files, one per source mzML file.

    Args:
        merged_idxml_path: Path to merged idXML file
        output_dir: Directory to write split idXML files

    Returns:
        Dict mapping mzML filename stem -> split idXML path
        Example: {"sample1": Path("output_dir/sample1.idXML"), ...}
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if file exists
    if not merged_idxml_path.exists():
        return {}

    # Load merged idXML using pyopenms container types
    protein_ids: List[ProteinIdentification] = []
    peptide_ids_container = PeptideIdentificationList()
    IdXMLFile().load(str(merged_idxml_path), protein_ids, peptide_ids_container)

    # Convert to Python list for easier manipulation
    peptide_ids = [peptide_ids_container.at(i) for i in range(peptide_ids_container.size())]

    if not protein_ids:
        return {}

    # Extract spectra_data (list of source mzML filenames)
    spectra_data_raw = protein_ids[0].getMetaValue("spectra_data")
    if not spectra_data_raw:
        return {}

    # Parse spectra_data to get filenames
    spectra_data: List[str] = []
    if isinstance(spectra_data_raw, (list, tuple)):
        spectra_data = [
            Path(s.decode("utf-8") if isinstance(s, bytes) else str(s)).name
            for s in spectra_data_raw
        ]
    elif isinstance(spectra_data_raw, (str, bytes)):
        sd = spectra_data_raw.decode("utf-8") if isinstance(spectra_data_raw, bytes) else spectra_data_raw
        sd = sd.strip("[]")
        spectra_data = [Path(s.strip()).name for s in sd.split(",")]

    if not spectra_data:
        return {}

    # Group PeptideIdentifications by id_merge_index
    peptide_groups: Dict[int, List[PeptideIdentification]] = defaultdict(list)
    for pep_id in peptide_ids:
        file_index = int(pep_id.getMetaValue("id_merge_index"))
        peptide_groups[file_index].append(pep_id)

    # Create output mapping
    result: Dict[str, Path] = {}

    # Write separate idXML for each source file
    for file_index, filename in enumerate(spectra_data):
        if file_index not in peptide_groups:
            continue  # No identifications for this file

        stem = Path(filename).stem
        output_path = output_dir / f"{stem}.idXML"

        # Skip if already cached
        if output_path.exists():
            result[stem] = output_path
            continue

        # Create new ProteinIdentification with single spectra_data entry
        new_protein_ids = []
        for prot_id in protein_ids:
            new_prot = ProteinIdentification(prot_id)  # Copy
            # Update spectra_data to only contain this file
            new_prot.setMetaValue("spectra_data", [filename.encode("utf-8")])
            new_protein_ids.append(new_prot)

        # Get peptide IDs for this file and reset id_merge_index to 0
        new_peptide_ids = PeptideIdentificationList()
        for pep_id in peptide_groups[file_index]:
            new_pep = PeptideIdentification(pep_id)  # Copy
            new_pep.setMetaValue("id_merge_index", 0)  # Reset to 0 (single file)
            new_peptide_ids.push_back(new_pep)

        # Write split idXML
        IdXMLFile().store(str(output_path), new_protein_ids, new_peptide_ids)
        result[stem] = output_path

    return result

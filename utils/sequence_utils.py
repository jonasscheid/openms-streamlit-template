"""Sequence utilities for fragment mass calculation and SequenceView.

Provides functions to parse peptide sequences and calculate theoretical
fragment masses for visualization in openms-insight SequenceView.
"""

from typing import List, Dict, Tuple, Optional, Any

import numpy as np
import polars as pl
from pyopenms import AASequence, TheoreticalSpectrumGenerator, MSSpectrum

# Proton mass for m/z calculation
PROTON_MASS = 1.007276


def parse_sequence(sequence_str: str) -> Tuple[List[str], List[Optional[float]]]:
    """Parse OpenMS sequence format to extract residues and modification mass shifts.

    Converts e.g. 'SHC(Carbamidomethyl)IAEVEK' to:
    - residues: ['S', 'H', 'C', 'I', 'A', 'E', 'V', 'E', 'K']
    - modifications: [None, None, 57.02, None, None, None, None, None, None]

    Args:
        sequence_str: Peptide sequence in OpenMS format with modifications

    Returns:
        Tuple of (residues list, modifications list where None means unmodified)
    """
    try:
        aa_seq = AASequence.fromString(sequence_str)
        residues = []
        modifications = []

        for i in range(aa_seq.size()):
            residue = aa_seq.getResidue(i)
            one_letter = residue.getOneLetterCode()
            residues.append(one_letter)

            mod = residue.getModification()
            if mod:
                diff_mono = mod.getDiffMonoMass()
                modifications.append(round(diff_mono, 2))
            else:
                modifications.append(None)

        return residues, modifications
    except Exception:
        # Fallback: extract single letters
        residues = []
        modifications = []
        i = 0
        while i < len(sequence_str):
            if sequence_str[i].isupper():
                residues.append(sequence_str[i])
                modifications.append(None)
                i += 1
            elif sequence_str[i] == '(':
                end = sequence_str.find(')', i)
                if end > i:
                    i = end + 1
                else:
                    i += 1
            else:
                i += 1
        return residues, modifications


def calculate_fragment_masses(
    sequence_str: str,
    ion_types: Optional[List[str]] = None,
) -> Dict[str, List[List[float]]]:
    """Calculate theoretical fragment masses using TheoreticalSpectrumGenerator.

    Args:
        sequence_str: Peptide sequence string (can include modifications)
        ion_types: List of ion types to generate. Default: ['a', 'b', 'c', 'x', 'y', 'z']

    Returns:
        Dict with fragment_masses_a, fragment_masses_b, etc.
        Each is a list of lists (one per position, supporting multiple masses).
    """
    if ion_types is None:
        ion_types = ['a', 'b', 'c', 'x', 'y', 'z']

    try:
        aa_seq = AASequence.fromString(sequence_str)
        n = aa_seq.size()

        # Configure TheoreticalSpectrumGenerator
        tsg = TheoreticalSpectrumGenerator()
        params = tsg.getParameters()

        params.setValue("add_a_ions", "true" if 'a' in ion_types else "false")
        params.setValue("add_b_ions", "true" if 'b' in ion_types else "false")
        params.setValue("add_c_ions", "true" if 'c' in ion_types else "false")
        params.setValue("add_x_ions", "true" if 'x' in ion_types else "false")
        params.setValue("add_y_ions", "true" if 'y' in ion_types else "false")
        params.setValue("add_z_ions", "true" if 'z' in ion_types else "false")
        params.setValue("add_metainfo", "true")

        tsg.setParameters(params)

        # Generate spectrum for charge 1 (neutral masses)
        spec = MSSpectrum()
        tsg.getSpectrum(spec, aa_seq, 1, 1)

        # Initialize result dict
        result = {f'fragment_masses_{ion}': [[] for _ in range(n)] for ion in ion_types}

        # Get ion names from StringDataArrays
        ion_names = []
        sdas = spec.getStringDataArrays()
        for sda in sdas:
            if sda.getName() == "IonNames":
                for i in range(sda.size()):
                    name = sda[i]
                    if isinstance(name, bytes):
                        name = name.decode('utf-8')
                    ion_names.append(name)
                break

        # Parse peaks and organize by ion type and position
        for i in range(spec.size()):
            peak = spec[i]
            mz = peak.getMZ()
            ion_name = ion_names[i] if i < len(ion_names) else ""

            if not ion_name:
                continue

            # Parse ion name (e.g., "b3+", "y5++")
            ion_type = None
            ion_number = None

            for t in ion_types:
                if ion_name.lower().startswith(t):
                    ion_type = t
                    try:
                        num_str = ""
                        for c in ion_name[1:]:
                            if c.isdigit():
                                num_str += c
                            else:
                                break
                        if num_str:
                            ion_number = int(num_str)
                    except (ValueError, IndexError):
                        pass
                    break

            if ion_type and ion_number and 1 <= ion_number <= n:
                idx = ion_number - 1
                key = f'fragment_masses_{ion_type}'
                if idx < len(result[key]):
                    result[key][idx].append(mz)

        return result

    except Exception as e:
        print(f"Error calculating fragments for {sequence_str}: {e}")
        return {f'fragment_masses_{ion}': [] for ion in (ion_types or ['a', 'b', 'c', 'x', 'y', 'z'])}


def get_theoretical_mass(sequence_str: str) -> float:
    """Calculate monoisotopic mass of a peptide sequence."""
    try:
        aa_seq = AASequence.fromString(sequence_str)
        return aa_seq.getMonoWeight()
    except Exception:
        return 0.0


def build_sequence_data(
    sequence_str: str,
    charge: int = 2,
    fragment_tolerance: float = 20.0,
    fragment_tolerance_ppm: bool = True,
) -> Dict:
    """Build sequence_data dict for SequenceView component.

    Args:
        sequence_str: Peptide sequence in OpenMS format
        charge: Precursor charge state
        fragment_tolerance: Fragment mass tolerance
        fragment_tolerance_ppm: True if tolerance is in ppm, False for Da

    Returns:
        Dict suitable for SequenceView._precomputed_sequence_data
    """
    residues, modifications = parse_sequence(sequence_str)
    fragment_masses = calculate_fragment_masses(sequence_str)
    theoretical_mass = get_theoretical_mass(sequence_str)

    return {
        "sequence": residues,
        "modifications": modifications,
        "theoretical_mass": theoretical_mass,
        "fixed_modifications": [],
        "fragment_masses_a": fragment_masses.get("fragment_masses_a", []),
        "fragment_masses_b": fragment_masses.get("fragment_masses_b", []),
        "fragment_masses_c": fragment_masses.get("fragment_masses_c", []),
        "fragment_masses_x": fragment_masses.get("fragment_masses_x", []),
        "fragment_masses_y": fragment_masses.get("fragment_masses_y", []),
        "fragment_masses_z": fragment_masses.get("fragment_masses_z", []),
        "fragment_tolerance": fragment_tolerance,
        "fragment_tolerance_ppm": fragment_tolerance_ppm,
    }


def compute_peak_annotations(
    peaks_df: pl.DataFrame,
    sequence_data: Dict,
    precursor_charge: int = 2,
    tolerance: float = 20.0,
    tolerance_ppm: bool = True,
    ion_types: Optional[List[str]] = None,
) -> pl.DataFrame:
    """Compute fragment ion annotations for peaks and return annotated DataFrame.

    Matches observed peaks against theoretical fragment masses considering
    charge states from 1 to precursor_charge.

    Args:
        peaks_df: DataFrame with peaks (must have peak_id, mass, intensity columns)
        sequence_data: Dict with fragment_masses_[b,y,...] from build_sequence_data()
        precursor_charge: Maximum charge state to consider
        tolerance: Mass tolerance for matching
        tolerance_ppm: If True, tolerance is in ppm; if False, in Daltons
        ion_types: List of ion types to consider (default: ['b', 'y'])

    Returns:
        DataFrame with added 'highlight' and 'annotation' columns
    """
    if ion_types is None:
        ion_types = ['b', 'y']

    if peaks_df.height == 0:
        return peaks_df.with_columns([
            pl.lit(False).alias("highlight"),
            pl.lit("").alias("annotation"),
        ])

    # Get observed m/z values and peak_ids
    observed_mz = peaks_df["mass"].to_numpy()
    peak_ids = peaks_df["peak_id"].to_list()

    # Superscript digits for charge display
    superscript = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
                   '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'}

    def to_superscript(n: int) -> str:
        return ''.join(superscript.get(d, d) for d in str(n))

    # Build annotations dict: peak_id -> annotation string
    annotations_dict: Dict[int, str] = {}

    # Build list of theoretical fragments
    for ion_type in ion_types:
        key = f"fragment_masses_{ion_type}"
        fragment_list = sequence_data.get(key, [])

        for ion_number, masses in enumerate(fragment_list, start=1):
            if not masses:
                continue

            for neutral_mass in masses:
                if neutral_mass <= 0:
                    continue

                # Generate m/z for each charge state
                for charge in range(1, precursor_charge + 1):
                    theoretical_mz = (neutral_mass + charge * PROTON_MASS) / charge

                    # Calculate tolerance in Da
                    if tolerance_ppm:
                        tol_da = theoretical_mz * tolerance / 1e6
                    else:
                        tol_da = tolerance

                    # Find matching peaks
                    for idx, obs_mz in enumerate(observed_mz):
                        mass_diff = abs(obs_mz - theoretical_mz)
                        if mass_diff <= tol_da:
                            peak_id = peak_ids[idx]
                            charge_str = to_superscript(charge) + "⁺"
                            label = f"{ion_type}{ion_number}{charge_str}"

                            if peak_id not in annotations_dict:
                                annotations_dict[peak_id] = label
                            else:
                                existing = annotations_dict[peak_id]
                                if label not in existing:
                                    annotations_dict[peak_id] = f"{existing}/{label}"

    # Create highlight and annotation columns
    highlights = [peak_id in annotations_dict for peak_id in peak_ids]
    annotations = [annotations_dict.get(peak_id, "") for peak_id in peak_ids]

    return peaks_df.with_columns([
        pl.Series("highlight", highlights),
        pl.Series("annotation", annotations),
    ])

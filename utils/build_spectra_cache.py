"""Build combined spectra parquet from mzML files.

Creates a single parquet file containing all MS2 peaks from all mzML files,
indexed by file_index and scan_id for efficient filtering.

Note: scan_id is extracted from the native ID (e.g., "scan=511"), NOT from
the spectrum index, to match how idXML references spectra.
"""

import re
from pathlib import Path
from typing import List, Dict, Tuple

import polars as pl
from pyopenms import MSExperiment, MzMLFile


def parse_scan_from_native_id(native_id: str) -> int:
    """Extract scan number from native ID string.

    Args:
        native_id: Native ID like "controllerType=0 controllerNumber=1 scan=511"

    Returns:
        Scan number as integer, or -1 if not found
    """
    match = re.search(r'scan=(\d+)', native_id)
    if match:
        return int(match.group(1))
    return -1


def build_spectra_cache(
    mzml_dir: Path,
    filename_to_index : dict
) -> Tuple[pl.DataFrame, Dict[int, str]]:
    """Build combined spectra parquet from all mzML files in a directory.

    Args:
        mzml_dir: Directory containing mzML files
        output_path: Path to write the combined parquet file

    Returns:
        Tuple of:
        - DataFrame with all peaks
        - Dict mapping file_index to filename
    """
    # Get mzML files
    mzml_files = sorted(mzml_dir.glob("*.mzML"))

    if not mzml_files:
        raise ValueError(f"No mzML files found in {mzml_dir}")

    # Extract peaks from all files
    all_peaks = []
    global_peak_id = 0

    for mzml_path in mzml_files:
        filename = mzml_path.name
        file_index = filename_to_index.get(filename)

        if file_index is None:
            continue

        exp = MSExperiment()
        MzMLFile().load(str(mzml_path), exp)

        for spec_index in range(exp.size()):
            spectrum = exp[spec_index]

            # Only extract MS2 spectra
            if spectrum.getMSLevel() != 2:
                continue

            # Extract scan number from native ID to match idXML references
            native_id = spectrum.getNativeID()
            scan_id = parse_scan_from_native_id(native_id)

            mz, intensity = spectrum.get_peaks()

            for m, inten in zip(mz, intensity):
                all_peaks.append({
                    "peak_id": global_peak_id,
                    "file_index": file_index,
                    "scan_id": scan_id,
                    "mass": float(m),
                    "intensity": float(inten),
                })
                global_peak_id += 1

    # Create DataFrame
    if all_peaks:
        df = pl.DataFrame(all_peaks).cast({
            "peak_id": pl.Int64,
            "file_index": pl.Int32,
            "scan_id": pl.Int32,
            "mass": pl.Float64,
            "intensity": pl.Float64,
        })
    else:
        df = pl.DataFrame({
            "peak_id": pl.Series([], dtype=pl.Int64),
            "file_index": pl.Series([], dtype=pl.Int32),
            "scan_id": pl.Series([], dtype=pl.Int32),
            "mass": pl.Series([], dtype=pl.Float64),
            "intensity": pl.Series([], dtype=pl.Float64),
        })

    return df, filename_to_index
"""Build combined spectra parquet from mzML files.

Creates a single parquet file containing all MS2 peaks from all mzML files,
indexed by file_index and scan_id for efficient filtering.

Note: scan_id is extracted from the native ID (e.g., "scan=511" or "index=4419"),
NOT from the spectrum index, to match how idXML references spectra.
"""

import re
from pathlib import Path
from typing import Dict, Set, Tuple

import polars as pl
from pyopenms import MSExperiment, MzMLFile, OnDiscMSExperiment


def parse_scan_from_native_id(native_id: str) -> int:
    """Extract scan/spectrum ID from native ID string.

    Handles both formats:
    - Thermo: "controllerType=0 controllerNumber=1 scan=394" -> 394
    - Index: "index=4419" -> 4419

    Args:
        native_id: Native ID string from mzML spectrum

    Returns:
        Scan/spectrum ID as integer, or -1 if not found
    """
    # Try scan= format first (Thermo)
    match = re.search(r'scan=(\d+)', native_id)
    if match:
        return int(match.group(1))

    # Try index= format
    match = re.search(r'index=(\d+)', native_id)
    if match:
        return int(match.group(1))

    return -1


def build_spectra_cache(
    mzml_dir: Path,
    filename_to_index: dict,
    required_scans: Set[Tuple[int, int]],
) -> Tuple[pl.DataFrame, Dict[int, str]]:
    """Build combined spectra parquet from all mzML files in a directory.

    Uses OnDiscMSExperiment for indexed mzML files to enable lazy loading,
    which dramatically reduces memory usage by only loading spectra that
    are actually needed (those in required_scans).

    Args:
        mzml_dir: Directory containing mzML files
        filename_to_index: Dict mapping filename to file_index
        required_scans: Set of (file_index, scan_id) tuples to extract.
            Only spectra matching these tuples will have peaks extracted.

    Returns:
        Tuple of:
        - DataFrame with all peaks for required spectra
        - Dict mapping file_index to filename (unchanged from input)
    """
    # Get mzML files
    mzml_files = sorted(mzml_dir.glob("*.mzML"))

    if not mzml_files:
        raise ValueError(f"No mzML files found in {mzml_dir}")

    # Extract peaks from all files
    all_peaks = []
    global_peak_id = 0

    for mzml_path in mzml_files:
        filename = mzml_path.name  # Use name (with extension) to match filename_to_index
        file_index = filename_to_index.get(filename)

        if file_index is None:
            continue

        # Try OnDiscMSExperiment for lazy loading (indexed mzML)
        od_exp = OnDiscMSExperiment()
        is_indexed = od_exp.openFile(str(mzml_path))

        if is_indexed:
            # INDEXED: Use lazy loading - only metadata in memory
            meta = od_exp.getMetaData()
            num_spectra = meta.size()
        else:
            # NON-INDEXED: Fallback to full load (still filter by required_scans)
            exp = MSExperiment()
            MzMLFile().load(str(mzml_path), exp)
            num_spectra = exp.size()

        for spec_index in range(num_spectra):
            if is_indexed:
                spec_meta = meta[spec_index]
            else:
                spec_meta = exp[spec_index]

            # Only extract MS2 spectra
            if spec_meta.getMSLevel() != 2:
                continue

            # Extract scan number from native ID to match idXML references
            native_id = spec_meta.getNativeID()
            scan_id = parse_scan_from_native_id(native_id)

            # FILTER: Skip if not in required scans
            if (file_index, scan_id) not in required_scans:
                continue

            # Load spectrum data (on-demand for indexed, already in memory for non-indexed)
            if is_indexed:
                spectrum = od_exp.getSpectrum(spec_index)
            else:
                spectrum = exp[spec_index]

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
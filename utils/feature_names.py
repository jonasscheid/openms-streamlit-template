"""Parse the MS²Rescore feature-name table for PSMFeatureExtractor."""
from pathlib import Path


def read_extra_features(feature_names_tsv: Path) -> list[str]:
    """Feature names for PSMFeatureExtractor -extra; search engine features (psm_file) are excluded."""
    lines = Path(feature_names_tsv).read_text().splitlines()[1:]
    rows = [line.split("\t") for line in lines if "\t" in line]
    return [row[1] for row in rows if "psm_file" not in row[0]]

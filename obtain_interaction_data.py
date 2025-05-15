#!/usr/bin/env python3
"""
ceRNA Data Acquisition Script

Downloads raw interaction tables from multiple databases,
parses and filters them for high-confidence edges, and
harmonizes gene identifiers.
"""

import os
import requests
import pandas as pd

# Configuration: data sources and parameters
DATA_DIR = "data/raw"
INTERACTION_SOURCES = {
    "starBase": {
        "url": "https://starbase.sysu.edu.cn/download/starBase_ceRNA_interactions.tsv.gz",
        "filename": "starBase_ceRNA_interactions.tsv.gz",
        "sep": "\t",
        "score_col": "clip_score",
        "threshold": 0.5
    },
    "LncBase": {
        "url": "https://diana.e-ce.uth.gr/lncbasev2/download/LncBase_interactions.tsv.gz",
        "filename": "LncBase_interactions.tsv.gz",
        "sep": "\t",
        "score_col": "lncbase_confidence",
        "threshold": 0.7
    },
    "miRTarBase": {
        "url": "https://mirtarbase.cuhk.edu.cn/cache/download/2023_MTI.tsv.gz",
        "filename": "miRTarBase_interactions.tsv.gz",
        "sep": "\t",
        "score_col": "SupportType",  # e.g., 'strong_evidence'
        "threshold": None  # assume all entries are experimentally validated
    },
    "miRcode": {
        "url": "http://www.mircode.org/download/mircode_v11.tsv",
        "filename": "miRcode_interactions.tsv",
        "sep": "\t",
        "score_col": None,  # no score column
        "threshold": None
    }
}

ID_MAPPING_FILE = "data/annotations/gene_id_mapping.tsv"  # Tab-delimited: raw_id → official_symbol


def ensure_data_dir():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(ID_MAPPING_FILE), exist_ok=True)


def download_file(url: str, dest_path: str):
    """Download a file from a URL to a local destination."""
    print(f"Downloading {url} → {dest_path}")
    resp = requests.get(url, stream=True)
    resp.raise_for_status()
    with open(dest_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"Downloaded {dest_path}")


def load_and_filter(filepath: str,
                    sep: str = "\t",
                    score_col: str = None,
                    threshold: float = None) -> pd.DataFrame:
    """
    Load an interaction table and filter by confidence score.

    Parameters
    ----------
    filepath : str
        Path to the raw interaction file (can be gzipped).
    sep : str
        Separator used in the file.
    score_col : str or None
        Name of the column containing confidence scores.
    threshold : float or None
        Minimum score to keep. If None, no filtering is applied.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['source', 'target', ...].
    """
    print(f"Loading {filepath}")
    df = pd.read_csv(filepath, sep=sep, compression='infer', low_memory=False)
    # ensure we have `source` and `target` columns (rename if needed)
    # Assuming raw columns named 'miRNA' and 'target_gene' or similar; adjust as needed:
    if 'miRNA' in df.columns and 'target_gene' in df.columns:
        df = df.rename(columns={'miRNA': 'source', 'target_gene': 'target'})
    else:
        # fallback generic
        df = df.rename(columns={df.columns[0]: 'source', df.columns[1]: 'target'})
    
    # filter by confidence score
    if score_col and threshold is not None and score_col in df.columns:
        initial = df.shape[0]
        df = df[df[score_col] >= threshold]
        print(f"  Filtered {initial - df.shape[0]} low-confidence edges (threshold={threshold})")
    
    # drop missing and duplicate edges
    df = df.dropna(subset=['source', 'target'])
    df = df.drop_duplicates(subset=['source', 'target'])
    print(f"  {df.shape[0]} interactions retained after filtering")
    return df[['source', 'target'] + [c for c in df.columns if c not in ['source', 'target']] ]


def harmonize_ids(df: pd.DataFrame, mapping_path: str) -> pd.DataFrame:
    """
    Harmonize gene/transcript IDs to official symbols using a mapping file.

    mapping_path : str
        TSV with columns ['raw_id', 'official_symbol'].
    """
    print(f"Loading ID mapping from {mapping_path}")
    mapping = pd.read_csv(mapping_path, sep="\t")
    mapping_dict = dict(zip(mapping['raw_id'], mapping['official_symbol']))
    
    df['source'] = df['source'].map(mapping_dict).fillna(df['source'])
    df['target'] = df['target'].map(mapping_dict).fillna(df['target'])
    return df


def main():
    ensure_data_dir()
    
    # 1. Download raw tables
    for name, info in INTERACTION_SOURCES.items():
        dest = os.path.join(DATA_DIR, info['filename'])
        if not os.path.exists(dest):
            download_file(info['url'], dest)
        else:
            print(f"{dest} already exists, skipping download.")
    
    # 2. Load, filter, and harmonize each source
    dfs = []
    for name, info in INTERACTION_SOURCES.items():
        path = os.path.join(DATA_DIR, info['filename'])
        df = load_and_filter(
            filepath=path,
            sep=info['sep'],
            score_col=info['score_col'],
            threshold=info['threshold']
        )
        df = harmonize_ids(df, ID_MAPPING_FILE)
        df['source_db'] = name
        dfs.append(df)
    
    # 3. Merge all into a single edge list
    merged = pd.concat(dfs, ignore_index=True)
    merged = merged.drop_duplicates(subset=['source', 'target'])
    print(f"Total unique interactions across all databases: {merged.shape[0]}")
    
    # 4. Save merged edge list
    out_path = "data/processed/ceRNA_interactions_merged.tsv"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"Merged interactions saved to {out_path}")


if __name__ == "__main__":
    main()

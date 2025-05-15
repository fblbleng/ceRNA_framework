#!/usr/bin/env python3
"""
scRNA-seq Data Preprocessing Script

Loads a normalized single-cell RNA-seq expression matrix,
filters genes based on expression and presence in the ceRNA network,
and saves the filtered expression matrix.
"""

import os
import argparse
import pandas as pd

def load_network_genes(network_path: str) -> set:
    """
    Load the ceRNA network edge list and return the set of all genes/lncRNAs/circRNAs.
    """
    df = pd.read_csv(network_path, sep='\t', usecols=['source', 'target'])
    genes = set(df['source']).union(df['target'])
    print(f"[Network] Loaded {len(genes)} unique nodes from {network_path}")
    return genes

def load_expression(expr_path: str) -> pd.DataFrame:
    """
    Load scRNA-seq expression matrix.
    Expects a tab- or comma-delimited file with genes as rows and cells as columns,
    and the first column named 'gene' or used as the row index.
    """
    # Try TSV first, then CSV
    sep = '\t'
    try:
        expr = pd.read_csv(expr_path, sep=sep, index_col=0)
    except Exception:
        sep = ','
        expr = pd.read_csv(expr_path, sep=sep, index_col=0)
    print(f"[Expression] Loaded expression matrix with {expr.shape[0]} genes and {expr.shape[1]} cells")
    return expr

def filter_expression(expr: pd.DataFrame, network_genes: set) -> pd.DataFrame:
    """
    Filter the expression matrix to:
      1) Keep only genes present in the ceRNA network.
      2) Remove genes with zero expression across all cells.
    """
    # Intersect with network genes
    genes_before = expr.shape[0]
    expr = expr.loc[expr.index.intersection(network_genes)]
    print(f"[Filter] {genes_before - expr.shape[0]} genes removed (not in network)")

    # Remove zero-expression genes
    nonzero = (expr.sum(axis=1) > 0)
    expr = expr.loc[nonzero]
    print(f"[Filter] {nonzero.size - nonzero.sum()} genes removed (zero expression)")

    print(f"[Filter] {expr.shape[0]} genes remain after filtering")
    return expr

def save_expression(expr: pd.DataFrame, out_path: str):
    """
    Save the filtered expression matrix to a file.
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    expr.to_csv(out_path, sep='\t')
    print(f"[Output] Filtered expression saved to {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Preprocess scRNA-seq expression for ceRNA integration")
    parser.add_argument('--expr', required=True,
                        help="Path to the normalized scRNA-seq expression matrix (genes x cells)")
    parser.add_argument('--network', required=True,
                        help="Path to the merged ceRNA network edge list (tab-delim, columns 'source','target')")
    parser.add_argument('--out', required=True,
                        help="Output path for the filtered expression matrix")

    args = parser.parse_args()

    network_genes = load_network_genes(args.network)
    expr = load_expression(args.expr)
    expr_filtered = filter_expression(expr, network_genes)
    save_expression(expr_filtered, args.out)

if __name__ == "__main__":
    main()

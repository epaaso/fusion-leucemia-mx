#!/usr/bin/env python3

import pandas as pd
import argparse
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Filter fusion calls for protocol-relevant genes.")
    parser.add_argument("--input", required=True, help="Path to fusion TSV (Arriba or STAR-Fusion)")
    parser.add_argument("--output", required=True, help="Path to output summary TSV")
    parser.add_argument(
        "--format",
        default="auto",
        choices=["auto", "arriba", "starfusion"],
        help="Input format (auto, arriba, starfusion)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file {input_path} not found.")
        sys.exit(1)

    # Load fusions
    try:
        fusions = pd.read_csv(input_path, sep="\t")
    except Exception as e:
        print(f"Error reading fusions file: {e}")
        sys.exit(1)

    if fusions.empty:
        print("No fusions found in input file.")
        # Create empty output with same columns plus protocol_flag
        fusions["protocol_flag"] = []
        fusions.to_csv(args.output, sep="\t", index=False)
        return

    # Define key genes
    key_genes = {
        "DUX4", "IGH", "ETV6", "RUNX1", "TCF3", "PBX1", "PBX3",
        "CRLF2", "JAK1", "JAK2", "EPOR", "ABL1", "ABL2", "PDGFRA", "PDGFRB", "CSF1R", "LYN",
    }

    def normalize_gene(value):
        if pd.isna(value):
            return ""
        gene = str(value)
        gene = gene.split("^")[0]
        gene = gene.split("|")[0]
        gene = gene.split(",")[0]
        return gene.strip()

    def detect_format(columns):
        if ("#gene1" in columns or "gene1" in columns) and "gene2" in columns:
            return "arriba"
        if "LeftGene" in columns and "RightGene" in columns:
            return "starfusion"
        return None

    input_format = args.format
    if input_format == "auto":
        input_format = detect_format(fusions.columns)

    if input_format == "arriba":
        gene1_col = "#gene1" if "#gene1" in fusions.columns else "gene1"
        gene2_col = "gene2"
    elif input_format == "starfusion":
        gene1_col = "LeftGene"
        gene2_col = "RightGene"
    else:
        print("Error: Unable to detect fusion format. Use --format arriba or starfusion.")
        sys.exit(1)

    if gene1_col not in fusions.columns or gene2_col not in fusions.columns:
        print(f"Error: Expected columns '{gene1_col}' and '{gene2_col}' not found.")
        sys.exit(1)

    fusions["protocol_flag"] = fusions.apply(
        lambda row: any(
            g in key_genes
            for g in (normalize_gene(row[gene1_col]), normalize_gene(row[gene2_col]))
        ),
        axis=1,
    )

    # Filter for protocol flag = True
    protocol_fusions = fusions[fusions["protocol_flag"]]
    
    print(f"Total fusions: {len(fusions)}")
    print(f"Protocol fusions: {len(protocol_fusions)}")

    # Write output
    protocol_fusions.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()

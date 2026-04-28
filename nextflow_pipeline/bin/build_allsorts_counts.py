#!/usr/bin/env python3

import argparse
import gzip
import re
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build an ALLSorts-compatible counts matrix (samples x genes) from "
            "STAR ReadsPerGene.out.tab files."
        )
    )
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        metavar="SAMPLE=PATH",
        help="Input pair in the form SAMPLE_ID=/path/to/ReadsPerGene.out.tab",
    )
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--gtf", required=True, help="Reference GTF used to map gene IDs to symbols")
    parser.add_argument(
        "--strand",
        choices=["no", "forward", "reverse"],
        default="reverse",
        help="STAR gene-count strand column to use (default: reverse)",
    )
    return parser.parse_args()


def parse_input_pair(raw_value: str):
    if "=" not in raw_value:
        raise ValueError(f"Invalid --input value '{raw_value}'. Expected SAMPLE=PATH")

    sample_id, path_str = raw_value.split("=", 1)
    sample_id = sample_id.strip()
    path_str = path_str.strip()

    if not sample_id:
        raise ValueError(f"Invalid --input value '{raw_value}'. Empty sample ID")

    input_path = Path(path_str)
    if not input_path.exists():
        raise ValueError(f"Reads file not found for sample '{sample_id}': {input_path}")

    return sample_id, input_path


def build_gene_lookup(gtf_path: Path):
    gene_id_pattern = re.compile(r'gene_id "([^"]+)"')
    gene_name_pattern = re.compile(r'gene_name "([^"]+)"')

    lookup = {}

    open_fn = gzip.open if gtf_path.suffix == ".gz" else open
    with open_fn(gtf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            attributes = fields[8]
            gene_id_match = gene_id_pattern.search(attributes)
            gene_name_match = gene_name_pattern.search(attributes)

            if not gene_id_match or not gene_name_match:
                continue

            gene_id = gene_id_match.group(1)
            gene_name = gene_name_match.group(1)

            if not gene_id or not gene_name:
                continue

            lookup.setdefault(gene_id, gene_name)
            if "." in gene_id:
                lookup.setdefault(gene_id.split(".", 1)[0], gene_name)

    if not lookup:
        raise ValueError(f"No gene mappings loaded from {gtf_path}")

    return lookup


def load_star_counts(reads_path: Path, strand: str):
    table = pd.read_csv(reads_path, sep="\t", header=None)
    if table.shape[1] < 4:
        raise ValueError(f"Unexpected STAR ReadsPerGene format in {reads_path}")

    table = table.iloc[:, :4]
    table.columns = ["gene_id", "no", "forward", "reverse"]
    table = table[~table["gene_id"].astype(str).str.startswith("N_")]

    counts = table[["gene_id", strand]].copy()
    counts[strand] = pd.to_numeric(counts[strand], errors="coerce").fillna(0)
    return counts


def map_gene_ids_to_symbols(counts_df: pd.DataFrame, gene_lookup: dict):
    counts_df = counts_df.copy()

    def resolve_symbol(gene_id: str):
        if gene_id in gene_lookup:
            return gene_lookup[gene_id]
        if "." in gene_id:
            return gene_lookup.get(gene_id.split(".", 1)[0], "")
        return ""

    counts_df["gene_symbol"] = counts_df["gene_id"].astype(str).map(resolve_symbol)
    counts_df = counts_df[counts_df["gene_symbol"] != ""]

    if counts_df.empty:
        return pd.Series(dtype=float)

    grouped = counts_df.groupby("gene_symbol", sort=True).sum(numeric_only=True)
    series = grouped.iloc[:, 0]
    return series


def main():
    args = parse_args()

    try:
        gtf_path = Path(args.gtf)
        if not gtf_path.exists():
            raise ValueError(f"GTF file not found: {gtf_path}")

        parsed_inputs = [parse_input_pair(raw) for raw in args.input]
        gene_lookup = build_gene_lookup(gtf_path)

        rows = []
        for sample_id, reads_path in parsed_inputs:
            counts = load_star_counts(reads_path, args.strand)
            gene_counts = map_gene_ids_to_symbols(counts, gene_lookup)
            gene_counts.name = sample_id
            rows.append(gene_counts)

        if not rows:
            raise ValueError("No valid inputs to build ALLSorts counts matrix")

        matrix = pd.DataFrame(rows)
        matrix.index.name = "sample_id"
        matrix = matrix.fillna(0).astype(int)
        matrix = matrix.reindex(sorted(matrix.columns), axis=1)

        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        matrix.to_csv(output_path)

        print(f"Wrote ALLSorts matrix with {matrix.shape[0]} samples x {matrix.shape[1]} genes to {output_path}")

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

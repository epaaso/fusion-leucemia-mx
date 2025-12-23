#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path


def parse_stats(path):
    input_reads = None
    matched_reads = None
    input_bases = None
    matched_bases = None

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            match_input = re.search(r"Input:\s+(\d+)\s+reads\s+(\d+)\s+bases", line)
            if match_input:
                input_reads = int(match_input.group(1))
                input_bases = int(match_input.group(2))
                continue

            match_matched = re.search(r"Matched:\s+(\d+)\s+reads\s+(\d+)\s+bases", line)
            if match_matched:
                matched_reads = int(match_matched.group(1))
                matched_bases = int(match_matched.group(2))
                continue

    return input_reads, matched_reads, input_bases, matched_bases


def sample_from_path(path):
    name = Path(path).stem
    if name.startswith("rrna_"):
        return name[len("rrna_") :]
    return name


def main():
    parser = argparse.ArgumentParser(description="Summarize rRNA filtering stats.")
    parser.add_argument("--inputs", nargs="+", required=True, help="BBDuk stats files.")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    args = parser.parse_args()

    rows = []
    for stats_path in args.inputs:
        path = Path(stats_path)
        if not path.exists():
            print(f"Warning: stats file not found: {path}", file=sys.stderr)
            continue

        input_reads, matched_reads, input_bases, matched_bases = parse_stats(path)
        if input_reads is None or matched_reads is None:
            print(f"Warning: could not parse stats from {path}", file=sys.stderr)
            continue

        rrna_fraction_reads = matched_reads / input_reads if input_reads else 0.0
        rrna_fraction_bases = matched_bases / input_bases if input_bases else 0.0

        rows.append(
            {
                "sample_id": sample_from_path(path),
                "input_reads": input_reads,
                "matched_reads": matched_reads,
                "rrna_fraction_reads": rrna_fraction_reads,
                "input_bases": input_bases if input_bases is not None else "",
                "matched_bases": matched_bases if matched_bases is not None else "",
                "rrna_fraction_bases": rrna_fraction_bases,
            }
        )

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(
            "sample_id\tinput_reads\tmatched_reads\trrna_fraction_reads\t"
            "input_bases\tmatched_bases\trrna_fraction_bases\n"
        )
        for row in rows:
            handle.write(
                f"{row['sample_id']}\t{row['input_reads']}\t{row['matched_reads']}\t"
                f"{row['rrna_fraction_reads']:.6f}\t{row['input_bases']}\t"
                f"{row['matched_bases']}\t{row['rrna_fraction_bases']:.6f}\n"
            )


if __name__ == "__main__":
    main()

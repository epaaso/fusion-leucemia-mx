#!/usr/bin/env python3

import argparse
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import quote_plus, urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

import pandas as pd

USER_AGENT = "fusion-annotator/1.0"


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
    if "Gene_1_symbol(5end_fusion_partner)" in columns and "Gene_2_symbol(3end_fusion_partner)" in columns:
        return "fusioncatcher"
    return None


def truncate(text, max_chars):
    if not text:
        return ""
    text = " ".join(str(text).split())
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 3].rstrip() + "..."


def first_sentence(text):
    if not text:
        return ""
    text = " ".join(str(text).split())
    match = re.search(r"(.+?[.!?])\s", text)
    if match:
        return match.group(1)
    return text


def fetch_url(url, timeout=15):
    try:
        req = Request(url, headers={"User-Agent": USER_AGENT})
        with urlopen(req, timeout=timeout) as response:
            return response.read()
    except Exception as exc:
        print(f"Warning: Failed to fetch {url}: {exc}")
        return None


def fetch_json(url, sleep_s):
    data = fetch_url(url)
    if sleep_s:
        time.sleep(sleep_s)
    if not data:
        return None
    try:
        return json.loads(data.decode("utf-8", errors="replace"))
    except Exception as exc:
        print(f"Warning: Failed to parse JSON from {url}: {exc}")
        return None


def fetch_text(url, sleep_s):
    data = fetch_url(url)
    if sleep_s:
        time.sleep(sleep_s)
    if not data:
        return ""
    return data.decode("utf-8", errors="replace")


def query_mygene_summary(gene, sleep_s, max_chars):
    query = f"symbol:{gene}"
    params = {
        "q": query,
        "species": "human",
        "fields": "summary,symbol,name",
        "size": 5,
    }
    url = "https://mygene.info/v3/query?" + urlencode(params, quote_via=quote_plus)
    data = fetch_json(url, sleep_s)
    if not data or "hits" not in data:
        return ""

    hits = data.get("hits", [])
    if not hits:
        return ""

    exact = next(
        (hit for hit in hits if str(hit.get("symbol", "")).upper() == gene.upper()),
        hits[0],
    )
    summary = exact.get("summary") or exact.get("name") or ""
    return truncate(summary, max_chars)


def build_pubmed_term(gene1, gene2):
    g1 = gene1.replace('"', "")
    g2 = gene2.replace('"', "")
    return (
        f'("{g1}-{g2}"[Title/Abstract] OR "{g2}-{g1}"[Title/Abstract] OR '
        f'({g1}[Title/Abstract] AND {g2}[Title/Abstract] AND fusion[Title/Abstract]))'
    )


def query_pubmed_summary(gene1, gene2, email, sleep_s, max_chars):
    term = build_pubmed_term(gene1, gene2)
    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "retmax": 1,
        "tool": "fusion-annotator",
    }
    if email:
        params["email"] = email
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        + urlencode(params, quote_via=quote_plus)
    )
    data = fetch_json(search_url, sleep_s)
    id_list = (
        data.get("esearchresult", {}).get("idlist", []) if data else []
    )
    if not id_list:
        return "No PubMed hits found."

    pmid = id_list[0]
    fetch_params = {"db": "pubmed", "id": pmid, "retmode": "xml", "tool": "fusion-annotator"}
    if email:
        fetch_params["email"] = email
    fetch_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        + urlencode(fetch_params, quote_via=quote_plus)
    )
    xml_text = fetch_text(fetch_url, sleep_s)
    if not xml_text:
        return f"PubMed record {pmid} retrieved without abstract."

    try:
        root = ET.fromstring(xml_text)
    except Exception as exc:
        print(f"Warning: Failed to parse PubMed XML for {pmid}: {exc}")
        return f"PubMed record {pmid} retrieved without abstract."

    title = root.findtext(".//ArticleTitle") or ""
    abstract_parts = []
    for elem in root.findall(".//AbstractText"):
        abstract_parts.append("".join(elem.itertext()))
    abstract = " ".join(part for part in abstract_parts if part).strip()

    summary_parts = []
    if title:
        summary_parts.append(title.strip().rstrip("."))
    if abstract:
        summary_parts.append(first_sentence(abstract).strip().rstrip("."))

    summary = ". ".join(part for part in summary_parts if part)
    summary = truncate(summary, max_chars)
    if summary:
        return f"{summary} (PMID:{pmid})"
    return f"PubMed record {pmid} retrieved without abstract."


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate fusion calls across samples and report common fusions."
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Fusion TSVs (Arriba or STAR-Fusion).",
    )
    parser.add_argument("--output", required=True, help="Path to output TSV")
    parser.add_argument(
        "--format",
        default="auto",
        choices=["auto", "arriba", "starfusion", "fusioncatcher"],
        help="Input format (auto, arriba, starfusion, fusioncatcher)",
    )
    parser.add_argument(
        "--annotate",
        dest="annotate",
        action="store_true",
        help="Annotate with fusion literature and gene summaries.",
    )
    parser.add_argument(
        "--no-annotate",
        dest="annotate",
        action="store_false",
        help="Disable literature and gene annotations.",
    )
    parser.add_argument(
        "--email",
        default="",
        help="Email for NCBI E-utilities requests (recommended).",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.34,
        help="Delay between external requests (seconds).",
    )
    parser.add_argument(
        "--max-abstract-chars",
        type=int,
        default=280,
        help="Max characters for fusion literature summaries.",
    )
    parser.add_argument(
        "--max-gene-summary-chars",
        type=int,
        default=220,
        help="Max characters for gene function summaries.",
    )
    parser.set_defaults(annotate=True)
    args = parser.parse_args()

    input_paths = [Path(p) for p in args.inputs]
    total_samples = len(input_paths)
    if total_samples == 0:
        print("Error: No inputs provided.")
        sys.exit(1)

    sample_counts = {}
    event_counts = {}

    for input_path in input_paths:
        if not input_path.exists():
            print(f"Warning: Input file {input_path} not found. Skipping.")
            continue

        try:
            fusions = pd.read_csv(input_path, sep="\t")
        except Exception as exc:
            print(f"Warning: Failed to read {input_path}: {exc}")
            continue

        if fusions.empty:
            continue

        input_format = args.format
        if input_format == "auto":
            input_format = detect_format(fusions.columns)

        if input_format == "arriba":
            gene1_col = "#gene1" if "#gene1" in fusions.columns else "gene1"
            gene2_col = "gene2"
        elif input_format == "starfusion":
            gene1_col = "LeftGene"
            gene2_col = "RightGene"
        elif input_format == "fusioncatcher":
            gene1_col = "Gene_1_symbol(5end_fusion_partner)"
            gene2_col = "Gene_2_symbol(3end_fusion_partner)"
        else:
            print("Error: Unable to detect fusion format. Use --format arriba or starfusion.")
            sys.exit(1)

        if gene1_col not in fusions.columns or gene2_col not in fusions.columns:
            print(
                f"Warning: Expected columns '{gene1_col}' and '{gene2_col}' not found in {input_path}."
            )
            continue

        fusions_in_sample = set()
        for _, row in fusions.iterrows():
            gene1 = normalize_gene(row[gene1_col])
            gene2 = normalize_gene(row[gene2_col])
            if not gene1 or not gene2:
                continue
            fusion_key = "--".join(sorted([gene1, gene2]))
            fusions_in_sample.add(fusion_key)
            event_counts[fusion_key] = event_counts.get(fusion_key, 0) + 1

        for fusion_key in fusions_in_sample:
            sample_counts[fusion_key] = sample_counts.get(fusion_key, 0) + 1

    rows = []
    for fusion_key, sample_count in sample_counts.items():
        rows.append(
            {
                "fusion": fusion_key,
                "sample_count": sample_count,
                "event_count": event_counts.get(fusion_key, 0),
                "sample_fraction": sample_count / total_samples,
                "total_samples": total_samples,
            }
        )

    if not rows:
        print("No fusions found across inputs.")
        pd.DataFrame(
            columns=[
                "fusion",
                "sample_count",
                "event_count",
                "sample_fraction",
                "total_samples",
                "fusion_literature_summary",
                "gene1_function",
                "gene2_function",
            ]
        ).to_csv(args.output, sep="\t", index=False)
        return

    df = pd.DataFrame(rows)
    df = df.sort_values(
        by=["sample_count", "event_count", "fusion"], ascending=[False, False, True]
    )
    if args.annotate:
        sleep_s = max(0.0, args.sleep)
        gene_pairs = df["fusion"].str.split("--", n=1, expand=True)
        gene1_list = gene_pairs[0].fillna("").tolist()
        gene2_list = gene_pairs[1].fillna("").tolist()

        gene_cache = {}
        genes = sorted({g for g in gene1_list + gene2_list if g})
        for gene in genes:
            summary = query_mygene_summary(gene, sleep_s, args.max_gene_summary_chars)
            gene_cache[gene] = summary or "No summary found."

        fusion_cache = {}
        fusion_summaries = []
        gene1_functions = []
        gene2_functions = []

        for gene1, gene2 in zip(gene1_list, gene2_list):
            if gene1 and gene2:
                fusion_key = f"{gene1}--{gene2}"
                if fusion_key not in fusion_cache:
                    fusion_cache[fusion_key] = query_pubmed_summary(
                        gene1, gene2, args.email, sleep_s, args.max_abstract_chars
                    )
                fusion_summaries.append(fusion_cache[fusion_key])
            else:
                fusion_summaries.append("")

            gene1_functions.append(gene_cache.get(gene1, ""))
            gene2_functions.append(gene_cache.get(gene2, ""))

        df["fusion_literature_summary"] = fusion_summaries
        df["gene1_function"] = gene1_functions
        df["gene2_function"] = gene2_functions
    else:
        df["fusion_literature_summary"] = ""
        df["gene1_function"] = ""
        df["gene2_function"] = ""
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()

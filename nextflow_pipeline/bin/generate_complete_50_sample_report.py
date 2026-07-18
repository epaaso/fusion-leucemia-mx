#!/usr/bin/env python3
"""Generate a complete 50-sample HTML report integrating Arriba, STAR-Fusion, FusionCatcher, and ALLSorts."""

import argparse
import html
import re
import sys
from pathlib import Path
import pandas as pd

RESULT_EXCLUDE_DIRS = {"summary", "multiqc", "allsorts"}
KNOWN_DRIVER_GROUPS = [
    "BCR--ABL1",
    "ABL-class",
    "TCF3--PBX1",
    "KMT2A",
    "DUX4/IGH",
    "ZNF384",
    "ETV6--RUNX1",
    "PAX5",
    "MEF2D",
    "NUTM1",
    "HLF",
    "STIL--TAL1",
    "CRLF2/P2RY8",
]

CONFIDENCE_RANK = {"high": 0, "medium": 1, "low": 2}


def esc(value) -> str:
    if pd.isna(value):
        return ""
    return html.escape(str(value), quote=True)


def format_float(value, digits=3) -> str:
    if value is None or pd.isna(value):
        return "-"
    return f"{float(value):.{digits}f}"


def split_gene_tokens(value) -> list[str]:
    if pd.isna(value):
        return []
    tokens = []
    for part in re.split(r"[,|;]", str(value)):
        token = re.sub(r"\([^)]*\)", "", part).strip()
        token = token.split("^")[0].strip()
        if token and token != ".":
            tokens.append(token)
    return tokens


def primary_gene(value) -> str:
    tokens = split_gene_tokens(value)
    return tokens[0] if tokens else ""


def normalize_fusion(gene1: str, gene2: str) -> str:
    g1 = primary_gene(gene1)
    g2 = primary_gene(gene2)
    genes = sorted([g1, g2])
    return "--".join(g for g in genes if g)


def classify_driver(gene1: str, gene2: str) -> tuple[str, str] | None:
    g1 = primary_gene(gene1)
    g2 = primary_gene(gene2)
    genes = {g1, g2}
    gene_blob = " ".join(sorted(genes)).upper()

    def has(*items: str) -> bool:
        return all(item in genes for item in items)

    if has("BCR", "ABL1"):
        return "BCR--ABL1", "BCR--ABL1"
    if {"ABL1", "ABL2", "PDGFRA", "PDGFRB", "CSF1R", "JAK2", "EPOR"} & genes:
        return "ABL-class", f"{g1}--{g2}"
    if has("TCF3", "PBX1"):
        return "TCF3--PBX1", "TCF3--PBX1"
    if "KMT2A" in genes:
        return "KMT2A", f"{g1}--{g2}"
    if "DUX4" in gene_blob and re.search(r"\bIGH|IGH[DJM]", gene_blob):
        return "DUX4/IGH", f"{g1}--{g2}"
    if "ZNF384" in genes:
        return "ZNF384", f"{g1}--{g2}"
    if has("ETV6", "RUNX1"):
        return "ETV6--RUNX1", "ETV6--RUNX1"
    if "PAX5" in genes:
        return "PAX5", f"{g1}--{g2}"
    if "MEF2D" in genes:
        return "MEF2D", f"{g1}--{g2}"
    if "NUTM1" in genes:
        return "NUTM1", f"{g1}--{g2}"
    if "HLF" in genes:
        return "HLF", f"{g1}--{g2}"
    if has("STIL", "TAL1"):
        return "STIL--TAL1", "STIL--TAL1"
    if {"CRLF2", "P2RY8"} & genes:
        return "CRLF2/P2RY8", f"{g1}--{g2}"
    return None


KNOWN_DRIVER_GENES = {
    "BCR", "ABL1", "ABL2", "PDGFRA", "PDGFRB", "CSF1R", "JAK2", "EPOR",
    "TCF3", "PBX1", "KMT2A", "DUX4", "ZNF384", "ETV6", "RUNX1", "PAX5",
    "MEF2D", "NUTM1", "HLF", "STIL", "TAL1", "CRLF2", "P2RY8"
}


def is_pseudogene(gene: str) -> bool:
    gene = gene.upper()
    if gene in KNOWN_DRIVER_GENES:
        return False
    if re.search(r'P\d+$|P$|-PS\d*$', gene):
        return True
    return False


def build_allsorts_audit(
    sample_ids: set[str],
    current_path: Path,
    historical_path: Path,
) -> pd.DataFrame:
    rows = []

    # Current probabilities
    if current_path.exists():
        try:
            df = pd.read_csv(current_path)
            for _, r in df.iterrows():
                sid = r.get("sample_id")
                st = r.get("Pred")
                prob = r.get(st) if pd.notna(st) else 0.0
                flag = r.get("B-ALL", True)
                if pd.notna(sid):
                    rows.append({
                        "sample_id": str(sid),
                        "subtype": str(st) if pd.notna(st) else "Unclassified",
                        "main_probability": float(prob) if pd.notna(prob) else 0.0,
                        "pred_flag": bool(flag),
                        "allsorts_source": "current_probabilities_csv",
                    })
        except Exception as e:
            print(f"Warning parsing {current_path}: {e}")

    # Historical probabilities from allsorts_50_sample_audit.tsv
    audit_path = Path("nextflow_pipeline/results/summary/allsorts_50_sample_audit.tsv")
    if audit_path.exists():
        try:
            df_audit = pd.read_csv(audit_path, sep="\t")
            for _, r in df_audit.iterrows():
                sid = r.get("sample_id")
                st = r.get("subtype")
                prob = r.get("main_probability")
                flag = r.get("pred_flag", True)
                if pd.notna(sid):
                    rows.append({
                        "sample_id": str(sid),
                        "subtype": str(st) if pd.notna(st) else "Unclassified",
                        "main_probability": float(prob) if pd.notna(prob) else 0.0,
                        "pred_flag": bool(flag),
                        "allsorts_source": "allsorts_50_sample_audit",
                    })
        except Exception as e:
            print(f"Warning parsing {audit_path}: {e}")

    # Historical probabilities from HTML
    if historical_path.exists():
        try:
            content = historical_path.read_text(encoding="utf-8")
            # Parse HTML table rows for ALLSorts
            # Sample row format: <td>103-MO-15</td><td><span class="badge badge-tcf3">TCF3-PBX1</span></td><td>1.000</td>
            matches = re.findall(
                r"<tr>\s*<td>([a-zA-Z0-9_-]+)</td>\s*<td><span class=\"badge [^\"]+\">([^<]+)</span></td>\s*<td>([\d.]+)</td>",
                content,
            )
            for sid, subtype, prob in matches:
                rows.append({
                    "sample_id": sid,
                    "subtype": subtype,
                    "main_probability": float(prob),
                    "pred_flag": True,
                    "allsorts_source": "historical_html",
                })
        except Exception as e:
            print(f"Warning parsing {historical_path}: {e}")

    combined = pd.DataFrame(rows)
    if combined.empty:
        # Create empty placeholder if no files found
        combined = pd.DataFrame(columns=["sample_id", "subtype", "main_probability", "pred_flag", "allsorts_source"])

    source_rank = {"current_probabilities_csv": 0, "allsorts_50_sample_audit": 1, "historical_html": 2}
    combined["source_rank"] = combined["allsorts_source"].map(source_rank).fillna(9)
    combined = combined.sort_values(["sample_id", "source_rank"]).drop_duplicates(
        "sample_id", keep="first"
    )

    # Ensure all sample_ids have a row, filling unclassified if missing
    for sid in sample_ids:
        if sid not in combined["sample_id"].values:
            combined = pd.concat(
                [combined, pd.DataFrame([{
                    "sample_id": sid,
                    "subtype": "Unclassified",
                    "main_probability": 0.0,
                    "pred_flag": False,
                    "allsorts_source": "missing",
                    "source_rank": 9,
                }])],
                ignore_index=True,
            )

    combined = combined[combined["sample_id"].isin(sample_ids)].copy()
    return combined.sort_values("sample_id").drop(columns=["source_rank"])


def badge_class(subtype: str) -> str:
    text = subtype.lower()
    if "ph-like" in text:
        return "badge-phlike"
    if text == "ph" or "ph group" in text:
        return "badge-ph"
    if "tcf3" in text:
        return "badge-tcf3"
    if "kmt2a" in text:
        return "badge-kmt2a"
    if "dux4" in text:
        return "badge-dux4"
    if "znf" in text:
        return "badge-znf"
    if "unclassified" in text:
        return "badge-unclass"
    return "badge-driver"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="nextflow_pipeline/results")
    parser.add_argument(
        "--output-html",
        default="nextflow_pipeline/results/reporte_fusiones_completo_50_muestras.html",
    )
    parser.add_argument(
        "--allsorts-current",
        default="nextflow_pipeline/results/summary/allsorts/probabilities.csv",
    )
    parser.add_argument(
        "--cohort",
        choices=["all", "50", "minerva"],
        default="all",
        help="Cohort to report: all (72), 50 (leukemia), or minerva (22 new)"
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_html = Path(args.output_html)

    # Find all sample directories
    sample_dirs = sorted([
        d for d in results_dir.iterdir()
        if d.is_dir() and d.name not in RESULT_EXCLUDE_DIRS
    ])

    if args.cohort == "50":
        sample_dirs = [d for d in sample_dirs if not d.name.startswith("H")]
    elif args.cohort == "minerva":
        sample_dirs = [d for d in sample_dirs if d.name.startswith("H")]

    sample_ids = {d.name for d in sample_dirs}
    total_samples = len(sample_ids)

    print(f"Analyzing {total_samples} samples...")

    # Load ALLSorts
    allsorts = build_allsorts_audit(
        sample_ids,
        Path(args.allsorts_current),
        Path(args.results_dir) / "reporte_fusiones_completo.html",
    )
    allsorts_dict = allsorts.set_index("sample_id").to_dict(orient="index")

    # Load fusions across all 3 tools
    cohort_fusions = []
    
    # Track statistics per tool
    arriba_call_counts = {}
    starfusion_call_counts = {}
    fusioncatcher_call_counts = {}
    sample_star_qc = {}

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        fusions_by_normalized = {} # normalized_fusion -> fusion_entry

        # STAR Alignment QC Parser
        star_qc = {"reads": "-", "map_rate": "-", "chimeric": "-"}
        log_files = list(sample_dir.glob("*.Log.final.out"))
        if log_files and log_files[0].exists():
            try:
                content = log_files[0].read_text(encoding="utf-8")
                reads_match = re.search(r"Number of input reads\s*\|\s*(\d+)", content)
                map_match = re.search(r"Uniquely mapped reads %\s*\|\s*([\d.]+)%", content)
                chim_match = re.search(r"% of chimeric reads\s*\|\s*([\d.]+)%", content)
                if reads_match:
                    star_qc["reads"] = f"{int(reads_match.group(1)) / 1_000_000:.1f} M"
                if map_match:
                    star_qc["map_rate"] = f"{map_match.group(1)}%"
                if chim_match:
                    star_qc["chimeric"] = f"{chim_match.group(1)}%"
            except Exception as e:
                print(f"Error parsing STAR log for {sample_id}: {e}")
        sample_star_qc[sample_id] = star_qc

        # --- 1. ARRIBA ---
        arriba_file = sample_dir / "arriba_fusions.tsv"
        arriba_call_counts[sample_id] = 0
        if arriba_file.exists():
            try:
                df = pd.read_csv(arriba_file, sep="\t")
                arriba_call_counts[sample_id] = len(df)
                g1_col = "#gene1" if "#gene1" in df.columns else "gene1"
                for _, row in df.iterrows():
                    g1, g2 = row[g1_col], row["gene2"]
                    norm = normalize_fusion(g1, g2)
                    if not norm:
                        continue
                    
                    sr1 = int(row.get("split_reads1", 0) or 0)
                    sr2 = int(row.get("split_reads2", 0) or 0)
                    dm = int(row.get("discordant_mates", 0) or 0)
                    conf = str(row.get("confidence", "low")).lower()
                    frame = str(row.get("reading_frame", "unknown")).lower()
                    ftype = str(row.get("type", "translocation")).lower()
                    ffilters = str(row.get("filters", "")).lower()

                    if norm not in fusions_by_normalized:
                        fusions_by_normalized[norm] = {
                            "sample_id": sample_id,
                            "normalized": norm,
                            "gene1": primary_gene(g1),
                            "gene2": primary_gene(g2),
                            "arriba_detected": True,
                            "arriba_support": f"{conf} ({sr1}+{sr2} SR, {dm} DM)",
                            "arriba_confidence": conf,
                            "arriba_filters": ffilters,
                            "starfusion_detected": False,
                            "starfusion_support": "-",
                            "fusioncatcher_detected": False,
                            "fusioncatcher_support": "-",
                            "reading_frame": frame,
                            "type": ftype,
                        }
                    else:
                        fusions_by_normalized[norm]["arriba_detected"] = True
                        fusions_by_normalized[norm]["arriba_support"] = f"{conf} ({sr1}+{sr2} SR, {dm} DM)"
                        fusions_by_normalized[norm]["arriba_confidence"] = conf
                        fusions_by_normalized[norm]["arriba_filters"] = ffilters
            except Exception as e:
                print(f"Error parsing Arriba for {sample_id}: {e}")

        # --- 2. STAR-FUSION ---
        sf_file = sample_dir / "star-fusion.fusion_predictions.tsv"
        starfusion_call_counts[sample_id] = 0
        if sf_file.exists():
            try:
                df = pd.read_csv(sf_file, sep="\t")
                starfusion_call_counts[sample_id] = len(df)
                for _, row in df.iterrows():
                    # Columns: #FusionName, JunctionReadCount, SpanningFragCount, LeftGene, RightGene
                    fname = str(row.get("#FusionName", ""))
                    g1, g2 = fname.split("--") if "--" in fname else (row.get("LeftGene", ""), row.get("RightGene", ""))
                    norm = normalize_fusion(g1, g2)
                    if not norm:
                        continue
                    
                    j_reads = int(row.get("JunctionReadCount", 0) or 0)
                    s_frags = int(row.get("SpanningFragCount", 0) or 0)

                    if norm not in fusions_by_normalized:
                        fusions_by_normalized[norm] = {
                            "sample_id": sample_id,
                            "normalized": norm,
                            "gene1": primary_gene(g1),
                            "gene2": primary_gene(g2),
                            "arriba_detected": False,
                            "arriba_support": "-",
                            "arriba_confidence": "",
                            "starfusion_detected": True,
                            "starfusion_support": f"{j_reads} J, {s_frags} S",
                            "fusioncatcher_detected": False,
                            "fusioncatcher_support": "-",
                            "reading_frame": "unknown",
                            "type": "unknown",
                        }
                    else:
                        fusions_by_normalized[norm]["starfusion_detected"] = True
                        fusions_by_normalized[norm]["starfusion_support"] = f"{j_reads} J, {s_frags} S"
            except Exception as e:
                print(f"Error parsing STAR-Fusion for {sample_id}: {e}")

        # --- 3. FUSIONCATCHER ---
        fc_file = sample_dir / "fusioncatcher_fusions.txt"
        fusioncatcher_call_counts[sample_id] = 0
        if fc_file.exists():
            try:
                df = pd.read_csv(fc_file, sep="\t")
                fusioncatcher_call_counts[sample_id] = len(df)
                for _, row in df.iterrows():
                    # Columns: Gene_1_symbol(5end_fusion_partner), Gene_2_symbol(3end_fusion_partner), Spanning_pairs_count, Spanning_unique_reads_count
                    g1 = row.get("Gene_1_symbol(5end_fusion_partner)")
                    g2 = row.get("Gene_2_symbol(3end_fusion_partner)")
                    norm = normalize_fusion(g1, g2)
                    if not norm:
                        continue

                    pairs = int(row.get("Spanning_pairs_count", 0) or 0)
                    uniq = int(row.get("Spanning_unique_reads_count", 0) or 0)

                    if norm not in fusions_by_normalized:
                        fusions_by_normalized[norm] = {
                            "sample_id": sample_id,
                            "normalized": norm,
                            "gene1": primary_gene(g1),
                            "gene2": primary_gene(g2),
                            "arriba_detected": False,
                            "arriba_support": "-",
                            "arriba_confidence": "",
                            "starfusion_detected": False,
                            "starfusion_support": "-",
                            "fusioncatcher_detected": True,
                            "fusioncatcher_support": f"{uniq}+{pairs} reads",
                            "reading_frame": "unknown",
                            "type": "unknown",
                        }
                    else:
                        fusions_by_normalized[norm]["fusioncatcher_detected"] = True
                        fusions_by_normalized[norm]["fusioncatcher_support"] = f"{uniq}+{pairs} reads"
            except Exception as e:
                print(f"Error parsing FusionCatcher for {sample_id}: {e}")

        for entry in fusions_by_normalized.values():
            # Check concordance score
            score = 0
            if entry["arriba_detected"]: score += 1
            if entry["starfusion_detected"]: score += 1
            if entry["fusioncatcher_detected"]: score += 1
            entry["concordance"] = score

            # Driver classification
            driver = classify_driver(entry["gene1"], entry["gene2"])
            if driver:
                entry["driver_group"] = driver[0]
                entry["driver_display"] = driver[1]
            else:
                is_high_conf = entry.get("arriba_confidence") == "high"
                ftype = str(entry.get("type", "")).lower()
                ffilters = str(entry.get("arriba_filters", "")).lower()
                is_read_through = "read-through" in ftype or "read_through" in ftype or "read-through" in ffilters or "read_through" in ffilters
                is_pseudo = is_pseudogene(entry["gene1"]) or is_pseudogene(entry["gene2"])
                
                if is_high_conf and not is_read_through and not is_pseudo:
                    entry["driver_group"] = "Non-Driver (High Confidence)"
                    entry["driver_display"] = f"{entry['gene1']}--{entry['gene2']} (High Confidence)"
                else:
                    entry["driver_group"] = ""
                    entry["driver_display"] = ""

            cohort_fusions.append(entry)

    df_fusions = pd.DataFrame(cohort_fusions)

    # If empty, create placeholder
    if df_fusions.empty:
        df_fusions = pd.DataFrame(columns=[
            "sample_id", "normalized", "gene1", "gene2", "arriba_detected", "arriba_support", "arriba_confidence",
            "starfusion_detected", "starfusion_support", "fusioncatcher_detected", "fusioncatcher_support",
            "reading_frame", "type", "concordance", "driver_group", "driver_display"
        ])

    # Extract driver fusions only
    driver_fusions = df_fusions[df_fusions["driver_group"] != ""].copy()

    # Concordance table
    concordance_summary = []
    if not driver_fusions.empty:
        for fusion_key, grp in driver_fusions.groupby("normalized"):
            samples = sorted(list(grp["sample_id"].unique()))
            # check what proportion of runs detected this fusion
            arriba_cnt = grp["arriba_detected"].sum()
            sf_cnt = grp["starfusion_detected"].sum()
            fc_cnt = grp["fusioncatcher_detected"].sum()
            total_events = len(grp)

            concordance_summary.append({
                "fusion": fusion_key,
                "samples": ", ".join(samples),
                "arriba_yes": arriba_cnt > 0,
                "starfusion_yes": sf_cnt > 0,
                "fusioncatcher_yes": fc_cnt > 0,
                "concordance_label": f"{int(arriba_cnt > 0) + int(sf_cnt > 0) + int(fc_cnt > 0)}/3",
                "sample_count": len(samples),
            })
    df_concordance = pd.DataFrame(concordance_summary)

    # Let's count some summary stats
    num_samples_with_driver = len(driver_fusions["sample_id"].unique()) if not driver_fusions.empty else 0
    num_samples_phlike = sum(1 for sid in sample_ids if "ph-like" in allsorts_dict.get(sid, {}).get("subtype", "").lower())
    num_samples_ph = sum(1 for sid in sample_ids if allsorts_dict.get(sid, {}).get("subtype", "").lower() == "ph")

    # Generate HTML content
    html_template = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Reporte Integral Multi-Algoritmo &mdash; Cohorte LLA-B Pediatrica</title>
<style>
  :root {{
    --bg: #0f172a;
    --card: #1e293b;
    --accent: #3b82f6;
    --accent2: #8b5cf6;
    --danger: #ef4444;
    --warn: #f59e0b;
    --success: #10b981;
    --text: #f8fafc;
    --muted: #94a3b8;
    --border: #334155;
    --ph-bg: linear-gradient(135deg, #7c2d12 0%, #451a03 100%);
    --ph-border: #f97316;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: 'Segoe UI', system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); line-height: 1.6; }}
  .container {{ max-width: 1500px; margin: 0 auto; padding: 2rem; }}

  /* Header */
  .header {{
    background: linear-gradient(135deg, #1e3a8a 0%, #3b82f6 50%, #8b5cf6 100%);
    color: white;
    padding: 3rem 2rem;
    border-radius: 16px;
    margin-bottom: 2rem;
    text-align: center;
    box-shadow: 0 10px 40px rgba(59,130,246,0.3);
  }}
  .header h1 {{ font-size: 2.5rem; font-weight: 800; margin-bottom: 0.5rem; letter-spacing: -0.5px; }}
  .header .subtitle {{ font-size: 1.2rem; opacity: 0.9; }}
  .header .date {{ font-size: 0.95rem; opacity: 0.7; margin-top: 0.5rem; }}

  /* Navigation */
  .nav {{
    display: flex;
    gap: 0.5rem;
    flex-wrap: wrap;
    margin-bottom: 2rem;
    position: sticky;
    top: 0;
    background: var(--bg);
    padding: 1rem 0;
    z-index: 100;
    border-bottom: 1px solid var(--border);
  }}
  .nav a {{
    padding: 0.6rem 1.4rem;
    border-radius: 20px;
    text-decoration: none;
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--accent);
    background: var(--card);
    border: 1px solid var(--border);
    transition: all 0.2s;
  }}
  .nav a:hover {{
    background: var(--accent);
    color: white;
    transform: translateY(-1px);
  }}

  /* Dashboard Stats Grid */
  .summary-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 1rem;
    margin-bottom: 2rem;
  }}
  .summary-box {{
    background: var(--card);
    border-radius: 12px;
    padding: 1.5rem;
    text-align: center;
    border: 1px solid var(--border);
    box-shadow: 0 4px 12px rgba(0,0,0,0.15);
  }}
  .summary-box .number {{ font-size: 2.8rem; font-weight: 800; color: var(--accent); }}
  .summary-box .label {{ font-size: 0.85rem; color: var(--muted); font-weight: 600; text-transform: uppercase; letter-spacing: 0.5px; }}
  .summary-box.ph {{ border-color: var(--ph-border); background: var(--ph-bg); }}
  .summary-box.ph .number {{ color: #fdba74; }}

  /* Cards */
  .card {{
    background: var(--card);
    border-radius: 12px;
    padding: 2rem;
    margin-bottom: 1.5rem;
    border: 1px solid var(--border);
    box-shadow: 0 4px 20px rgba(0,0,0,0.1);
  }}
  .card h2 {{ font-size: 1.5rem; color: var(--accent); margin-bottom: 1.2rem; display: flex; align-items: center; gap: 0.6rem; }}
  .card h3 {{ font-size: 1.15rem; color: var(--text); margin: 1.5rem 0 0.8rem; }}

  /* Search Box */
  .search-container {{
    display: flex;
    gap: 1rem;
    margin-bottom: 1.5rem;
    flex-wrap: wrap;
  }}
  .search-input {{
    background: var(--bg);
    border: 1px solid var(--border);
    color: var(--text);
    padding: 0.6rem 1rem;
    border-radius: 8px;
    font-size: 0.9rem;
    flex: 1;
    min-width: 250px;
  }}
  .search-input:focus {{
    outline: none;
    border-color: var(--accent);
  }}

  /* Tables */
  table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; }}
  th {{
    background: #1e293b;
    font-weight: 700;
    text-transform: uppercase;
    font-size: 0.75rem;
    letter-spacing: 0.5px;
    color: var(--muted);
    border-bottom: 2px solid var(--border);
  }}
  th, td {{ padding: 0.8rem 1rem; text-align: left; border-bottom: 1px solid var(--border); }}
  tr:hover {{ background: #334155; }}
  .table-scroll {{ overflow-x: auto; max-height: 500px; overflow-y: auto; }}

  /* Badges */
  .badge {{ display: inline-block; padding: 0.25rem 0.7rem; border-radius: 12px; font-size: 0.75rem; font-weight: 700; text-align: center; }}
  .badge-high {{ background: #065f46; color: #34d399; }}
  .badge-medium {{ background: #78350f; color: #fbbf24; }}
  .badge-low {{ background: #7f1d1d; color: #fca5a5; }}
  .badge-ph {{ background: #7c2d12; color: #fdba74; border: 1px solid #f97316; }}
  .badge-phlike {{ background: #431407; color: #ff9d00; border: 1px solid #ff7b00; }}
  .badge-driver {{ background: #1e3a8a; color: #93c5fd; border: 1px solid #3b82f6; }}
  .badge-tcf3 {{ background: #312e81; color: #c7d2fe; }}
  .badge-kmt2a {{ background: #5c0632; color: #fbcfe8; }}
  .badge-dux4 {{ background: #064e3b; color: #a7f3d0; }}
  .badge-znf {{ background: #4c1d95; color: #ddd6fe; }}
  .badge-unclass {{ background: #334155; color: #cbd5e1; }}
  .badge-pass {{ background: #064e3b; color: #a7f3d0; }}
  .badge-warn {{ background: #78350f; color: #fde68a; }}
  .badge-fail {{ background: #7f1d1d; color: #fecaca; }}

  /* Ph-like highlight block */
  .ph-highlight {{
    background: linear-gradient(90deg, #382400, #1e293b);
    border-left: 4px solid var(--warn);
    padding: 1.2rem 1.8rem;
    border-radius: 0 8px 8px 0;
    margin: 1rem 0;
  }}
  .ph-highlight strong {{ color: #fbbf24; }}

  /* Footer */
  .footer {{ text-align: center; padding: 3rem; color: var(--muted); font-size: 0.85rem; border-top: 1px solid var(--border); margin-top: 3rem; }}

  /* Print */
  @media print {{ .nav {{ display: none; }} .card {{ break-inside: avoid; }} }}
</style>
</head>
<body>
<div class="container">

<!-- ============ HEADER ============ -->
<div class="header">
  <h1>Reporte Integral de Fusiones Genomicas Multi-Herramienta</h1>
  <div class="subtitle">Cohorte LLA-B Pediatrica &mdash; {total_samples} muestras</div>
  <div class="date">Pipeline: Nextflow (Arriba + STAR-Fusion + FusionCatcher + ALLSorts)</div>
</div>

<!-- ============ NAV ============ -->
<div class="nav">
  <a href="#resumen">Resumen</a>
  <a href="#allsorts">ALLSorts</a>
  <a href="#fusiones">Driver por Muestra</a>
  <a href="#concordancia">Concordancia</a>
  <a href="#calidad">Calidad &amp; Llamadas</a>
  <a href="#metodos">Metodos</a>
</div>

<!-- ============ RESUMEN ============ -->
<div class="card" id="resumen">
  <h2>&#x1F4CA; Resumen de la Cohorte</h2>
  <div class="summary-grid">
    <div class="summary-box">
      <div class="number">{total_samples}</div>
      <div class="label">Muestras Totales</div>
    </div>
    <div class="summary-box">
      <div class="number">{num_samples_with_driver}</div>
      <div class="label">Muestras con Driver Detectado</div>
    </div>
    <div class="summary-box ph">
      <div class="number">{num_samples_phlike}</div>
      <div class="label">Subtipo Ph-like (ALLSorts)</div>
    </div>
    <div class="summary-box">
      <div class="number">{num_samples_ph}</div>
      <div class="label">Subtipo Ph (BCR-ABL1)</div>
    </div>
  </div>
</div>

<!-- ============ ALLSORTS ============ -->
<div class="card" id="allsorts">
  <h2>&#x1F9E0; Clasificacion por Expresion (ALLSorts)</h2>
  <p style="color:var(--muted); margin-bottom:1rem;">Predicciones de subtipo basadas en firmas de expresion genica.</p>
  <div class="search-container">
    <input type="text" id="allsorts-search" class="search-input" placeholder="Buscar por muestra o subtipo..." onkeyup="filterTable('allsorts-table', 'allsorts-search')">
  </div>
  <div class="table-scroll">
    <table id="allsorts-table">
      <thead>
        <tr><th>Muestra</th><th>Subtipo Predicho</th><th>Probabilidad</th><th>Prediccion Confiable</th></tr>
      </thead>
      <tbody>
"""
    for _, row in allsorts.sort_values("sample_id").iterrows():
        sid = row["sample_id"]
        subtype = row["subtype"]
        prob = format_float(row["main_probability"])
        flag = "PASS" if row["pred_flag"] else "WARN"
        flag_class = "badge-pass" if row["pred_flag"] else "badge-warn"

        html_template += f"""        <tr>
          <td><strong>{esc(sid)}</strong></td>
          <td><span class="badge {badge_class(subtype)}">{esc(subtype)}</span></td>
          <td>{prob}</td>
          <td><span class="badge {flag_class}">{flag}</span></td>
        </tr>
"""

    html_template += f"""      </tbody>
    </table>
  </div>
</div>

<!-- ============ FUSIONES POR MUESTRA ============ -->
<div class="card" id="fusiones">
  <h2>&#x1F52C; Fusiones Driver por Muestra</h2>
  <p style="color:var(--muted); margin-bottom:1rem;">Detecciones de fusiones prioritarias con detalles de soporte por Arriba, STAR-Fusion y FusionCatcher.</p>
  <div class="search-container">
    <input type="text" id="fusiones-search" class="search-input" placeholder="Buscar por muestra, fusion o subtipo..." onkeyup="filterTable('fusiones-table', 'fusiones-search')">
  </div>
  <div class="table-scroll">
    <table id="fusiones-table">
      <thead>
        <tr>
          <th>Muestra</th>
          <th>Subtipo</th>
          <th>Fusion Driver</th>
          <th>Arriba</th>
          <th>STAR-Fusion</th>
          <th>FusionCatcher</th>
          <th>Concordancia</th>
          <th>Frame</th>
        </tr>
      </thead>
      <tbody>
"""
    # Sort by sample_id, then driver group
    if not driver_fusions.empty:
        ordered_drivers = driver_fusions.sort_values(["sample_id", "driver_group"])
        for _, row in ordered_drivers.iterrows():
            sid = row["sample_id"]
            sub = allsorts_dict.get(sid, {}).get("subtype", "Unclassified")
            fusion = row["gene1"] + "&ndash;" + row["gene2"]
            
            arriba_det = row["arriba_support"] if row["arriba_detected"] else '<span style="color:var(--danger)">&#10007;</span>'
            sf_det = row["starfusion_support"] if row["starfusion_detected"] else '<span style="color:var(--danger)">&#10007;</span>'
            fc_det = row["fusioncatcher_support"] if row["fusioncatcher_detected"] else '<span style="color:var(--danger)">&#10007;</span>'
            
            conc = f"{row['concordance']}/3"
            conc_class = "badge-high" if row["concordance"] == 3 else "badge-medium" if row["concordance"] == 2 else "badge-low"
            frame = row["reading_frame"]
            
            # Row styling based on subtype
            bg_style = ""
            if "tcf3" in sub.lower(): bg_style = ' style="background:rgba(99,102,241,0.08);"'
            elif "ph" == sub.lower() or "bcr" in fusion.lower(): bg_style = ' style="background:rgba(239,68,68,0.08);"'
            elif "kmt2a" in sub.lower(): bg_style = ' style="background:rgba(236,72,153,0.08);"'
            elif "dux4" in sub.lower(): bg_style = ' style="background:rgba(16,185,129,0.08);"'
            elif "ph-like" in sub.lower(): bg_style = ' style="background:rgba(245,158,11,0.08);"'

            html_template += f"""        <tr{bg_style}>
          <td><strong>{esc(sid)}</strong></td>
          <td><span class="badge {badge_class(sub)}">{esc(sub)}</span></td>
          <td><strong>{fusion}</strong></td>
          <td>{arriba_det}</td>
          <td>{sf_det}</td>
          <td>{fc_det}</td>
          <td><span class="badge {conc_class}">{conc}</span></td>
          <td>{esc(frame)}</td>
        </tr>
"""
    else:
        html_template += """        <tr><td colspan="8" style="text-align:center;">No se detectaron fusiones driver en la cohorte.</td></tr>"""

    html_template += f"""      </tbody>
    </table>
  </div>
</div>

<!-- ============ CONCORDANCIA ============ -->
<div class="card" id="concordancia">
  <h2>&#x1F91D; Concordancia de Llamadores para Drivers</h2>
  <p style="color:var(--muted); margin-bottom:1rem;">Resumen de deteccion cruzada entre herramientas para fusions driver identificadas.</p>
  <div class="table-scroll">
    <table>
      <thead>
        <tr><th>Fusion Driver</th><th>Muestras</th><th>Arriba</th><th>STAR-Fusion</th><th>FusionCatcher</th><th>Concordancia</th></tr>
      </thead>
      <tbody>
"""
    if not df_concordance.empty:
        for _, row in df_concordance.sort_values(by="sample_count", ascending=False).iterrows():
            fus = row["fusion"]
            samples = row["samples"]
            arr = '<span style="color:var(--success);font-weight:bold;">&#10003;</span>' if row["arriba_yes"] else '<span style="color:var(--danger)">&#10007;</span>'
            sf = '<span style="color:var(--success);font-weight:bold;">&#10003;</span>' if row["starfusion_yes"] else '<span style="color:var(--danger)">&#10007;</span>'
            fc = '<span style="color:var(--success);font-weight:bold;">&#10003;</span>' if row["fusioncatcher_yes"] else '<span style="color:var(--danger)">&#10007;</span>'
            
            c_label = row["concordance_label"]
            c_class = "badge-high" if "3" in c_label else "badge-medium" if "2" in c_label else "badge-low"

            html_template += f"""        <tr>
          <td><strong>{esc(fus)}</strong></td>
          <td>{esc(samples)}</td>
          <td style="text-align:center;">{arr}</td>
          <td style="text-align:center;">{sf}</td>
          <td style="text-align:center;">{fc}</td>
          <td><span class="badge {c_class}">{c_label}</span></td>
        </tr>
"""
    else:
        html_template += """        <tr><td colspan="6" style="text-align:center;">Sin datos de concordancia disponibles.</td></tr>"""

    html_template += f"""      </tbody>
    </table>
  </div>
</div>

<!-- ============ CALIDAD & LLAMADAS ============ -->
<div class="card" id="calidad">
  <h2>&#x2705; Metricas de Calidad de Secuenciacion y Mapeo</h2>
  <p style="color:var(--muted); margin-bottom:1rem;">Medidas de calidad de secuenciacion y mapeo (STAR QC) junto con el numero de llamadas totales por herramienta.</p>
  <div class="table-scroll">
    <table>
      <thead>
        <tr><th>Muestra</th><th>Reads Secuenciadas</th><th>Mapeo Unico (STAR)</th><th>% Quimericas (STAR)</th><th>Llamadas Arriba</th><th>Llamadas STAR-Fusion</th><th>Llamadas FusionCatcher</th></tr>
      </thead>
      <tbody>
"""
    for sid in sorted(list(sample_ids)):
        arr_cnt = arriba_call_counts.get(sid, 0)
        sf_cnt = starfusion_call_counts.get(sid, 0)
        fc_cnt = fusioncatcher_call_counts.get(sid, 0)
        
        qc = sample_star_qc.get(sid, {"reads": "-", "map_rate": "-", "chimeric": "-"})
        reads_str = qc["reads"]
        map_str = qc["map_rate"]
        chim_str = qc["chimeric"]
        
        # Color coding for reads count
        reads_display = reads_str
        if reads_str != "-":
            try:
                reads_val = float(reads_str.replace(" M", ""))
                if reads_val < 15.0:
                    reads_display = f"<span class='badge badge-low'>{reads_str}</span>"
                elif reads_val < 30.0:
                    reads_display = f"<span class='badge badge-medium'>{reads_str}</span>"
                else:
                    reads_display = f"<span class='badge badge-pass'>{reads_str}</span>"
            except:
                pass
                
        # Color coding for mapping rate
        map_display = map_str
        if map_str != "-":
            try:
                map_val = float(map_str.replace("%", ""))
                if map_val < 70.0:
                    map_display = f"<span class='badge badge-low'>{map_str}</span>"
                elif map_val < 80.0:
                    map_display = f"<span class='badge badge-medium'>{map_str}</span>"
                else:
                    map_display = f"<span class='badge badge-pass'>{map_str}</span>"
            except:
                pass
                
        html_template += f"""        <tr>
          <td><strong>{esc(sid)}</strong></td>
          <td>{reads_display}</td>
          <td>{map_display}</td>
          <td>{chim_str}</td>
          <td>{arr_cnt}</td>
          <td>{sf_cnt}</td>
          <td>{fc_cnt}</td>
        </tr>
"""

    html_template += f"""      </tbody>
    </table>
  </div>
</div>

<!-- ============ METODOS ============ -->
<div class="card" id="metodos">
  <h2>&#x1F4D6; Metodos y Parametros</h2>
  <h3>Ambiente de Analisis</h3>
  <p>Todos los alineamientos y llamadas se realizaron utilizando el genoma de referencia GRCh38 y anotaciones de GENCODE v50. El pipeline de Nextflow se ejecuto con integracion de Conda.</p>
  
  <h3>Criterios de Deteccion</h3>
  <ul>
    <li><strong>Arriba:</strong> Streaming de STAR, filtrado de ruido con blacklist de hg38.</li>
    <li><strong>STAR-Fusion:</strong> Procesamiento en base a libreria genómica CTAT.</li>
    <li><strong>FusionCatcher:</strong> Filtrado exhaustivo y mapeo multi-aligner (Bowtie/Bowtie2/STAR).</li>
  </ul>
</div>

<div class="footer">
  Reporte Generado Automaticamente el {pd.Timestamp.now().strftime('%d de %B de %Y')} &bull; Cohorte LLA-B 50 muestras.
</div>

</div>

<script>
function filterTable(tableId, inputId) {{
  var input, filter, table, tr, td, i, txtValue;
  input = document.getElementById(inputId);
  filter = input.value.toUpperCase();
  table = document.getElementById(tableId);
  tr = table.getElementsByTagName("tr");
  
  for (i = 1; i < tr.length; i++) {{
    tr[i].style.display = "none";
    td = tr[i].getElementsByTagName("td");
    for (var j = 0; j < td.length; j++) {{
      if (td[j]) {{
        txtValue = td[j].textContent || td[j].innerText;
        if (txtValue.toUpperCase().indexOf(filter) > -1) {{
          tr[i].style.display = "";
          break;
        }}
      }}
    }}
  }}
}}
</script>

</body>
</html>
"""
    output_html.write_text(html_template, encoding="utf-8")
    print(f"Unified report successfully written to {output_html}")


if __name__ == "__main__":
    main()

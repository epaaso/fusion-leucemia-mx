#!/usr/bin/env python3
"""Generate a 50-sample Arriba + ALLSorts HTML report."""

from __future__ import annotations

import argparse
import html
import re
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


def parse_float(value) -> float | None:
    if pd.isna(value):
        return None
    match = re.search(r"\d+(?:\.\d+)?", str(value))
    if not match:
        return None
    return float(match.group(0))


def format_float(value, digits=3) -> str:
    if value is None or pd.isna(value):
        return "-"
    return f"{float(value):.{digits}f}"


def truthy(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


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


def normalize_fusion(row) -> str:
    genes = sorted([primary_gene(row["gene1"]), primary_gene(row["gene2"])])
    return "--".join(g for g in genes if g)


def display_fusion(row) -> str:
    return f"{primary_gene(row['gene1'])}--{primary_gene(row['gene2'])}"


def all_genes(row) -> set[str]:
    return set(split_gene_tokens(row["gene1"]) + split_gene_tokens(row["gene2"]))


def support_text(row) -> str:
    sr1 = int(row.get("split_reads1", 0) or 0)
    sr2 = int(row.get("split_reads2", 0) or 0)
    dm = int(row.get("discordant_mates", 0) or 0)
    return f"{sr1}+{sr2} SR, {dm} DM"


def classify_driver(row) -> tuple[str, str] | None:
    genes = all_genes(row)
    gene_blob = " ".join(sorted(genes)).upper()

    def has(*items: str) -> bool:
        return all(item in genes for item in items)

    if has("BCR", "ABL1"):
        return "BCR--ABL1", "BCR--ABL1"
    if {"ABL1", "ABL2", "PDGFRA", "PDGFRB", "CSF1R", "JAK2", "EPOR"} & genes:
        return "ABL-class", display_fusion(row)
    if has("TCF3", "PBX1"):
        return "TCF3--PBX1", "TCF3--PBX1"
    if "KMT2A" in genes:
        return "KMT2A", display_fusion(row)
    if "DUX4" in gene_blob and re.search(r"\bIGH|IGH[DJM]", gene_blob):
        return "DUX4/IGH", display_fusion(row)
    if "ZNF384" in genes:
        return "ZNF384", display_fusion(row)
    if has("ETV6", "RUNX1"):
        return "ETV6--RUNX1", "ETV6--RUNX1"
    if "PAX5" in genes:
        return "PAX5", display_fusion(row)
    if "MEF2D" in genes:
        return "MEF2D", display_fusion(row)
    if "NUTM1" in genes:
        return "NUTM1", display_fusion(row)
    if "HLF" in genes:
        return "HLF", display_fusion(row)
    if has("STIL", "TAL1"):
        return "STIL--TAL1", "STIL--TAL1"
    if {"CRLF2", "P2RY8"} & genes:
        return "CRLF2/P2RY8", display_fusion(row)
    return None


def find_arriba_files(results_dir: Path) -> list[Path]:
    files = []
    for sample_dir in sorted(p for p in results_dir.iterdir() if p.is_dir()):
        if sample_dir.name in RESULT_EXCLUDE_DIRS:
            continue
        fusions = sample_dir / "arriba_fusions.tsv"
        if fusions.exists():
            files.append(fusions)
    return files


def load_arriba(results_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    files = find_arriba_files(results_dir)
    if len(files) != 50:
        raise SystemExit(f"Expected 50 Arriba fusion files, found {len(files)}")

    rows = []
    sample_rows = []
    for path in files:
        sample = path.parent.name
        df = pd.read_csv(path, sep="\t")
        gene1_col = "#gene1" if "#gene1" in df.columns else "gene1"
        df = df.rename(columns={gene1_col: "gene1"})
        df["sample_id"] = sample
        rows.append(df)
        sample_rows.append(
            {
                "sample_id": sample,
                "arriba_calls": len(df),
                "high_confidence": int((df["confidence"] == "high").sum()),
                "medium_confidence": int((df["confidence"] == "medium").sum()),
            }
        )

    fusions = pd.concat(rows, ignore_index=True)
    fusions["fusion_key"] = fusions.apply(normalize_fusion, axis=1)
    fusions["fusion_display"] = fusions.apply(display_fusion, axis=1)
    fusions["driver_classification"] = fusions.apply(classify_driver, axis=1)
    fusions["driver_group"] = fusions["driver_classification"].apply(
        lambda value: value[0] if value else ""
    )
    fusions["driver_fusion"] = fusions["driver_classification"].apply(
        lambda value: value[1] if value else ""
    )

    sample_qc = pd.DataFrame(sample_rows).sort_values("sample_id")
    driver_counts = (
        fusions[fusions["driver_group"] != ""]
        .groupby("sample_id")
        .size()
        .rename("driver_like_calls")
    )
    sample_qc = sample_qc.merge(driver_counts, on="sample_id", how="left")
    sample_qc["driver_like_calls"] = sample_qc["driver_like_calls"].fillna(0).astype(int)

    return fusions, sample_qc, pd.DataFrame({"sample_id": sorted(sample_qc["sample_id"])})


def build_common_fusions(fusions: pd.DataFrame, total_samples: int) -> pd.DataFrame:
    sample_counts = (
        fusions.dropna(subset=["fusion_key"])
        .groupby("fusion_key")["sample_id"]
        .nunique()
        .rename("sample_count")
    )
    event_counts = fusions.groupby("fusion_key").size().rename("event_count")
    common = pd.concat([sample_counts, event_counts], axis=1).reset_index()
    common = common.rename(columns={"fusion_key": "fusion"})
    common["sample_fraction"] = common["sample_count"] / total_samples
    common["total_samples"] = total_samples
    return common.sort_values(
        ["sample_count", "event_count", "fusion"], ascending=[False, False, True]
    )


def load_current_allsorts(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    subtype_cols = [
        c for c in df.columns if c not in {"sample_id", "B-ALL", "Pred"}
    ]
    rows = []
    for _, row in df.iterrows():
        subtype = str(row.get("Pred", "")).strip()
        probs = {col: parse_float(row[col]) for col in subtype_cols}
        top_label = max(probs, key=lambda col: probs[col] if probs[col] is not None else -1)
        main_prob = probs.get(subtype)
        if main_prob is None:
            main_prob = probs[top_label]
        rows.append(
            {
                "sample_id": row["sample_id"],
                "subtype": subtype,
                "main_probability": main_prob,
                "main_probability_label": top_label if top_label != subtype else "",
                "ph_like_probability": probs.get("Ph-like"),
                "ph_probability": probs.get("Ph"),
                "pred_flag": truthy(row.get("B-ALL")),
                "allsorts_source": "current_probabilities_csv",
            }
        )
    return pd.DataFrame(rows)


def load_historical_allsorts(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        tables = pd.read_html(path)
    except Exception:
        return pd.DataFrame()
    if not tables:
        return pd.DataFrame()
    table = tables[0]
    required = {"Muestra", "Subtipo", "Prob. Principal", "Ph-like Prob.", "Ph Prob.", "Pred"}
    if not required.issubset(set(table.columns)):
        return pd.DataFrame()
    rows = []
    for _, row in table.iterrows():
        sample_id = str(row.get("Muestra", "")).strip()
        if not sample_id or sample_id == "nan" or not re.search(r"MO|LAL", sample_id):
            continue
        prob_text = str(row.get("Prob. Principal", ""))
        label_match = re.search(r"\(([^)]+)\)", prob_text)
        rows.append(
            {
                "sample_id": sample_id,
                "subtype": str(row.get("Subtipo", "")).strip(),
                "main_probability": parse_float(prob_text),
                "main_probability_label": label_match.group(1) if label_match else "",
                "ph_like_probability": parse_float(row.get("Ph-like Prob.")),
                "ph_probability": parse_float(row.get("Ph Prob.")),
                "pred_flag": truthy(row.get("Pred")),
                "allsorts_source": "historical_html",
            }
        )
    return pd.DataFrame(rows)


def build_allsorts_audit(
    sample_ids: set[str],
    current_path: Path,
    historical_path: Path,
) -> pd.DataFrame:
    current = load_current_allsorts(current_path)
    historical = load_historical_allsorts(historical_path)
    combined = pd.concat([current, historical], ignore_index=True)
    if combined.empty:
        raise SystemExit("No ALLSorts rows were found.")

    source_rank = {"current_probabilities_csv": 0, "historical_html": 1}
    combined["source_rank"] = combined["allsorts_source"].map(source_rank).fillna(9)
    combined = combined.sort_values(["sample_id", "source_rank"]).drop_duplicates(
        "sample_id", keep="first"
    )
    combined = combined[combined["sample_id"].isin(sample_ids)].copy()

    missing = sorted(sample_ids - set(combined["sample_id"]))
    extra = sorted(set(combined["sample_id"]) - sample_ids)
    if missing or extra or len(combined) != 50:
        raise SystemExit(
            "ALLSorts coverage validation failed. "
            f"rows={len(combined)}, missing={missing}, extra={extra}"
        )
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


def allsorts_driver_map(driver_fusions: pd.DataFrame) -> dict[str, str]:
    rows = []
    for sample, group in driver_fusions.groupby("sample_id"):
        entries = []
        for _, row in group.sort_values(
            by=["driver_group", "confidence"],
            key=lambda col: col.map(CONFIDENCE_RANK).fillna(9) if col.name == "confidence" else col,
        ).iterrows():
            label = row["driver_fusion"] or row["fusion_display"]
            if label not in entries:
                entries.append(label)
        rows.append((sample, "; ".join(entries[:4])))
    return dict(rows)


def render_allsorts_rows(allsorts: pd.DataFrame, driver_map: dict[str, str]) -> str:
    ordered = allsorts.sort_values(
        ["subtype", "sample_id"], ascending=[True, True]
    )
    parts = []
    for _, row in ordered.iterrows():
        subtype = row["subtype"]
        main = format_float(row["main_probability"])
        if row.get("main_probability_label"):
            main = f"{main} ({esc(row['main_probability_label'])})"
        pred_class = "badge-pass" if row["pred_flag"] else "badge-fail"
        parts.append(
            "<tr>"
            f"<td>{esc(row['sample_id'])}</td>"
            f"<td><span class=\"badge {badge_class(subtype)}\">{esc(subtype)}</span></td>"
            f"<td>{main}</td>"
            f"<td>{format_float(row['ph_like_probability'])}</td>"
            f"<td>{format_float(row['ph_probability'])}</td>"
            f"<td><span class=\"badge {pred_class}\">{esc(row['pred_flag'])}</span></td>"
            f"<td>{esc(driver_map.get(row['sample_id'], '-'))}</td>"
            f"<td>{esc(row['allsorts_source'])}</td>"
            "</tr>"
        )
    return "\n".join(parts)


def render_subtype_bar(allsorts: pd.DataFrame) -> str:
    colors = [
        "#ef4444", "#f97316", "#6366f1", "#22c55e", "#ec4899",
        "#8b5cf6", "#06b6d4", "#3b82f6", "#14b8a6", "#94a3b8",
    ]
    counts = allsorts["subtype"].value_counts()
    parts = []
    for i, (subtype, count) in enumerate(counts.items()):
        color = colors[i % len(colors)]
        text_color = " color:#333;" if color in {"#94a3b8"} else ""
        parts.append(
            f"<div style=\"flex:{count}; background:{color};{text_color}\" "
            f"title=\"{esc(subtype)}\">{esc(subtype)} ({count})</div>"
        )
    return "\n".join(parts)


def render_driver_rows(driver_fusions: pd.DataFrame, allsorts: pd.DataFrame) -> str:
    subtype_map = dict(zip(allsorts["sample_id"], allsorts["subtype"]))
    if driver_fusions.empty:
        return "<tr><td colspan=\"9\">No se detectaron fusiones driver por Arriba.</td></tr>"
    rows = []
    sortable = driver_fusions.copy()
    sortable["confidence_rank"] = sortable["confidence"].map(CONFIDENCE_RANK).fillna(9)
    sortable = sortable.sort_values(["driver_group", "sample_id", "confidence_rank"])
    for _, row in sortable.iterrows():
        subtype = subtype_map.get(row["sample_id"], "-")
        rows.append(
            "<tr>"
            f"<td>{esc(row['sample_id'])}</td>"
            f"<td><span class=\"badge {badge_class(subtype)}\">{esc(subtype)}</span></td>"
            f"<td>{esc(row['driver_group'])}</td>"
            f"<td><strong>{esc(row['driver_fusion'] or row['fusion_display'])}</strong></td>"
            f"<td>{esc(row.get('type', '-'))}</td>"
            f"<td><span class=\"badge badge-{esc(row.get('confidence', 'low'))}\">{esc(row.get('confidence', '-'))}</span></td>"
            f"<td>{esc(support_text(row))}</td>"
            f"<td>{esc(row.get('reading_frame', '-'))}</td>"
            f"<td>{esc(row.get('tags', '-'))}</td>"
            "</tr>"
        )
    return "\n".join(rows)


def render_phlike_rows(phlike: pd.DataFrame, driver_fusions: pd.DataFrame) -> str:
    rows = []
    driver_by_sample = {sample: group for sample, group in driver_fusions.groupby("sample_id")}
    for _, subtype_row in phlike.sort_values(
        ["ph_like_probability", "ph_probability", "sample_id"], ascending=[False, False, True]
    ).iterrows():
        sample = subtype_row["sample_id"]
        drivers = driver_by_sample.get(sample)
        if drivers is None or drivers.empty:
            rows.append(
                "<tr>"
                f"<td>{esc(sample)}</td>"
                f"<td><span class=\"badge {badge_class(subtype_row['subtype'])}\">{esc(subtype_row['subtype'])}</span></td>"
                f"<td>{format_float(subtype_row['ph_like_probability'])}</td>"
                f"<td>{format_float(subtype_row['ph_probability'])}</td>"
                "<td colspan=\"5\">Sin fusion driver detectada por Arriba en las categorias priorizadas.</td>"
                "</tr>"
            )
            continue
        for _, row in drivers.iterrows():
            rows.append(
                "<tr>"
                f"<td>{esc(sample)}</td>"
                f"<td><span class=\"badge {badge_class(subtype_row['subtype'])}\">{esc(subtype_row['subtype'])}</span></td>"
                f"<td>{format_float(subtype_row['ph_like_probability'])}</td>"
                f"<td>{format_float(subtype_row['ph_probability'])}</td>"
                f"<td><strong>{esc(row['driver_fusion'] or row['fusion_display'])}</strong></td>"
                f"<td>{esc(row.get('confidence', '-'))}</td>"
                f"<td>{esc(support_text(row))}</td>"
                f"<td>{esc(row.get('reading_frame', '-'))}</td>"
                f"<td>{esc(row.get('type', '-'))}</td>"
                "</tr>"
            )
    return "\n".join(rows)


def render_recurrent_rows(recurrent: pd.DataFrame) -> str:
    if recurrent.empty:
        return "<tr><td colspan=\"5\">No hay fusiones recurrentes no-driver con el umbral seleccionado.</td></tr>"
    rows = []
    for _, row in recurrent.head(30).iterrows():
        rows.append(
            "<tr>"
            f"<td>{esc(row['fusion'])}</td>"
            f"<td>{int(row['sample_count'])}</td>"
            f"<td>{int(row['event_count'])}</td>"
            f"<td>{row['sample_fraction']:.2f}</td>"
            "<td>Recurrente en cohorte; interpretar con cautela si corresponde a read-through, pseudogen o region repetitiva.</td>"
            "</tr>"
        )
    return "\n".join(rows)


def render_sample_qc_rows(sample_qc: pd.DataFrame) -> str:
    rows = []
    for _, row in sample_qc.iterrows():
        rows.append(
            "<tr>"
            f"<td>{esc(row['sample_id'])}</td>"
            f"<td>{int(row['arriba_calls'])}</td>"
            f"<td>{int(row['high_confidence'])}</td>"
            f"<td>{int(row['medium_confidence'])}</td>"
            f"<td>{int(row['driver_like_calls'])}</td>"
            "<td><span class=\"badge badge-pass\">OK</span></td>"
            "</tr>"
        )
    return "\n".join(rows)


def build_html(
    allsorts: pd.DataFrame,
    sample_qc: pd.DataFrame,
    common: pd.DataFrame,
    driver_fusions: pd.DataFrame,
    recurrent: pd.DataFrame,
    output_stem: str,
) -> str:
    driver_map = allsorts_driver_map(driver_fusions)
    phlike = allsorts[
        allsorts["subtype"].str.contains("Ph", case=False, na=False)
        | (allsorts["ph_like_probability"].fillna(0) >= 0.5)
        | (allsorts["ph_probability"].fillna(0) >= 0.5)
    ].copy()
    subtype_counts = allsorts["subtype"].value_counts()
    phlike_count = len(phlike)
    driver_sample_count = driver_fusions["sample_id"].nunique()
    high_calls = int(sample_qc["high_confidence"].sum())
    medium_calls = int(sample_qc["medium_confidence"].sum())
    total_calls = int(sample_qc["arriba_calls"].sum())

    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Reporte Arriba + ALLSorts - 50 muestras</title>
<style>
  :root {{
    --bg:#f8f9fa; --card:#ffffff; --accent:#2563eb; --accent2:#7c3aed;
    --danger:#dc2626; --warn:#f59e0b; --success:#16a34a;
    --text:#1e293b; --muted:#64748b; --border:#e2e8f0;
    --ph-bg:linear-gradient(135deg,#fef3c7 0%,#fde68a 100%);
    --ph-border:#f59e0b;
  }}
  * {{ box-sizing:border-box; margin:0; padding:0; }}
  body {{ font-family:'Segoe UI',system-ui,-apple-system,sans-serif; background:var(--bg); color:var(--text); line-height:1.6; }}
  .container {{ max-width:1400px; margin:0 auto; padding:2rem; }}
  .header {{ background:linear-gradient(135deg,#1e3a5f 0%,#2563eb 50%,#7c3aed 100%); color:white; padding:3rem 2rem; border-radius:16px; margin-bottom:2rem; text-align:center; box-shadow:0 10px 40px rgba(37,99,235,0.3); }}
  .header h1 {{ font-size:2.2rem; font-weight:800; margin-bottom:0.5rem; }}
  .header .subtitle {{ font-size:1.1rem; opacity:0.92; }}
  .header .date {{ font-size:0.9rem; opacity:0.78; margin-top:0.5rem; }}
  .nav {{ display:flex; gap:0.5rem; flex-wrap:wrap; margin-bottom:2rem; position:sticky; top:0; background:var(--bg); padding:1rem 0; z-index:100; border-bottom:1px solid var(--border); }}
  .nav a {{ padding:0.5rem 1.2rem; border-radius:20px; text-decoration:none; font-size:0.85rem; font-weight:600; color:var(--accent); background:white; border:1px solid var(--border); }}
  .card {{ background:var(--card); border-radius:12px; padding:1.5rem 2rem; margin-bottom:1.5rem; box-shadow:0 1px 3px rgba(0,0,0,0.08); border:1px solid var(--border); }}
  .card h2 {{ font-size:1.4rem; color:var(--accent); margin-bottom:1rem; }}
  .card h3 {{ font-size:1.1rem; color:var(--text); margin:1.2rem 0 0.6rem; }}
  .summary-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr)); gap:1rem; margin-bottom:2rem; }}
  .summary-box {{ background:white; border-radius:12px; padding:1.5rem; text-align:center; border:1px solid var(--border); box-shadow:0 2px 8px rgba(0,0,0,0.04); }}
  .summary-box .number {{ font-size:2.35rem; font-weight:800; }}
  .summary-box .label {{ font-size:0.82rem; color:var(--muted); font-weight:600; text-transform:uppercase; letter-spacing:0.5px; }}
  .summary-box.ph {{ border-color:var(--ph-border); background:var(--ph-bg); }}
  .summary-box.ph .number {{ color:#b45309; }}
  table {{ width:100%; border-collapse:collapse; font-size:0.84rem; }}
  th {{ background:#f1f5f9; font-weight:700; text-transform:uppercase; font-size:0.74rem; letter-spacing:0.4px; color:var(--muted); }}
  th,td {{ padding:0.58rem 0.75rem; text-align:left; border-bottom:1px solid var(--border); vertical-align:top; }}
  tr:hover {{ background:#f8fafc; }}
  .table-scroll {{ overflow-x:auto; }}
  .badge {{ display:inline-block; padding:0.2rem 0.6rem; border-radius:12px; font-size:0.75rem; font-weight:700; white-space:nowrap; }}
  .badge-high {{ background:#dcfce7; color:#166534; }}
  .badge-medium {{ background:#fef9c3; color:#854d0e; }}
  .badge-low {{ background:#fee2e2; color:#991b1b; }}
  .badge-ph {{ background:#fef3c7; color:#92400e; border:1px solid #f59e0b; }}
  .badge-phlike {{ background:#fff7ed; color:#c2410c; border:1px solid #fb923c; }}
  .badge-driver {{ background:#dbeafe; color:#1e40af; border:1px solid #93c5fd; }}
  .badge-tcf3 {{ background:#e0e7ff; color:#3730a3; }}
  .badge-kmt2a {{ background:#fce7f3; color:#9d174d; }}
  .badge-dux4 {{ background:#d1fae5; color:#065f46; }}
  .badge-znf {{ background:#ede9fe; color:#5b21b6; }}
  .badge-unclass {{ background:#f1f5f9; color:#475569; }}
  .badge-pass {{ background:#dcfce7; color:#166534; }}
  .badge-fail {{ background:#fee2e2; color:#991b1b; }}
  .subtype-bar {{ display:flex; border-radius:8px; overflow:hidden; min-height:34px; margin:1rem 0; }}
  .subtype-bar div {{ display:flex; align-items:center; justify-content:center; font-size:0.7rem; font-weight:700; color:white; min-width:34px; padding:0 0.35rem; }}
  .note {{ background:#eff6ff; border-left:4px solid var(--accent); padding:1rem 1.5rem; border-radius:0 8px 8px 0; margin:1rem 0; }}
  .ph-highlight {{ background:linear-gradient(90deg,#fffbeb,#fef3c7); border-left:4px solid #f59e0b; padding:1rem 1.5rem; border-radius:0 8px 8px 0; margin:1rem 0; }}
  .footer {{ text-align:center; padding:2rem; color:var(--muted); font-size:0.8rem; }}
  @media print {{ .nav {{ display:none; }} .card {{ break-inside:avoid; }} }}
</style>
</head>
<body>
<div class="container">
<div class="header">
  <h1>Reporte Integral de Fusiones Genomicas</h1>
  <div class="subtitle">Cohorte LLA-B Pediatrica - 50 muestras</div>
  <div class="date">Pipeline: Nextflow (Arriba + ALLSorts). STAR-Fusion y FusionCatcher no disponibles para toda la cohorte; no se usa concordancia multi-herramienta.</div>
</div>

<div class="nav">
  <a href="#resumen">Resumen</a>
  <a href="#allsorts">ALLSorts</a>
  <a href="#phlike">Ph-like</a>
  <a href="#fusiones">Fusiones</a>
  <a href="#recurrentes">Recurrentes</a>
  <a href="#evidencia">Evidencia</a>
  <a href="#qc">QC</a>
  <a href="#metodos">Metodos</a>
</div>

<div class="card" id="resumen">
  <h2>&#x1F4CA; Resumen Ejecutivo</h2>
  <div class="summary-grid">
    <div class="summary-box"><div class="number" style="color:var(--accent);">50</div><div class="label">Muestras Analizadas</div></div>
    <div class="summary-box"><div class="number" style="color:var(--accent2);">{total_calls}</div><div class="label">Llamadas Arriba</div></div>
    <div class="summary-box"><div class="number" style="color:var(--success);">{high_calls}</div><div class="label">Alta Confianza</div></div>
    <div class="summary-box"><div class="number" style="color:var(--warn);">{medium_calls}</div><div class="label">Confianza Media</div></div>
    <div class="summary-box"><div class="number" style="color:var(--danger);">{driver_sample_count}</div><div class="label">Muestras con Driver</div></div>
    <div class="summary-box ph"><div class="number">{phlike_count}</div><div class="label">Ph / Ph-like / Ph Group</div></div>
  </div>
  <h3>Distribucion de Subtipos (ALLSorts)</h3>
  <div class="subtype-bar">{render_subtype_bar(allsorts)}</div>
  <div class="note">Este reporte integra evidencia completa de Arriba para 50/50 muestras y clasificacion ALLSorts para 50/50 muestras. Las conclusiones de fusiones se basan exclusivamente en Arriba.</div>
</div>

<div class="card" id="allsorts">
  <h2>&#x1F9EC; Clasificacion ALLSorts</h2>
  <p style="color:var(--muted); margin-bottom:1rem;">Clasificacion de subtipos de LLA-B basada en expresion genica. La columna Fuente indica si la fila proviene del CSV actual o fue recuperada del reporte historico.</p>
  <div class="table-scroll"><table>
    <thead><tr><th>Muestra</th><th>Subtipo</th><th>Prob. Principal</th><th>Ph-like Prob.</th><th>Ph Prob.</th><th>Pred</th><th>Fusion Driver Arriba</th><th>Fuente</th></tr></thead>
    <tbody>{render_allsorts_rows(allsorts, driver_map)}</tbody>
  </table></div>
</div>

<div class="card" id="phlike">
  <h2>&#x26A0;&#xFE0F; Fenotipo Ph-like: Analisis Detallado</h2>
  <div class="ph-highlight"><strong>{phlike_count} de 50 muestras</strong> tienen subtipo ALLSorts que contiene Ph/Ph-like o probabilidad Ph/Ph-like relevante. La evidencia de fusiones en esta tabla corresponde solo a Arriba.</div>
  <div class="table-scroll"><table>
    <thead><tr><th>Muestra</th><th>Subtipo</th><th>Ph-like Prob.</th><th>Ph Prob.</th><th>Fusion Arriba</th><th>Confianza</th><th>Soporte</th><th>Frame</th><th>Tipo</th></tr></thead>
    <tbody>{render_phlike_rows(phlike, driver_fusions)}</tbody>
  </table></div>
</div>

<div class="card" id="fusiones">
  <h2>&#x1F52C; Fusiones Clinicamente Relevantes por Muestra</h2>
  <p style="color:var(--muted); margin-bottom:1rem;">Se priorizan categorias driver de LLA-B detectadas por Arriba: {esc(', '.join(KNOWN_DRIVER_GROUPS))}.</p>
  <div class="table-scroll"><table>
    <thead><tr><th>Muestra</th><th>Subtipo</th><th>Categoria</th><th>Fusion</th><th>Tipo</th><th>Confianza</th><th>Soporte</th><th>Frame</th><th>Tags</th></tr></thead>
    <tbody>{render_driver_rows(driver_fusions, allsorts)}</tbody>
  </table></div>
</div>

<div class="card" id="recurrentes">
  <h2>&#x1F501; Fusiones Recurrentes de Cohorte</h2>
  <p style="color:var(--muted); margin-bottom:1rem;">Fusiones no-driver presentes en al menos 5 muestras. Esta tabla ayuda a separar eventos recurrentes de cohorte de fusiones driver clinicamente interpretables.</p>
  <div class="table-scroll"><table>
    <thead><tr><th>Fusion</th><th>Muestras</th><th>Eventos</th><th>Fraccion</th><th>Nota</th></tr></thead>
    <tbody>{render_recurrent_rows(recurrent)}</tbody>
  </table></div>
</div>

<div class="card" id="evidencia">
  <h2>&#x1F4CC; Disponibilidad de Evidencia</h2>
  <div class="table-scroll"><table>
    <thead><tr><th>Fuente</th><th>Cobertura</th><th>Uso en este reporte</th><th>Limitacion</th></tr></thead>
    <tbody>
      <tr><td>Arriba</td><td><span class="badge badge-pass">50/50</span></td><td>Base de todas las conclusiones de fusiones</td><td>Un solo llamador de fusiones para la cohorte completa</td></tr>
      <tr><td>ALLSorts</td><td><span class="badge badge-pass">50/50</span></td><td>Clasificacion de subtipo por expresion</td><td>No es un llamador de fusiones</td></tr>
      <tr><td>STAR-Fusion</td><td><span class="badge badge-fail">No cohorte completa</span></td><td>No usado para conclusiones 50-sample</td><td>Disponible solo para una parte historica de la cohorte</td></tr>
      <tr><td>FusionCatcher</td><td><span class="badge badge-fail">No cohorte completa</span></td><td>No usado para conclusiones 50-sample</td><td>Disponible solo para una parte historica de la cohorte</td></tr>
    </tbody>
  </table></div>
</div>

<div class="card" id="qc">
  <h2>&#x2705; Control de Calidad de Salidas</h2>
  <h3>Resumen por Muestra</h3>
  <div class="table-scroll"><table>
    <thead><tr><th>Muestra</th><th>Arriba Total</th><th>High</th><th>Medium</th><th>Driver-like</th><th>Estado</th></tr></thead>
    <tbody>{render_sample_qc_rows(sample_qc)}</tbody>
  </table></div>
</div>

<div class="card" id="metodos">
  <h2>&#x1F4D6; Metodos y Limitaciones</h2>
  <h3>Deteccion de Fusiones</h3>
  <p>Las fusiones se llamaron con Arriba sobre alineamiento STAR en modo streaming. Para este reporte, la tabla de fusiones comunes se regenero desde los 50 archivos <code>arriba_fusions.tsv</code> publicados en <code>results/</code>.</p>
  <h3>Clasificacion de Subtipos</h3>
  <p>ALLSorts se uso como clasificador de expresion y se reporta separado de la llamada de fusiones. Las filas actuales provienen de <code>summary/allsorts/probabilities.csv</code>; las filas historicas fueron recuperadas de <code>reporte_fusiones_completo.html</code> cuando no existia CSV machine-readable en la carpeta actual.</p>
  <h3>Limitaciones</h3>
  <p>STAR-Fusion y FusionCatcher no estan disponibles para las 50 muestras, por lo que no se hacen conclusiones basadas en soporte multi-herramienta. Las fusiones frecuentes tipo read-through, pseudogen o regiones repetitivas deben revisarse con cautela.</p>
</div>

<div class="footer">Archivo base: {esc(output_stem)}. Auditorias TSV generadas en <code>results/summary/</code>.</div>
</div>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="nextflow_pipeline/results")
    parser.add_argument(
        "--output-html",
        default="nextflow_pipeline/results/reporte_fusiones_50_muestras_arriba_allsorts.html",
    )
    parser.add_argument(
        "--allsorts-current",
        default="nextflow_pipeline/results/summary/allsorts/probabilities.csv",
    )
    parser.add_argument(
        "--historical-report",
        default="nextflow_pipeline/results/reporte_fusiones_completo.html",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_html = Path(args.output_html)
    summary_dir = results_dir / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    fusions, sample_qc, sample_table = load_arriba(results_dir)
    sample_ids = set(sample_table["sample_id"])
    allsorts = build_allsorts_audit(
        sample_ids,
        Path(args.allsorts_current),
        Path(args.historical_report),
    )

    common = build_common_fusions(fusions, total_samples=50)
    if common["total_samples"].nunique() != 1 or int(common["total_samples"].iloc[0]) != 50:
        raise SystemExit("Common fusion table validation failed: total_samples is not 50.")

    driver_fusions = fusions[fusions["driver_group"] != ""].copy()
    recurrent = common[
        (common["sample_count"] >= 5)
        & ~common["fusion"].isin(set(driver_fusions["fusion_key"]))
    ].copy()

    common_path = summary_dir / "arriba_common_fusions_50_samples.tsv"
    sample_qc_path = summary_dir / "arriba_sample_qc_50_samples.tsv"
    allsorts_path = summary_dir / "allsorts_50_sample_audit.tsv"
    driver_path = summary_dir / "clinically_relevant_arriba_50_samples.tsv"
    recurrent_path = summary_dir / "recurrent_arriba_50_samples.tsv"

    common.to_csv(common_path, sep="\t", index=False)
    sample_qc.to_csv(sample_qc_path, sep="\t", index=False)
    allsorts.to_csv(allsorts_path, sep="\t", index=False)
    driver_fusions.to_csv(driver_path, sep="\t", index=False)
    recurrent.to_csv(recurrent_path, sep="\t", index=False)

    html_text = build_html(
        allsorts=allsorts,
        sample_qc=sample_qc,
        common=common,
        driver_fusions=driver_fusions,
        recurrent=recurrent,
        output_stem=output_html.name,
    )
    output_html.write_text(html_text, encoding="utf-8")

    print(f"Wrote {output_html}")
    print(f"Wrote {common_path}")
    print(f"Wrote {allsorts_path}")


if __name__ == "__main__":
    main()

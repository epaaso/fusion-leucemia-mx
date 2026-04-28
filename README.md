# Fusion Detection + ALLSorts (Nextflow)

Nextflow DSL2 pipeline to call gene fusions with Arriba, STAR-Fusion, or FusionCatcher, with optional ALLSorts subtype classification. It performs read QC with FastQC, runs the selected fusion caller, summarizes fusions of interest, aggregates QC metrics with MultiQC, and can run ALLSorts from current or precomputed STAR gene-count outputs.

## Requirements
- Nextflow (DSL2) and Java 11+
- Conda if using the provided profile/environment (`-profile conda`)
- Arriba/STAR reference bundle present on disk (for Arriba/STAR-Fusion/STAR counts)
- FusionCatcher database (if using FusionCatcher method)

## Inputs
- Paired-end gzipped FASTQs under `FASTQ_DIR/SAMPLE_ID/*_R{1,2}_*.fastq.gz`
- Defaults are set in `nextflow_pipeline/nextflow.config`; override via `--param value` flags.
- Optional ALLSorts inputs for reuse mode:
  - `--allsorts_counts_csv` (prebuilt `samples x genes` matrix), or
  - `--allsorts_reads_glob` (glob of STAR `ReadsPerGene.out.tab` files from `work/` or `results/`)

## Quickstart
```bash
# Arriba
cd /home/epaaso/fusions
nextflow run nextflow_pipeline/main.nf \
  -profile conda \
  --fastq_dir /path/to/fastqs \
  --ref_dir /path/to/arriba_refs \
  --outdir /path/to/results \
  --threads 16 \
  --method arriba

# FusionCatcher (requires database)
nextflow run nextflow_pipeline/main.nf \
  -profile conda \
  --fastq_dir /path/to/fastqs \
  --fusioncatcher_dir /path/to/fusioncatcher_data \
  --outdir /path/to/results \
  --skip_fastqc true \
  --method fusioncatcher

# Run ALLSorts together with Arriba and reuse pipeline cache with -resume
nextflow run nextflow_pipeline/main.nf \
  -profile conda \
  --fastq_dir /path/to/fastqs \
  --ref_dir /path/to/arriba_refs \
  --outdir /path/to/results \
  --method arriba \
  --run_allsorts true \
  -resume

# ALLSorts-only from existing STAR-Fusion work files (no fusion caller rerun)
nextflow run nextflow_pipeline/main.nf \
  -profile conda \
  --method none \
  --run_allsorts true \
  --allsorts_reads_glob '/path/to/work/*/*/starfusion_out/ReadsPerGene.out.tab' \
  --gtf /path/to/annotation.gtf \
  --outdir /path/to/results
```
Add `--sample_id SAMPLE123` to run a single sample, and `-resume` to reuse prior work.

## Key Parameters (defaults in `nextflow_pipeline/nextflow.config`)
- `--fastq_dir`: `/datos/migccl/exposoma_fusion/fastqs`
- `--outdir`: `${launchDir}/results`
- `--sample_id`: run a single sample (default: all samples found)
- `--threads`: 16; `--star_ram`: `90.GB`
- `--method`: `arriba`, `starfusion`, `fusioncatcher`, or `none` (`none` is ALLSorts-only mode)
- `--annotate_common_fusions`: add literature/gene annotations to common fusion report
- `--literature_email`: email for NCBI E-utilities (recommended for PubMed queries)
- `--filter_rrna`: run optional rRNA filtering before fusion calling
- `--rrna_qc`: write rRNA fraction summary when filtering
- `--skip_fastqc`: skip both FastQC and MultiQC steps (default: `false`)
- `--rrna_ref`: rRNA reference FASTA for filtering (required if `--filter_rrna`)
- Reference paths: `--ref_dir` (default `${fastq_dir}/arriba_refs`), `--genome_fasta`, `--gtf`, `--star_index`
- Arriba resources: `--blacklist`, `--known_fusions`, `--protein_domains`, `--cytobands`
- STAR-Fusion resources: `--starfusion_genome_lib`, `--starfusion_extra_args`
- FusionCatcher resources: `--fusioncatcher_dir`
- ALLSorts: `--run_allsorts`, `--allsorts_counts_csv`, `--allsorts_reads_glob`, `--allsorts_strand`, `--allsorts_model_dir`, `--allsorts_parents`, `--allsorts_ball`

## Workflow Steps
- FASTQC: QC for each pair of reads
- STAR + Arriba: stream alignments to Arriba; uses STAR `Chimeric.out.junction` for fusion calling
- STAR-Fusion: runs STAR-Fusion directly on FASTQs using a CTAT genome lib
- FusionCatcher: runs FusionCatcher on FASTQs (automatically interleaves paired-end reads)
- STAR counts (for ALLSorts when needed): creates per-sample `ReadsPerGene.out.tab`
- ALLSorts: builds a counts matrix and predicts B-ALL subtype probabilities
- Summarize: `bin/summarize_fusions.py` flags fusions involving protocol-relevant genes
- Common fusions: aggregates fusion calls across all samples
- rRNA filter: optional BBDuk filter against rRNA reference with summary table
- MultiQC: aggregates QC reports

## Outputs
- Arriba: `${outdir}/${sample_id}/arriba_fusions.tsv` and `arriba_fusions.discarded.tsv`
- STAR-Fusion: `${outdir}/${sample_id}/star-fusion.fusion_predictions.tsv`
- FusionCatcher: `${outdir}/${sample_id}/fusioncatcher_fusions.txt` and `${outdir}/${sample_id}/fusioncatcher_summary.txt`
- STAR gene counts: `${outdir}/${sample_id}/*.ReadsPerGene.out.tab`
- Protocol summary: `${outdir}/${sample_id}/${method}_fusions_protocol.tsv` (filtered genes)
- Common fusions: `${outdir}/summary/${method}_common_fusions.tsv` (includes `fusion_literature_summary`, `gene1_function`, `gene2_function`)
- rRNA fraction: `${outdir}/summary/rrna_fraction.tsv`
- ALLSorts: `${outdir}/summary/allsorts/predictions.csv` and `${outdir}/summary/allsorts/probabilities.csv`
- `${outdir}/${sample_id}/qc/` from FastQC
- `${outdir}/multiqc/multiqc_report.html`
- STAR/Arriba logs for each sample alongside outputs

## Notes
- Ensure the reference bundle exists before running; `check_refs()` currently emits warnings instead of stopping the run.
- Update `nextflow_pipeline/nextflow.config` if your paths or resource requirements differ.
- Common fusion annotations are pulled from PubMed abstracts and MyGene.info summaries; set `--literature_email` to comply with NCBI E-utilities usage.
- For maximum reuse of existing preprocessing, keep `nextflow_pipeline/work/` and run with `-resume`; for historical STAR-Fusion work directories, use `--method none --run_allsorts true --allsorts_reads_glob ...`.

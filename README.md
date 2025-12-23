# Fusion Detection (Arriba/STAR-Fusion) (Nextflow)

Minimal Nextflow DSL2 pipeline to call gene fusions with Arriba or STAR-Fusion. It performs read QC with FastQC, runs the selected fusion caller, summarizes fusions of interest, and aggregates QC metrics with MultiQC.

## Requirements
- Nextflow (DSL2) and Java 11+
- Conda if using the provided profile/environment (`-profile conda` and `nextflow_pipeline/environment.yml`)
- Arriba/STAR reference bundle present on disk (the pipeline only warns when references are missing)

## Inputs
- Paired-end gzipped FASTQs under `FASTQ_DIR/SAMPLE_ID/*_R{1,2}_*.fastq.gz`
- Defaults are set in `nextflow_pipeline/nextflow.config`; override via `--param value` flags.

## Quickstart
```bash
cd /home/epaaso/fusions
nextflow run nextflow_pipeline/main.nf \
  -profile conda \
  --fastq_dir /path/to/fastqs \
  --ref_dir /path/to/arriba_refs \
  --outdir /path/to/results \
  --threads 16 \
  --method arriba
```
Add `--sample_id SAMPLE123` to run a single sample, and `-resume` to reuse prior work.

## Key Parameters (defaults in `nextflow_pipeline/nextflow.config`)
- `--fastq_dir`: `/datos/migccl/exposoma_fusion/fastqs`
- `--outdir`: `${launchDir}/results`
- `--sample_id`: run a single sample (default: all samples found)
- `--threads`: 16; `--star_ram`: `90.GB`
- `--method`: `arriba` or `starfusion`
- `--annotate_common_fusions`: add literature/gene annotations to common fusion report
- `--literature_email`: email for NCBI E-utilities (recommended for PubMed queries)
- `--filter_rrna`: run optional rRNA filtering before fusion calling
- `--rrna_qc`: write rRNA fraction summary when filtering
- `--rrna_ref`: rRNA reference FASTA for filtering (required if `--filter_rrna`)
- Reference paths: `--ref_dir` (default `${fastq_dir}/arriba_refs`), `--genome_fasta`, `--gtf`, `--star_index`
- Arriba resources: `--blacklist`, `--known_fusions`, `--protein_domains`, `--cytobands`
- STAR-Fusion resources: `--starfusion_genome_lib`, `--starfusion_extra_args`

## Workflow Steps
- FASTQC: QC for each pair of reads
- STAR + Arriba: stream alignments to Arriba; uses STAR `Chimeric.out.junction` for fusion calling
- STAR-Fusion: runs STAR-Fusion directly on FASTQs using a CTAT genome lib
- Summarize: `bin/summarize_fusions.py` flags fusions involving protocol-relevant genes
- Common fusions: aggregates fusion calls across all samples
- rRNA filter: optional BBDuk filter against rRNA reference with summary table
- MultiQC: aggregates QC reports

## Outputs
- Arriba: `${outdir}/${sample_id}/arriba_fusions.tsv` and `arriba_fusions.discarded.tsv`
- STAR-Fusion: `${outdir}/${sample_id}/star-fusion.fusion_predictions.tsv`
- Protocol summary: `${outdir}/${sample_id}/${method}_fusions_protocol.tsv` (filtered genes)
- Common fusions: `${outdir}/summary/${method}_common_fusions.tsv` (includes `fusion_literature_summary`, `gene1_function`, `gene2_function`)
- rRNA fraction: `${outdir}/summary/rrna_fraction.tsv`
- `${outdir}/${sample_id}/qc/` from FastQC
- `${outdir}/multiqc/multiqc_report.html`
- STAR/Arriba logs for each sample alongside outputs

## Notes
- Ensure the reference bundle exists before running; `check_refs()` currently emits warnings instead of stopping the run.
- Update `nextflow_pipeline/nextflow.config` if your paths or resource requirements differ.
- Common fusion annotations are pulled from PubMed abstracts and MyGene.info summaries; set `--literature_email` to comply with NCBI E-utilities usage.

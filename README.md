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
- STAR + Arriba sensitivity: `--star_out_filter_multimap_nmax 50`, `--star_chim_out_type 'WithinBAM HardClip'`, `--star_chim_multimap_nmax 50`
- `--method`: `arriba`, `starfusion`, `fusioncatcher`, or `none` (`none` is ALLSorts-only mode)
- `--annotate_common_fusions`: add literature/gene annotations to common fusion report
- `--literature_email`: email for NCBI E-utilities (recommended for PubMed queries)
- `--filter_rrna`: run optional rRNA filtering before fusion calling
- `--rrna_qc`: write rRNA fraction summary when filtering
- `--skip_fastqc`: skip both FastQC and MultiQC steps (default: `false`)
- `--rrna_ref`: rRNA reference FASTA for filtering (required if `--filter_rrna`)
- Reference paths: `--ref_dir` (default `${fastq_dir}/arriba_refs`), `--genome_fasta`, `--gtf`, `--star_index`
- Arriba resources: `--blacklist`, `--known_fusions`, `--protein_domains`, `--cytobands`
- STAR-Fusion resources: `--starfusion_genome_lib`, `--starfusion_extra_args`, `--starfusion_threads 19`, `--starfusion_ram 100.GB`, `--starfusion_max_forks 1`. A single run may require 80-100 GB, so this 125-GiB host processes one STAR-Fusion sample at a time.
- FusionCatcher resources: `--fusioncatcher_dir`
- ALLSorts: `--run_allsorts`, `--allsorts_counts_csv`, `--allsorts_reads_glob`, `--allsorts_strand`, `--allsorts_model_dir`, `--allsorts_parents`, `--allsorts_ball`

## Workflow Steps
- FASTQC: QC for each pair of reads
- STAR + Arriba: streams an unsorted BAM containing supplementary chimeric alignments (`--chimOutType WithinBAM HardClip`) from STAR to Arriba through stdin. Multimap limits are set to 50 to improve sensitivity for repetitive and complex rearrangements such as `DUX4::IGH`, which can provide evidence consistent with enhancer hijacking; the fusion call alone does not prove the regulatory mechanism.
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

## Shared STAR Alignment (arriba_starfusion)

The combined `--method arriba_starfusion` mode runs STAR alignment once per sample using a CTAT pseudo-masked reference, streams the unsorted BAM directly to Arriba, and retains `Chimeric.out.junction` in the Nextflow work cache for concurrent STAR-Fusion post-processing. This reduces alignment overhead and optimizes resource usage.

> [!IMPORTANT]
> **FusionCatcher Alignment Exemption:** FusionCatcher cannot participate in the shared STAR alignment workflow because it does not accept pre-aligned BAM or chimeric junction files as input. It requires raw FASTQ reads and will always perform its own separate, duplicate alignment from scratch under the hood (using its internal multi-stage Bowtie/Bowtie2/STAR alignment pipeline).


### Quickstart Example
```bash
nextflow run nextflow_pipeline/main.nf \
  -profile conda \
  --fastq_dir /path/to/fastqs \
  --ref_dir /path/to/arriba_refs \
  --outdir /path/to/results \
  --method arriba_starfusion \
  --star_index /path/to/arriba_refs/star_grch38_v50_unmasked \
  --starfusion_genome_lib /path/to/arriba_refs/GRCh38_gencode_v50_CTAT_lib/ctat_genome_lib_build_dir \
  --gtf /path/to/arriba_refs/gencode.v50.basic.annotation.gtf
```

### Reference Locking and Masking Behavior
- **Primary Reference:** Built on GRCh38 assembly and GENCODE v50.
- **CTAT Pseudo-Masked Reference:** A pseudo-masked reference bundle is used for STAR-Fusion and the shared alignment.
- **Unmasked Reference:** A matching unmasked STAR index is built for fallback Arriba alignments.
- **Note:** Reference versions must never be mixed within a production cohort.

### Efficient Resource Scheduling
Process resources are configured to maximize throughput on a 20-CPU/125-GiB host:
- `STAR_SHARED_ARRIBA`: memory-heavy task allocated `16 cpus`, `100.GB` RAM, and restricted to `maxForks: 1` to ensure only one alignment runs at a time.
- `STARFUSION_CALL`: lightweight junction-processing task allocated `1 cpu`, `12.GB` RAM, and running concurrently with `maxForks: 4`.

### Validation & Fallback Workflow
Before enabling the combined mode for the full cohort, validate using:
1. `LAL1276_S9` (known `DUX4::IGH` candidate)
2. `103-MO-15-563486477` (strong `TCF3::PBX1` control)

Compare shared-workflow results against caller-specific outputs. The validation fails if:
- Expected `DUX4::IGH` or `TCF3::PBX1` evidence disappears.
- Clinically relevant fusions are lost.
- Support decreases materially or reference masking causes unacceptable Arriba behavior.

**Fallback Workflow:** If validation fails, production uses separate alignments:
- Unmasked GRCh38/GENCODE v50 STAR index for Arriba streaming (`--method arriba`).
- CTAT pseudo-masked GRCh38/GENCODE v50 reference library for STAR-Fusion (`--method starfusion`).

## Notes
- Ensure the reference bundle exists before running; `check_refs()` currently emits warnings instead of stopping the run.
- Update `nextflow_pipeline/nextflow.config` if your paths or resource requirements differ.
- Common fusion annotations are pulled from PubMed abstracts and MyGene.info summaries; set `--literature_email` to comply with NCBI E-utilities usage.
- For maximum reuse of existing preprocessing, keep `nextflow_pipeline/work/` and run with `-resume`; for historical STAR-Fusion work directories, use `--method none --run_allsorts true --allsorts_reads_glob ...`.

## Generating HTML Reports with Dispersed Results

> [!IMPORTANT]
> **Clinical Metadata Requirement:** To build consolidated HTML reports and map patients in a clinical context, a clinical metadata table (e.g. `Datos_clinicos_con_fusiones.csv` containing patient IDs, diagnostics, and clinical features) is required. For privacy and data protection compliance, this file must remain strictly local on your machine and must never be committed to Git or pushed to GitHub (it is ignored by default in `.gitignore`).

If your fusion caller outputs (Arriba, STAR-Fusion, FusionCatcher) and your ALLSorts subtype classification results are located in different directories, you can generate a consolidated HTML report by following these methods:

### Method A: Using Command-Line Arguments (Direct)
You can directly pass different paths to the report generator utility script:

```bash
python3 nextflow_pipeline/bin/generate_complete_50_sample_report.py \
  --results-dir /path/to/fusion_caller_results \
  --allsorts-current /path/to/allsorts_summary/allsorts/probabilities.csv \
  --output-html /path/to/reporte_fusiones.html
```

* `--results-dir`: Path to the directory containing individual sample directories (where each subdirectory contains `arriba_fusions.tsv`, etc.).
* `--allsorts-current`: Path to the ALLSorts `probabilities.csv` file.

### Method B: Consolidating via Symbolic Links (Recommended for Multiple Cohorts)
If you want to merge multiple cohorts (e.g. 50 original samples and 22 new samples) that were run in different places:

1. Create a central results directory and symlink the sample folders from both locations:
   ```bash
   mkdir -p consolidated_results
   # Link cohort 1 samples
   for dir in /path/to/cohort1/*; do ln -s "$dir" consolidated_results/$(basename "$dir"); done
   # Link cohort 2 samples
   for dir in /path/to/cohort2/*; do ln -s "$dir" consolidated_results/$(basename "$dir"); done
   ```

2. Merge the ALLSorts probabilities into a single CSV using a Python snippet:
   ```python
   import pandas as pd
   p1 = pd.read_csv('/path/to/cohort1/summary/allsorts/probabilities.csv')
   p2 = pd.read_csv('/path/to/cohort2/summary/allsorts/probabilities.csv')
   pd.concat([p1, p2]).drop_duplicates(subset=['sample_id']).to_csv('consolidated_probabilities.csv', index=False)
   ```

3. Run the report generator:
   ```bash
   python3 nextflow_pipeline/bin/generate_complete_50_sample_report.py \
     --results-dir consolidated_results \
     --allsorts-current consolidated_probabilities.csv \
     --output-html reporte_combinado.html
   ```


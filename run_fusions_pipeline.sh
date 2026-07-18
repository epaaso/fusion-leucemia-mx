#!/bin/bash
set -e

echo "=== INICIANDO PIPELINE DE FUSIONES COMPLETO ==="
date

# 1. Ejecutar Arriba + STAR-Fusion
echo "--> Paso 1: Ejecutando arriba_starfusion..."
nextflow run nextflow_pipeline/main.nf \
  -profile conda \
  --fastq_dir /datos/migccl/sined/fastq573 \
  --outdir /datos/migccl/sined/qc_results \
  --method arriba_starfusion \
  --filter_rrna true \
  --rrna_ref /datos/migccl/rRNA_refs/human_rRNA.fasta \
  --skip_fastqc true \
  --threads 9 \
  --shared_star_threads 9 \
  --starfusion_threads 9 \
  --rrna_threads 9 \
  --star_ram 60.GB \
  --shared_star_ram 60.GB \
  --starfusion_ram 60.GB \
  --rrna_ram 60.GB \
  -resume

# 2. Ejecutar FusionCatcher
echo "--> Paso 2: Ejecutando fusioncatcher..."
nextflow run nextflow_pipeline/main.nf \
  -profile conda \
  --fastq_dir /datos/migccl/sined/fastq573 \
  --outdir /datos/migccl/sined/qc_results \
  --method fusioncatcher \
  --fusioncatcher_dir /datos/migccl/fusioncatcher_data \
  --filter_rrna true \
  --rrna_ref /datos/migccl/rRNA_refs/human_rRNA.fasta \
  --skip_fastqc true \
  --threads 9 \
  --fusioncatcher_threads 9 \
  --rrna_threads 9 \
  --fusioncatcher_ram 60.GB \
  --rrna_ram 60.GB \
  -resume

echo "=== PIPELINE FINALIZADO CON ÉXITO ==="
date

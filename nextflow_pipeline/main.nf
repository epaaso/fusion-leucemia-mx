#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def method = (params.method ?: 'arriba').toString().toLowerCase()

def check_refs(String method) {
    def missing = []

    if (method == 'arriba') {
        if (!file(params.genome_fasta).exists()) missing.add(params.genome_fasta)
        if (!file(params.gtf).exists()) missing.add(params.gtf)
        if (!file(params.star_index).exists()) missing.add(params.star_index)
        if (!file(params.blacklist).exists()) missing.add(params.blacklist)
        if (!file(params.known_fusions).exists()) missing.add(params.known_fusions)
        if (!file(params.protein_domains).exists()) missing.add(params.protein_domains)
    } else if (method == 'starfusion') {
        if (!file(params.starfusion_genome_lib).exists()) missing.add(params.starfusion_genome_lib)
    } else {
        error "Unknown method '${params.method}'. Use 'arriba' or 'starfusion'."
    }

    if (params.filter_rrna) {
        if (!params.rrna_ref) {
            missing.add("params.rrna_ref (not set)")
        } else if (!file(params.rrna_ref).exists()) {
            missing.add(params.rrna_ref)
        }
    }

    if (missing) {
        log.warn "Missing reference files: ${missing.join(', ')}"
        log.warn "Please ensure references exist under ${params.ref_dir} (or override with params)."
    }
}

def make_samples_channel() {
    // If --fastq_glob is provided, use it as-is. Otherwise use fastq_dir + sample_id/all.
    def input_pattern
    if (params.fastq_glob) {
        input_pattern = params.fastq_glob
    } else if (params.sample_id) {
        input_pattern = "${params.fastq_dir}/${params.sample_id}/*_R{1,2}_*.fastq.gz"
    } else {
        input_pattern = "${params.fastq_dir}/*/*_R{1,2}_*.fastq.gz"
    }

    // Group all lane FASTQs by sample directory name, then split into R1/R2 lists.
    Channel
        .fromPath(input_pattern, checkIfExists: true)
        .map { f -> tuple(f.getParent().getName(), f) }
        .groupTuple()
        .map { sample_id, files ->
            def r1 = files.findAll { it.getName().contains('_R1_') }.sort()
            def r2 = files.findAll { it.getName().contains('_R2_') }.sort()
            if (!r1 || !r2) error("No R1/R2 FASTQs found for sample '${sample_id}'")
            if (r1.size() != r2.size()) error("Mismatched R1/R2 counts for sample '${sample_id}' (${r1.size()} vs ${r2.size()})")
            tuple(sample_id, r1, r2)
        }
}

workflow {
    if (!(method in ['arriba', 'starfusion'])) {
        error "Unknown method '${params.method}'. Use 'arriba' or 'starfusion'."
    }

    check_refs(method)

    ch_samples_raw = make_samples_channel()

    if (params.filter_rrna) {
        RRNA_FILTER(ch_samples_raw)
        ch_samples = RRNA_FILTER.out.reads
        if (params.rrna_qc) {
            RRNA_QC(RRNA_FILTER.out.stats.collect())
        }
    } else {
        ch_samples = ch_samples_raw
    }

    FASTQC(ch_samples)

    if (method == 'arriba') {
        STAR_ARRIBA(ch_samples)
        SUMMARIZE(STAR_ARRIBA.out.fusions)
        COMMON_FUSIONS(STAR_ARRIBA.out.fusions.map { it[1] }.collect())
    }
    if (method == 'starfusion') {
        STARFUSION(ch_samples)
        SUMMARIZE(STARFUSION.out.fusions)
        COMMON_FUSIONS(STARFUSION.out.fusions.map { it[1] }.collect())
    }

    MULTIQC(FASTQC.out.zip.collect().ifEmpty([]))
}

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(r1_reads), path(r2_reads)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    def fastqs = (r1_reads + r2_reads).join(' ')
    """
    fastqc -t ${task.cpus} ${fastqs}
    """
}

process RRNA_FILTER {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/rrna_filter", mode: 'copy'

    input:
    tuple val(sample_id), path(r1_reads), path(r2_reads)

    output:
    tuple val(sample_id),
        path("${sample_id}_R1.filtered.fastq.gz"),
        path("${sample_id}_R2.filtered.fastq.gz"),
        emit: reads
    path "rrna_${sample_id}.stats", emit: stats

    script:
    def r1_csv = r1_reads.join(',')
    def r2_csv = r2_reads.join(',')
    """
    set -euo pipefail
    bbduk.sh \\
        in1=${r1_csv} in2=${r2_csv} \\
        out1=${sample_id}_R1.filtered.fastq.gz out2=${sample_id}_R2.filtered.fastq.gz \\
        ref=${params.rrna_ref} k=${params.rrna_k} hdist=${params.rrna_hdist} \\
        stats=rrna_${sample_id}.stats \\
        overwrite=t threads=${task.cpus}
    """
}

process STAR_ARRIBA {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1_reads), path(r2_reads)

    output:
    tuple val(sample_id), path("arriba_fusions.tsv"), emit: fusions
    tuple val(sample_id), path("arriba_fusions.discarded.tsv"), emit: discarded
    path "*.{out,tab,log,junction}", emit: logs, optional: true

    script:
    def r1_csv = r1_reads.join(',')
    def r2_csv = r2_reads.join(',')

    """
    set -euo pipefail
    mkdir -p star_out

    # Streaming mode: STAR writes BAM to stdout; Arriba reads from stdin.
    # Do NOT use '-c Chimeric.out.junction' with this mode.
    STAR --runThreadN ${task.cpus} \\
        --genomeDir ${params.star_index} --genomeLoad NoSharedMemory \\
        --readFilesIn ${r1_csv} ${r2_csv} --readFilesCommand zcat \\
        --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \\
        --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 \\
        --outFileNamePrefix star_out/${sample_id}. \\
    | arriba \\
        -x /dev/stdin \\
        -o arriba_fusions.tsv -O arriba_fusions.discarded.tsv \\
        -a ${params.genome_fasta} -g ${params.gtf} \\
        -b ${params.blacklist} -k ${params.known_fusions} -t ${params.known_fusions} -p ${params.protein_domains}

    cp star_out/*.out star_out/*.tab star_out/*.log star_out/*.junction . 2>/dev/null || true
    """
}

process STARFUSION {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1_reads), path(r2_reads)

    output:
    tuple val(sample_id), path("star-fusion.fusion_predictions.tsv"), emit: fusions
    path "star-fusion.fusion_predictions.abridged.tsv", emit: abridged, optional: true
    path "*.log", emit: logs, optional: true

    script:
    def r1_csv = r1_reads.join(',')
    def r2_csv = r2_reads.join(',')
    def extra_args = params.starfusion_extra_args ?: ""

    """
    set -euo pipefail
    mkdir -p starfusion_out

    STAR-Fusion \\
        --genome_lib_dir ${params.starfusion_genome_lib} \\
        --left_fq ${r1_csv} --right_fq ${r2_csv} \\
        --CPU ${task.cpus} \\
        ${extra_args} \\
        --output_dir starfusion_out

    cp starfusion_out/star-fusion.fusion_predictions.tsv .
    cp starfusion_out/star-fusion.fusion_predictions.abridged.tsv . 2>/dev/null || true
    cp starfusion_out/*.log . 2>/dev/null || true
    """
}

process SUMMARIZE {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fusions)

    output:
    path "${method}_fusions_protocol.tsv"

    script:
    """
    summarize_fusions.py --input ${fusions} --output ${method}_fusions_protocol.tsv --format ${method}
    """
}

process COMMON_FUSIONS {
    tag "${method}"
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    val fusion_files

    output:
    path "${method}_common_fusions.tsv"

    script:
    def annotate_flag = params.annotate_common_fusions ? "--annotate" : "--no-annotate"
    def email_arg = params.literature_email ? "--email ${params.literature_email}" : ""
    def inputs_arg = fusion_files.collect { it.toString() }.join(' ')
    """
    aggregate_fusions.py --inputs ${inputs_arg} --output ${method}_common_fusions.tsv --format ${method} ${annotate_flag} ${email_arg}
    """
}

process RRNA_QC {
    tag "rrna_qc"
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    val stats_files

    output:
    path "rrna_fraction.tsv"

    script:
    def inputs_arg = stats_files.collect { it.toString() }.join(' ')
    """
    rrna_qc.py --inputs ${inputs_arg} --output rrna_fraction.tsv
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path fastqc_files

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

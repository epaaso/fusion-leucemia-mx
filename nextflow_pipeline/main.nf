#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def method = (params.method ?: 'arriba').toString().toLowerCase()
def skip_fastqc = (params.skip_fastqc?.toString()?.toBoolean()) ?: false
def run_allsorts = (params.run_allsorts?.toString()?.toBoolean()) ?: false
def has_external_allsorts_counts = params.allsorts_counts_csv ? true : false
def has_external_allsorts_reads = params.allsorts_reads_glob ? true : false

def resolve_fusioncatcher_db_dir() {
    def configured = file(params.fusioncatcher_dir)
    if (!configured.exists()) return configured

    if (file("${configured}/organism.txt").exists()) return configured

    def candidates = configured
        .toFile()
        .listFiles()
        ?.findAll { it.isDirectory() && new File(it, 'organism.txt').exists() } ?: []

    if (candidates.size() == 1) {
        def detected = file(candidates[0].toString())
        log.warn "FusionCatcher DB auto-detected at '${detected}' because '${configured}' has no organism.txt"
        return detected
    }

    return configured
}

def fusioncatcher_db_dir = method == 'fusioncatcher' ? resolve_fusioncatcher_db_dir() : null

def check_refs(
    String method,
    def fusioncatcher_db_dir,
    boolean run_allsorts,
    boolean has_external_allsorts_counts,
    boolean has_external_allsorts_reads
) {
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
    } else if (method == 'fusioncatcher') {
        if (!fusioncatcher_db_dir.exists()) {
            missing.add(params.fusioncatcher_dir)
        } else if (!file("${fusioncatcher_db_dir}/organism.txt").exists()) {
            missing.add("${fusioncatcher_db_dir}/organism.txt")
        }
    } else if (method == 'none') {
        // Allows ALLSorts-only runs when precomputed inputs are provided.
    } else {
        error "Unknown method '${params.method}'. Use 'arriba', 'starfusion', 'fusioncatcher', or 'none'."
    }

    if (params.filter_rrna) {
        if (!params.rrna_ref) {
            missing.add("params.rrna_ref (not set)")
        } else if (!file(params.rrna_ref).exists()) {
            missing.add(params.rrna_ref)
        }
    }

    if (run_allsorts) {
        if (has_external_allsorts_counts) {
            if (!file(params.allsorts_counts_csv).exists()) {
                missing.add(params.allsorts_counts_csv)
            }
        } else {
            if (!file(params.gtf).exists()) {
                missing.add(params.gtf)
            }
            if (!has_external_allsorts_reads && !file(params.star_index).exists()) {
                missing.add(params.star_index)
            }
        }
        if (params.allsorts_model_dir && !file(params.allsorts_model_dir).exists()) {
            missing.add(params.allsorts_model_dir)
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

def infer_sample_id_from_reads_per_gene(def reads_path) {
    def file_name = reads_path.getFileName().toString()

    if (file_name.endsWith('.ReadsPerGene.out.tab') && file_name != 'ReadsPerGene.out.tab') {
        return file_name.substring(0, file_name.length() - '.ReadsPerGene.out.tab'.length())
    }

    if (file_name == 'ReadsPerGene.out.tab') {
        def parent = reads_path.getParent()
        if (parent && parent.getFileName().toString() == 'starfusion_out') {
            def work_dir = parent.getParent()
            if (work_dir) {
                try {
                    def links = work_dir.toFile().listFiles()?.findAll {
                        java.nio.file.Files.isSymbolicLink(it.toPath()) && it.getName().contains('_R1_')
                    } ?: []
                    if (links) {
                        def resolved = links[0].toPath().toRealPath()
                        return resolved.getParent().getFileName().toString()
                    }
                } catch (Exception ignored) {
                    // Fall through to script parsing and parent directory fallback below.
                }

                // Fallback: parse staged input declarations if symlinks are unavailable.
                def command_run = work_dir.resolve('.command.run').toFile()
                if (command_run.exists()) {
                    def match = (command_run.text =~ /\/fastqs\/([^\/\s]+)\/[^\/\s]*_R[12]_[^\/\s]*\.fastq\.gz/)
                    if (match.find()) {
                        return match.group(1)
                    }
                }
            }
        }
    }

    def parent = reads_path.getParent()
    return parent ? parent.getFileName().toString() : file_name
}

def make_allsorts_external_reads_channel() {
    Channel
        .fromPath(params.allsorts_reads_glob, checkIfExists: true)
        .map { p -> tuple(infer_sample_id_from_reads_per_gene(p), p) }
        .groupTuple()
        .map { sample_id, files ->
            if (!files || files.isEmpty()) {
                error("No ReadsPerGene file found for sample '${sample_id}' from --allsorts_reads_glob")
            }
            def selected = files.size() == 1
                ? files[0]
                : files.sort { a, b -> a.toFile().lastModified() <=> b.toFile().lastModified() }[-1]
            if (files.size() > 1) {
                log.warn "Multiple ReadsPerGene files found for '${sample_id}'. Using latest: ${selected}"
            }
            tuple(sample_id, selected)
        }
}

workflow {
    if (!(method in ['arriba', 'starfusion', 'fusioncatcher', 'none'])) {
        error "Unknown method '${params.method}'. Use 'arriba', 'starfusion', 'fusioncatcher', or 'none'."
    }

    if (method == 'none' && !run_allsorts) {
        error "method=none requires --run_allsorts true"
    }

    if (run_allsorts && has_external_allsorts_counts && has_external_allsorts_reads) {
        error "Use only one ALLSorts input source: --allsorts_counts_csv OR --allsorts_reads_glob"
    }

    if (method == 'none' && !has_external_allsorts_counts && !has_external_allsorts_reads) {
        error "method=none requires --allsorts_counts_csv or --allsorts_reads_glob"
    }

    check_refs(
        method,
        fusioncatcher_db_dir,
        run_allsorts,
        has_external_allsorts_counts,
        has_external_allsorts_reads
    )

    def ch_samples = null
    def ch_allsorts_reads = null
    def ch_allsorts_counts = null

    if (run_allsorts && has_external_allsorts_counts) {
        ch_allsorts_counts = Channel.value(file(params.allsorts_counts_csv))
    }

    if (run_allsorts && has_external_allsorts_reads) {
        ch_allsorts_reads = make_allsorts_external_reads_channel()
    }

    if (method != 'none') {
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

        if (!skip_fastqc) {
            FASTQC(ch_samples)
        } else {
            log.info "Skipping FASTQC and MULTIQC because --skip_fastqc is true"
        }

        if (method == 'arriba') {
            STAR_ARRIBA(ch_samples)
            SUMMARIZE(STAR_ARRIBA.out.fusions)
            COMMON_FUSIONS(STAR_ARRIBA.out.fusions.map { it[1] }.collect())
            if (run_allsorts && !has_external_allsorts_counts && !has_external_allsorts_reads) {
                ch_allsorts_reads = STAR_ARRIBA.out.reads_per_gene
            }
        }
        if (method == 'starfusion') {
            STARFUSION(ch_samples)
            SUMMARIZE(STARFUSION.out.fusions)
            COMMON_FUSIONS(STARFUSION.out.fusions.map { it[1] }.collect())
            if (run_allsorts && !has_external_allsorts_counts && !has_external_allsorts_reads) {
                ch_allsorts_reads = STARFUSION.out.reads_per_gene
            }
        }
        if (method == 'fusioncatcher') {
            FUSIONCATCHER(ch_samples)
            SUMMARIZE(FUSIONCATCHER.out.fusions)
            COMMON_FUSIONS(FUSIONCATCHER.out.fusions.map { it[1] }.collect())
        }

        if (!skip_fastqc) {
            MULTIQC(FASTQC.out.zip.collect().ifEmpty([]))
        }
    }

    if (run_allsorts && !has_external_allsorts_counts) {
        if (ch_allsorts_reads == null) {
            if (method == 'none') {
                error "ALLSorts reads source could not be resolved for method=none"
            }
            STAR_COUNTS(ch_samples)
            ch_allsorts_reads = STAR_COUNTS.out.reads_per_gene
        }
        BUILD_ALLSORTS_COUNTS(ch_allsorts_reads.collect())
        ch_allsorts_counts = BUILD_ALLSORTS_COUNTS.out
    }

    if (run_allsorts) {
        ALLSORTS(ch_allsorts_counts)
    }
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
    tuple val(sample_id), path("${sample_id}.ReadsPerGene.out.tab"), emit: reads_per_gene
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
        --quantMode GeneCounts \\
        --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \\
        --outFilterMultimapNmax ${params.star_out_filter_multimap_nmax} --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentMin 10 --chimOutType ${params.star_chim_out_type} --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax ${params.star_chim_multimap_nmax} \\
        --outFileNamePrefix star_out/${sample_id}. \\
    | arriba \\
        -x /dev/stdin \\
        -o arriba_fusions.tsv -O arriba_fusions.discarded.tsv \\
        -a ${params.genome_fasta} -g ${params.gtf} \\
        -b ${params.blacklist} -k ${params.known_fusions} -t ${params.known_fusions} -p ${params.protein_domains}

    cp star_out/*.out star_out/*.tab star_out/*.log star_out/*.junction . 2>/dev/null || true
    cp star_out/${sample_id}.ReadsPerGene.out.tab .
    """
}

process STARFUSION {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1_reads), path(r2_reads)

    output:
    tuple val(sample_id), path("star-fusion.fusion_predictions.tsv"), emit: fusions
    tuple val(sample_id), path("ReadsPerGene.out.tab"), emit: reads_per_gene
    path "star-fusion.fusion_predictions.abridged.tsv", emit: abridged, optional: true
    path "*.log", emit: logs, optional: true

    script:
    def extra_args = params.starfusion_extra_args ?: ""

    """
    set -euo pipefail
    mkdir -p starfusion_out

    # Merge FASTQs to avoid issues with comma-separated lists in FusionInspector
    cat ${r1_reads} > combined_R1.fastq.gz
    cat ${r2_reads} > combined_R2.fastq.gz

    STAR-Fusion \\
        --genome_lib_dir ${params.starfusion_genome_lib} \\
        --left_fq combined_R1.fastq.gz --right_fq combined_R2.fastq.gz \\
        --CPU ${task.cpus} \\
        ${extra_args} \\
        --output_dir starfusion_out

    cp starfusion_out/star-fusion.fusion_predictions.tsv .
    cp starfusion_out/star-fusion.fusion_predictions.abridged.tsv . 2>/dev/null || true
    cp starfusion_out/ReadsPerGene.out.tab . 2>/dev/null || true
    cp starfusion_out/*.log . 2>/dev/null || true
    """
}

process STAR_COUNTS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(r1_reads), path(r2_reads)

    output:
    tuple val(sample_id), path("${sample_id}.ReadsPerGene.out.tab"), emit: reads_per_gene
    path "${sample_id}.Log.out", optional: true
    path "${sample_id}.Log.final.out", optional: true
    path "${sample_id}.Log.progress.out", optional: true

    script:
    def r1_csv = r1_reads.join(',')
    def r2_csv = r2_reads.join(',')

    """
    set -euo pipefail
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${params.star_index} --genomeLoad NoSharedMemory \\
        --readFilesIn ${r1_csv} ${r2_csv} --readFilesCommand zcat \\
        --quantMode GeneCounts \\
        --outSAMtype None \\
        --outFileNamePrefix ${sample_id}.
    """
}

process BUILD_ALLSORTS_COUNTS {
    publishDir "${params.outdir}/summary/allsorts", mode: 'copy'

    input:
    val reads_entries

    output:
    path "allsorts_counts.csv"

    script:
    def pairs
    if (reads_entries && reads_entries[0] instanceof List) {
        // Shape: [[sample_id, path], [sample_id, path], ...]
        pairs = reads_entries.collect { it.take(2) }
    } else {
        // Shape after collect() over tuple channel: [sample_id, path, sample_id, path, ...]
        if ((reads_entries?.size() ?: 0) % 2 != 0) {
            throw new IllegalArgumentException("Unexpected ReadsPerGene entry list size: ${reads_entries?.size()}")
        }
        pairs = reads_entries.collate(2)
    }

    def entries_arg = pairs.collect { pair ->
        def sample_id = pair[0].toString()
        def reads_file = pair[1].toString()
        "--input \"${sample_id}=${reads_file}\""
    }.join(' ')

    """
    set -euo pipefail
    build_allsorts_counts.py \\
        ${entries_arg} \\
        --output allsorts_counts.csv \\
        --gtf ${params.gtf} \\
        --strand ${params.allsorts_strand}
    """
}

process ALLSORTS {
    publishDir "${params.outdir}/summary/allsorts", mode: 'copy'
    conda "${projectDir}/allsorts.yml"

    input:
    path counts_csv

    output:
    path "predictions.csv"
    path "probabilities.csv"
    path "processed_counts.csv", optional: true
    path "*.png", optional: true

    script:
    def parents_flag = (params.allsorts_parents?.toString()?.toBoolean()) ? "-p" : ""
    def ball_flag = params.allsorts_ball ? "-b ${params.allsorts_ball}" : ""
    def save_counts_flag = (params.allsorts_save_processed_counts?.toString()?.toBoolean()) ? "-counts" : ""
    def model_flag = params.allsorts_model_dir ? "-model_dir ${params.allsorts_model_dir}" : ""

    """
    set -euo pipefail

    # Patch umap layouts.py: numba int32 vs int64 incompatibility on 64-bit
    UMAP_DIR=\$(python -c "import importlib.util; spec=importlib.util.find_spec('umap'); print(spec.submodule_search_locations[0])" 2>/dev/null || true)
    if [ -n "\$UMAP_DIR" ] && grep -q 'numba.types.int32' "\$UMAP_DIR/layouts.py" 2>/dev/null; then
        sed -i 's/numba\\.types\\.int32/numba.types.intp/g' "\$UMAP_DIR/layouts.py"
        # Clear numba cache so it recompiles with the patched types
        find "\$UMAP_DIR" -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
        echo "Patched umap int32->intp in \$UMAP_DIR/layouts.py"
    fi

    mkdir -p allsorts_out
    ALLSorts \\
        -s ${counts_csv} \\
        -d allsorts_out \\
        ${parents_flag} \\
        ${ball_flag} \\
        ${save_counts_flag} \\
        ${model_flag}

    cp allsorts_out/predictions.csv .
    cp allsorts_out/probabilities.csv .
    cp allsorts_out/processed_counts.csv . 2>/dev/null || true
    cp allsorts_out/*.png . 2>/dev/null || true
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

process FUSIONCATCHER {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    conda "${projectDir}/fusioncatcher.yml"

    input:
    tuple val(sample_id), path(r1_reads), path(r2_reads)

    output:
    tuple val(sample_id), path("fusioncatcher_fusions.txt"), emit: fusions
    path "fusioncatcher_summary.txt", optional: true
    path "*.log", optional: true
    path "final-list_candidate-fusion-genes.txt", optional: true

    script:
    // Interleave R1/R2 for FusionCatcher input
    // r1_reads and r2_reads are lists of files
    def input_files = []
    r1_reads.eachWithIndex { r1, idx ->
        input_files.add(r1)
        if (idx < r2_reads.size()) {
            input_files.add(r2_reads[idx])
        }
    }
    def input_csv = input_files.join(',')
    def fusioncatcher_db_abs = fusioncatcher_db_dir.toAbsolutePath().toString()

    """
    set -euo pipefail
    
    # Run FusionCatcher
    fusioncatcher \\
        -d "${fusioncatcher_db_abs}" \\
        -i ${input_csv} \\
        -o fusioncatcher_out \\
        -p ${task.cpus}

    # Rename output for consistency with pipeline
    cp fusioncatcher_out/final-list_candidate-fusion-genes.txt fusioncatcher_fusions.txt
    cp fusioncatcher_out/summary_candidate_fusions.txt fusioncatcher_summary.txt 2>/dev/null || true
    cp fusioncatcher_out/*.log . 2>/dev/null || true
    """
}

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline parameters
params.samplesheet = null
params.genome = null
params.gtf = null
params.star_index = null
params.outdir = "results"
params.stranded = 2
params.help = false

// Help message
if (params.help) {
    log.info"""
    RNA-seq Count Pipeline
    ======================
    Usage:
        nextflow run main.nf --samplesheet <csv> --genome <fasta> --gtf <gtf>
    
    Required:
        --samplesheet   CSV file with columns: sample_id,read1,read2,cell_line,replicate
        --genome        Reference genome FASTA file
        --gtf           Gene annotation GTF file
    
    Optional:
        --star_index    Pre-built STAR index directory
        --outdir        Output directory (default: results)
        --stranded      Strandedness: 0 (unstranded), 1 (stranded), 2 (reverse) (default: 0)
    """.stripIndent()
    exit 0
}

// Validate required parameters
if (!params.samplesheet || !params.genome || !params.gtf) {
    error "ERROR: --samplesheet, --genome, and --gtf are required!"
}

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{zip,html}"
    
    script:
    """
    fastqc -q -t ${task.cpus} ${reads}
    """
}

process STAR_INDEX {
    publishDir "${params.outdir}/star_index", mode: 'copy'
    
    input:
    path genome
    path gtf
    
    output:
    path "star_index"
    
    script:
    """
    mkdir star_index
    STAR --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles ${genome} \
         --sjdbGTFfile ${gtf} \
         --sjdbOverhang 99 \
         --runThreadN ${task.cpus}
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy', pattern: "*.bam*"
    publishDir "${params.outdir}/star_logs", mode: 'copy', pattern: "*.Log.final.out"
    publishDir "${params.outdir}/star_logs", mode: 'copy', pattern: "*.SJ.out.tab"
    
    input:
    tuple val(sample_id), path(reads)
    path index
    
    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"), emit: bam
    path "${sample_id}.Log.final.out", emit: log
    path "${sample_id}.SJ.out.tab", emit: splice_junctions
    
    script:
    """
    STAR --genomeDir ${index} \
         --readFilesIn ${reads[0]} ${reads[1]} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}. \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --runThreadN ${task.cpus} \
         --limitBAMsortRAM ${task.memory.toBytes()}
    
    samtools index ${sample_id}.Aligned.sortedByCoord.out.bam
    """
}

process FEATURECOUNTS {
    publishDir "${params.outdir}/counts", mode: 'copy'
    
    input:
    path gtf
    path bams
    
    output:
    path "counts.txt", emit: counts
    path "counts.txt.summary", emit: summary
    
    script:
    def strandedness = params.stranded
    """
    featureCounts -T ${task.cpus} \
                  -p -B -C \
                  -s ${strandedness} \
                  -a ${gtf} \
                  -o counts.txt \
                  ${bams}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path '*'
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"
    
    script:
    """
    multiqc -f .
    """
}

workflow {
    // Input validation
    log.info """
    ============================================
    RNA-seq Count Pipeline
    ============================================
    samplesheet  : ${params.samplesheet}
    genome       : ${params.genome}
    gtf          : ${params.gtf}
    star_index   : ${params.star_index ?: 'Will be generated'}
    outdir       : ${params.outdir}
    stranded     : ${params.stranded}
    ============================================
    """.stripIndent()
    

    def reads_fastqc = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample_id
            def read1 = file(row.read1, checkIfExists: true)
            def read2 = file(row.read2, checkIfExists: true)
            return tuple(sample_id, [read1, read2])
        }
    
    def reads_star = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample_id
            def read1 = file(row.read1, checkIfExists: true)
            def read2 = file(row.read2, checkIfExists: true)
            return tuple(sample_id, [read1, read2])
        }
    
    // Create genome and GTF channels
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)
    
    reads_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { tuple(it.sample_id, [file(it.read1), file(it.read2)]) }

    
    // STAR index
    if (params.star_index && file(params.star_index).exists()) {
        index_ch = Channel.fromPath(params.star_index, type: 'dir').first()
    } else {
        index_ch = STAR_INDEX(genome_ch, gtf_ch)
    }
    FASTQC(reads_ch)
    // Alignment
    STAR_ALIGN(reads_ch, index_ch)
    
    // Count features
    bams_ch = STAR_ALIGN.out.bam
        .map { sample_id, bam, bai -> bam }
        .collect()
    FEATURECOUNTS(gtf_ch, bams_ch)
    
    // MultiQC
    qc_ch = FASTQC.out
        .mix(STAR_ALIGN.out.log)
        .mix(FEATURECOUNTS.out.summary)
        .collect()
    MULTIQC(qc_ch)
}

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline parameters
params.reads = "data/*_{1,2}.fastq.gz"
params.genome = null
params.gtf = null
params.star_index = null
params.outdir = "results"
params.stranded = 0
params.help = false

// Help message
if (params.help) {
    log.info"""
    RNA-seq Count Pipeline
    ======================
    Usage:
        nextflow run main.nf --reads '<pattern>' --genome <fasta> --gtf <gtf>
    
    Required:
        --reads         Path to paired-end reads (e.g., 'data/*_{1,2}.fastq.gz')
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
if (!params.genome || !params.gtf) {
    error "ERROR: --genome and --gtf are required!"
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
    
    input:
    tuple val(sample_id), path(reads)
    path index
    
    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), emit: bam
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

// Change parameter
params.reads = null  // Remove this
params.samplesheet = "samplesheet.csv"  // Add this

// Validate
if (!params.samplesheet) {
    error "ERROR: --samplesheet is required!"
}


workflow {
    // Input validation
    log.info """
    ============================================
    RNA-seq Count Pipeline
    ============================================
    reads        : ${params.reads}
    genome       : ${params.genome}
    gtf          : ${params.gtf}
    outdir       : ${params.outdir}
    stranded     : ${params.stranded}
    ============================================
    """.stripIndent()
    
        Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample_id
            def read1 = file(row.read1, checkIfExists: true)
            def read2 = file(row.read2, checkIfExists: true)
            return tuple(sample_id, [read1, read2])
        }
        .set { reads_ch }
    
    // FastQC
    FASTQC(reads_ch)
    
    // STAR index
    if (params.star_index && file(params.star_index).exists()) {
        index_ch = Channel.fromPath(params.star_index, checkIfExists: true)
    } else {
        index_ch = STAR_INDEX(genome_ch, gtf_ch)
    }
    
    // Alignment
    STAR_ALIGN(reads_ch, index_ch)
    
    // Count features
    bams_ch = STAR_ALIGN.out.bam.map { it[1] }.collect()
    FEATURECOUNTS(gtf_ch, bams_ch)
    
    // MultiQC
    qc_ch = FASTQC.out
        .mix(STAR_ALIGN.out.log)
        .mix(FEATURECOUNTS.out.summary)
        .collect()
    MULTIQC(qc_ch)
}
#!/usr/bin/env nextflow

params.genome             = "tmp/ref.fa"
params.reads_se           = "$baseDir/input/se_reads/*.fq.gz"
params.reads_pe           = "$baseDir/input/pe_reads/*-R{1,2}.fq.gz"
params.bam                = "results/alignments"

// params.threads            = 1000
// params.cov_out            = "results/coverage"
params.env                = "env.sh"

log.info """\
          VARIANT CALLING
         =============================
         genome:                 ${params.genome}
         env:                    ${params.env}
         Single-ended reads:     ${params.reads_se}
         Pair-ended reads:       ${params.reads_pe}
         results directory:      ${params.bam}
          """
          .stripIndent()

// Output directory
bam       = file(params.bam)
mkdir_bam = bam.mkdirs()
println mkdir_bam ? "OK" : "Cannot create directory: $mkdir_bam"

// Inputs
Channel
    .fromPath(params.genome)
    .set { reference_assembly }

Channel
    .fromPath( params.reads_se )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads_se}" }
    .map { file -> tuple(file.baseName - ~/\.fq/, file) }
    .into { se_reads; se_reads_check }

se_reads_check
  .subscribe {it -> println "Foudn the following SE reads: ${it}."}

Channel
    .fromFilePairs( params.reads_pe )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads_pe}" }
    .map { file_name, file -> [ file_name, file[0], file[1] ] }
    .into { pe_reads; pe_reads_check }

pe_reads_check
  .subscribe {it -> println "Found the following PE reads: ${it}."}

ENV=file(params.env)

// Group the single-ended and pair-ended
se_reads
    .join(pe_reads, by: 0)
    .into { reads_channel; reads_channel_check }

reads_channel_check
  .subscribe {it -> println "Aligning the following samples: ${it}."}

// Align the samples

process makeAlignmentIndex {

    errorStrategy { 'retry' }
    maxRetries 3

    executor 'sge'
    penv 'smp'
    maxForks 1
    cpus 1
    time { 1.hour * task.attempt }
    clusterOptions "-l h_vmem=${12 * task.attempt}G"

    input:
    file env from ENV
    file genome from reference_assembly

    output:
    file "genome_index*" into genome_index_channel

    """
    source ${env}

    bowtie2-build $genome genome_index

    """

}

// Combine the genome index with the reads channel
reads_channel
    .combine(genome_index_channel)
    .map { it -> [ it[0], it[1], it[2], it[3], it[ 4..-1 ] ]}
    .into { index_and_reads_channel ; index_and_reads_channel_check }

index_and_reads_channel_check.println()

process alignSamplesToReference {

    publishDir params.bam, mode: 'move'

    maxForks 1
    //
    // errorStrategy { 'retry' }
    // maxRetries 3
    //
    // executor 'sge'
    // penv 'smp'
    // cpus 12
    // time { 4.hour * task.attempt }
    // clusterOptions "-l h_vmem=${1 * task.attempt}G"

    input:
    file env from ENV
    set val(id),
        file(se),
        file(pe1),
        file(pe2),
        file(_index) from index_and_reads_channel

    output:
    file "${id}.bam" into bam_channel

    script:
    """
    source ${env}

    bowtie2 \
      -p 12 \
      -x genome_index \
      -U ${se} \
      -1 ${pe1} \
      -2 ${pe2} \
      | samtools sort - \
      > ${id}.bam

    """

}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

#!/usr/bin/env nextflow

// Inputs
params.vcf                      = "results/filtered_vcf/*.vcf.gz"
params.env                      = "env.sh"
// Outputs
params.results                  = "results/filtered_vcf_decomposed"

log.info """\
          VARIANT CALLING
         ================================================
         VCF:                     ${params.vcf}
         File loading software:   ${params.env}
         Directory for output:    ${params.results}
          """
          .stripIndent()


// Output directory
mkdir_res = file(params.results).mkdirs()
println mkdir_res ? "OK" : "Cannot create directory: ${params.results}"

// --------------------------------------------------------------------------- #
// Input: VCF file
// --------------------------------------------------------------------------- #

Channel
  .fromPath(params.vcf)
  .map { file ->
    if ( !(file =~ /gz/ )) { exit 1, "You have to supply a compressed VCF file: ${file}"}
    if ( !file.exists() )  { exit 1, "VCF file does not exist: ${file}"}
    tuple(file.baseName - ~/\.vcf/, file) }
  .set { vcf_files }

// --------------------------------------------------------------------------- #
// Complex  variants: divide to "primitive" genotypes:
//      e.g. single variant ATA/CTG becomes two variants, A/C and A/G
// --------------------------------------------------------------------------- #

process decomposeComplexVariantsIntoSnps {

  publishDir params.results, mode: 'copy'

  executor 'sge'
  penv 'smp'
  maxForks 25
  cpus 1
  time {30.minute * ((2 * task.attempt) -  1)}
  clusterOptions "-l h_vmem=4G"
  maxRetries 3

  input:
  file env from file(params.env)
  set val(vcf_name), file('vcf') from vcf_files

  output:
  set file("${vcf_name}.decomposed.vcf.gz"),
    file("${vcf_name}.decomposed.vcf.gz.csi") into decomposed_vcf

  """
  source ${env}

  zcat ${vcf} | vcfallelicprimitives --keep-geno --keep-info | bgzip -c > ${vcf_name}.decomposed.vcf.gz
  bcftools index ${vcf_name}.decomposed.vcf.gz

  """
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

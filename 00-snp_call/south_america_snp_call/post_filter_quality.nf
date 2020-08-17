#!/usr/bin/env nextflow

// Inputs
params.vcf          = "results/filtered_vcf/*.vcf.gz"
params.env          = "env.sh"
// Outputs
params.report         = "results/post_filter_quality_reports"
params.used_scaffolds = "results/scaffold_subset/subset_genome_index"

log.info """\
          VARIANT CALLING
         =============================
         VCF:                                             ${params.vcf}
         File loading software:                           ${params.env}
         File with scaffolds used during calling:         ${params.used_scaffolds}
         Directory for report       :                     ${params.report}
          """
          .stripIndent()


// Output directory
mkdir_res = file(params.report).mkdirs()
println mkdir_res ? "OK" : "Cannot create directory: ${params.report}"

// --------------------------------------------------------------------------- #
// Input: VCF file
// --------------------------------------------------------------------------- #

Channel
  .fromPath(params.vcf)
  .map { file ->
    file_name = file.baseName - ~/\.vcf/
    tuple(file_name, file) }
  .into { vcf_files; vcf_files_check }

vcf_files_check
  .subscribe {it -> println "Reporting on the following VCF files: ${it}."}

// Decompress the VCF File

process decompressVcfFile {

  input:
  file env from file(params.env)
  set vcf_name, file('vcf') from vcf_files

  output:
  set vcf_name, file("${vcf_name}") into vcf_decompressed

  """
  zcat $vcf > ${vcf_name}
  """
}



// --------------------------------------------------------------------------- #
// Get specific fields from VCF file
// --------------------------------------------------------------------------- #

vcf_decompressed.into { vcf_gt ; vcf_ao ; vcf_ro ; vcf_qual ; vcf_combination}


process getGT {

  input:
  file env from file(params.env)
  set vcf_name, file('vcf') from vcf_gt

  output:
  set vcf_name, file('meas') into gt

  """
  bcftools query -f '[%GT ]\\n' ${vcf} \
    > meas
  """
}

process getAO {

  input:
  file env from file(params.env)
  set vcf_name, file('vcf') from vcf_ao

  output:
  set vcf_name, file('meas') into ao

  """
  source ${env}
  bcftools query -f '[%AO ]\\n' ${vcf} \
    > meas
  """
}

process getRO {

  input:
  file env from file(params.env)
  set vcf_name, file('vcf') from vcf_ro

  output:
  set vcf_name, file('meas') into ro

  """
  source ${env}
  bcftools query -f '[%RO ]\\n' ${vcf} \
    > meas
  """
}

process getQual {

  input:
  file env from file(params.env)
  set vcf_name, file('vcf') from vcf_qual

  output:
  set vcf_name, file('meas') into qual

  """
  source ${env}
  bcftools query -f '%QUAL\\n' ${vcf} \
    > meas
  """
}

// --------------------------------------------------------------------------- #
// Collect all files into a single channel...
// --------------------------------------------------------------------------- #

combined_files = vcf_combination
  .combine(gt, by: 0)
  .combine(ao, by: 0)
  .combine(ro, by: 0)
  .combine(qual, by: 0)
  //.subscribe { println it }

combined_files.into { combined_files; combined_files_print}

combined_files_print.subscribe { it ->
  println "The following combination will be sent to rmarkdown:\n ${it}" }

// --------------------------------------------------------------------------- #
// Report: Rmarkdown script called using ezknitr
// --------------------------------------------------------------------------- #

process generateReport {
  publishDir params.report, mode: 'copy'

  executor 'sge'
  penv 'smp'
  maxForks 25
  cpus 1
  time {2.hour * ((2 * task.attempt) -  1)}
  clusterOptions "-l h_vmem=56G"
  maxRetries 3

  input:
  file env from file(params.env)
  file 'used_scaffolds' from file(params.used_scaffolds)
  set val(vcf_name),
      file('vcf'),
      file('gt'),
      file('ao'),
      file('ro'),
      file('qual') from combined_files

  output:
  set vcf_name, file("${vcf_name}/") into report_files

  """
  source ${env}
  # get sample names
  bcftools query -l ${vcf} > samples

  # run rmarkdown script
  Rscript '$baseDir/run_postfilter_qual_report.R' \
    $baseDir/postfilter_qual_report.Rmd \
    ${gt} \
    ${ao} \
    ${ro} \
    ${qual} \
    ${used_scaffolds} \
    samples \
    ${vcf_name}/

  """
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

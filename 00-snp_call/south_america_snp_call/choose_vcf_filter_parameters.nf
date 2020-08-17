#!/usr/bin/env nextflow

// Inputs
params.vcf          = "results/prefilter/q25_1_1_0.2.vcf.gz"
params.q_pre        = 25
params.sample       = "tmp/choose_filter_samples/*.samples"
params.env          = "env.sh"
// Outputs
params.report         = "results/choose_filter_report"
params.used_scaffolds = "results/scaffold_subset/subset_genome_index"

log.info """\
          VARIANT CALLING
         =============================
         VCF:                                             ${params.vcf}
         QUAL pre-filter:                                 ${params.q_pre}
         Files listing samples:                           ${params.sample}
         File loading software:                           ${params.env}
         File with scaffolds used during calling:         ${params.used_scaffolds}
         Directory for report       :                     ${params.report}
          """
          .stripIndent()


// Output directory
mkdir_res = file(params.report).mkdirs()
println mkdir_res ? "OK" : "Cannot create directory: ${params.report}"

// Each subset of samples
Channel
  .fromPath(params.sample)
  .map { file -> tuple(file.baseName, file) }
  .into { sample_files; sample_files_check }

sample_files_check
  .subscribe {it -> println "Following sample files used: ${it}."}

// --------------------------------------------------------------------------- #
// Prefilter the VCF file
// --------------------------------------------------------------------------- #

Channel
  .fromPath(params.vcf)
  .set { vcf_files }

Channel
  .from(params.q_pre)
  .combine(vcf_files)
  .set { q_vcf_combination }

process preFilterVcfByQual {

  maxRetries 2

  executor 'sge'
  penv 'smp'
  maxForks 25
  cpus 1
  time {6.hour * ((2 * task.attempt) -  1)}
  clusterOptions "-l h_vmem=${32 * task.attempt}G"

  input:
  file env from file(params.env)
  set val(q), file('vcf') from q_vcf_combination

  output:
  set q, file("out.vcf") into q_prefilteredvcf_combination

  """
  source ${env}
  bcftools index ${vcf}
  bcftools view --min-ac 1:minor --exclude 'QUAL<${q}' ${vcf} > out.vcf
  """
}

// Make different VCFs, each with a different subset of samples
//      e.g. all_samples, bb, lb

q_prefilteredvcf_combination
  .combine(sample_files)
  .set { q_prefilteredvcf_sample_combination }

process vcfWithSubsetOfSamples {

  maxRetries 2

  maxForks 25

  executor 'sge'
  penv 'smp'
  cpus 1
  time {1.hour * task.attempt}
  clusterOptions "-l h_vmem=${2 * task.attempt}G"

  module 'bcftools/1.4:htslib/1.6'

  input:
  file env from file(params.env)
  set q, file("in_vcf"), val(sample_set), file('include_samples') from q_prefilteredvcf_sample_combination


  output:
  set sample_set, q, file("out.vcf") into subsampled_q_vcf_combination

  """
  source ${env}
  bcftools view ${in_vcf} --samples-file ${include_samples} \
     | bcftools view --min-ac 1:minor - \
     | bcftools view -M2 -m2 -v snps - \
     > out.vcf
  """
}

// --------------------------------------------------------------------------- #
// Get specific fields from VCF file
// --------------------------------------------------------------------------- #

subsampled_q_vcf_combination.into { vcf_gt ; vcf_ao ; vcf_ro ; vcf_qual ; subsampled_q_vcf_combination}


process getGT {

  input:
  file env from file(params.env)
  set sample, q, file('vcf') from vcf_gt

  output:
  set sample, q, file('meas') into gt

  """
  bcftools query -f '[%GT ]\\n' ${vcf} \
    > meas
  """
}

process getAO {

  input:
  file env from file(params.env)
  set sample, q, file('vcf') from vcf_ao

  output:
  set sample, q, file('meas') into ao

  """
  source ${env}
  bcftools query -f '[%AO ]\\n' ${vcf} \
    > meas
  """
}

process getRO {

  input:
  file env from file(params.env)
  set sample, q, file('vcf') from vcf_ro

  output:
  set sample, q, file('meas') into ro

  """
  source ${env}
  bcftools query -f '[%RO ]\\n' ${vcf} \
    > meas
  """
}

process getQual {

  input:
  file env from file(params.env)
  set sample, q, file('vcf') from vcf_qual

  output:
  set sample, q, file('meas') into qual

  """
  source ${env}
  bcftools query -f '%QUAL\\n' ${vcf} \
    > meas
  """
}

// --------------------------------------------------------------------------- #
// Collect all files into a single channel...
// --------------------------------------------------------------------------- #

combined_files = subsampled_q_vcf_combination
  .combine(gt, by: [0, 1])
  .combine(ao, by: [0, 1])
  .combine(ro, by: [0, 1])
  .combine(qual, by: [0, 1])
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
  set val(sample),
      val(q),
      file('vcf'),
      file('gt'),
      file('ao'),
      file('ro'),
      file('qual') from combined_files

  output:
  set sample, q, file("${sample}_${q}/") into report_files

  """
  source ${env}
  # get sample names
  bcftools query -l ${vcf} > samples

  # run rmarkdown script
  Rscript '$baseDir/run_report.R' \
    $baseDir/report.Rmd \
    ${gt} \
    ${ao} \
    ${ro} \
    ${qual} \
    ${used_scaffolds} \
    ${q} \
    samples \
    ${sample}_${q}/

  """
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

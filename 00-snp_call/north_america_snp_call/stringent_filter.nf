#!/usr/bin/env nextflow

// Inputs
params.vcf          = "results/prefilter/q25_1_1_0.2.vcf.gz"
params.q_pre        = 25
params.sample       = "tmp/choose_filter_samples/*.samples"
params.env          = "env.sh"
params.mean_dp      = 16
params.min_cov      = 1
params.called_prop  = 0.6
// Outputs
params.report       = "results/summary_table"
params.filtered_out = "results/filtered_vcf"

log.info """\
          VARIANT CALLING
         =============================
         VCF:                                             ${params.vcf}
         QUAL pre-filter:                                 ${params.q_pre}
         Files listing samples:                           ${params.sample}
         Maximum mean coverage:                           ${params.mean_dp}
         Minimum coverage:                                ${params.min_cov}
         Minimum proportion of reads covering called
           allele versus non-called allele:               ${params.called_prop}
         File loading software:                           ${params.env}
         Directory for filtered VCF:                      ${params.filtered_out}
         Directory for summary table:                     ${params.report}
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

process getVcfFilterColumn {
  publishDir params.report, pattern: "*_summary.txt", mode: 'copy'

  executor 'sge'
  penv 'smp'
  maxForks 25
  cpus 1
  time {2.hour * ((2 * task.attempt) -  1)}
  clusterOptions "-l h_vmem=${32 * task.attempt}G"
  maxRetries 3

  input:
  file env from file(params.env)
  set val(sample),
      val(q),
      file('vcf'),
      file('gt'),
      file('ao'),
      file('ro'),
      file('qual') from combined_files
  val mean_dp from params.mean_dp
  val min_cov from params.min_cov
  val called_prop from params.called_prop

  output:
  set sample, q, file("${sample}_${q}_summary.txt") into summary_table
  set sample, q, mean_dp, min_cov, called_prop,
    file("${sample}_${q}_filter"), file("${vcf}") into vcf_filter_column

  """
  source ${env}

  # run rmarkdown script
  Rscript '$baseDir/make_filter_column.R' \
    --gt ${gt} \
    --ao ${ao} \
    --ro ${ro} \
    --mean_dp ${mean_dp} \
    --min_cov ${min_cov} \
    --called_prop ${called_prop} \
    --outsummary ${sample}_${q}_summary.txt \
    --outfile ${sample}_${q}_filter
  """
}


process addFilterColumnToVcf {

  publishDir params.filtered_out, mode: 'copy'

  executor 'sge'
  penv 'smp'
  maxForks 25
  cpus 1
  time {2.hour * ((2 * task.attempt) -  1)}
  clusterOptions "-l h_vmem=${2 * task.attempt}G"
  maxRetries 2

  input:
  file env from file(params.env)
  set val(sample),
      val(q),
      val(mean_dp),
      val(min_cov),
      val(called_prop),
      file('filter_column'),
      file('vcf') from vcf_filter_column

  output:
  set val(sample),
      val(q),
      file("${sample}_q${q}_${mean_dp}_${min_cov}_${called_prop}_filtered.vcf.gz.tbi"),
      file("${sample}_q${q}_${mean_dp}_${min_cov}_${called_prop}_filtered.vcf.gz") into filtered_vcf

  """
  # Load right software
  source ${env}

  # Remove the header from the VCF file
  grep -v "#" ${vcf} > no_header

  # Make sure the column and the VCF have the same number of lines
  wc -l no_header        | cut -f 1 -d ' ' > vcf.lc
  wc -l ${filter_column} | cut -f 1 -d ' ' > filter_column.lc


  if diff vcf.lc filter_column.lc >/dev/null ; then
    echo Filtering ${vcf}
  else
    echo "The vcf and the filter column have different sizes."
    exit 1
  fi

  # Paste the VCF and the filter column together
  paste ${filter_column} no_header > column_vcf

  # Make a new VCF, using only the lines that say PASS, concatenated onto the
    # header, and making sure the additional column with the filter is removed.
  bcftools view -h ${vcf}            > header
  grep -a "PASS" column_vcf | cut -f 2- >> header

  # If the previous runs out of memory, it doesn't throw an error...
  if grep -q "Binary file column_vcf matches" header; then
    exit 2
  fi

  # Zip up and index, with an informative name
  cat header | bgzip -c \
    > ${sample}_q${q}_${mean_dp}_${min_cov}_${called_prop}_filtered.vcf.gz
  tabix -p vcf ${sample}_q${q}_${mean_dp}_${min_cov}_${called_prop}_filtered.vcf.gz

  """

}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

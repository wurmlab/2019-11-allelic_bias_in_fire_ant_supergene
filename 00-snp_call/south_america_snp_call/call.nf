#!/usr/bin/env nextflow

params.genome             = "input/ref.fa"
params.ploidy             = 1
params.alternate_count    = 1
params.alternate_fraction = 0.2
params.joined_bam         = "tmp/joined.bam"
params.joined_bam_index   = "tmp/joined.bam.bai"
params.q_prefilter        = 25
params.out_raw            = "results/raw"
params.out_prefilter      = "results/prefilter"
params.used_scaffolds     = "results/scaffold_subset"
params.scaffs_per_chunk   = 10
params.min_scaffold_size  = 5000
// params.threads            = 1000
// params.cov_out            = "results/coverage"
params.env                = "env.sh"
params.sample             = "argentina_samples"

log.info """\
          VARIANT CALLING
         =============================
         genome:                 ${params.genome}
         ploidy:                 ${params.ploidy}
         min scaffold size:      ${params.min_scaffold_size}
         min alternate count:    ${params.alternate_count}
         min alternate fraction: ${params.alternate_fraction}
         bam file:               ${params.joined_bam}
         prefilter:              ${params.q_prefilter}
         raw vcf dir:            ${params.out_raw}
         prefilter vcf dir:      ${params.out_prefilter}
          """
          .stripIndent()



// Output directory
prefilter = file(params.out_prefilter)
raw       = file(params.out_raw)
mkdir_raw = raw.mkdirs()
println mkdir_raw ? "OK" : "Cannot create directory: $mkdir_raw"
mkdir_res = prefilter.mkdirs()
println mkdir_res ? "OK" : "Cannot create directory: $prefilter"
mkdir_sca = file(params.used_scaffolds).mkdirs()
println mkdir_sca ? "OK" : "Cannot create directory: ${file(params.used_scaffolds)}"

// Inputs
   // Each channel can be used only once
bam_file_split     = Channel.fromPath(params.joined_bam)
bam_file_freebayes = Channel.fromPath(params.joined_bam)

bam_index          = Channel.fromPath(params.joined_bam_index)

genome_file_index      = Channel.fromPath(params.genome)
genome_file_freebayes  = Channel.fromPath(params.genome)

// Index fasta
process indexFasta {

  errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }

  maxRetries 3

  executor 'sge'
  penv 'smp'
  maxForks 20
  cpus 1
  time {20.minute * ((2 * task.attempt) -  1)}
  clusterOptions "-l h_vmem=${2 * task.attempt}G"

  input:
  file env from file(params.env)
  file genome from genome_file_index

  output:
  file "${genome}.fai" into fai_original

  """
  source ${env}
  samtools faidx ${genome}
  """
}

fai_original.into { fai_to_remove_tiny_scaffolds; genome_index }

process removeTinyScaffoldsFromReferenceIndex {

  publishDir params.used_scaffolds, mode: 'copy'

  executor 'sge'
  penv 'smp'
  maxForks 20
  cpus 1
  time {1.minute}
  clusterOptions "-l h_vmem=1G"


  input:
  file genome_index from fai_to_remove_tiny_scaffolds
  val min_scaffold_size from params.min_scaffold_size

  output:
  file 'subset_genome_index' into fai_tiny_scaffs_removed

  """
  cat ${genome_index} \
    | awk '\$2 > ${min_scaffold_size} {print}' \
    > subset_genome_index
  """
}

// Split the fasta into chunks of X scaffolds to run freebayes in parallel
fai_chunks = fai_tiny_scaffs_removed.splitText(by: params.scaffs_per_chunk, file: true)

// Count number of chunks
fai_chunks
  .into { fai_chunks ; fai_chunks_count }

fai_chunks_count
  .count()
  .subscribe { it ->
    println "After removing scaffolds smaller than ${params.min_scaffold_size}, we end up with ${it} chunks of ${params.scaffs_per_chunk} scaffolds." }

process faiToBed {
  input:
  each file(region) from fai_chunks

  output:
  file 'region_*' into bed_chunks

  """
  cat $region | awk '{print \$1 " 0 " \$2}' > region_
  """
}
// // Divide the reference into chunks of similar coverage
// process totalCoverage {
//   executor 'sge'
//   clusterOptions='-pe smp 1 -l h_vmem=32G,h_rt=48:0:0'
//   publishDir params.cov_out
//   module 'bamtools/2.4.1'
//
//   input:
//   file bam from bam_file_split
//
//   output:
//   file 'regions.bed' into regions_bed
//
//   """
//   bamtools coverage -in $bam > regions.bed
//   """
// }
//
// // Divide the reference into chunks of similar coverage
// process regionSplit {
//   executor 'sge'
//   clusterOptions='-pe smp 1 -l h_vmem=128G,h_rt=24:0:0'
//
//   input:
//   file 'regions.bed' from regions_bed
//   val threads from params.threads
//   file genome_fai from genome_index
//
//   output:
//   file 'region.split' into region_split
//
//   """
//   cat regions.bed \
//     | /data/home/btw749/software/freebayes/scripts/coverage_to_regions.py $genome_fai $threads \
//     > region.split
//   """
// }

// // Split region file into lines, giving each a task number
// region_split.into{ regions_freebayes; regions_concat }
//
// regions_freebayes_split = regions_freebayes.splitText{line ->
//  lsplit = line.split(":")[1]
//  coords = lsplit.split("-").collect { it as int }
//  if (coords[1] > (coords[0] + 1)) {
//    return line - '\n' // removes trailing '\n'
//   }
//  }.filter{ it != null}
//
// regions_concat_split = regions_concat.splitText{line ->
//   lsplit = line.split(":")[1]
//   coords = lsplit.split("-").collect { it as int }
//   if (coords[1] > (coords[0] + 1)) {
//     return line - '\n' // removes trailing '\n'
//    }
//   }.filter{ it != null}

// runs freebayes
process runFreebayes {

  errorStrategy { 'retry' }
  maxRetries 3

  executor 'sge'
  penv 'smp'
  maxForks 85
  cpus 1
  time { 1.hour * ((task.attempt * 48) - 24) }
  clusterOptions "-l h_vmem=${16 * task.attempt}G"

  input:
  file env from file(params.env)
  file genome from genome_file_freebayes
  file fai from genome_index
  file bam from bam_file_freebayes
  file bai from bam_index
  file samples from file(params.sample)
  each ploidy from params.ploidy
  each count from params.alternate_count
  each fraction from params.alternate_fraction
  each file('region') from bed_chunks

  output:
  set ploidy, count, fraction, region, 'vcf' into vcf_files

  """
  source ${env}

  freebayes \
       --ploidy ${ploidy} \
       --min-alternate-count ${count} \
       --min-alternate-fraction ${fraction} \
       --targets ${region} \
       --fasta-reference ${genome} \
       --samples ${samples} \
       ${bam} \
       > vcf
  """
}

vcf_files_for_concat = vcf_files
       .groupTuple(by: [0,1,2])
       .map {field ->
         param_set = field[0,1,2].join("_")
         return [param_set, field[4]]}

// ploidy   = Channel.from(params.ploidy)
// count    = Channel.from(params.alternate_count)
// fraction = Channel.from(params.alternate_fraction)
//
// // Parameter set
// parameter_set = ploidy
//     .combine(count)
//     .combine(fraction)
//     .map { set -> set.join("_") }

// Creates tuples grouping all the regions for a given set of parameters
//     [${ploidy}_${count}_${fraction}, reg1_${ploidy}_${count}_${fraction}.vcf, reg1_${ploidy}_${count}_${fraction}.vcf, ...]
// region_vcf_files = vcf_files.groupTuple()

// Concatenates all regions for each set of parameters
process ConcatenateVCFs {

  errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
  maxRetries 3

  executor 'sge'
  penv 'smp'
  maxForks 20
  cpus 1
  time { 1.hour * task.attempt }
  clusterOptions "-l h_vmem=${1 * task.attempt}G"

  publishDir params.out_raw, mode: 'copy'

  input:
  file env from file(params.env)
  set val(parameter_set), file('*.vcf') from vcf_files_for_concat

  output:
  file "${parameter_set}.vcf"  into concat_vcf

  """
  source ${env}

  vcf-concat *.vcf > ${parameter_set}.vcf
  """
}

// Pre filter VCF file
process preFilterVcf {
  publishDir params.out_prefilter, mode: 'copy'

  errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
  maxRetries 2

  executor 'sge'
  penv 'smp'
   25
  cpus 1
  time {12.hour * ((2 * task.attempt) -  1)}
  clusterOptions "-l h_vmem=${1 * task.attempt}G"

  input:
  file env from file(params.env)
  file vcf from concat_vcf
  val q from params.q_prefilter

  output:
  file "q${q}_${vcf}" into prefiltered_vcf

  """
  source ${env}

  bcftools view --min-ac 1:minor --exclude 'QUAL<${q}' ${vcf} > q${q}_${vcf}
  """

}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

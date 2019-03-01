#!/usr/bin/env nextflow

input_dir = params.in //input directory containing all fastq files
input_data = "${input_dir}/*.fastq.gz"
fastq_files = Channel.fromPath(input_data).map{file -> tuple(file.simpleName, file)}

output_dir = params.out //output directory containing results for all files
threads = params.threads

gdc_ref_gen_idx_dir = params.ref_gen_idx_dir
gdc_ref_gen = params.ref_gen
gdc_gtf = params.gtf
rnaseqc_gtf = params.rnaseqc_gtf

path_to_data = "\$HOME/.nextflow/assets/Reda94/GDC-RNAseq-pipeline/data/"

//Star alignment:

process starAlignment {

  label 'parallel'

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(fq) from fastq_files

  output:
  set sample_name, file("${sample_name}_star_alignment/*Aligned.sortedByCoord.out.bam") into star_alignment_results
  set sample_name, file("${sample_name}_star_alignment/*Aligned.sortedByCoord.out.bam"), file("${sample_name}_star_alignment/*Aligned.sortedByCoord.out.bam.bai") into rseqc_input

  """
  module load STAR/2.5.2a-foss-2016b
  module load SAMtools/1.3.1-foss-2016b
  mkdir ${sample_name}_star_alignment

  STAR \\
  --runThreadN $threads \\
  --genomeDir $gdc_ref_gen_idx_dir \\
  --readFilesIn $fq \\
  --readFilesCommand zcat \\
  --outFilterType BySJout \\
  --outFilterMultimapNmax 20 \\
  --alignSJoverhangMin 8 \\
  --alignSJDBoverhangMin 1 \\
  --outFilterMismatchNmax 999 \\
  --outFilterMismatchNoverLmax 0.6 \\
  --alignIntronMin 20 \\
  --alignIntronMax 1000000 \\
  --alignMatesGapMax 1000000 \\
  --outSAMattributes NH HI NM MD \\
  --outSAMtype BAM SortedByCoordinate \\
  --outFileNamePrefix ./${sample_name}_star_alignment/${sample_name}

  samtools index ./${sample_name}_star_alignment/${sample_name}Aligned.sortedByCoord.out.bam
  """
}

//Raw read counting:

process rawReadCount {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(bam) from star_alignment_results

  output:
  set sample_name, file("${sample_name}_raw_read_counts/*") into raw_counts_results

  """
  module load HTSeq/0.6.1p1-foss-2016b-Python-2.7.12
  mkdir ${sample_name}_raw_read_counts

  htseq-count \\
  -m intersection-nonempty \\
  -s yes \\
  -f bam \\
  -r pos \\
  -i gene_id \\
  $bam \\
  $gdc_gtf \\
  > ./${sample_name}_raw_read_counts/${sample_name}raw_counts.txt
  """
}

process FPKM_TPM {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(raw_counts) from raw_counts_results

  output:
  set sample_name, file("${sample_name}_FPKM_TPM/*.FPKM.txt"), file("${sample_name}_FPKM_TPM/*.TPM.txt") into FPKM_TPM_results

  """
  module load Python/3.6.6-foss-2018b
  mkdir ${sample_name}_FPKM_TPM

  FPKM_script.py $raw_counts ./${sample_name}_FPKM_TPM/${sample_name} $path_to_data
  """
}

process RNASeQC {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(bam), file(bai) from rseqc_input

  output:
  set sample_name, file("${sample_name}_QC/*") into qc_results

  """
  module load RNA-SeQC/1.1.8-Java-1.7.0_80
  mkdir ${sample_name}_QC

  java -jar \${EBROOTRNAMINSEQC}/RNA-SeQC_v1.1.8.jar -o ${sample_name}_QC -r $gdc_ref_gen -s \"${sample_name}|$bam|notes\" -t $rnaseqc_gtf
  """
}
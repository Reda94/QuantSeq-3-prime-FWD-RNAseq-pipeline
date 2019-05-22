#!/usr/bin/env nextflow

input_dir = params.in //input directory containing all fastq files
input_data = "${input_dir}/*.fastq.gz"
fastq_files = Channel.fromPath(input_data).map{file -> tuple(file.simpleName, file)}
fastq_files_copy = Channel.fromPath(input_data).map{file -> tuple(file.simpleName, file)}

output_dir = params.out //output directory containing results for all files
threads = params.threads

gdc_ref_gen_idx_dir = params.ref_gen_idx_dir
gdc_ref_gen = params.ref_gen
gdc_gtf = params.gtf
rnaseqc_gtf = params.rnaseqc_gtf
refs_trimming = params.refs_trimming

path_to_data = "\$HOME/.nextflow/assets/Reda94/GDC-RNAseq-pipeline/data/"

//Trimming

process trimmingReads {

  label 'parallel'

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(fq) from fastq_files

  output:
  set sample_name, file("${sample_name}_trimming/*.clean.fastq.gz") into trimming_results
  """
  module load BBMap/36.20-foss-2016b-Java-1.8.0_92
  mkdir ${sample_name}_trimming

  bbduk.sh \\
  in=$fq \\
  out=./${sample_name}_trimming/${sample_name}.clean.fastq.gz \\
  ref=$refs_trimming \\
  k=13 \\
  ktrim=r \\
  useshortkmers=t \\
  mink=5 \\
  qtrim=r \\
  trimq=10 \\
  minlength=20 \\
  threads=$threads
  """
}

//Star alignment:

process starAlignment {

  label 'parallel'

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(fq_trimmed) from trimming_results

  output:
  set sample_name, file("${sample_name}_star_alignment/*Aligned.sortedByCoord.out.bam") into star_alignment_results, sar_for_trm_reads, sar_for_algnd_reads, sar_for_ualgnd_reads
  set sample_name, file("${sample_name}_star_alignment/*Log.final.out") into fb_log_final_out
  set sample_name, file("${sample_name}_star_alignment/*.out") into fb_out
  set sample_name, file("${sample_name}_star_alignment/*Log.out") into fb_log_out

  """
  module load STAR/2.5.2a-foss-2016b
  module load SAMtools/1.3.1-foss-2016b
  mkdir ${sample_name}_star_alignment

  STAR \\
  --runThreadN $threads \\
  --genomeDir $gdc_ref_gen_idx_dir \\
  --readFilesIn $fq_trimmed \\
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
  --outSAMunmapped Within \\
  --outFileNamePrefix ./${sample_name}_star_alignment/${sample_name} \\
  --outSAMattrRGline ID:${sample_name} SM:${sample_name}

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
  module load Pysam/0.9.1.4-foss-2016b-Python-2.7.12
  mkdir ${sample_name}_raw_read_counts

  htseq-count \\
  -m intersection-nonempty \\
  -s yes \\
  -f bam \\
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
  file("${sample_name}_FPKM_TPM/*.FPKM.txt") into FPKM_results
  file("${sample_name}_FPKM_TPM/*.TPM.txt") into TPM_results

  """
  module load Python/3.6.6-foss-2018b
  mkdir ${sample_name}_FPKM_TPM

  FPKM_script.py $raw_counts ./${sample_name}_FPKM_TPM/${sample_name} $path_to_data
  """
}

process stats_total_reads {

  label 'stats'

  input:
  set sample_name, file(fq) from fastq_files_copy

  output:
  file("${sample_name}_stats_total_reads/*.txt") into stats_total_reads_results
  """
  mkdir ${sample_name}_stats_total_reads

  a=\"${sample_name}\t\"
  b=\$((\$(zcat $fq|wc -l)/4))
  c="\$a\$b"
  echo \$c > ./${sample_name}_stats_total_reads/${sample_name}.txt
  """
}

process stats_trm_reads {

  label 'stats'

  input:
  set sample_name, file(bam) from sar_for_trm_reads

  output:
  file("${sample_name}_stats_trm_reads/*.txt") into stats_trm_reads_results
  """
  module load SAMtools/1.3.1-foss-2016b
  mkdir ${sample_name}_stats_trm_reads

  a=\"${sample_name}\t\"
  b=\$(samtools view -c -F 256 $bam)
  c="\$a\$b"
  echo \$c > ./${sample_name}_stats_trm_reads/${sample_name}.txt
  """
}

process stats_algnd_reads {

  label 'stats'

  input:
  set sample_name, file(bam) from sar_for_algnd_reads

  output:
  file("${sample_name}_stats_algnd_reads/*.txt") into stats_algnd_reads_results
  """
  module load SAMtools/1.3.1-foss-2016b
  mkdir ${sample_name}_stats_algnd_reads

  a=\"${sample_name}\t\"
  b=\$(samtools view -c -F 260 $bam)
  c="\$a\$b"
  echo \$c > ./${sample_name}_stats_algnd_reads/${sample_name}.txt
  """
}

process stats_ualgnd_reads {

  label 'stats'

  input:
  set sample_name, file(bam) from sar_for_ualgnd_reads

  output:
  file("${sample_name}_stats_ualgnd_reads/*.txt") into stats_ualgnd_reads_results
  """
  module load SAMtools/1.3.1-foss-2016b
  mkdir ${sample_name}_stats_ualgnd_reads

  a=\"${sample_name}\t\"
  b=\$(samtools view -c -q 255 -F 260 $bam)
  c="\$a\$b"
  echo \$c > ./${sample_name}_stats_ualgnd_reads/${sample_name}.txt
  """
}

stats_total_reads_results.collectFile(name: "${output_dir}/total_reads.txt")
stats_trm_reads_results.collectFile(name: "${output_dir}/trimmed_reads.txt")
stats_algnd_reads_results.collectFile(name: "${output_dir}/algnd_reads.txt")
stats_ualgnd_reads_results.collectFile(name: "${output_dir}/ualgnd_reads.txt")
FPKM_results.collectFile(name: "${output_dir}/FPKM_results.txt")
TPM_results.collectFile(name: "${output_dir}/TPM_results.txt")

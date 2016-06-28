#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --array=1-4%4
#SBATCH --job-name=bismark

if [ $# -ne 4 ]; then
    echo "usage: sbatch bismark.sbatch <genome folder> <read1.fq> <min_insert> <max_insert>"
    echo "Align pe reads with bismark. Bowtie version is 2."
    exit -1; 
fi

genome_dir=$1 # Index directory containing *.fa reference, Bisulfite_Genome/* index made with chosen bowtie version.
read1_fq=$2   # fastq of of paired-end read1
min_insert=$3 # Minimum insert parameter (bismark -I nnn)
max_insert=$4 # Maximum insert parameter (bismark -X nnn)

# SLURM ARRAY 4 lanes per one job

read1_fq=${read1_fq%_L001_R1_001.fastq.gz}
read1_fq=${read1_fq}_L00${SLURM_ARRAY_TASK_ID}_R1_001.fastq.gz

# Note bismark uses 2 per multicore so --multi ncpus/2
if [ $SLURM_CPUS_PER_TASK -gt 1 ]; then
    ncpus=`expr $SLURM_CPUS_PER_TASK / 2`
else
    ncpus=$SLURM_CPUS_PER_TASK
fi

ml trimgalore
ml bismark

read2_fq=${read1_fq/R1_001/R2_001}
base=${read1_fq##*/}
base=${base%_R1_001.fastq.gz}

echo "-- Expect to create '${base}_bismark_pe.bam' and '${base}_bismark_se.bam'"
echo "-- SLURM JOB ID: $SLURM_JOB_ID"
echo ""

# Trimming reads
# Note first 6 bp of EpiGnome libraries should be trimmed.
echo "-- Trimming reads..."
set -x
mkdir -p $base
trim_galore -q 20 --phred33 --fastqc --fastqc_args "-o $base --noextract --nogroup" \
	    --illumina --stringency 1 -e 0.1 --dont_gzip --length 20 \
	    -o $base --clip_R1 6 --clip_R2 6 \
	    --paired --retain_unpaired -r1 21 -r2 21 $read1_fq $read2_fq &>${base}/trim_galore.log
rename "s/_001_val_\d//g" ${base}/${base}*_val_*.fq*
rename "s/_001_val_\d/_trimmed/g" ${base}/${base}*_fastqc.* 
rename "s/001.fastq.gz_//g" ${base}/${base}*_trimming_report.txt
set +x

read1_fq=${base}/${base}_R1.fq
read2_fq=${base}/${base}_R2.fq
echo ""

# Mapping with bismark/bowtie2
# Note --bowtie2 and -p $nthreads are both SLOWER than single threaded bowtie1
echo "-- Mapping to reference with bismark/bowtie2..."
set -x
srun bismark -q -I $min_insert -X $max_insert --multicore $ncpus \
	     --un --ambiguous -o $base --temp_dir $base --nucleotide_coverage \
	     --bowtie2 -N 1 -L 28 --dovetail --score_min L,0,-0.6 \
             $genome_dir -1 $read1_fq -2 $read2_fq &>${base}/bismark_pe.log
samtools sort -O sam -T ${base}/${base} -@ $ncpus ${base}/${base}_R1_bismark_bt2_pe.bam |
  sed "s/_1:N:0:[ACGTN]*//g" - | samtools view -b -@ $ncpus - >${base}/${base}_bismark_bt2_pe.bam
rm ${base}/${base}_R1_bismark_bt2_pe.bam
rename "s/R1_//g" ${base}/${base}_R1_bismark_bt2_*
zcat ${base}/${base}_R1.fq_unmapped_reads_1.fq.gz | 
  cat ${base}/${base}_R1_001_unpaired_1.fq - > ${base}/${base}_R1.fq
zcat ${base}/${base}_R2.fq_unmapped_reads_2.fq.gz |
  cat ${base}/${base}_R2_001_unpaired_2.fq - > ${base}/${base}_R2.fq
set +x
echo ""

# Mapping unpaired and unmapped reads with single end mode
echo "-- Mapping unpaired and unmapped reads with single end mode..."
set -x
srun bismark -q --multicore $ncpus \
             -o $base --temp_dir $base --nucleotide_coverage \
             --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 \
             $genome_dir $read1_fq &>${base}/bismark_se1.log
srun bismark -q --multicore $ncpus \
             --pbat -o $base --temp_dir $base --nucleotide_coverage \
             --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 \
             $genome_dir $read2_fq &>>${base}/bismark_se2.log
srun samtools sort -O bam -T ${base}/${base} -@ $ncpus ${base}/${base}_R1_bismark_bt2.bam >${base}/${base}_R1_bismark_bt2_se.bam
srun samtools sort -O bam -T ${base}/${base} -@ $ncpus ${base}/${base}_R2_bismark_bt2.bam >${base}/${base}_R2_bismark_bt2_se.bam
samtools merge ${base}/${base}_bismark_bt2_se.bam ${base}/${base}_R*_bismark_bt2_se.bam
set +x
echo ""

# Deleting temporary fastq and bam files
echo "-- Deleting temporary fastq and bam files..."
set -x
rm ${base}/${base}_R*_bismark_bt2*.bam
rm ${base}/${base}_R*.fq
rm ${base}/${base}_R*.fq.gz
set +x
echo ""

# Print the files generated
echo "-- The results..."
ls -l ${base}/

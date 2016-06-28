#!/usr/bin/env bash
# Author: Hyung Joo Lee

#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --job-name=methXtractor

if [ $# -ne 1 ]; then
    echo "usage v1: methXtractor.sh <bismark_pe.bam>"
    echo "Extracts methylation from paired-end bismark bam."
    exit -1; 
fi

bismark_bam=$1 # bam input aligned with bismark and bowtie version that matches genome_dir.
ncpus=$SLURM_CPUS_PER_TASK

ml bismark

base=${bismark_bam##*/}
base=${base%_bismark_bt2_pe.bam}

se_bam=${bismark_bam/_pe/_se}

sorted_pe_bam=${base}/${base}_pe.bam
sorted_se_bam=${base}/${base}_se.bam
lambda_pe_bam=lambda/${base}_pe.bam
lambda_se_bam=lambda/${base}_se.bam

echo "-- Expect to create CpG resolution methylation level file"
echo "-- Started on $(date) with SLURM JOB ID: $SLURM_JOB_ID"
echo ""

echo "-- Separating alignment to Lambda DNA..."
set -x
mkdir -p lambda/
mkdir -p ${base}/
srun samtools view -h -@ $ncpus $bismark_bam |
    tee >(awk '($3!="chrLambda" && $2!="SN:chrLambda"){OFS="\t"; print}' - |
            samtools sort -n -o $sorted_pe_bam -T ${base} -@ $ncpus - )  |
          awk '($3=="chrLambda" || $1~/^@[HP]/ || $2=="SN:chrLambda"){OFS="\t"; print}' - |
            samtools sort -n -o $lambda_pe_bam -T lambda/${base} -@ $ncpus -
srun samtools view -h -@ $ncpus $se_bam |
    tee >(awk '($3!="chrLambda" && $2!="SN:chrLambda"){OFS="\t"; print}' - |
            samtools sort -n -o $sorted_se_bam -T ${base} -@ $ncpus - )  |
          awk '($3=="chrLambda" || $1~/^@[HP]/ || $2=="SN:chrLambda"){OFS="\t"; print}' - |
            samtools sort -n -o $lambda_se_bam -T lambda/${base} -@ $ncpus -
set +x
echo ""

log_pe=lambda/${base}_pe.bismark_methylation_extractor.log
log_se=lambda/${base}_se.bismark_methylation_extractor.log

echo "-- Analyse conversion rate in Lambda DNA using $ncpus threads..."
set -x
srun bismark_methylation_extractor --paired-end --no_overlap --report \
                              -o lambda/ --mbias_only --multicore $ncpus \
                              $lambda_pe_bam &>$log_pe
srun bismark_methylation_extractor --single-end --report \
                              -o lambda/ --mbias_only --multicore $ncpus \
                              $lambda_se_bam &>$log_se
set +x
echo ""

log_pe=${base}/${base}_pe.bismark_methylation_extractor.log
log_se=${base}/${base}_se.bismark_methylation_extractor.log

echo "-- Analyse methylation in $bismark_bam and $se_bam using $ncpus threads..."
set -x
srun bismark_methylation_extractor --paired-end --no_overlap --comprehensive --report \
                              -o ${base}/ --gzip --multicore $ncpus \
                              $sorted_pe_bam &>$log_pe
srun bismark_methylation_extractor --single-end --no_overlap --comprehensive --report --no_header \
                              -o ${base}/ --gzip --multicore $ncpus \
                              $sorted_se_bam &>$log_se
set +x
echo ""

cpg_pe=CpG_context_${base}_pe.txt.gz
cpg_se=CpG_context_${base}_se.txt.gz
cpg=CpG_context_${base}.txt.gz
bedGraph=${base}.bedGraph.gz
log=${base}.bismark2bedGraph.log

echo "-- Generate bedGraph file..."
set -x
cd ${base}/
zcat $cpg_pe $cpg_se | gzip >$cpg
bismark2bedGraph --dir ./ --cutoff 1 --buffer_size=20G --scaffolds --zero_based \
                 -o $bedGraph $cpg &>$log
mv ${bedGraph}.bismark.zero.cov ${base}.bismark.zero.cov
cd ../
set +x
echo ""

echo "-- Delete temporary files..."
set -x
rm ${base}/${base}_*.bam
rm lambda/${base}_*.bam
rm ${base}/${base}.bismark.cov.gz
rm ${base}/$cpg
set +x
echo ""

echo "-- The results..."
ls -l ${base}/ lambda/*${base}*
echo ""
echo "-- Finished on $(date)"

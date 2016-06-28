#!/usr/bin/env bash
# Author: Hyung Joo Lee

#SBATCH --mem=75G
#SBATCH --array=1-8
#SBATCH --job-name=markDup

module load picard-tools
module load samtools

dir=/scratch/twlab/hyungjoo.lee/MethylC/bismark/aligned
run1=${dir}/160301_NS500153_0157/WangT-MethylC-0${SLURM_ARRAY_TASK_ID}-*_S${SLURM_ARRAY_TASK_ID}
run2=${dir}/160418_NS500153_0165/WangT_MethylC_0${SLURM_ARRAY_TASK_ID}_*_S${SLURM_ARRAY_TASK_ID}
run3=${dir}/160519_NS500153_0176/WangT_MethylC_0${SLURM_ARRAY_TASK_ID}_*_S${SLURM_ARRAY_TASK_ID}
run4=${dir}/160527_NS500153_0177/WangT_MethylC_0${SLURM_ARRAY_TASK_ID}_*_S${SLURM_ARRAY_TASK_ID}

if [ $SLURM_ARRAY_TASK_ID -ge 2 ] && [ $SLURM_ARRAY_TASK_ID -le 4 ]
then
    r1l1=$(ls ${run1}_L001/*_L001_bismark_bt2_se.bam)
    r1l2=$(ls ${run1}_L002/*_L002_bismark_bt2_se.bam)
    r1l3=$(ls ${run1}_L003/*_L003_bismark_bt2_se.bam)
    r1l4=$(ls ${run1}_L004/*_L004_bismark_bt2_se.bam)
fi

r2l1=$(ls ${run2}_L001/*_L001_bismark_bt2_se.bam)
r2l2=$(ls ${run2}_L002/*_L002_bismark_bt2_se.bam)
r2l3=$(ls ${run2}_L003/*_L003_bismark_bt2_se.bam)
r2l4=$(ls ${run2}_L004/*_L004_bismark_bt2_se.bam)

r3l1=$(ls ${run3}_L001/*_L001_bismark_bt2_se.bam)
r3l2=$(ls ${run3}_L002/*_L002_bismark_bt2_se.bam)
r3l3=$(ls ${run3}_L003/*_L003_bismark_bt2_se.bam)
r3l4=$(ls ${run3}_L004/*_L004_bismark_bt2_se.bam)

r4l1=$(ls ${run4}_L001/*_L001_bismark_bt2_se.bam)
r4l2=$(ls ${run4}_L002/*_L002_bismark_bt2_se.bam)
r4l3=$(ls ${run4}_L003/*_L003_bismark_bt2_se.bam)
r4l4=$(ls ${run4}_L004/*_L004_bismark_bt2_se.bam)

base=${r4l1##*/}
base=${base%_GTAC*_L001_bismark_bt2_se.bam}_bismark_bt2
bam_rmDup=${base}_se.bam
txt_metrics=${base}_se.markDuplicates.txt
log=${base}_se.markDuplicates.log

set -x
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
    java -Xmx64G -jar $PICARD_HOME/picard.jar MarkDuplicates \
	I=$r2l1 I=$r2l2 I=$r2l3 I=$r2l4 \
	I=$r3l1 I=$r3l2 I=$r3l3 I=$r3l4 \
	I=$r4l1 I=$r4l2 I=$r4l3 I=$r4l4 \
	O=${bam_rmDup/_bismark/_run234_bismark} M=$txt_metrics REMOVE_DUPLICATES=true &> $log
    samtools merge $bam_rmDup ../merged_run1/WangT-MethylC-01-4dpaGFPpo-rep1-GTAC1-TGAGGTT_S1_bismark_bt2_se.rmdup.bam ${bam_rmDup/_bismark/_run234_bismark}
elif [ $SLURM_ARRAY_TASK_ID -ge 2 ] && [ $SLURM_ARRAY_TASK_ID -le 4 ]
then
    java -Xmx64G -jar $PICARD_HOME/picard.jar MarkDuplicates \
	I=$r1l1 I=$r1l2 I=$r1l3 I=$r1l4 \
	I=$r2l1 I=$r2l2 I=$r2l3 I=$r2l4 \
	I=$r3l1 I=$r3l2 I=$r3l3 I=$r3l4 \
	I=$r4l1 I=$r4l2 I=$r4l3 I=$r4l4 \
	O=$bam_rmDup M=$txt_metrics REMOVE_DUPLICATES=true &> $log
elif [ $SLURM_ARRAY_TASK_ID -ge 5 ]
then
    java -Xmx64G -jar $PICARD_HOME/picard.jar MarkDuplicates \
	I=$r2l1 I=$r2l2 I=$r2l3 I=$r2l4 \
	I=$r3l1 I=$r3l2 I=$r3l3 I=$r3l4 \
	I=$r4l1 I=$r4l2 I=$r4l3 I=$r4l4 \
	O=$bam_rmDup M=$txt_metrics REMOVE_DUPLICATES=true &> $log
fi
set +x

r1l1=$(ls ${run1}_L001/*_L001_bismark_bt2_pe.bam)
r1l2=$(ls ${run1}_L002/*_L002_bismark_bt2_pe.bam)
r1l3=$(ls ${run1}_L003/*_L003_bismark_bt2_pe.bam)
r1l4=$(ls ${run1}_L004/*_L004_bismark_bt2_pe.bam)

r2l1=$(ls ${run2}_L001/*_L001_bismark_bt2_pe.bam)
r2l2=$(ls ${run2}_L002/*_L002_bismark_bt2_pe.bam)
r2l3=$(ls ${run2}_L003/*_L003_bismark_bt2_pe.bam)
r2l4=$(ls ${run2}_L004/*_L004_bismark_bt2_pe.bam)

r3l1=$(ls ${run3}_L001/*_L001_bismark_bt2_pe.bam)
r3l2=$(ls ${run3}_L002/*_L002_bismark_bt2_pe.bam)
r3l3=$(ls ${run3}_L003/*_L003_bismark_bt2_pe.bam)
r3l4=$(ls ${run3}_L004/*_L004_bismark_bt2_pe.bam)

r4l1=$(ls ${run4}_L001/*_L001_bismark_bt2_pe.bam)
r4l2=$(ls ${run4}_L002/*_L002_bismark_bt2_pe.bam)
r4l3=$(ls ${run4}_L003/*_L003_bismark_bt2_pe.bam)
r4l4=$(ls ${run4}_L004/*_L004_bismark_bt2_pe.bam)

bam_rmDup=${base}_pe.bam
txt_metrics=${base}_pe.markDuplicates.txt
log=${base}_pe.markDuplicates.log

set -x
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
    java -Xmx64G -jar $PICARD_HOME/picard.jar MarkDuplicates \
	I=$r2l1 I=$r2l2 I=$r2l3 I=$r2l4 \
	I=$r3l1 I=$r3l2 I=$r3l3 I=$r3l4 \
	I=$r4l1 I=$r4l2 I=$r4l3 I=$r4l4 \
	O=${bam_rmDup/_bismark/_run234_bismark} M=$txt_metrics REMOVE_DUPLICATES=true &> $log
    samtools merge $bam_rmDup ../merged_run1/WangT-MethylC-01-4dpaGFPpo-rep1-GTAC1-TGAGGTT_S1_bismark_bt2_pe.rmdup.bam ${bam_rmDup/_bismark/_run234_bismark}
elif [ $SLURM_ARRAY_TASK_ID -ge 2 ] && [ $SLURM_ARRAY_TASK_ID -le 4 ]
then
    java -Xmx64G -jar $PICARD_HOME/picard.jar MarkDuplicates \
	I=$r1l1 I=$r1l2 I=$r1l3 I=$r1l4 \
	I=$r2l1 I=$r2l2 I=$r2l3 I=$r2l4 \
	I=$r3l1 I=$r3l2 I=$r3l3 I=$r3l4 \
	I=$r4l1 I=$r4l2 I=$r4l3 I=$r4l4 \
	O=$bam_rmDup M=$txt_metrics REMOVE_DUPLICATES=true &> $log
elif [ $SLURM_ARRAY_TASK_ID -ge 5 ]
then
    java -Xmx64G -jar $PICARD_HOME/picard.jar MarkDuplicates \
	I=$r2l1 I=$r2l2 I=$r2l3 I=$r2l4 \
	I=$r3l1 I=$r3l2 I=$r3l3 I=$r3l4 \
	I=$r4l1 I=$r4l2 I=$r4l3 I=$r4l4 \
	O=$bam_rmDup M=$txt_metrics REMOVE_DUPLICATES=true &> $log
fi
set +x

#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=scratch:200G
/scratch/$SLURM_JOB_USER/$SLURM_JOB_ID
TMPDIR=/scratch/$USER
mkdir -p "$TMPDIR"
cd "$TMPDIR"
cp <Lacie/pomerantzj/201123_A00269_0437_AHVH7MDRXX_fastqs_analysis/HVH7MDRXX/500_activated> $TMPDIR         
module load CBI
module load cellranger/5.0.1


cellranger-atac count --id=HVH7MDRXX --reference=/opt/refdata-cellranger-atac-GRCh38-1.2.0 \

                      --fastqs=$TMPDIR/500_activated --sample=500_activated \

                      --localcores=8 \

                      --localmem=64 
                      mv output.bam ~

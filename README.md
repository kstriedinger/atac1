#!/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=scratch:200G
/scratch/$SLURM_JOB_USER/$SLURM_JOB_ID
TMPDIR=/scratch/$USER
mkdir -p "$TMPDIR"
cd "$TMPDIR"


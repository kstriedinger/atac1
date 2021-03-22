 #!/bin/env bash
#SBATCH --ntasks=2
#SBATCH --mem=8gb 
#SBATCH --gres=scratch:200G
#SBATCH --time=3-00:00:00
cd $TMPDIR
cp <Lacie/pomerantzj/201123_A00269_0437_AHVH7MDRXX_fastqs_analysis/HVH7MDRXX/500_activated> $TMPDIR        
module load CBI
module load cellranger/5.0.1
mv output.bam ~

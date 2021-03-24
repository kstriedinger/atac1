#SBATCH --job-name=atac1
#SBATCH --mail-type=END,FAIL     
#SBATCH --time=12:01:00
#SBATCH --export=NONE
#SBATCH --output=atac1_%j.log
#SBATCH --ntasks=8
#SBATCH --mem=64gb
cp $HOME/500_activated $TMPDIR
module load CBI
module load cellranger/5.0.1
cellranger-atac count --id=HVH7MDRXX --reference=/opt/refdata-cellranger-atac-GRCh38-1.2.0 \

                      --fastqs=$TMPDIR/500_activated --sample=500_activated \

                      --localcores=8 \

                      --localmem=64 
                      mv output.bam ~

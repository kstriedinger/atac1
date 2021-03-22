cp /Volumes/LaCie/Complete\ ATAC\ files\ 2020/pomerantzj/201123_A00269_0437_AHVH7MDRXX_fastqs_analysis/HVH7MDRXX/500_activated $TMPDIR         
module load CBI
module load cellranger/5.0.1
module load cellranger-atac
$ cd /Volumes/LaCie/Complete\ ATAC\ files\ 2020/pomerantzj/201123_A00269_0437_AHVH7MDRXX_fastqs_analysis/HVH7MDRXX/500_activated
$ cellranger-atac count --id= \ HVH7MDRXX
                   --reference=/opt/refdata-cellranger-atac-GRCh38-1.2.0 \
                   --fastqs= /Volumes/LaCie/Complete\ ATAC\ files\ 2020/pomerantzj/201123_A00269_0437_AHVH7MDRXX_fastqs_analysis/HVH7MDRXX/500_activated
                   --sample=500_activated
                   --localcores=8 \
                   --localmem=64 
mv output.bam 
  

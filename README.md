## install and load SeuratDisk for saving Seurat objects

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")


library(SeuratDisk)
##

## install/update BiocManager (using v3.14 since this allows updates to a lot of packages)
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::version()
##

## install/update all packages loaded below (force re-install for Bioconductor packages)
bioc_pkgs <- c("Signac", "Seurat", "GenomeInfoDb", "EnsDb.Hsapiens.v86", "biovizBase")
BiocManager::install(bioc_pkgs)
install.packages(c("ggplot2", "patchwork", "hdf5r"))
## KJL: add packages for pathway analysis 
BiocManager::install(c("ReactomeGSA", "progeny"))

install_miniconda(path = miniconda_path, update = TRUE, force = FALSE)

install_CondaTools(tools="macs2", env="PeakCalling_analysis", pathToMiniConda="~/miniconda3")

install.packages("tidyverse")


#Call all libraries-------------
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(hdf5r)
library(biovizBase)
library(ReactomeGSA)
library(progeny)
library(tidyverse)
library(fs)

## make sure R is v4.1.1 at least and GenomeInfoDb is and v1.30.0 at least
sessionInfo()
## 

getwd() 

#Aggregate analysis HUMUSC ACTIVATED and NOT ACTIVATED---------


## creating directory reference for files 

counts500Ag <- Read10X_h5 (filename = "~/ATAC_500/outs/filtered_peak_bc_matrix.h5")
metadata500Ag <- read.csv (
file = "~/ATAC_500/outs/singlecell.csv", header = TRUE, row.names = 1)
metadata500Ag$sample_id = sapply( row.names(metadata500Ag), function(x){ unlist(strsplit(x, "-"))[2] } )


head(metadata500Ag) ## KJL 
metadata500Ag$sample_id = sapply( row.names(metadata500Ag), function(x){ unlist(strsplit(x, "-"))[2] } )
head(metadata500Ag) ## KJL

chrom_assay500Ag <- CreateChromatinAssay (
  counts = counts500Ag,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '~/ATAC_500/outs/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200,
  validate.fragments = TRUE
)

## 
chrom_assay500Ag 
slotNames(chrom_assay500Ag)
head(chrom_assay500Ag@fragments)
assaydat <- GetAssayData(chrom_assay500Ag, slot = "data")
##

husc500Ag <- CreateSeuratObject (
  counts = chrom_assay500Ag,
  assay = 'peaks500Ag',
  meta.data = metadata500Ag
)

## 
##select samples 1 and 2 (first create new objects to make sure this worked for the large objects)

## husc500Ag
colnames(husc500Ag)
length(colnames(husc500Ag)) ## 10432 samples
husc500Ag_12 <- husc500Ag[, grep("1|2", colnames(husc500Ag))]
colnames(husc500Ag_12)
length(colnames(husc500Ag_12)) ## 9232 samples
husc500Ag <- husc500Ag_12 ## assign to husc500Ag
rm(husc500Ag_12) ## get rid of this new object

## metadata500Ag
table(metadata500Ag$sample_id, useNA="ifany")
metadata500Ag <- metadata500Ag[grep("1|2", metadata500Ag$sample_id), ]
table(metadata500Ag$sample_id, useNA="ifany") ## only 1 and 2 samples

## assaydat
dim(assaydat)
head(assaydat@Dimnames[2])
vec123 <- unlist(assaydat@Dimnames[2])
length(vec123)
keep12 <- vec123[grep("1|2", vec123)]
length(keep12)
head(colnames(assaydat))
assaydat_12 <- assaydat[, colnames(assaydat) %in% keep12] 
dim(assaydat_12) ## 9232 samples
assaydat <- assaydat_12
rm(assaydat_12)

## counts500Ag
dim(counts500Ag)
counts500Ag_12 <- counts500Ag[, colnames(counts500Ag) %in% keep12] 
dim(counts500Ag_12) ## 9232 samples
counts500Ag <- counts500Ag_12
rm(counts500Ag_12)

## chrom_assay500Ag
dim(chrom_assay500Ag)
chrom_assay500Ag_12 <- CreateChromatinAssay (
  counts = counts500Ag, ## now only contains the 9232 samples 
  sep = c(":", "-"),
  fragments = '~/ATAC_500/outs/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200,
  validate.fragments = TRUE
)
dim(chrom_assay500Ag_12) ## 9232 samples
chrom_assay500Ag <- chrom_assay500Ag_12
rm(chrom_assay500Ag_12)

husc500Ag
##An object of class Seurat 
##96610 features across 9232 samples within 1 assay 
##Active assay: peaks500Ag (96610 features, 0 variable features)


husc500Ag[['peaks500Ag']]

granges(husc500Ag)

## :
slotNames(husc500Ag)
##

# extract gene annotations from ensdb---------


##: changing annotations to hg38 (v75 is for hg19) 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
annotations ## confirming GRCh38
seqlevels(annotations)
##


# change to UCSC style since the data was mapped to hg19-------
seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg19" ## KJL changing to hg38
genome(annotations) <- "hg38"
annotations ## KJL confirming hg38 (UCSC)

# add the gene information to the object-------
Annotation(husc500Ag) <- annotations

##: 
Annotation(husc500Ag)
Fragments(husc500Ag)
granges(husc500Ag)
##

#Computing QC Metrics----

# compute nucleosome signal score per cell
husc500Ag <- NucleosomeSignal(object = husc500Ag)

# compute TSS enrichment score per cell
husc500Ag <- TSSEnrichment(object = husc500Ag, fast = FALSE)

## :
head(husc500Ag@meta.data) 
FragmentHistogram(husc500Ag) ## enrichments around 50 and 200 (nuc-free cutoff should be between)
FragmentHistogram(husc500Ag, log.scale=TRUE) 
##


# add blacklist ratio and fraction of reads in peaks
husc500Ag$pct_reads_in_peaks <- husc500Ag$peak_region_fragments / husc500Ag$passed_filters * 100
husc500Ag$blacklist_ratio <- husc500Ag$blacklist_region_fragments / husc500Ag$peak_region_fragments

## :
summary(husc500Ag$pct_reads_in_peaks)
hist(husc500Ag$pct_reads_in_peaks)
##

husc500Ag$high.tss <- ifelse(husc500Ag$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(husc500Ag, group.by = 'high.tss') + NoLegend()


husc500Ag$nucleosome_group <- ifelse(husc500Ag$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = husc500Ag, group.by = 'nucleosome_group')

VlnPlot(
  object = husc500Ag,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

## 
olddim <- dim(husc500Ag@meta.data)
##

husc500Ag <- subset(
  x = husc500Ag,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
husc500Ag

## :
VlnPlot(
  object = husc500Ag,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
##


##:
newdim <- dim(husc500Ag@meta.data)
print(paste("Percent of peaks retained:", round(newdim[1]/olddim[1]*100, 3)))
##


#Normalization and linear dimensional reduction------

husc500Ag <- RunTFIDF(husc500Ag)
husc500Ag <- FindTopFeatures(husc500Ag, min.cutoff = 'q0')
husc500Ag <- RunSVD(husc500Ag)

DepthCor(husc500Ag)
## KJL: there is a strong correlation between the first component and the counts/cell, 
## so this is why it is removed below 

## :
husc500Ag
# An object of class Seurat 
# 96610 features across 2325 samples within 1 assay 
# Active assay: peaks500Ag (96610 features, 96610 variable features)
# 1 dimensional reduction calculated: lsi
##

#Non-linear dimension reduction and clustering------

husc500Ag <- RunUMAP(object = husc500Ag, reduction = 'lsi', dims = 2:30)
husc500Ag <- FindNeighbors(object = husc500Ag, reduction = 'lsi', dims = 2:30)

husc500Ag <- FindClusters(object = husc500Ag, verbose = FALSE, algorithm = 3)
DimPlot(object = husc500Ag, label = TRUE) + NoLegend()

## KJL: look at some different values of resolution in clustering
husc500Ag_copy <- husc500Ag

## resolution=0.8 (default and ultimately used for clustering and analysis)
husc500Ag <- FindClusters(object = husc500Ag, verbose = FALSE, algorithm = 3, resolution=0.8)
DimPlot(object = husc500Ag_copy, label = TRUE) + NoLegend() ## 5 clusters
husc500Ag_copy_df <- as.data.frame(husc500Ag_copy@meta.data)
table(husc500Ag_copy_df$seurat_clusters, husc500Ag_copy_df$sample_id)
prop.table(table(husc500Ag_copy_df$seurat_clusters, husc500Ag_copy_df$sample_id), margin=1)

## resolution=1
husc500Ag_copy <- FindClusters(object = husc500Ag, verbose = FALSE, algorithm = 3, resolution=1)
DimPlot(object = husc500Ag_copy, label = TRUE) + NoLegend() ## 6 clusters
husc500Ag_copy_df <- as.data.frame(husc500Ag_copy@meta.data)
table(husc500Ag_copy_df$seurat_clusters, husc500Ag_copy_df$sample_id)
prop.table(table(husc500Ag_copy_df$seurat_clusters, husc500Ag_copy_df$sample_id), margin=1)

## resolution=1.2
husc500Ag_copy <- FindClusters(object = husc500Ag, verbose = FALSE, algorithm = 3, resolution=1.2)
DimPlot(object = husc500Ag_copy, label = TRUE) + NoLegend() ## 7 clusters
husc500Ag_copy_df <- as.data.frame(husc500Ag_copy@meta.data)
table(husc500Ag_copy_df$seurat_clusters, husc500Ag_copy_df$sample_id)
prop.table(table(husc500Ag_copy_df$seurat_clusters, husc500Ag_copy_df$sample_id), margin=1)

# VizDimLoadings(husc500Ag, dims = 2:3, reduction = "lsi")
husc500Ag@reductions[["lsi"]]
## write a fragment bed files for each sample type
SplitFragments(object = husc500Ag, group.by="sample_id")
file.rename("1.bed", "activated.bed")
file.rename("2.bed", "not_activated.bed")
file.rename("3.bed", "otherfrags.bed")
## 

## Save
SaveH5Seurat(husc500Ag, overwrite=TRUE) 

##

#Create a gene activity matrix------

gene.activities <- GeneActivity(husc500Ag)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
husc500Ag[['RNA']] <- CreateAssayObject(counts = gene.activities)

## :
husc500Ag@assays[["RNA"]]@data
sumct <- rowSums(husc500Ag@assays[["RNA"]]@data)
summary(sumct)
length(which(sumct==0)) ## 255 genes have no expression
length(which(sumct==0)) / nrow(husc500Ag@assays[["RNA"]]@data)
## 



husc500Ag <- NormalizeData(
  object = husc500Ag,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(husc500Ag$nCount_RNA)
)

DefaultAssay(husc500Ag) <- 'RNA'




DimPlot(husc500Ag, label = TRUE) + NoLegend()


## umap by gene
FeaturePlot(husc500Ag, features = c("PAX7", "MYOD1", "MYOG","CCL20", "CCL2", "CXCL14", "TNFRSF12A", "IL33", "CXCL8", "IL1RL1", "CDKN1C", "MKI67"), 
            pt.size = 0.2, ncol = 3)

FeaturePlot(husc500Ag, features = c("CDKN2A", "CDKN1C", "CDKN1B", "PAX3", "CD9", "CD44", "CD34", "MYOD1", "CAV1", "SDC1"), 
            pt.size = 0.2, ncol = 3)

FeaturePlot(husc500Ag, features = c("GAS1", "HES1", "CALCR", "MYF5", "HEY1", "MYC", "FOS", "JUN"), 
            pt.size = 0.2, ncol = 3)

FeaturePlot(husc500Ag, features = c("CXCL8","CDKN1C", "CDKN2A", "CDKN1B", "CAV1", "YAP1", "MKI67", "PAX3", "FOX"), 
            pt.size = 0.2, ncol = 3)

FeaturePlot(husc500Ag, features = c("LIF","CDKN1C", "CDKN2A", "CDKN1B", "HEY1",  "FOS", "FOX1"), 
            pt.size = 0.2, ncol = 3)


## violins by gene

VlnPlot(husc500Ag, features = c("PAX7", "PAX3", "MYF5", "MKI67", "IL6", "CXCL8", "LIF", "CCL20", "CXCL14", "CCL2", "TNFRSF12A", "IL33", "CxCL8"),
        pt.size = 0.2, ncol = 3)

VlnPlot(husc500Ag, features = c("CXCL14", "CCL2", "TNFRSF12A", "IL33", "CxCL8"),
        pt.size = 0.2, ncol = 3)

VlnPlot(husc500Ag, features = c("FOX", "PAX3", "CXCL8", "IL33", "CCL2",
                                "CD9", "CD34", "CD44", "CDK6", "SPRY1"),
        pt.size = 0.2, ncol = 3)

VlnPlot(husc500Ag, features = c("IL6", "LIF", "CCL20", "CXCL14", "CCL2", "TNFRSF12A"),
        pt.size = 0.2, ncol = 3)


ExpressionPlot(
  object = husc500Ag,
  features = c("PAX7", "MYF5", "MYOG", "MYOD1", "MKI67"),
  group.by = "sample_id",
  assay = "RNA"
)
##

ExpressionPlot(
  object = husc500Ag,
  features = c("CCL20", "CCL2", "CXCL14", "LIF", "TNFRSF12A", "IL33", "CXCL8",
               "CDKN1C", "CDKN2A", "YAP1"),
  group.by = "sample_id",
  assay = "RNA"
)

ExpressionPlot(
  object = husc500Ag,
  features = c("CDKN1B","CDKN1C", "CDKN2A", "YAP1"),
  group.by = "sample_id",
  assay = "RNA"
)

ExpressionPlot(
  object = husc500Ag,
  features = c("PAX7", "PAX3", "MYOG", "MYF5", "MYOD1","MYC", "MKI67"),
  group.by = "sample_id",
  assay = "RNA"
)

ExpressionPlot(
  object = husc500Ag,
  features = c("CDKN2A", "CDKN1C", "CDKN1B", "MKI67"),
  group.by = "sample_id",
  assay = "RNA"
)

ExpressionPlot(
  object = husc500Ag,
  features = c("GAS1", "HES1", "CALCR", "MYF5", "HEY1", "MYC", "FOS", "JUN", "MKI67","CD9", "CD44", "CD34", "MYOD1", "CAV1", "SDC1" ),
  group.by = "sample_id",
  assay = "RNA"
)

egenes <- c("PAX7", "MYOD1", "MYOG","CCL20", "CCL2", "CXCL14", "TNFRSF12A", "IL33", "CXCL8", "IL1RL1", "MKI67")
egenes %in% row.names(husc500Ag@assays[["RNA"]]@data)


VlnPlot(husc500Ag, features = c("PAX7", "MYOD1", "MYOG","MYC", "CCL20", "CCL2", "CXCL14", "TNFRSF12A", "IL33", "CXCL8", "IL1RL1", "CDKN1C", "MKI67"),
        pt.size = 0.2, ncol = 3)

VlnPlot(husc500Ag, features = c("PAX7", "PAX3", "MYOD1", "MYF5", "MKI67"),
        pt.size = 0.2, ncol = 3)

VlnPlot(husc500Ag, features = c("CCL20", "CCL2", "CXCL14", "TNFRSF12A", "IL33"),
        pt.size = 0.2, ncol = 3)

ExpressionPlot(
  object = husc500Ag,
  features = c("CCL20", "CCL2", "CXCL14", "TNFRSF12A", "IL33", "CXCL8", "IL1RL1", "CD9", "CD44", "CD34"),
  group.by = "sample_id",
  assay = "RNA"
)


ExpressionPlot(
  object = husc500Ag,
  features = c("GAS1", "HES1", "CALCR", "HEY1", "MYC", "FOS", "JUN"),
  group.by = "sample_id",
  assay = "RNA"
)


#Find differentially accessible peaks------

# change back to working with peaks instead of gene activities

DefaultAssay(husc500Ag) <- 'peaks500Ag'

da_peaks <- FindMarkers(
  object = husc500Ag,
  ident.1 = "1",
  ident.2 = "2",
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)

open_1 <- rownames(da_peaks[da_peaks$avg_log2FC > 0.5, ])
open_2 <- rownames(da_peaks[da_peaks$avg_log2FC < -0.5, ])

closest_genes_1 <- ClosestFeature(husc500Ag, regions = open_1)
closest_genes_2 <- ClosestFeature(husc500Ag, regions = open_2)

head(open_2)

list(open_1)
list(open_2)

## looks like none meet the +/-0.5 log2fc threshold criteria:
summary(da_peaks$avg_log2FC)
##

## look at adjusted p-values and use this as a threshold instead? 
## e.g. you could use 0.05 as a threshold as below (or you can use 0.01 or 0.1 if you
## want to be more/less stringent)

## let's redo the above but use a 0.05 adjusted p-value to indicate "significant"
## up/down log2fc
summary(da_peaks$p_val_adj)
plot(da_peaks$p_val_adj, da_peaks$avg_log2FC) 

open_1 <- rownames(da_peaks[da_peaks$p_val_adj<0.05 & da_peaks$avg_log2FC>0, ])
length(open_1) ## 106 are significantly up
open_2 <- rownames(da_peaks[da_peaks$p_val_adj<0.05 & da_peaks$avg_log2FC<0, ])
length(open_2) ## 52 are significantly down

closest_genes_1 <- ClosestFeature(husc500Ag, regions = open_1)
closest_genes_2 <- ClosestFeature(husc500Ag, regions = open_2)

list(open_1)
list(open_2)

## let's save these genes in a data frame with just the significant ones
sigpeaks <- c(open_1, open_2)
length(sigpeaks)
sig_da_peaks <- da_peaks[which(rownames(da_peaks) %in% c(open_1, open_2)),]
dim(sig_da_peaks)
head(sig_da_peaks)
summary(sig_da_peaks$p_val_adj) ## all p<0.05
write.csv(sig_da_peaks, file="sig_da_peaks.csv", quote=FALSE)
##


## violin plot code 

plot1 <- VlnPlot(
  object = husc500Ag,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("1","2")
)

plot2 <- FeaturePlot(
  object = husc500Ag,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

plot3 <- VlnPlot(
  object = husc500Ag,
  features = rownames(da_peaks)[2],
  pt.size = 0.1,
  idents = c("1","2")
)
plot4 <- FeaturePlot(
  object = husc500Ag,
  features = rownames(da_peaks)[2],
  pt.size = 0.1
)

plot3 | plot4





##levels(husc500Ag) <- c("PAX7", "MYOD1", "MYOG","CCL20", "CCL2", "CXCL14", "TNFRSF12A", "IL33", "CXCL8", "IL1RL1", "CDKN1C", "MKI67")

CoveragePlot(
object = husc500Ag,
region = rownames(da_peaks),
group.by = "sample_id",
extend.upstream = 40000,
extend.downstream = 20000
)


-----
  #fold change 
  
fc <- FoldChange(husc500Ag, ident.1 = "1", ident.2 = "2")
head(fc)
list(fc)


# Calling peaks----


DimPlot(husc500Ag)

peaks500Ag <- CallPeaks(
  object = husc500Ag,
  group.by = "sample_id")
  
  ## Coverage Plots below

CoveragePlot(
  object = husc500Ag,
  region = "LIF",
  ranges = peaks500Ag,
  group.by = "sample_id",
  ranges.title = "MACS2")

  CoveragePlot(
    object = husc500Ag,
    region = "IL6",
    ranges = peaks500Ag,
    group.by = "sample_id",
    ranges.title = "MACS2")
  
  
  CoveragePlot(
    object = husc500Ag,
    region = "CCl2",
    ranges = peaks500Ag,
    extend.upstream = 1000,
    extend.downstream = 1000,
    group.by = "sample_id",
    ranges.title = "MACS2")
  
  
  CoveragePlot(
    object = husc500Ag,
    region = "CCl20",
    ranges = peaks500Ag,
    extend.upstream = 1000,
    extend.downstream = 1000,
    group.by = "sample_id",
    ranges.title = "MACS2")
  
  CoveragePlot(
    object = husc500Ag,
    region = "CXCL14",
    ranges = peaks500Ag,
    group.by = "sample_id",
    ranges.title = "MACS2")
  
  CoveragePlot(
    object = husc500Ag,
    region = "CCL8",
    ranges = peaks500Ag,
    group.by = "sample_id",
    ranges.title = "MACS2")
  
  CoveragePlot(
    object = husc500Ag,
    region = "IL33",
    ranges = peaks500Ag,
    group.by = "sample_id",
    ranges.title = "MACS2")
  
  CoveragePlot(
    object = husc500Ag,
    region = "CDKN1C",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  
  CoveragePlot(
    object = husc500Ag,
    region = "CD14",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "CDKN2A",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  
  CoveragePlot(
    object = husc500Ag,
    region = "NR3C1",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "CCL20",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "IL6",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  CoveragePlot(
    object = husc500Ag,
    region = "CCL2",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  CoveragePlot(
    object = husc500Ag,
    region = "LIF",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "CD44",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "CD9",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "CXCL8",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "PAX7",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "MYOD1",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "IL33",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "MKI67",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "MYF5",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "TNFRSF12A",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "CD34",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "IL1R1",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "CXCL14",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )
  
  CoveragePlot(
    object = husc500Ag,
    region = "PAX3",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    extend.upstream = 1000,
    extend.downstream = 1000,
    ncol = 1,
    group.by = "sample_id",
  )

  #Quantification of these genes-----
  
  
  ##
  
  FeatureMatrix(
    fragments = fragments,
    features = granges(peaks500Ag)
  )
  # how to represent the quantification of the feature matrix
  gene_plot <- AnnotationPlot(
    object = husc500Ag,
    region = "PAX7",
  )
  gene_plot
  
  
  # )
  # ________
  # Visualization of genomic regions : https://satijalab.org/signac/articles/visualization.html
  
  
  tile_plot <- TilePlot (
    object = husc500Ag,
    region = c("PAX7", "MYOD1"),
    idents = c("1", "2"),
  )
  tile_plot
  
  
  ExpressionPlot(
    object = husc500Ag,
    features = c("PAX7", "MYOD1", 'CXCL14', 'TNFRSF12A', 'CCL2', 'IL33',  
                 'CDKN1C', 'CDKN2A', 'YAP1'), 
    assay = "RNA")
  
  CoveragePlot(
    object = husc500Ag,
    region = "chr2-87011729-87035519",
    features = "CD8A",
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE,
    group.by = "sample_id",
  )
  
  FeaturePlot(
    object = husc500Ag,
    features = c("PAX7", "MYF5", 'CXCL14', 'TNFRSF12A', 'CCL20', 'IL33', 'IL8', "LIF", 'CDKN1C', 'CDKN2A', 'YAP1'),
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3,
  )
  
  
  CombineTracks(
    plotlist = list(FeaturePlot, gene_plot),
    expression.plot = expr_plot,
    heights = c(10, 6, 1, 2, 3),
    widths = c(10, 1)
  )
  
#Pathway analysis---------
 
## using ReactomeGSA
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require(ReactomeGSA))
  BiocManager::install("ReactomeGSA")
#> Loading required package: ReactomeGSA

# install the ReactomeGSA.data package for the example data
if (!require(ReactomeGSA.data))
  BiocManager::install("ReactomeGSA.data")
#> Loading required package: ReactomeGSA.data
#> Loading required package: limma
#> Loading required package: edgeR
#> Loading required package: Seurat
#> Attaching SeuratObject
#> Attaching sp


## see https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html

gsva_result <- analyse_sc_clusters(husc500Ag, verbose = TRUE)
gsva_result

#ReactomeAnalysisResult object
#ReactomeAnalysisResult object
#Reactome Release: 82
#Results:
  #- Seurat:
  #1730 pathways
#13871 fold changes for genes
#No Reactome visualizations available

## pathway-level expression values per cell cluster
pathway_expression <- pathways(gsva_result)

## simplify the column names by removing the default dataset identifier
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
head(pathway_expression)

## find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min

## sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

## heatmap
dev.off() 
plot_gsva_heatmap(gsva_result, max_pathways = 10, group.by = "sample_id")

#graphs for individual pathways
plot_gsva_pathway(gsva_result, pathway_id = "R-HSA-9014843")

plot_gsva_pathway(gsva_result, pathway_id = "R-HSA-1059683")

plot_gsva_pathway(gsva_result, pathway_id = "R-HSA-196299")

plot_gsva_pathway(gsva_result, pathway_id = "R-HSA-9013700")

plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[8])

#graphs for relevant pathways

relevant_pathways <- c("R-HSA-9014843", "R-HSA-1059683", "R-HSA-109606", "R-HSA-109703", "R-HSA-109704", "R-HSA-110056","R-HSA-1299361", "R-HSA-112122", "R-HSA-141334", "R-HSA-1483196", "R-HSA-166020", "R-HSA-391908")
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,20), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left

relevant_pathways <- c("R-HSA-9014843", "R-HSA-1059683", "R-HSA-196299", "R-HSA-9013700")
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,20), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left

## plot 8 most enriched (based on heatmap above)

head(max_difference, n=20) ## names of first 20


relevant_pathways <- c("R-HSA-112122", "R-HSA-1299287", "R-HSA-2408552", "R-HSA-451307", "R-HSA-5263617","R-HSA-211994", "R-HSA-1299503", "R-HSA-9014843", "R-HSA-141334", "R-HSA-1483196", "R-HSA-166020", "R-HSA-391908")
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,20), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left


relevant_pathways <- c("R-HSA-1483076", "R-HSA-1296061", "R-HSA-3248023", "R-HSA-1299361", "R-HSA-70688","R-HSA-391908", "R-HSA-916853", "R-HSA-1296025", "R-HSA-166020", "R-HSA-1855231")
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,20), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left

##Most relevant
relevant_pathways <- c("R-HSA-9014843", "R-HSA-109704", "R-HSA-177929", "R-HSA-1059683", "R-HSA-422356")
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,20), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left

relevant_pathways <- c("R-HSA-110056", "R-HSA-6783589", "R-HSA-196299", "R-HSA-512988", "R-HSA-1059683", "R-HSA-8983432", "R-HSA-525793", "R-HSA-451927", "R-HSA-9008059")
plot_gsva_heatmap(gsva_result, 
                  pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,20), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left




plot_gsva_pca(gsva_result)

png(file=paste0("pathway1_dapeaks_bycluster.png"))
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])
dev.off()
png(file=paste0("pathway2_dapeaks_bycluster.png"))
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[2])
dev.off()
png(file=paste0("pathway3_dapeaks_bycluster.png"))
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[3])
dev.off()
png(file=paste0("pathway4_dapeaks_bycluster.png"))
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[4])
dev.off()
png(file=paste0("pathway5_dapeaks_bycluster.png"))
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[5])
dev.off()
png(file=paste0("pathway6_dapeaks_bycluster.png"))
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[6])
dev.off()
png(file=paste0("pathway7_dapeaks_bycluster.png"))
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[7])
dev.off()
png(file=paste0("pathway8_dapeaks_bycluster.png"))
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[8])
dev.off()



#save rds objects----

saveRDS(metadata500Ag, 'metadata500Ag.rds')
saveRDS(husc500Ag, 'husc500Ag.rds')
saveRDS(plot1, 'plot1.rds')
saveRDS(plot2, 'plot2.rds')
saveRDS(gene.activities, 'gene.activities.rds')
saveRDS(fragments, 'fragments.rds')
saveRDS(fc, 'fc.rds')
saveRDS(counts500Ag, 'counts500Ag.rds')
saveRDS(closest_genes_2, 'closest_genes_2.rds')
saveRDS(chrom_assay500Ag,'chrom_assay500Ag.rds')
saveRDS(annotations, 'annotations.rds')
save.image(file = 'Differential_Analysis.Rdata') 
# 
readRDS(file = 'metadata500Ag.rds')



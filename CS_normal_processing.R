library(Seurat)               
library(ggplot2)              
library(dplyr)                
library(Matrix)               
library(DoubletFinder) 

setwd("~/scRNAseq_SilviaG/Endo")

filename <- file.path("F:/singleCell/Carrera4/20240422/A2_cr/outs/per_sample_outs/3/count/sample_filtered_feature_bc_matrix.h5")
hdf5_obj <- Read10X_h5(filename = filename,use.names = TRUE,
                       unique.features = TRUE)
CS1 <- CreateSeuratObject(counts = hdf5_obj, min.cells = 0, min.features = 0)
CS1$orig.ident <- "CS1"
filename <- file.path("F:/singleCell/Carrera4/20240422/A2_cr/outs/per_sample_outs/4/count/sample_filtered_feature_bc_matrix.h5")
hdf5_obj <- Read10X_h5(filename = filename,use.names = TRUE,
                       unique.features = TRUE)
CS2 <- CreateSeuratObject(counts = hdf5_obj, min.cells = 0, min.features = 0)
CS2$orig.ident <- "CS2"
filename <- file.path("F:/singleCell/Carrera4/20240422/D2_cr/outs/per_sample_outs/15/count/sample_filtered_feature_bc_matrix.h5")
hdf5_obj <- Read10X_h5(filename = filename,use.names = TRUE,
                       unique.features = TRUE)
CS3 <- CreateSeuratObject(counts = hdf5_obj, min.cells = 0, min.features = 0)
CS3$orig.ident <- "CS3"
filename <- file.path("F:/singleCell/Carrera4/20240422/B2_cr/outs/per_sample_outs/5/count/sample_filtered_feature_bc_matrix.h5")
hdf5_obj <- Read10X_h5(filename = filename,use.names = TRUE,
                       unique.features = TRUE)
CS4 <- CreateSeuratObject(counts = hdf5_obj, min.cells = 0, min.features = 0)
CS4$orig.ident <- "CS4"
filename <- file.path("F:/singleCell/Carrera4/20240422/A2_cr/outs/per_sample_outs/2/count/sample_filtered_feature_bc_matrix.h5")
hdf5_obj <- Read10X_h5(filename = filename,use.names = TRUE,
                       unique.features = TRUE)
CS5 <- CreateSeuratObject(counts = hdf5_obj, min.cells = 0, min.features = 0)
CS5$orig.ident <- "CS5"
filename <- file.path("F:/singleCell/Carrera4/20240422/A2_cr/outs/per_sample_outs/1/count/sample_filtered_feature_bc_matrix.h5")
hdf5_obj <- Read10X_h5(filename = filename,use.names = TRUE,
                       unique.features = TRUE)
CS6 <- CreateSeuratObject(counts = hdf5_obj, min.cells = 0, min.features = 0)
CS6$orig.ident <- "CS6"
filename <- file.path("F:/singleCell/Carrera4/20240422/C2_cr/outs/per_sample_outs/11/count/sample_filtered_feature_bc_matrix.h5")
hdf5_obj <- Read10X_h5(filename = filename,use.names = TRUE,
                       unique.features = TRUE)
N1 <- CreateSeuratObject(counts = hdf5_obj, min.cells = 0, min.features = 0)
N1$orig.ident <- "N1"
filename <- file.path("F:/singleCell/Carrera4/20240422/C2_cr/outs/per_sample_outs/12/count/sample_filtered_feature_bc_matrix.h5")
hdf5_obj <- Read10X_h5(filename = filename,use.names = TRUE,
                       unique.features = TRUE)
N2 <- CreateSeuratObject(counts = hdf5_obj, min.cells = 0, min.features = 0)
N2$orig.ident <- "N2"

samples_v <- c("CS1", "CS2", "CS3", "CS4", "CS5", "CS6", "N1", "N2")

if (length(samples_v) > 1){
  rawData_combined <- merge(x = get(samples_v[1]), y = sapply(samples_v[-1], get), add.cell.ids = samples_v)
}

## QC and Filtering ----------

rawData_combined[["percent_mt"]] <- PercentageFeatureSet(rawData_combined, pattern = "^MT-")
rawData_combined$log10GenesPerUMI <- log10(rawData_combined$nFeature_RNA) / log10(rawData_combined$nCount_RNA)
min_nFeature <- 200
min_log10GenesPerUMI <- 0.8
min_percent_mt <- 20
data_combined <- subset(rawData_combined, subset = nFeature_RNA > min_nFeature & log10GenesPerUMI > min_log10GenesPerUMI & percent_mt < min_percent_mt)
VlnPlot(data_combined, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 4, pt.size=0)

# pre-process standard workflow-------
data_combined <- NormalizeData(object = data_combined)
data_combined <- FindVariableFeatures(object = data_combined)
data_combined <- ScaleData(object = data_combined)
data_combined <- RunPCA(object = data_combined)
ElbowPlot(data_combined)
data_combined <- FindNeighbors(object = data_combined, dims = 1:20)
data_combined <- FindClusters(object = data_combined, resolution =0.3)
data_combined <- RunUMAP(object = data_combined, dims = 1:20)

# Doublets detection---------
## pK Identification (no ground-truth) 
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
sweep.res.list_ <- paramSweep(data_combined, PCs = 1:15, sct = TRUE)
sweep.stats_ <- summarizeSweep(sweep.res.list_, GT = FALSE)
bcmvn_ <- find.pK(sweep.stats_)

ggplot(bcmvn_, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_ %>%   filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 
print(pK)
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- data_combined@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*nrow(data_combined@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
data_combined_doublets_singlets <- doubletFinder(data_combined, 
                                                    PCs = 1:20, 
                                                    pN = 0.25, 
                                                    pK = pK, 
                                                    nExp = nExp_poi.adj,
                                                    reuse.pANN = FALSE, sct = FALSE)


colnames(data_combined_doublets_singlets@meta.data)[colnames(data_combined_doublets_singlets@meta.data) == "DF.classifications_0.25_0.02_7265"] <- "Singlets_doublets"

# Singlets pre-process---------
data_combined_singlets <- subset(data_combined_doublets_singlets, subset = Singlets_doublets == "Singlet")

split_data <- SplitObject(data_combined_singlets, split.by = "orig.ident")

options(future.globals.maxSize = 2 * 1024^3)  

split_data <- lapply(X = split_data, FUN = function(x) {
  x@meta.data$orig.ident[1]
  x <- SCTransform(x, vars.to.regress = c("percent_mt"), vst.flavor = "v2", verbose = FALSE)
  return(x)
})

data_combined_singlets <- RunPCA(data_combined_singlets, assay = "SCT")
ElbowPlot(data_combined_singlets, ndims = 50)
data_combined_singlets <- RunUMAP(data_combined_singlets, dims = 1:30)
data_combined_singlets <- FindNeighbors(data_combined_singlets, dims = 1:30)
data_combined_singlets <- FindClusters(data_combined_singlets, verbose = FALSE, resolution = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4))
resolution_find_clusters <- "SCT_snn_res.0.6"
Idents(object = data_combined_singlets) <- resolution_find_clusters
data_combined_singlets$seurat_clusters <- data_combined_singlets$SCT_snn_res.0.6

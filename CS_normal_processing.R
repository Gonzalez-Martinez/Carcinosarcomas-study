library(Seurat)               
library(ggplot2)              
library(dplyr)                
library(Matrix)               
library(DoubletFinder) 

setwd("~/scRNAseq_SilviaG/Endo")

filename <- file.path("E:/singleCell/Carrera4/20240422/Agregados/outs/count/filtered_feature_bc_matrix.h5")
hdf5_obj <- Read10X_h5(filename = filename,use.names = TRUE,
                       unique.features = TRUE)
Endometrial_matrix <- CreateSeuratObject(counts = hdf5_obj)

## QC and Filtering
Endometrial_matrix$mitoPercent <- PercentageFeatureSet(Endometrial_matrix, pattern = '^MT-')
VlnPlot(Endometrial_matrix, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3, pt.size=0)
Endometrial_matrix <- subset(Endometrial_matrix, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & nCount_RNA > 250 & nCount_RNA < 50000 & mitoPercent < 15)

# pre-process standard workflow
Endometrial_matrix <- NormalizeData(object = Endometrial_matrix)
Endometrial_matrix <- FindVariableFeatures(object = Endometrial_matrix)
Endometrial_matrix <- ScaleData(object = Endometrial_matrix)
Endometrial_matrix <- RunPCA(object = Endometrial_matrix)
ElbowPlot(Endometrial_matrix)
Endometrial_matrix <- FindNeighbors(object = Endometrial_matrix, dims = 1:15)
Endometrial_matrix <- FindClusters(object = Endometrial_matrix, resolution =0.3)
Endometrial_matrix <- RunUMAP(object = Endometrial_matrix, dims = 1:15)

#SUBSET POOLES

Endometrial_matrix_pool_1 <- subset(Endometrial_matrix, subset = Pool %in% c("1"))
Endometrial_matrix_pool_2 <- subset(Endometrial_matrix, subset = Pool %in% c("2"))
Endometrial_matrix_pool_3 <- subset(Endometrial_matrix, subset = Pool %in% c("3"))
Endometrial_matrix_pool_4 <- subset(Endometrial_matrix, subset = Pool %in% c("4"))
save(Endometrial_matrix_pool_1, file="Endometrial_matrix_pool_1.RData")
save(Endometrial_matrix_pool_2, file="Endometrial_matrix_pool_2.RData")
save(Endometrial_matrix_pool_3, file="Endometrial_matrix_pool_3.RData")
save(Endometrial_matrix_pool_4, file="Endometrial_matrix_pool_4.RData")

#DOUBLETS

###Pool 1
## pK Identification
sweep.res.list_ <- paramSweep(Endometrial_matrix_pool_1, PCs = 1:15, sct = FALSE)
sweep.stats_ <- summarizeSweep(sweep.res.list_, GT = FALSE)
bcmvn_ <- find.pK(sweep.stats_)

ggplot(bcmvn_, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_ %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 
print(pK)
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Endometrial_matrix_pool_1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.076*nrow(Endometrial_matrix_pool_1@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Endometrial_matrix_pool_1_doublets <- doubletFinder(Endometrial_matrix_pool_1, 
                                                    PCs = 1:15, 
                                                    pN = 0.25, 
                                                    pK = pK, 
                                                    nExp = nExp_poi.adj,
                                                    reuse.pANN = FALSE, sct = FALSE)

colnames(Endometrial_matrix_pool_1_doublets@meta.data)[colnames(Endometrial_matrix_pool_1_doublets@meta.data) == "DF.classifications_0.25_0.04_2620"] <- "Singlets_doublets"
table(Endometrial_matrix_pool_1_doublets@meta.data$Singlets_doublets)
table(Endometrial_matrix_pool_1_doublets@meta.data$Sample_ID, 
      Endometrial_matrix_pool_1_doublets@meta.data$Singlets_doublets)
save(Endometrial_matrix_pool_1_doublets, file="Endometrial_matrix_pool_1_doublets.RData")

###Pool 2

sweep.res.list_ <- paramSweep(Endometrial_matrix_pool_2, PCs = 1:15, sct = FALSE)
sweep.stats_ <- summarizeSweep(sweep.res.list_, GT = FALSE)
bcmvn_ <- find.pK(sweep.stats_)

ggplot(bcmvn_, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- sweep.stats_ %>%
  filter(BCreal == max(BCreal)) %>%
  dplyr::select(pK)
print(pK)

pK <- as.numeric(as.character(pK[[1]]))

#Probar esto
pK <- bcmvn_ %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 
print(pK)
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Endometrial_matrix_pool_2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.076*nrow(Endometrial_matrix_pool_2@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Endometrial_matrix_pool_2_doublets <- doubletFinder(Endometrial_matrix_pool_2, 
                                                    PCs = 1:15, 
                                                    pN = 0.25, 
                                                    pK = pK, 
                                                    nExp = nExp_poi.adj,
                                                    reuse.pANN = FALSE, sct = FALSE)

colnames(Endometrial_matrix_pool_2_doublets@meta.data)[colnames(Endometrial_matrix_pool_2_doublets@meta.data) == "DF.classifications_0.25_0.25_3227"] <- "Singlets_doublets"
table(Endometrial_matrix_pool_2_doublets@meta.data$Singlets_doublets)
table(Endometrial_matrix_pool_2_doublets@meta.data$Sample_ID, 
      Endometrial_matrix_pool_2_doublets@meta.data$Singlets_doublets) 

save(Endometrial_matrix_pool_2_doublets, file="Endometrial_matrix_pool_2_doublets.RData")

###Pool 3

sweep.res.list_ <- paramSweep(Endometrial_matrix_pool_3, PCs = 1:15, sct = FALSE)
sweep.stats_ <- summarizeSweep(sweep.res.list_, GT = FALSE)
bcmvn_ <- find.pK(sweep.stats_)

ggplot(bcmvn_, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_ %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 
print(pK)
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Endometrial_matrix_pool_3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.076*nrow(Endometrial_matrix_pool_3@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Endometrial_matrix_pool_3_doublets <- doubletFinder(Endometrial_matrix_pool_3, 
                                                    PCs = 1:15, 
                                                    pN = 0.25, 
                                                    pK = pK, 
                                                    nExp = nExp_poi.adj,
                                                    reuse.pANN = FALSE, sct = FALSE)

colnames(Endometrial_matrix_pool_3_doublets@meta.data)[colnames(Endometrial_matrix_pool_3_doublets@meta.data) == "DF.classifications_0.25_0.3_2701"] <- "Singlets_doublets"
table(Endometrial_matrix_pool_3_doublets@meta.data$Singlets_doublets)
table(Endometrial_matrix_pool_3_doublets@meta.data$Sample_ID, 
      Endometrial_matrix_pool_3_doublets@meta.data$Singlets_doublets)  
save(Endometrial_matrix_pool_3_doublets, file="Endometrial_matrix_pool_3_doublets.RData")

### Pool 4
sweep.res.list_pbmc <- paramSweep(Endometrial_matrix_pool_4, PCs = 1:15, sct = FALSE)
sweep.stats_ <- summarizeSweep(sweep.res.list_, GT = FALSE)
bcmvn_ <- find.pK(sweep.stats_)

ggplot(bcmvn_, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_ %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 
print(pK)
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- Endometrial_matrix_pool_4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)          
nExp_poi <- round(0.076*nrow(Endometrial_matrix_pool_4@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder 
Endometrial_matrix_pool_4_doublets <- doubletFinder(Endometrial_matrix_pool_4, 
                                                    PCs = 1:15, 
                                                    pN = 0.25, 
                                                    pK = pK, 
                                                    nExp = nExp_poi.adj,
                                                    reuse.pANN = FALSE, sct = FALSE)

colnames(Endometrial_matrix_pool_4_doublets@meta.data)[colnames(Endometrial_matrix_pool_4_doublets@meta.data) == "DF.classifications_0.25_0.24_2591"] <- "Singlets_doublets"
table(Endometrial_matrix_pool_4_doublets@meta.data$Singlets_doublets)
table(Endometrial_matrix_pool_4_doublets@meta.data$Sample_ID, 
      Endometrial_matrix_pool_4_doublets@meta.data$Singlets_doublets) 
save(Endometrial_matrix_pool_4_doublets, file="Endometrial_matrix_pool_4_doublets.RData")

#Merge

Endometrial_matrix_doublets <- merge(Endometrial_matrix_pool_1_doublets, y = list(Endometrial_matrix_pool_2_doublets, Endometrial_matrix_pool_3_doublets, Endometrial_matrix_pool_4_doublets))

# pre-process standard workflow
Endometrial_matrix_doublets <- NormalizeData(object = Endometrial_matrix_doublets)
Endometrial_matrix_doublets <- FindVariableFeatures(object = Endometrial_matrix_doublets)
Endometrial_matrix_doublets <- ScaleData(object = Endometrial_matrix_doublets)
Endometrial_matrix_doublets <- RunPCA(object = Endometrial_matrix_doublets)
ElbowPlot(Endometrial_matrix_doublets)
Endometrial_matrix_doublets <- FindNeighbors(object = Endometrial_matrix_doublets, dims = 1:15)
Endometrial_matrix_doublets <- FindClusters(object = Endometrial_matrix_doublets, resolution =0.3)
Endometrial_matrix_doublets <- RunUMAP(object = Endometrial_matrix_doublets, dims = 1:15)

Endometrial_matrix_singlets <- subset(Endometrial_matrix_doublets, subset = Singlets_doublets %in% c("Singlet"))

#Carcinosarcomas + normal analysis-----------

CS_normal <- subset(Endometrial_matrix_singlets, subset = Histological_type %in% c("Carcinosarcoma", "Normal"))

CS_normal <- NormalizeData(object = CS_normal)
CS_normal <- FindVariableFeatures(object = CS_normal)
CS_normal <- ScaleData(object = CS_normal)
CS_normal <- RunPCA(object = CS_normal)
ElbowPlot(CS_normal)
CS_normal <- FindNeighbors(object = CS_normal, dims = 1:15)
CS_normal <- FindClusters(object = CS_normal, resolution =0.4)
CS_normal <- RunUMAP(object = CS_normal, dims = 1:15)

save(CS_normal, file="CS_normal.RData")

markers_CS_normal <- FindAllMarkers(CS_normal, only.pos = TRUE)
save(markers_CS_normal, file="markers_CS_normal.RData")

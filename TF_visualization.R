#Expression matrix generation to pySCENIC-----
library(Seurat)
library(Matrix)

setwd("~/Documentos/Carcinosarcomas")

# CS 1

load("~/Documentos/Carcinosarcomas/CS1_epi_mes_CS.RData")

exprMat_CS1 <- GetAssayData(CS1_epi_mes_CS, assay = 'RNA', layer = 'counts')
exprMat_CS1 <- as.matrix(exprMat_CS1)

# save CSV
write.csv(exprMat_CS1, file = "exprMat_CS1.csv", row.names = TRUE)

# CS 2

load("~/Documentos/Carcinosarcomas/CS2_epi_mes_CS.RData")

exprMat_CS2 <- GetAssayData(CS2_epi_mes_CS, assay = 'RNA', layer = 'counts')
exprMat_CS2 <- as.matrix(exprMat_CS2)

# save CSV
write.csv(exprMat_CS2, file = "exprMat_CS2.csv", row.names = TRUE)

# CS 3

load("~/Documentos/Carcinosarcomas/CS1_epi_mes_CS.RData")

exprMat_CS3 <- GetAssayData(CS3_epi_mes_CS, assay = 'RNA', layer = 'counts')
exprMat_CS3 <- as.matrix(exprMat_CS3)

# save CSV
write.csv(exprMat_CS3, file = "exprMat_CS3.csv", row.names = TRUE)

# CS 4

load("~/Documentos/Carcinosarcomas/CS4_epi_mes_CS.RData")

exprMat_CS4 <- GetAssayData(CS4_epi_mes_CS, assay = 'RNA', layer = 'counts')
exprMat_CS4 <- as.matrix(exprMat_CS4)

# save CSV
write.csv(exprMat_CS4, file = "exprMat_CS4.csv", row.names = TRUE)

# CS 5

load("~/Documentos/Carcinosarcomas/CS5_epi_mes_CS.RData")

exprMat_CS5 <- GetAssayData(CS5_epi_mes_CS, assay = 'RNA', layer = 'counts')
exprMat_CS5 <- as.matrix(exprMat_CS5)

# save CSV
write.csv(exprMat_CS5, file = "exprMat_CS5.csv", row.names = TRUE)

# CS 6

load("~/Documentos/Carcinosarcomas/CS6_epi_mes_CS.RData")

exprMat_CS6 <- GetAssayData(CS6_epi_mes_CS, assay = 'RNA', layer = 'counts')
exprMat_CS6 <- as.matrix(exprMat_CS6)

# save CSV
write.csv(exprMat_CS6, file = "exprMat_CS6.csv", row.names = TRUE)

#Seurat Extend SCENIC----------

#R

library(Seurat)
library(SeuratExtend)
library(loomR)
library(hdf5r)
library(SCopeLoomR)
library(SCENIC)

#CS 1

scenic_loom_path <- "/home/silvia/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS1/CS1_scenic_integrated-output.loom"

# loom cells filter
loom_cells_data <- h5read("~/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS1/CS1_scenic_integrated-output.loom", "col_attrs/CellID")
head(loom_cells_data)
seurat_cells <- colnames(CS1_epi_mes_CS)

common_cells <- intersect(seurat_cells, loom_cells_data)

#Common cells
CS1_epi_mes_CS_filtered <- subset(CS1_epi_mes_CS, cells = common_cells)
dim(CS1_epi_mes_CS_filtered)

# Importing SCENIC Loom Files into Seurat
scenic_output <- ImportPyscenicLoom(scenic_loom_path, seu = CS1_epi_mes_CS_filtered)

#Visualizing-scenic-results

#Viewing the outputs
tf_auc <- scenic_output@misc$SCENIC$RegulonsAUC
head(tf_auc, 4:5)

tf_gene_list <- scenic_output@misc$SCENIC$Regulons
head(tf_gene_list, 40)


# Cell_type check
head(scenic_output$Cell_type)

# Heatmap
tf_zscore <- CalcStats(tf_auc, f = scenic_output$Cell_type, order = "p", n = 6, t = TRUE)
Heatmap(tf_zscore, lab_fill = "zscore")


# Setting the default assay to "TF" for easier access to regulon activity
DefaultAssay(scenic_output) <- "TF"

scenic_output <- SetIdent(scenic_output, value="Epi_mes")

# Creating a waterfall plot to compare regulon activity
WaterfallPlot(
  scenic_output,
  features = rownames(scenic_output),  # Using all available TFs in the "TF" assay
  ident.1 = "Epi",     
  ident.2 = "Mes",    
  exp.transform = FALSE,      
  top.n = 20                 
)

#CS 2

scenic_loom_path <- "/home/silvia/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS2/CS2_scenic_integrated-output.loom"

loom_cells_data <- h5read("~/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS2/CS2_scenic_integrated-output.loom", "col_attrs/CellID")
head(loom_cells_data)
seurat_cells <- colnames(CS2_epi_mes_CS)
common_cells <- intersect(seurat_cells, loom_cells_data)
CS2_epi_mes_CS_filtered <- subset(CS2_epi_mes_CS, cells = common_cells)
dim(CS2_epi_mes_CS_filtered)

# Importing SCENIC Loom Files into Seurat
scenic_output <- ImportPyscenicLoom(scenic_loom_path, seu = CS2_epi_mes_CS_filtered)

#Visualizing-scenic-results

#Viewing the outputs
tf_auc <- scenic_output@misc$SCENIC$RegulonsAUC
head(tf_auc, 4:5)

tf_gene_list <- scenic_output@misc$SCENIC$Regulons
head(tf_gene_list, 40)

# Chequeo que esten los Cell_type
head(scenic_output$Cell_type)

# Heatmap
tf_zscore <- CalcStats(tf_auc, f = scenic_output$Cell_type, order = "p", n = 6, t = TRUE)
Heatmap(tf_zscore, lab_fill = "zscore")

# Setting the default assay to "TF" for easier access to regulon activity
DefaultAssay(scenic_output) <- "TF"

scenic_output <- SetIdent(scenic_output, value="Epi_mes")

# Creating a waterfall plot to compare regulon activity 
WaterfallPlot(
  scenic_output,
  features = rownames(scenic_output),  # Using all available TFs in the "TF" assay
  ident.1 = "Epi",      
  ident.2 = "Mes",     
  exp.transform = FALSE,     
  top.n = 20                 
)

#CS 3

scenic_loom_path <- "/home/silvia/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS3/CS3_scenic_integrated_output.loom"

loom_cells_data <- h5read("~/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS3/CS3_scenic_integrated_output.loom", "col_attrs/CellID")
head(loom_cells_data)
seurat_cells <- colnames(CS3_epi_mes_CS)
common_cells <- intersect(seurat_cells, loom_cells_data)
CS3_epi_mes_CS_filtered <- subset(CS3_epi_mes_CS, cells = common_cells)
dim(CS3_epi_mes_CS_filtered)

# Importing SCENIC Loom Files into Seurat
scenic_output <- ImportPyscenicLoom(scenic_loom_path, seu = CS3_epi_mes_CS_filtered)

#Visualizing-scenic-results

#Viewing the outputs
tf_auc <- scenic_output@misc$SCENIC$RegulonsAUC
head(tf_auc, 4:5)

tf_gene_list <- scenic_output@misc$SCENIC$Regulons
head(tf_gene_list, 40)

# Heatmap
tf_zscore <- CalcStats(tf_auc, f = scenic_output$Cell_type, order = "p", n = 6, t = TRUE)
Heatmap(tf_zscore, lab_fill = "zscore")

# Setting the default assay to "TF" for easier access to regulon activity
DefaultAssay(scenic_output) <- "TF"

scenic_output <- SetIdent(scenic_output, value="Epi_mes")

# Creating a waterfall plot to compare regulon activity 
WaterfallPlot(
  scenic_output,
  features = rownames(scenic_output),  # Using all available TFs in the "TF" assay
  ident.1 = "Epi",      
  ident.2 = "Mes",     
  exp.transform = FALSE,     
  top.n = 20                 
)                  

#CS 4

scenic_loom_path <- "/home/silvia/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS4/CS4_scenic_integrated-output.loom"

loom_cells_data <- h5read("~/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS4/CS4_scenic_integrated-output.loom", "col_attrs/CellID")
head(loom_cells_data)
seurat_cells <- colnames(CS4_epi_mes_CS)
common_cells <- intersect(seurat_cells, loom_cells_data)
CS4_epi_mes_CS_filtered <- subset(CS4_epi_mes_CS, cells = common_cells)
dim(CS4_epi_mes_CS_filtered)

# Importing SCENIC Loom Files into Seurat
scenic_output <- ImportPyscenicLoom(scenic_loom_path, seu = CS4_epi_mes_CS_filtered)

#Visualizing-scenic-results

#Viewing the outputs
tf_auc <- scenic_output@misc$SCENIC$RegulonsAUC
head(tf_auc, 4:5)

tf_gene_list <- scenic_output@misc$SCENIC$Regulons
head(tf_gene_list, 40)

# Heatmap
tf_zscore <- CalcStats(tf_auc, f = scenic_output$Cell_type, order = "p", n = 6, t = TRUE)
Heatmap(tf_zscore, lab_fill = "zscore")

# Setting the default assay to "TF" for easier access to regulon activity
DefaultAssay(scenic_output) <- "TF"

scenic_output <- SetIdent(scenic_output, value="Epi_mes")

# Creating a waterfall plot to compare regulon activity 
WaterfallPlot(
  scenic_output,
  features = rownames(scenic_output),  # Using all available TFs in the "TF" assay
  ident.1 = "Epi",      
  ident.2 = "Mes",     
  exp.transform = FALSE,     
  top.n = 20                 
)

#CS 5

scenic_loom_path <- "/home/silvia/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS5/CS5_scenic_integrated-output.loom"

loom_cells_data <- h5read("~/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS5/CS5_scenic_integrated-output.loom", "col_attrs/CellID")
head(loom_cells_data)
seurat_cells <- colnames(CS5_epi_mes_CS)
common_cells <- intersect(seurat_cells, loom_cells_data)
CS5_epi_mes_CS_filtered <- subset(CS5_epi_mes_CS, cells = common_cells)
dim(CS5_epi_mes_CS_filtered)

scenic_output <- ImportPyscenicLoom(scenic_loom_path, seu = CS5_epi_mes_CS_filtered)

#Visualizing-scenic-results

#Viewing the outputs
tf_auc <- scenic_output@misc$SCENIC$RegulonsAUC
head(tf_auc, 4:5)

tf_gene_list <- scenic_output@misc$SCENIC$Regulons
head(tf_gene_list, 40)

# Heatmap
tf_zscore <- CalcStats(tf_auc, f = scenic_output$Cell_type, order = "p", n = 6, t = TRUE)
Heatmap(tf_zscore, lab_fill = "zscore")

# Setting the default assay to "TF" for easier access to regulon activity
DefaultAssay(scenic_output) <- "TF"

scenic_output <- SetIdent(scenic_output, value="Epi_mes")

# Creating a waterfall plot to compare regulon activity 
WaterfallPlot(
  scenic_output,
  features = rownames(scenic_output),  # Using all available TFs in the "TF" assay
  ident.1 = "Epi",      
  ident.2 = "Mes",     
  exp.transform = FALSE,     
  top.n = 20                 
)

# CS 6

scenic_loom_path <- "/home/silvia/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS6/CS6_scenic_integrated-output.loom"

loom_cells_data <- h5read("~/Documentos/Carcinosarcomas/pyscenic/pyscenic_CS6/CS6_scenic_integrated-output.loom", "col_attrs/CellID")
head(loom_cells_data)

seurat_cells <- colnames(CS6_epi_mes_CS)
common_cells <- intersect(seurat_cells, loom_cells_data)
CS6_epi_mes_CS_filtered <- subset(CS6_epi_mes_CS, cells = common_cells)
dim(CS6_epi_mes_CS_filtered)

# Importing SCENIC Loom Files into Seurat
scenic_output <- ImportPyscenicLoom(scenic_loom_path, seu = CS6_epi_mes_CS_filtered)

#Visualizing-scenic-results

#Viewing the outputs
tf_auc <- scenic_output@misc$SCENIC$RegulonsAUC
head(tf_auc, 4:5)

tf_gene_list <- scenic_output@misc$SCENIC$Regulons
head(tf_gene_list, 40)

# Heatmap
tf_zscore <- CalcStats(tf_auc, f = scenic_output$Cell_type, order = "p", n = 6, t = TRUE)
Heatmap(tf_zscore, lab_fill = "zscore")

# Setting the default assay to "TF" for easier access to regulon activity
DefaultAssay(scenic_output) <- "TF"

scenic_output <- SetIdent(scenic_output, value="Epi_mes")

# Creating a waterfall plot to compare regulon activity 
WaterfallPlot(
  scenic_output,
  features = rownames(scenic_output),  # Using all available TFs in the "TF" assay
  ident.1 = "Epi",      
  ident.2 = "Mes",     
  exp.transform = FALSE,     
  top.n = 20                 
)

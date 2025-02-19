#Seurat Extend SCENIC

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

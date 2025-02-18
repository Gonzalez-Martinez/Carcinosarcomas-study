#InferCNVs carcinosarcomas

library(infercnv)

#Normal epithelial cells as reference

load("~/scRNAseq_SilviaG/Endometrioid/Individual_analysis_endometrium/Normal_epithelial.RData")

#InferCNVs_CS1

##Merge CS1 y Reference cells

Object_to_CNVs <- merge(CS1_epi_mes_CS, y = list(Normal_epithelial))
Object_to_CNVs <- JoinLayers(Object_to_CNVs)
Object_to_CNVs <- NormalizeData(object = Object_to_CNVs)
Object_to_CNVs <- FindVariableFeatures(object = Object_to_CNVs)
Object_to_CNVs <- ScaleData(object = Object_to_CNVs)
Object_to_CNVs <- RunPCA(object = Object_to_CNVs)
ElbowPlot(Object_to_CNVs)
Object_to_CNVs <- FindNeighbors(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- FindClusters(object = Object_to_CNVs, resolution =0.20)
Object_to_CNVs <- RunUMAP(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- SetIdent(Object_to_CNVs, value = "Cell_type")

gene_order_file = "~/Documentos/Scripts_scRNAseq/order"
out_dir_name <- "~/Documentos/Carcinosarcomas/InferCNVs_CS1"

DefaultAssay(Object_to_CNVs) = 'RNA'
counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
sample_annotation = Object_to_CNVs@meta.data[,"Cell_type", drop=FALSE] 

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=sample_annotation,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names=c("Normal Epi"))


infercnv_obj = infercnv::run(infercnv_obj
                             ,cutoff=0.1 
                             ,out_dir=out_dir_name
                             ,denoise=TRUE
                             ,HMM=TRUE
)


#InferCNVs_CS2

##Merge CS2 y Reference cells

Object_to_CNVs <- merge(CS2_epi_mes_CS, y = list(Normal_epithelial))
Object_to_CNVs <- JoinLayers(Object_to_CNVs)
Object_to_CNVs <- NormalizeData(object = Object_to_CNVs)
Object_to_CNVs <- FindVariableFeatures(object = Object_to_CNVs)
Object_to_CNVs <- ScaleData(object = Object_to_CNVs)
Object_to_CNVs <- RunPCA(object = Object_to_CNVs)
ElbowPlot(Object_to_CNVs)
Object_to_CNVs <- FindNeighbors(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- FindClusters(object = Object_to_CNVs, resolution =0.20)
Object_to_CNVs <- RunUMAP(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- SetIdent(Object_to_CNVs, value = "Cell_type")

out_dir_name <- "~/Documentos/Carcinosarcomas/InferCNVs_CS2"

DefaultAssay(Object_to_CNVs) = 'RNA'
counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
sample_annotation = Object_to_CNVs@meta.data[,"Cell_type", drop=FALSE] 

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=sample_annotation,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names=c("Normal Epi"))


infercnv_obj = infercnv::run(infercnv_obj
                             ,cutoff=0.1 
                             ,out_dir=out_dir_name
                             ,denoise=TRUE
                             ,HMM=TRUE
)

#InferCNVs_CS3

##Merge CS3 y Reference cells

Object_to_CNVs <- merge(CS3_epi_mes_CS, y = list(Normal_epithelial))
Object_to_CNVs <- JoinLayers(Object_to_CNVs)
Object_to_CNVs <- NormalizeData(object = Object_to_CNVs)
Object_to_CNVs <- FindVariableFeatures(object = Object_to_CNVs)
Object_to_CNVs <- ScaleData(object = Object_to_CNVs)
Object_to_CNVs <- RunPCA(object = Object_to_CNVs)
ElbowPlot(Object_to_CNVs)
Object_to_CNVs <- FindNeighbors(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- FindClusters(object = Object_to_CNVs, resolution =0.20)
Object_to_CNVs <- RunUMAP(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- SetIdent(Object_to_CNVs, value = "Cell_type")

out_dir_name <- "~/Documentos/Carcinosarcomas/InferCNVs_CS3"

DefaultAssay(Object_to_CNVs) = 'RNA'
counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
sample_annotation = Object_to_CNVs@meta.data[,"Cell_type", drop=FALSE] 

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=sample_annotation,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names=c("Normal Epi"))


infercnv_obj = infercnv::run(infercnv_obj
                             ,cutoff=0.1 
                             ,out_dir=out_dir_name
                             ,denoise=TRUE
                             ,HMM=TRUE
)

#InferCNVs_CS4

##Merge CS4 y Reference cells

Object_to_CNVs <- merge(CS4_epi_mes_CS, y = list(Normal_epithelial))
Object_to_CNVs <- JoinLayers(Object_to_CNVs)
Object_to_CNVs <- NormalizeData(object = Object_to_CNVs)
Object_to_CNVs <- FindVariableFeatures(object = Object_to_CNVs)
Object_to_CNVs <- ScaleData(object = Object_to_CNVs)
Object_to_CNVs <- RunPCA(object = Object_to_CNVs)
ElbowPlot(Object_to_CNVs)
Object_to_CNVs <- FindNeighbors(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- FindClusters(object = Object_to_CNVs, resolution =0.20)
Object_to_CNVs <- RunUMAP(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- SetIdent(Object_to_CNVs, value = "Cell_type")

out_dir_name <- "~/Documentos/Carcinosarcomas/InferCNVs_CS4"

DefaultAssay(Object_to_CNVs) = 'RNA'
counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
sample_annotation = Object_to_CNVs@meta.data[,"Cell_type", drop=FALSE] 

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=sample_annotation,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names=c("Normal Epi"))


infercnv_obj = infercnv::run(infercnv_obj
                             ,cutoff=0.1 
                             ,out_dir=out_dir_name
                             ,denoise=TRUE
                             ,HMM=TRUE
)

#InferCNVs_CS5

##Merge CS5 y Reference cells

Object_to_CNVs <- merge(CS5_epi_mes_CS, y = list(Normal_epithelial))
Object_to_CNVs <- JoinLayers(Object_to_CNVs)
Object_to_CNVs <- NormalizeData(object = Object_to_CNVs)
Object_to_CNVs <- FindVariableFeatures(object = Object_to_CNVs)
Object_to_CNVs <- ScaleData(object = Object_to_CNVs)
Object_to_CNVs <- RunPCA(object = Object_to_CNVs)
ElbowPlot(Object_to_CNVs)
Object_to_CNVs <- FindNeighbors(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- FindClusters(object = Object_to_CNVs, resolution =0.20)
Object_to_CNVs <- RunUMAP(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- SetIdent(Object_to_CNVs, value = "Cell_type")

out_dir_name <- "~/Documentos/Carcinosarcomas/InferCNVs_CS5"

DefaultAssay(Object_to_CNVs) = 'RNA'
counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
sample_annotation = Object_to_CNVs@meta.data[,"Cell_type", drop=FALSE] 

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=sample_annotation,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names=c("Normal Epi"))


infercnv_obj = infercnv::run(infercnv_obj
                             ,cutoff=0.1 
                             ,out_dir=out_dir_name
                             ,denoise=TRUE
                             ,HMM=TRUE
)


#InferCNVs_CS6

##Merge CS6 y Reference cells

Object_to_CNVs <- merge(CS6_epi_mes_CS, y = list(Normal_epithelial))
Object_to_CNVs <- JoinLayers(Object_to_CNVs)
Object_to_CNVs <- NormalizeData(object = Object_to_CNVs)
Object_to_CNVs <- FindVariableFeatures(object = Object_to_CNVs)
Object_to_CNVs <- ScaleData(object = Object_to_CNVs)
Object_to_CNVs <- RunPCA(object = Object_to_CNVs)
ElbowPlot(Object_to_CNVs)
Object_to_CNVs <- FindNeighbors(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- FindClusters(object = Object_to_CNVs, resolution =0.20)
Object_to_CNVs <- RunUMAP(object = Object_to_CNVs, dims = 1:20)
Object_to_CNVs <- SetIdent(Object_to_CNVs, value = "Cell_type")

out_dir_name <- "~/Documentos/Carcinosarcomas/InferCNVs_CS6"

DefaultAssay(Object_to_CNVs) = 'RNA'
counts_matrix <- GetAssayData(Object_to_CNVs, layer = "counts")
sample_annotation = Object_to_CNVs@meta.data[,"Cell_type", drop=FALSE] 

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=sample_annotation,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names=c("Normal Epi"))


infercnv_obj = infercnv::run(infercnv_obj
                             ,cutoff=0.1 
                             ,out_dir=out_dir_name
                             ,denoise=TRUE
                             ,HMM=TRUE
)

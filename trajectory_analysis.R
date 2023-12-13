

setwd("/work/LAS/geetu-lab/sarwar/bcb570_project")

# install monocle
#devtools::install_github("cysouw/qlcMatrix")
#BiocManager::install("monocle")

library(monocle)
# ======================

obj_caf <- readRDS(file = "objCAF.caf_cells.rds")
obj_caf 

mkgene <- read.csv("CAF_int_DE.csv",
         row.names = 1)
#View(mkgene)

mkgene %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> mktop5

#View(mktop5)

data <- as(as.matrix(LayerData(obj_caf, assay="integrated", layer='scale.data')), 'sparseMatrix')
dim(x = data)
data[1:2,1:3]



pd <- new('AnnotatedDataFrame', data = obj_caf@meta.data)
pd[1:5,]

#View(pd@data)

fData <- data.frame(gene_short_name = row.names(data), 
                    row.names = row.names(data))
#fData

fd <- new('AnnotatedDataFrame', data = fData)
fd[1:2,]
fd@varMetadata

#Construct monocle cds
HSMM2 <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        #lowerDetectionLimit = 0.5,
                        expressionFamily =uninormal())


# gene selection
ordering_genes <-  mktop5$gene

HSMM2 <- setOrderingFilter(HSMM2, ordering_genes)
print(dim(exprs(HSMM2)))

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM2 <- reduceDimension(HSMM2, norm_method="none",
                         reduction_method="DDRTree",
                         max_components=2,
                         scaling=TRUE,
                         verbose=TRUE,
                         pseudo_expr=0)
#View data
pData(HSMM2)$celltype = obj_caf@meta.data$seurat_clusters

pData(HSMM2)
fData(HSMM2) # gene attribute

HSMM2 <- orderCells(HSMM2)
HSMM2


# .....trajectory plot

plot_cell_trajectory(HSMM2, 
                     color_by = "celltype",
                     theta = -15,
                     show_branch_points = T,
                     show_tree = TRUE,
                     cell_size = 1)

plot_cell_trajectory(HSMM2, 
                     color_by = "Pseudotime",
                     theta = -15,
                     show_branch_points = T,
                     show_tree = TRUE,
                     cell_size = 1.3)

plot_cell_trajectory(HSMM2, 
                     color_by = "State",
                     theta = -15,
                     show_branch_points = F,
                     show_tree = TRUE,
                     cell_size = 1.3)


VlnPlot(obj_caf, features = c("COL1A2", 
                               "FAP", 
                              "PDPN", 
                               "DCN", 
                              "COL3A1", 
                               "COL6A1"))

# https://biomarkerres.biomedcentral.com/articles/10.1186/s40364-022-00406-z/tables/1
# icc mCAFs marker: cluste 1

VlnPlot(obj_caf, features = c("COL5A1", 
                              "COL5A2", 
                              "COL6A3", 
                              "POSTN", 
                              "FN1", 
                              "LUM", 
                              "DCN",
                              "VCAN"))




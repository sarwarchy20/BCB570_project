
library(clustree)
library(ggplot2)

# CAF cell analysis
objCAF <- readRDS(file = "obj_integrated_all.rds")
objCAF

caf_cells <- WhichCells(object = objCAF, idents = "CAF")
length(caf_cells)


# *================ this "set1_sort.mtx" data from main_analysis.R, run before this
caf_cells.set1 <- set1_sort.mtx[, intersect(caf_cells,
                                        colnames(set1_sort.mtx))]
dim(caf_cells.set1)

# *================ this "set1_sort.mtx" data from main_analysis.R, run before this
caf_cells.set2 <- set2_sort.mtx[,intersect(caf_cells,
                                        colnames(set2_sort.mtx))]

dim(caf_cells.set2)

# Seurat 
Obj.set1.caf_cells<- CreateSeuratObject(counts = caf_cells.set1, 
                              project = "set1",
                              min.cells = 5, 
                              min.features = 500)
Obj.set1.caf_cells

Obj.set1.caf_cells <- NormalizeData(Obj.set1.caf_cells, normalization.method = "LogNormalize", scale.factor = 10000)

Obj.set1.caf_cells<- FindVariableFeatures(Obj.set1.caf_cells, selection.method = 'mean.var.plot', 
                                       mean.cutoff = c(0.0125, 3), 
                                       dispersion.cutoff = c(0.5, Inf),
                                       nfeatures = Inf)

Obj.set1.caf_cells

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Obj.set1.caf_cells), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Obj.set1.caf_cells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
plot2 <- plot2 + 
  xlab("Average Expression") +
  ylab("Dispersion")
plot2

# set2 T cells 
Obj.set2.caf_cells<- CreateSeuratObject(counts = caf_cells.set2, 
                                     project = "set2",
                                     min.cells = 5, 
                                     min.features = 500)
Obj.set2.caf_cells


Obj.set2.caf_cells <- NormalizeData(Obj.set2.caf_cells, normalization.method = "LogNormalize", scale.factor = 10000)

Obj.set2.caf_cells<- FindVariableFeatures(Obj.set2.caf_cells, selection.method = 'mean.var.plot', 
                                       mean.cutoff = c(0.0125, 3), 
                                       dispersion.cutoff = c(0.5, Inf),
                                       nfeatures = Inf)

Obj.set2.caf_cells

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(Obj.set2.caf_cells), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Obj.set2.caf_cells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2 <- plot2 + 
  xlab("Average Expression") +
  ylab("Dispersion")
plot2


obj.list = c(Obj.set1.caf_cells, Obj.set2.caf_cells)


anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                  normalization.method = "LogNormalize",
                                  anchor.features = 8000, dims = 1:30,
                                  reduction = "cca")

objCAF.caf_cells<- IntegrateData(anchorset = anchors,
                          normalization.method = "LogNormalize")

DefaultAssay(object = objCAF.caf_cells) <- "integrated"


objCAF.caf_cells <- ScaleData(object = objCAF.caf_cells)

# Dimensionality Reduction

length(VariableFeatures(object = objCAF.caf_cells))

objCAF.caf_cells = RunPCA(objCAF.caf_cells, 
                        features = VariableFeatures(object = objCAF.caf_cells),
                  npcs = 30)


DimPlot(objCAF.caf_cells, group.by = "orig.ident", reduction = "pca") +
  ggtitle("")

DimPlot(objCAF.caf_cells, reduction = "pca", split.by = "orig.ident")+
  ggtitle("")


# ===== determine significant PCS
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html

ElbowPlot(objCAF.caf_cells, ndims = 30)


pct <- objCAF.caf_cells[["pca"]]@stdev / sum(objCAF.caf_cells[["pca"]]@stdev) * 100

cumu <- cumsum(pct)

co1 <- which(cumu > 90 & pct < 5)[1]

co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
pcs <- min(co1, co2)

pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

#View(plot_df)

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw() +
  ylab("standard deviation (%)") +
  xlab("cumulative standard deviation (%)") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 5)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5))


DimHeatmap(objCAF.caf_cells, dims = 1:9, cells = 500, balanced = TRUE)
DimHeatmap(objCAF.caf_cells, dims = 10:18, cells = 500, balanced = TRUE)

objCAF.caf_cells <- FindNeighbors(objCAF.caf_cells, dims = 1:9)

objCAF.caf_cells <- FindClusters(objCAF.caf_cells, 
                                 resolution = 0.2)

# ======================= cluster tree 

#install.packages("clustree")

# https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
#jpeg("tree.jpg", width = 800, height = 2000)

clustree(objCAF.caf_cells, prefix = "integrated_snn_res.")

#dev.off()

# ========================
# Resolution= 0.08, min.dist = 0.8, n.neighbors = 20L
objCAF.caf_cells<- RunUMAP(objCAF.caf_cells, reduction = "pca", 
                    dims = 1:9,verbose=T,
                    min.dist = 0.8,
                    n.neighbors = 20L)

DimPlot(objCAF.caf_cells, label=T,
        reduction = "umap",pt.size=1.5,
        label.size = 8)


#=================== bar plot

bar.cells <- objCAF.caf_cells@meta.data
#View(bar.cells)
bar.cells <- table(select(bar.cells,seurat_clusters))

tab_anno <- data.frame(bar.cells)
colnames(tab_anno)<- c("Types","Freq")

#View(tab_anno)

colnames(tab_anno)<- c("Types","Freq")

ggplot(tab_anno,aes(x=Types,y=Freq, fill=Types)) + 
  geom_col(position="dodge") +
  geom_text(aes(label=Freq), 
            position=position_dodge(width=0.9),
            vjust=-0.25) + 
  xlab("Clusters") + ylab("Number of Cells")+
  theme(legend.position="none")

# ======================
#install.packages('ape') for BuildClusterTree

objCAF.caf_cells@active.ident

objCAF.caf_cells@active.assay

#saveRDS(objCAF.caf_cells, file = "obj_integrated_CAF.rds")

all.markers.caf_cells = FindAllMarkers(objCAF.caf_cells,
                             test.use = "wilcox", logfc.threshold = 0.25,
                             only.pos = FALSE)
dim(all.markers.caf_cells)

caf.all.de <- all.markers.caf_cells

#write.csv(all.markers.caf_cells, "CAF_int_DE.csv",
   #     row.names = T)

all.markers.caf_cells %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

View(top10)

# heatmap

DoHeatmap(objCAF.caf_cells, features = top10$gene)

# cluster 0
clst0.marker.up <- dplyr::filter(all.markers.caf_cells, 
                              cluster=="2",avg_log2FC >1)
#View(clst0.marker.up)

FeaturePlot(objCAF.caf_cells, features = clst0.marker.up$gene[1:6])

VlnPlot(objCAF.caf_cells, features = clst0.marker.up$gene[1:6])

saveRDS(objCAF.caf_cells, file = "objCAF.caf_cells.rds")

# =============================

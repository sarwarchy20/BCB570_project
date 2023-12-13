
library(dplyr)
#install.packages("Rcpp") 
library(Seurat)
#remove.packages("Seurat")
#install.packages("Seurat")

library(ggplot2)

library(patchwork)
library(tibble)
library("stringr")

# ======= Set1 pre-processing =================================

set1 <- Read10X(data.dir = "GSE125449/set1/") # read Set1 data
dim(set1) # 20124 X 5115
colnames(set1)
GSE125449.set1.mtx = as.data.frame(set1)
#View(GSE125449.set1.mtx[1:2,1:3])

set1.samples<-read.table("GSE125449/set1/GSE125449_Set1_samples.txt",
                         header=TRUE,
                         sep ='\t')
dim(set1.samples)
#View(set1.samples)

#levels(as.factor(set1.samples$Sample))

# Assigne sample id
set1.samples[set1.samples == "S16_P10_LCP18"] <- "S1H18"
set1.samples[set1.samples == "S02_P01_LCP21"] <- "S1H21"
set1.samples[set1.samples == "S07_P02_LCP28"] <- "S1H28"
set1.samples[set1.samples == "S08_P03_LCP26"] <- "S1C26"
set1.samples[set1.samples == "S09_P04_LCP25"] <- "S1C25"
set1.samples[set1.samples == "S10_P05_LCP23"] <- "S1H23"
set1.samples[set1.samples == "S11_P06_LCP29"] <- "S1C29"
set1.samples[set1.samples == "S12_P07_LCP30"] <- "S1H30"
set1.samples[set1.samples == "S15_P09_LCP38"] <- "S1H38"
set1.samples[set1.samples == "S21_P13_LCP37"] <- "S1H37"
set1.samples[set1.samples == "S20_P12_LCP35"] <- "S1C35"
set1.samples[set1.samples == "S19_P11_LCP39"] <- "S1C39"

set1.samples<- dplyr:: filter(set1.samples, Type != "unclassified")
#View(set1.samples)
# sort columns 
set1_sort.mtx <- GSE125449.set1.mtx[,match(set1.samples$Cell_Barcode,
                           colnames(GSE125449.set1.mtx))]
dim(set1_sort.mtx)
#View(head(set1_sort.mtx))

# Seurat 
Obj.set1<- CreateSeuratObject(counts = set1_sort.mtx, 
                              project = "set1",
                              meta.data = set1.samples,
                              min.cells = 5, 
                              min.features = 500)
Obj.set1

Obj.set1[["percent.mt"]] <- PercentageFeatureSet(Obj.set1, 
                                                 pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Obj.set1, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

set1.plot1 <- FeatureScatter(Obj.set1, feature1 = "nCount_RNA", feature2 = "percent.mt")
set1.plot2 <- FeatureScatter(Obj.set1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

jpeg("all-pre-2.jpeg", width=620)
set1.plot1 + set1.plot2
dev.off()

# filtering
Obj.set1 <- subset(Obj.set1,
                   subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20) 
Obj.set1

# after filtering 
set1.plot1 <- FeatureScatter(Obj.set1, feature1 = "nCount_RNA", feature2 = "percent.mt")
set1.plot2 <- FeatureScatter(Obj.set1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

jpeg("all-pre-2_flt.jpeg", width=620)
set1.plot1 + set1.plot2
dev.off()


Obj.set1 <- NormalizeData(Obj.set1, normalization.method = "LogNormalize", scale.factor = 10000)

Obj.set1<- FindVariableFeatures(Obj.set1, selection.method = "vst", nfeatures = 2244)

Obj.set1

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Obj.set1), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Obj.set1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = F)

jpeg("all-hvg-set1.jpeg", width=800)
plot1 + plot2
dev.off()

#all.genes <- rownames(Obj.set1)
#Obj.set1 <- ScaleData(Obj.set1, features = all.genes)


# ================================================
#================= Set2 pre-processing ===========

set2 <- Read10X(data.dir = "GSE125449/set2/") # read set2 data
dim(set2) # 19572 X 4831
colnames(set2)
GSE125449.set2.mtx = as.data.frame(set2)
#View(GSE125449.set2.mtx[1:2,1:3])

set2.samples <-read.table("GSE125449/set2/GSE125449_Set2_samples.txt",
                         header=TRUE,
                         sep ='\t')
dim(set2.samples)
#View(set2.samples)

#levels(as.factor(set2.samples$Sample))

# Assigne sample id
set2.samples[set2.samples == "S300_P02_LCP60"] <- "S2C60"
set2.samples[set2.samples == "S305_P06_LCP56"] <- "S2C56"
set2.samples[set2.samples == "S351_P10_LCP34"] <- "S2H34"
set2.samples[set2.samples == "S355_P13_LCP42"] <- "S2C42"
set2.samples[set2.samples == "S358_P16_LCP46"] <- "S2C46"
set2.samples[set2.samples == "S364_P21_LCP65"] <- "S2H65"
set2.samples[set2.samples == "S365_P22_LCP66"] <- "S2C66"

set2.samples<- dplyr:: filter(set2.samples, Type != "unclassified")
#View(set2.samples)
# sort columns 
set2_sort.mtx <- GSE125449.set2.mtx[,match(set2.samples$Cell.Barcode,
                                                       colnames(GSE125449.set2.mtx))]
dim(set2_sort.mtx)
#View(head(set2_sort.mtx))

# Seurat 
Obj.set2<- CreateSeuratObject(counts = set2_sort.mtx, 
                              project = "set2",
                              meta.data = set2.samples,
                              min.cells = 5, 
                              min.features = 500)
Obj.set2

Obj.set2[["percent.mt"]] <- PercentageFeatureSet(Obj.set2, 
                                                 pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Obj.set2, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

set2.plot1 <- FeatureScatter(Obj.set2, feature1 = "nCount_RNA", feature2 = "percent.mt")
set2.plot2 <- FeatureScatter(Obj.set2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

jpeg("all-pre-3_set2.jpeg", width=620)
set2.plot1 + set2.plot2
dev.off()


Obj.set2 <- subset(Obj.set2,
                   subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20) 
Obj.set2

set2.plot1 <- FeatureScatter(Obj.set2, feature1 = "nCount_RNA", feature2 = "percent.mt")
set2.plot2 <- FeatureScatter(Obj.set2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

jpeg("all-pre-4_set2.jpeg", width=620)
set2.plot1 + set2.plot2
dev.off()

Obj.set2 <- NormalizeData(Obj.set2, normalization.method = "LogNormalize", scale.factor = 10000)

Obj.set2<- FindVariableFeatures(Obj.set2, selection.method = "vst", nfeatures = 2244)

Obj.set2

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(Obj.set2), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Obj.set2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = F)

jpeg("all-hvg-set2.jpeg", width=800)
plot1 + plot2
dev.off()

#all.genes <- rownames(Obj.set2)
#Obj.set2 <- ScaleData(Obj.set2, features = all.genes)


#....................INTEGRATION: Set1 and Set2

obj.list = c(Obj.set1, Obj.set2)


anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                  normalization.method = "LogNormalize",
                                  anchor.features = 8000, dims = 1:30,
                                  reduction = "cca")

obj.intg <- IntegrateData(anchorset = anchors,
                          normalization.method = "LogNormalize")

DefaultAssay(object = obj.intg) <- "integrated"

length(VariableFeatures(object = obj.intg))

obj.intg <- ScaleData(object = obj.intg)

set2.plot1 <- FeatureScatter(obj.intg, feature1 = "nCount_RNA", feature2 = "percent.mt")
set2.plot2 <- FeatureScatter(obj.intg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

jpeg("int_1.jpeg", width=620)
set2.plot1 + set2.plot2
dev.off()

# Dimensionality Reduction

obj.intg = RunPCA(obj.intg, features = VariableFeatures(object = obj.intg),
                  npcs = 30)

jpeg("int_pca1.jpeg", width=600)
DimPlot(obj.intg, group.by = "orig.ident", 
        reduction = "pca") 
dev.off()

jpeg("int_pca2.jpeg", width=800)
DimPlot(obj.intg, reduction = "pca", split.by = "orig.ident")
dev.off()

jpeg("int_pca3.jpeg", width=600)
ElbowPlot(obj.intg, ndims = 30)
dev.off()

# ============= PC variance explained 
#mat <- LayerData(obj.intg, assay="integrated", layer='scale.data')
#dim(mat) # 2903 9365

pca.val <- obj.intg[["pca"]]

# Get the total variance:
#total_variance <- obj.intg[["pca"]]@misc$total.variance

#eigValues = (pca.val@stdev)^2  ## EigenValues
#varExplained = eigValues / total_variance

#======================

# ===== determine significant PCS
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- obj.intg[["pca"]]@stdev / sum(obj.intg[["pca"]]@stdev) * 100

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
jpeg("int_pca4.jpeg", width=820)

ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw() +
  ylab("standard deviation (%)") +
  xlab("cumulative standard deviation (%)") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 5)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5))

dev.off()

# prcomp
#View(head(obj.intg@assays$integrated@scale.data))

#pcaRes <-  prcomp(t(x = obj.intg@assays$integrated@scale.data), center = F, scale. = F)

jpeg("int_pca5.jpeg", width=820)
DimHeatmap(obj.intg, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

jpeg("int_pca6.jpeg", width=820)
DimHeatmap(obj.intg, dims = 10:18, cells = 500, balanced = TRUE)
dev.off()

jpeg("int_pca7.jpeg", width=820)
DimHeatmap(obj.intg, dims = 19:27, cells = 500, balanced = TRUE)
dev.off()

jpeg("int_pca8.jpeg", width=820)
DimHeatmap(obj.intg, dims = 28:30, cells = 500, balanced = TRUE)
dev.off()

# ============== UMAP
# change n.neighbors: 5 to 50 suggested
obj.intg <- RunUMAP(obj.intg, reduction = "pca", 
                    dims = 1:20,verbose=T,
                    min.dist = 0.3,
                    n.neighbors = 30L)

jpeg("int_umap1.jpeg", width=600)
DimPlot(obj.intg, label=F, group.by = "Type",
        reduction = "umap") + 
  ggtitle("min.dist = 0.3;neighbors = 30")
dev.off()


obj.intg <- RunUMAP(obj.intg, reduction = "pca", 
                    dims = 1:20,verbose=T,
                    min.dist = 0.3,
                    n.neighbors = 20L)

jpeg("int_umap2.jpeg", width=600)
DimPlot(obj.intg, label=F, group.by = "Type",
        reduction = "umap")+ 
  ggtitle("min.dist = 0.3;neighbors = 20")
dev.off()



obj.intg <- RunUMAP(obj.intg, reduction = "pca", 
                    dims = 1:20,verbose=T,
                    min.dist = 0.3,
                    n.neighbors = 10L)

jpeg("int_umap3.jpeg", width=600)
DimPlot(obj.intg, label=F, group.by = "Type",
        reduction = "umap")+ 
  ggtitle("min.dist = 0.3;neighbors = 10")
dev.off()

# change min.dist: 0.001 to 0.5 suggested

obj.intg <- RunUMAP(obj.intg, reduction = "pca", 
                    dims = 1:20,verbose=T,
                    min.dist = 0.2,
                    n.neighbors = 30L)

jpeg("int_umap4.jpeg", width=600)
DimPlot(obj.intg, label=F, group.by = "Type",
        reduction = "umap")+ 
  ggtitle("min.dist = 0.2;neighbors = 30")
dev.off()


obj.intg <- RunUMAP(obj.intg, reduction = "pca", 
                    dims = 1:20,verbose=T,
                    min.dist = 0.4,
                    n.neighbors = 30L)

jpeg("int_umap5.jpeg", width=600)
DimPlot(obj.intg, label=F, group.by = "Type",
        reduction = "umap")+ 
  ggtitle("min.dist = 0.4;neighbors = 30")
dev.off()



obj.intg <- RunUMAP(obj.intg, reduction = "pca", 
                    dims = 1:20,verbose=T,
                    min.dist = 0.5,
                    n.neighbors = 30L)
jpeg("int_umap6.jpeg", width=600)
DimPlot(obj.intg, label=F, group.by = "Type",
        reduction = "umap")+ 
  ggtitle("min.dist = 0.5;neighbors = 30")
dev.off()



# final UMAP plots
obj.intg <- RunUMAP(obj.intg, reduction = "pca", 
                    dims = 1:20,verbose=T,
                    min.dist = 0.3,
                    n.neighbors = 20L)

#jpeg("int_umap_f1.jpeg", width=600)
DimPlot(obj.intg, label=F, group.by = "Type",
        reduction = "umap") + ggtitle("")
#dev.off()


DimPlot(obj.intg, label=F, group.by = "orig.ident",
        reduction = "umap") + ggtitle("")

DimPlot(obj.intg, label=F, group.by = "Type",reduction = "pca")

DimPlot(obj.intg, reduction = "umap", split.by = "orig.ident")

DimPlot(obj.intg, label=T, group.by = "Type",
        reduction = "umap") + ggtitle("")

# percentage bar plot==============

meta.int <- obj.intg@meta.data

meta.int<-  select(meta.int, Sample, Type)

#View(meta.int)

hcc_barplot_Dat <- meta.int %>% group_by(Sample, Type) %>% 
  summarise(frequency = n()) %>%
  mutate(C = sum(frequency)) %>%
  mutate(percent = frequency/C*100)

#View(hcc_barplot_Dat)

ggplot(hcc_barplot_Dat, aes(x = Sample, y = percent, fill = Type))+
  geom_bar(stat = "identity")+scale_y_continuous(labels=scales::percent_format(scale = 1))+
  theme(axis.text = element_text(color="black",size=12),
        axis.text.x=element_text(angle=60),
        axis.title = element_blank(),
        legend.title = element_blank())

#=================== bar plot

bar.cells <- obj.intg@meta.data
#View(bar.cells)
bar.cells <- table(select(bar.cells, Type))

tab_anno <- data.frame(bar.cells)
colnames(tab_anno)<- c("Types","Freq")

#View(tab_anno)

colnames(tab_anno)<- c("Types","Freq")

ggplot(tab_anno,aes(x=Types,y=Freq, fill=Types)) + 
  geom_col(position="dodge") +
  geom_text(aes(label=Freq), 
            position=position_dodge(width=0.9),
            vjust=-0.25) + 
  xlab("Cell types") + ylab("Number of Cells")+
  theme(legend.position="none")

# ===============

obj.intg@assays$RNA

obj.intg@assays$integrated

obj.intg@active.assay

table(obj.intg$orig.ident)

#DefaultAssay(obj.intg) <- "RNA"

DefaultAssay(obj.intg) <- "integrated"


obj.intg = SetIdent(obj.intg, value = obj.intg$Type)

obj.intg@active.ident

 saveRDS(obj.intg, file = "obj_integrated_all.rds")

caf.markers <- FindMarkers(obj.intg, ident.1 = "CAF", 
                           ident.2 =c("B cell",
                                      "HPC-like",
                                      "T cell",
                                      "TAM",
                                      "TEC"))
dim(caf.markers)

#write.csv(caf.markers, "CAF_all_DE.csv",
 #         row.names = T)

caf.markers.flt<- dplyr::filter(caf.markers,
                                avg_log2FC > 1, 
                                p_val_adj < 0.05)
dim(caf.markers.flt)

#write.csv(caf.markers.flt, "CAF_all_DE_flt.csv",
  #        row.names = T)

#View(caf.markers.flt)

#png("int_hmap.png", width = 500, height = 800)
DoHeatmap(obj.intg, 
          features = rownames(caf.markers.flt)[1:50])
#dev.off()

# =================

#all.markers = FindAllMarkers(obj.intg,
   #                          test.use = "wilcox", logfc.threshold = 0.25,
   #                                only.pos = FALSE)
#View(all.markers)

#all.markers %>%
#  group_by(cluster) %>%
 # dplyr::filter(avg_log2FC > 1) %>%
 # slice_head(n = 5) %>%
 # ungroup() -> top10

#DoHeatmap(obj.intg, features = top10$gene) + NoLegend()

# =============
─ Session info ────────────────────────────────────────────────────
 setting  value
 version  R version 4.1.0 (2021-05-18)
 os       Ubuntu 20.04.3 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Etc/UTC
 date     2023-12-13
 rstudio  1.4.1717 Juliet Rose (server)
 pandoc   2.11.4 @ /usr/local/bin/pandoc

─ Packages ────────────────────────────────────────────────────────
 package          * version date (UTC) lib source
 abind              1.4-5   2016-07-21 [2] RSPM (R 4.1.0)
 Biobase          * 2.54.0  2021-10-26 [1] Bioconductor
 BiocGenerics     * 0.40.0  2021-10-26 [1] Bioconductor
 cachem             1.0.5   2021-05-15 [2] RSPM (R 4.1.0)
 callr              3.7.3   2022-11-02 [1] CRAN (R 4.1.0)
 cli              * 3.6.1   2023-03-23 [1] CRAN (R 4.1.0)
 cluster            2.1.2   2021-04-17 [3] CRAN (R 4.1.0)
 clustree         * 0.5.1   2023-11-05 [1] CRAN (R 4.1.0)
 codetools          0.2-18  2020-11-04 [3] CRAN (R 4.1.0)
 colorspace         2.0-2   2021-06-24 [2] RSPM (R 4.1.0)
 combinat           0.0-8   2012-10-29 [1] CRAN (R 4.1.0)
 cowplot            1.1.1   2020-12-30 [1] CRAN (R 4.1.0)
 crayon             1.4.1   2021-02-08 [2] RSPM (R 4.1.0)
 data.table       * 1.14.0  2021-02-21 [2] RSPM (R 4.1.0)
 DDRTree          * 0.1.5   2017-04-30 [1] CRAN (R 4.1.0)
 deldir             1.0-9   2023-05-17 [1] CRAN (R 4.1.0)
 densityClust       0.3.2   2022-03-06 [1] CRAN (R 4.1.0)
 devtools           2.4.5   2022-10-11 [1] CRAN (R 4.1.0)
 digest             0.6.33  2023-07-07 [1] CRAN (R 4.1.0)
 docopt             0.7.1   2020-06-24 [2] RSPM (R 4.1.0)
 dotCall64          1.1-0   2023-10-17 [1] CRAN (R 4.1.0)
 dplyr            * 1.1.3   2023-09-03 [1] CRAN (R 4.1.0)
 ellipsis           0.3.2   2021-04-29 [2] RSPM (R 4.1.0)
 fansi              0.5.0   2021-05-25 [2] RSPM (R 4.1.0)
 farver             2.1.0   2021-02-28 [2] RSPM (R 4.1.0)
 fastDummies        1.7.3   2023-07-06 [1] CRAN (R 4.1.0)
 fastICA            1.2-4   2023-11-27 [1] CRAN (R 4.1.0)
 fastmap            1.1.0   2021-01-25 [2] RSPM (R 4.1.0)
 fitdistrplus       1.1-11  2023-04-25 [1] CRAN (R 4.1.0)
 FNN                1.1.3   2019-02-15 [2] RSPM (R 4.1.0)
 fs                 1.6.3   2023-07-20 [1] CRAN (R 4.1.0)
 future             1.33.0  2023-07-01 [1] CRAN (R 4.1.0)
 future.apply       1.11.0  2023-05-21 [1] CRAN (R 4.1.0)
 generics           0.1.3   2022-07-05 [1] CRAN (R 4.1.0)
 ggforce            0.4.1   2022-10-04 [1] CRAN (R 4.1.0)
 ggplot2          * 3.4.4   2023-10-12 [1] CRAN (R 4.1.0)
 ggraph           * 2.1.0   2022-10-09 [1] CRAN (R 4.1.0)
 ggrepel            0.9.4   2023-10-13 [1] CRAN (R 4.1.0)
 ggridges           0.5.4   2022-09-26 [1] CRAN (R 4.1.0)
 globals            0.16.2  2022-11-21 [1] CRAN (R 4.1.0)
 glue               1.6.2   2022-02-24 [1] CRAN (R 4.1.0)
 goftest            1.2-2   2019-12-02 [2] RSPM (R 4.1.0)
 graphlayouts       1.0.2   2023-11-03 [1] CRAN (R 4.1.0)
 gridExtra          2.3     2017-09-09 [2] RSPM (R 4.1.0)
 gtable             0.3.0   2019-03-25 [2] RSPM (R 4.1.0)
 HSMMSingleCell     1.14.0  2021-10-30 [1] Bioconductor
 htmltools          0.5.6.1 2023-10-06 [1] CRAN (R 4.1.0)
 htmlwidgets        1.5.3   2020-12-10 [2] RSPM (R 4.1.0)
 httpuv             1.6.1   2021-05-07 [2] RSPM (R 4.1.0)
 httr               1.4.2   2020-07-20 [2] RSPM (R 4.1.0)
 ica                1.0-3   2022-07-08 [1] CRAN (R 4.1.0)
 igraph             1.5.1   2023-08-10 [1] CRAN (R 4.1.0)
 irlba            * 2.3.5.1 2022-10-03 [1] CRAN (R 4.1.0)
 jsonlite           1.8.7   2023-06-29 [1] CRAN (R 4.1.0)
 KernSmooth         2.23-20 2021-05-03 [3] CRAN (R 4.1.0)
 later              1.2.0   2021-04-23 [2] RSPM (R 4.1.0)
 lattice            0.20-44 2021-05-02 [3] CRAN (R 4.1.0)
 lazyeval           0.2.2   2019-03-15 [2] RSPM (R 4.1.0)
 leiden             0.4.3   2022-09-10 [1] CRAN (R 4.1.0)
 lifecycle          1.0.3   2022-10-07 [1] CRAN (R 4.1.0)
 limma              3.50.3  2022-04-07 [1] Bioconductor
 listenv            0.8.0   2019-12-05 [2] RSPM (R 4.1.0)
 lmtest             0.9-40  2022-03-21 [1] CRAN (R 4.1.0)
 magrittr         * 2.0.3   2022-03-30 [1] CRAN (R 4.1.0)
 MASS               7.3-54  2021-05-03 [3] CRAN (R 4.1.0)
 Matrix           * 1.6-1.1 2023-09-18 [1] CRAN (R 4.1.0)
 matrixStats        1.0.0   2023-06-02 [1] CRAN (R 4.1.0)
 memoise            2.0.1   2021-11-26 [1] CRAN (R 4.1.0)
 MERINGUE         * 1.0     2023-10-21 [1] Github (JEFworks-Lab/MERINGUE@ca9e2cc)
 mgcv               1.8-35  2021-04-18 [3] CRAN (R 4.1.0)
 mime               0.11    2021-06-23 [2] RSPM (R 4.1.0)
 miniUI             0.1.1.1 2018-05-18 [2] RSPM (R 4.1.0)
 monocle          * 2.22.0  2021-10-26 [1] Bioconductor
 munsell            0.5.0   2018-06-12 [2] RSPM (R 4.1.0)
 nlme               3.1-152 2021-02-04 [3] CRAN (R 4.1.0)
 parallelly         1.36.0  2023-05-26 [1] CRAN (R 4.1.0)
 patchwork        * 1.1.3   2023-08-14 [1] CRAN (R 4.1.0)
 pbapply            1.7-2   2023-06-27 [1] CRAN (R 4.1.0)
 pheatmap           1.0.12  2019-01-04 [1] CRAN (R 4.1.0)
 pillar             1.9.0   2023-03-22 [1] CRAN (R 4.1.0)
 pkgbuild           1.4.2   2023-06-26 [1] CRAN (R 4.1.0)
 pkgconfig          2.0.3   2019-09-22 [2] RSPM (R 4.1.0)
 pkgload            1.3.3   2023-09-22 [1] CRAN (R 4.1.0)
 plotly             4.10.3  2023-10-21 [1] CRAN (R 4.1.0)
 plyr               1.8.6   2020-03-03 [2] RSPM (R 4.1.0)
 png                0.1-7   2013-12-03 [2] RSPM (R 4.1.0)
 polyclip           1.10-0  2019-03-14 [2] RSPM (R 4.1.0)
 prettyunits        1.1.1   2020-01-24 [2] RSPM (R 4.1.0)
 processx           3.8.2   2023-06-30 [1] CRAN (R 4.1.0)
 profvis            0.3.8   2023-05-02 [1] CRAN (R 4.1.0)
 progressr          0.14.0  2023-08-10 [1] CRAN (R 4.1.0)
 promises           1.2.0.1 2021-02-11 [2] RSPM (R 4.1.0)
 ps                 1.7.5   2023-04-18 [1] CRAN (R 4.1.0)
 purrr              1.0.2   2023-08-10 [1] CRAN (R 4.1.0)
 qlcMatrix          0.9.7   2023-12-12 [1] Github (cysouw/qlcMatrix@d919872)
 R6                 2.5.1   2021-08-19 [1] CRAN (R 4.1.0)
 RANN               2.6.1   2019-01-08 [1] CRAN (R 4.1.0)
 RColorBrewer       1.1-2   2014-12-07 [2] RSPM (R 4.1.0)
 Rcpp               1.0.11  2023-07-06 [1] CRAN (R 4.1.0)
 RcppAnnoy          0.0.21  2023-07-02 [1] CRAN (R 4.1.0)
 RcppHNSW           0.5.0   2023-09-19 [1] CRAN (R 4.1.0)
 remotes            2.4.2.1 2023-07-18 [1] CRAN (R 4.1.0)
 reshape2           1.4.4   2020-04-09 [1] CRAN (R 4.1.0)
 reticulate         1.34.0  2023-10-12 [1] CRAN (R 4.1.0)
 rlang              1.1.1   2023-04-28 [1] CRAN (R 4.1.0)
 ROCR               1.0-11  2020-05-02 [1] CRAN (R 4.1.0)
 RSpectra           0.16-1  2022-04-24 [1] CRAN (R 4.1.0)
 rstudioapi         0.13    2020-11-12 [2] RSPM (R 4.1.0)
 Rtsne              0.16    2022-04-17 [1] CRAN (R 4.1.0)
 scales             1.2.1   2022-08-20 [1] CRAN (R 4.1.0)
 scattermore        1.2     2023-06-12 [1] CRAN (R 4.1.0)
 sctransform        0.4.1   2023-10-19 [1] CRAN (R 4.1.0)
 sessioninfo        1.2.2   2021-12-06 [1] CRAN (R 4.1.0)
 Seurat           * 5.0.0   2023-11-04 [1] CRAN (R 4.1.0)
 SeuratObject     * 5.0.0   2023-10-26 [1] CRAN (R 4.1.0)
 shiny              1.6.0   2021-01-25 [2] RSPM (R 4.1.0)
 slam               0.1-50  2022-01-08 [1] CRAN (R 4.1.0)
 sp               * 2.1-1   2023-10-16 [1] CRAN (R 4.1.0)
 spam               2.10-0  2023-10-23 [1] CRAN (R 4.1.0)
 sparsesvd          0.2-2   2023-01-14 [1] CRAN (R 4.1.0)
 spatstat.data      3.0-3   2023-10-24 [1] CRAN (R 4.1.0)
 spatstat.explore   3.2-5   2023-10-22 [1] CRAN (R 4.1.0)
 spatstat.geom      3.2-7   2023-10-20 [1] CRAN (R 4.1.0)
 spatstat.random    3.2-1   2023-10-21 [1] CRAN (R 4.1.0)
 spatstat.sparse    3.0-3   2023-10-24 [1] CRAN (R 4.1.0)
 spatstat.utils     3.0-4   2023-10-24 [1] CRAN (R 4.1.0)
 stringi            1.7.3   2021-07-16 [2] RSPM (R 4.1.0)
 stringr          * 1.5.0   2022-12-02 [1] CRAN (R 4.1.0)
 survival           3.2-11  2021-04-26 [3] CRAN (R 4.1.0)
 tensor             1.5     2012-05-05 [2] RSPM (R 4.1.0)
 tibble           * 3.2.1   2023-03-20 [1] CRAN (R 4.1.0)
 tidygraph          1.2.3   2023-02-01 [1] CRAN (R 4.1.0)
 tidyr              1.3.0   2023-01-24 [1] CRAN (R 4.1.0)
 tidyselect         1.2.0   2022-10-10 [1] CRAN (R 4.1.0)
 tweenr             2.0.2   2022-09-06 [1] CRAN (R 4.1.0)
 urlchecker         1.0.1   2021-11-30 [1] CRAN (R 4.1.0)
 usethis            2.2.2   2023-07-06 [1] CRAN (R 4.1.0)
 utf8               1.2.2   2021-07-24 [2] RSPM (R 4.1.0)
 uwot               0.1.16  2023-06-29 [1] CRAN (R 4.1.0)
 vctrs              0.6.4   2023-10-12 [1] CRAN (R 4.1.0)
 VGAM             * 1.1-9   2023-09-19 [1] CRAN (R 4.1.0)
 viridis            0.6.1   2021-05-11 [2] RSPM (R 4.1.0)
 viridisLite        0.4.0   2021-04-13 [2] RSPM (R 4.1.0)
 withr              2.5.1   2023-09-26 [1] CRAN (R 4.1.0)
 xtable             1.8-4   2019-04-21 [2] RSPM (R 4.1.0)
 zoo                1.8-9   2021-03-09 [2] RSPM (R 4.1.0)



───────────────────────────────────────────────────────────────────


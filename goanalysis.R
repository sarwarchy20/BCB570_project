
library(reactome.db)
require(viridis)
library(org.Hs.eg.db)
library(data.table)
require(DOSE)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(gplots)

caf.de <-  read.csv("CAF_all_DE.csv", row.names = 1)
dim(caf.de)

caf.de<- caf.de[order(caf.de$avg_log2FC,decreasing=TRUE),]
#View(caf.de)

# get gene symbols; human = 'Hs'
SYMBOL2EG <-
  eval(parse(text = sprintf(
    'org.%s.egSYMBOL2EG', 'Hs'
  )))


log2FC_score <- caf.de$avg_log2FC

de_genes<- rownames(caf.de)

names(log2FC_score) <- de_genes


genes <- intersect(de_genes, mappedkeys(SYMBOL2EG)) #  Get the gene symbol that are mapped to an entrez gene identifiers
length(x = genes) 


log2FC_score<- log2FC_score[genes] # access log2FC score 
length(x = log2FC_score) 

gene_entrez <-genes %>% SYMBOL2EG[.] %>% as.list %>% map( ~ .[1]) %>% simplify
length(x = gene_entrez)

names(log2FC_score) <- gene_entrez

log2FC_score.flt <- log2FC_score[abs(log2FC_score) > 1]
#log2FC_score.flt <- log2FC_score

length(log2FC_score.flt)

# === GO 
gse <- gseGO(geneList= log2FC_score.flt, 
             ont ="CC",
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             minGSSize = 10,
             maxGSSize = 500)

length(gse@result$p.adjust)

#gse@result$p.adjust <- (-(log10(gse@result$p.adjust)))

#jpeg("Path_CAF_GO_BP.jpeg", 
#     width = 700, height = 800)

#dotplot(gse, showCategory=20) + 
#  ggtitle("CAFs:GO-MF") + 
 # scale_color_gradient(low ="red",high ="blue")

#dev.off()

# using ggplot
# https://bioinformatics.ccr.cancer.gov/docs/btep-coding-club/FunctionalEnrich_clusterProfiler/
# ==================
jpeg("Path_CAF_GO_CC.jpeg", 
     width = 700, 
     height = 600)
gse %>% filter(p.adjust < 0.05) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  #geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradient(low ="red",high ="blue") +
  scale_size_continuous(range=c(3, 8)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("CAFs:GO-CC") +
  theme(text = element_text(size=15, face="bold"))

dev.off()

# ================== KEGG pathway

kk2 <- gseKEGG(geneList     = log2FC_score.flt,
               organism     = 'hsa',
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05,
               verbose      = T)

length(kk2@result$p.adjust)


#jpeg("Path_CAF_KEGG.jpeg", 
  #   width = 700, 
  #   height = 600)

#dotplot(kk2, showCategory=20, 
 #       font.size = 12) + 
 # ggtitle("CAFs:KEGG") +
 # scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  #theme_minimal()

#dev.off()

jpeg("Path_CAF_KEGG_updated.jpeg", 
     width = 700, 
     height = 600)
kk2 %>% filter(p.adjust < 0.05) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  #geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradient(low ="red",high ="blue") +
  scale_size_continuous(range=c(3, 8)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("CAFs:KEGG") +
  theme(text = element_text(size=15, face="bold"))

dev.off()

# devtools::session_info()%>% capture.output(file="session_info.txt")
# ===================================
─ Session info ─────────────────────────────────────────────
 setting  value
 version  R version 4.3.1 (2023-06-16 ucrt)
 os       Windows 10 x64 (build 19045)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  English_United States.utf8
 ctype    English_United States.utf8
 tz       America/Chicago
 date     2023-12-13
 rstudio  2022.07.1+554 Spotted Wakerobin (desktop)
 pandoc   NA

─ Packages ─────────────────────────────────────────────────
 package          * version   date (UTC) lib source
 abind              1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
 AnnotationDbi    * 1.62.2    2023-07-02 [1] Bioconductor
 ape                5.7-1     2023-03-13 [1] CRAN (R 4.3.2)
 aplot              0.2.2     2023-10-06 [1] CRAN (R 4.3.2)
 askpass            1.2.0     2023-09-03 [1] CRAN (R 4.3.1)
 Biobase          * 2.60.0    2023-04-25 [1] Bioconductor
 BiocGenerics     * 0.46.0    2023-04-25 [1] Bioconductor
 BiocManager        1.30.22   2023-08-08 [1] CRAN (R 4.3.2)
 BiocParallel       1.34.2    2023-05-22 [1] Bioconductor
 Biostrings         2.68.1    2023-05-16 [1] Bioconductor
 bit                4.0.5     2022-11-15 [1] CRAN (R 4.3.1)
 bit64              4.0.5     2020-08-30 [1] CRAN (R 4.3.1)
 bitops             1.0-7     2021-04-24 [1] CRAN (R 4.3.1)
 blob               1.2.4     2023-03-17 [1] CRAN (R 4.3.1)
 cachem             1.0.8     2023-05-01 [1] CRAN (R 4.3.1)
 callr              3.7.3     2022-11-02 [1] CRAN (R 4.3.1)
 caTools            1.18.2    2021-03-28 [1] CRAN (R 4.3.1)
 cli                3.6.1     2023-03-23 [1] CRAN (R 4.3.1)
 cluster          * 2.1.4     2022-08-22 [2] CRAN (R 4.3.1)
 clusterProfiler  * 4.8.3     2023-08-31 [1] Bioconductor
 codetools          0.2-19    2023-02-01 [2] CRAN (R 4.3.1)
 colorspace         2.1-0     2023-01-23 [1] CRAN (R 4.3.1)
 corpcor            1.6.10    2021-09-16 [1] CRAN (R 4.3.1)
 cowplot            1.1.1     2020-12-30 [1] CRAN (R 4.3.1)
 crayon             1.5.2     2022-09-29 [1] CRAN (R 4.3.1)
 data.table       * 1.14.8    2023-02-17 [1] CRAN (R 4.3.1)
 DBI                1.1.3     2022-06-18 [1] CRAN (R 4.3.1)
 deldir             1.0-9     2023-05-17 [1] CRAN (R 4.3.1)
 devtools           2.4.5     2022-10-11 [1] CRAN (R 4.3.1)
 digest             0.6.33    2023-07-07 [1] CRAN (R 4.3.1)
 doParallel         1.0.17    2022-02-07 [1] CRAN (R 4.3.1)
 DOSE             * 3.26.2    2023-10-09 [1] Bioconductor
 doSNOW             1.0.20    2022-02-04 [1] CRAN (R 4.3.2)
 dotCall64          1.1-0     2023-10-17 [1] CRAN (R 4.3.2)
 downloader         0.4       2015-07-09 [1] CRAN (R 4.3.2)
 dplyr            * 1.1.3     2023-09-03 [1] CRAN (R 4.3.1)
 ellipsis           0.3.2     2021-04-29 [1] CRAN (R 4.3.1)
 enrichplot       * 1.20.3    2023-09-15 [1] Bioconductor
 evaluate           0.23      2023-11-01 [1] CRAN (R 4.3.2)
 factoextra       * 1.0.7     2020-04-01 [1] CRAN (R 4.3.2)
 fansi              1.0.5     2023-10-08 [1] CRAN (R 4.3.1)
 farver             2.1.1     2022-07-06 [1] CRAN (R 4.3.1)
 fastDummies        1.7.3     2023-07-06 [1] CRAN (R 4.3.2)
 fastmap            1.1.1     2023-02-24 [1] CRAN (R 4.3.1)
 fastmatch          1.1-4     2023-08-18 [1] CRAN (R 4.3.1)
 fgsea            * 1.26.0    2023-04-25 [1] Bioconductor
 fitdistrplus       1.1-11    2023-04-25 [1] CRAN (R 4.3.2)
 forcats            1.0.0     2023-01-29 [1] CRAN (R 4.3.1)
 foreach            1.5.2     2022-02-02 [1] CRAN (R 4.3.1)
 fs                 1.6.3     2023-07-20 [1] CRAN (R 4.3.1)
 future             1.33.0    2023-07-01 [1] CRAN (R 4.3.2)
 future.apply       1.11.0    2023-05-21 [1] CRAN (R 4.3.2)
 generics           0.1.3     2022-07-05 [1] CRAN (R 4.3.1)
 GenomeInfoDb       1.36.4    2023-10-02 [1] Bioconductor
 GenomeInfoDbData   1.2.10    2023-10-21 [1] Bioconductor
 ggcorrplot       * 0.1.4.1   2023-09-05 [1] CRAN (R 4.3.2)
 ggforce            0.4.1     2022-10-04 [1] CRAN (R 4.3.2)
 ggfortify        * 0.4.16    2023-03-20 [1] CRAN (R 4.3.2)
 ggfun              0.1.3     2023-09-15 [1] CRAN (R 4.3.2)
 ggplot2          * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
 ggplotify          0.1.2     2023-08-09 [1] CRAN (R 4.3.1)
 ggraph             2.1.0     2022-10-09 [1] CRAN (R 4.3.2)
 ggrepel          * 0.9.4     2023-10-13 [1] CRAN (R 4.3.1)
 ggridges           0.5.4     2022-09-26 [1] CRAN (R 4.3.2)
 ggsci            * 3.0.0     2023-03-08 [1] CRAN (R 4.3.1)
 ggtree             3.8.2     2023-07-24 [1] Bioconductor
 globals            0.16.2    2022-11-21 [1] CRAN (R 4.3.1)
 glue               1.6.2     2022-02-24 [1] CRAN (R 4.3.1)
 GO.db              3.17.0    2023-12-10 [1] Bioconductor
 goftest            1.2-3     2021-10-07 [1] CRAN (R 4.3.1)
 GOSemSim           2.26.1    2023-07-10 [1] Bioconductor
 gplots           * 3.1.3     2022-04-25 [1] CRAN (R 4.3.1)
 graphlayouts       1.0.2     2023-11-03 [1] CRAN (R 4.3.2)
 gridBase           0.4-7     2014-02-24 [1] CRAN (R 4.3.2)
 gridExtra          2.3       2017-09-09 [1] CRAN (R 4.3.1)
 gridGraphics       0.5-1     2020-12-13 [1] CRAN (R 4.3.1)
 gson               0.1.0     2023-03-07 [1] CRAN (R 4.3.2)
 gtable             0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
 gtools             3.9.4     2022-11-27 [1] CRAN (R 4.3.1)
 HDO.db             0.99.1    2023-12-10 [1] Bioconductor
 highr              0.10      2022-12-22 [1] CRAN (R 4.3.1)
 htmltools          0.5.6.1   2023-10-06 [1] CRAN (R 4.3.1)
 htmlwidgets        1.6.2     2023-03-17 [1] CRAN (R 4.3.1)
 httpuv             1.6.11    2023-05-11 [1] CRAN (R 4.3.1)
 httr               1.4.7     2023-08-15 [1] CRAN (R 4.3.1)
 ica                1.0-3     2022-07-08 [1] CRAN (R 4.3.1)
 igraph             1.5.1     2023-08-10 [1] CRAN (R 4.3.1)
 IRanges          * 2.34.1    2023-06-22 [1] Bioconductor
 irlba              2.3.5.1   2022-10-03 [1] CRAN (R 4.3.1)
 iterators          1.0.14    2022-02-05 [1] CRAN (R 4.3.1)
 jsonlite           1.8.7     2023-06-29 [1] CRAN (R 4.3.1)
 KEGGREST           1.40.1    2023-09-29 [1] Bioconductor
 KernSmooth         2.23-22   2023-07-10 [1] CRAN (R 4.3.1)
 knitr              1.45      2023-10-30 [1] CRAN (R 4.3.2)
 labeling           0.4.3     2023-08-29 [1] CRAN (R 4.3.1)
 later              1.3.1     2023-05-02 [1] CRAN (R 4.3.1)
 lattice            0.21-9    2023-10-01 [1] CRAN (R 4.3.1)
 lazyeval           0.2.2     2019-03-15 [1] CRAN (R 4.3.1)
 leiden             0.4.3     2022-09-10 [1] CRAN (R 4.3.2)
 lifecycle          1.0.4     2023-11-07 [1] CRAN (R 4.3.2)
 listenv            0.9.0     2022-12-16 [1] CRAN (R 4.3.2)
 lmtest             0.9-40    2022-03-21 [1] CRAN (R 4.3.2)
 M3C              * 1.22.0    2023-04-25 [1] Bioconductor
 magrittr           2.0.3     2022-03-30 [1] CRAN (R 4.3.1)
 MASS               7.3-60    2023-05-04 [2] CRAN (R 4.3.1)
 Matrix             1.6-1.1   2023-09-18 [1] CRAN (R 4.3.1)
 matrixcalc         1.0-6     2022-09-14 [1] CRAN (R 4.3.1)
 matrixStats        1.0.0     2023-06-02 [1] CRAN (R 4.3.1)
 memoise            2.0.1     2021-11-26 [1] CRAN (R 4.3.1)
 mime               0.12      2021-09-28 [1] CRAN (R 4.3.0)
 miniUI             0.1.1.1   2018-05-18 [1] CRAN (R 4.3.1)
 munsell            0.5.0     2018-06-12 [1] CRAN (R 4.3.1)
 nlme               3.1-163   2023-08-09 [1] CRAN (R 4.3.1)
 NMF              * 0.26      2023-03-20 [1] CRAN (R 4.3.2)
 openssl            2.1.1     2023-09-25 [1] CRAN (R 4.3.1)
 org.Hs.eg.db     * 3.17.0    2023-12-10 [1] Bioconductor
 parallelly         1.36.0    2023-05-26 [1] CRAN (R 4.3.1)
 patchwork        * 1.1.3     2023-08-14 [1] CRAN (R 4.3.2)
 pbapply            1.7-2     2023-06-27 [1] CRAN (R 4.3.2)
 pillar             1.9.0     2023-03-22 [1] CRAN (R 4.3.1)
 pkgbuild           1.4.2     2023-06-26 [1] CRAN (R 4.3.1)
 pkgconfig          2.0.3     2019-09-22 [1] CRAN (R 4.3.1)
 pkgload            1.3.3     2023-09-22 [1] CRAN (R 4.3.1)
 plotly             4.10.3    2023-10-21 [1] CRAN (R 4.3.2)
 plyr               1.8.9     2023-10-02 [1] CRAN (R 4.3.2)
 png                0.1-8     2022-11-29 [1] CRAN (R 4.3.1)
 polyclip           1.10-6    2023-09-27 [1] CRAN (R 4.3.1)
 prettyunits        1.2.0     2023-09-24 [1] CRAN (R 4.3.1)
 processx           3.8.2     2023-06-30 [1] CRAN (R 4.3.1)
 profvis            0.3.8     2023-05-02 [1] CRAN (R 4.3.1)
 progressr          0.14.0    2023-08-10 [1] CRAN (R 4.3.2)
 promises           1.2.1     2023-08-10 [1] CRAN (R 4.3.1)
 ps                 1.7.5     2023-04-18 [1] CRAN (R 4.3.1)
 purrr            * 1.0.2     2023-08-10 [1] CRAN (R 4.3.1)
 qvalue             2.32.0    2023-04-25 [1] Bioconductor
 R6                 2.5.1     2021-08-19 [1] CRAN (R 4.3.1)
 RANN               2.6.1     2019-01-08 [1] CRAN (R 4.3.1)
 RColorBrewer     * 1.1-3     2022-04-03 [1] CRAN (R 4.3.1)
 Rcpp               1.0.11    2023-07-06 [1] CRAN (R 4.3.1)
 RcppAnnoy          0.0.21    2023-07-02 [1] CRAN (R 4.3.1)
 RcppHNSW           0.5.0     2023-09-19 [1] CRAN (R 4.3.2)
 RCurl              1.98-1.12 2023-03-27 [1] CRAN (R 4.3.1)
 reactome.db      * 1.84.0    2023-12-10 [1] Bioconductor
 registry         * 0.5-1     2019-03-05 [1] CRAN (R 4.3.1)
 remotes            2.4.2.1   2023-07-18 [1] CRAN (R 4.3.1)
 reshape2           1.4.4     2020-04-09 [1] CRAN (R 4.3.2)
 reticulate         1.34.0    2023-10-12 [1] CRAN (R 4.3.2)
 rlang              1.1.1     2023-04-28 [1] CRAN (R 4.3.1)
 rngtools         * 1.5.2     2021-09-20 [1] CRAN (R 4.3.2)
 ROCR               1.0-11    2020-05-02 [1] CRAN (R 4.3.1)
 RSpectra           0.16-1    2022-04-24 [1] CRAN (R 4.3.2)
 RSQLite            2.3.3     2023-11-04 [1] CRAN (R 4.3.2)
 rstudioapi         0.15.0    2023-07-07 [1] CRAN (R 4.3.1)
 Rtsne            * 0.16      2022-04-17 [1] CRAN (R 4.3.2)
 S4Vectors        * 0.38.2    2023-09-22 [1] Bioconductor
 scales             1.2.1     2022-08-20 [1] CRAN (R 4.3.1)
 scattermore        1.2       2023-06-12 [1] CRAN (R 4.3.2)
 scatterpie         0.2.1     2023-06-07 [1] CRAN (R 4.3.2)
 sctransform        0.4.1     2023-10-19 [1] CRAN (R 4.3.2)
 sessioninfo        1.2.2     2021-12-06 [1] CRAN (R 4.3.1)
 Seurat           * 5.0.0     2023-11-04 [1] CRAN (R 4.3.2)
 SeuratObject     * 5.0.0     2023-10-26 [1] CRAN (R 4.3.2)
 shadowtext         0.1.2     2022-04-22 [1] CRAN (R 4.3.2)
 shiny              1.7.5.1   2023-10-14 [1] CRAN (R 4.3.1)
 snow               0.4-4     2021-10-27 [1] CRAN (R 4.3.1)
 sp               * 2.1-1     2023-10-16 [1] CRAN (R 4.3.1)
 spam               2.10-0    2023-10-23 [1] CRAN (R 4.3.2)
 spatstat.data      3.0-3     2023-10-24 [1] CRAN (R 4.3.2)
 spatstat.explore   3.2-5     2023-10-22 [1] CRAN (R 4.3.2)
 spatstat.geom      3.2-7     2023-10-20 [1] CRAN (R 4.3.2)
 spatstat.random    3.2-1     2023-10-21 [1] CRAN (R 4.3.2)
 spatstat.sparse    3.0-3     2023-10-24 [1] CRAN (R 4.3.2)
 spatstat.utils     3.0-4     2023-10-24 [1] CRAN (R 4.3.2)
 stringi            1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
 stringr            1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
 survival           3.5-7     2023-08-14 [1] CRAN (R 4.3.1)
 tensor             1.5       2012-05-05 [1] CRAN (R 4.3.1)
 tibble             3.2.1     2023-03-20 [1] CRAN (R 4.3.1)
 tidygraph          1.2.3     2023-02-01 [1] CRAN (R 4.3.2)
 tidyr              1.3.0     2023-01-24 [1] CRAN (R 4.3.1)
 tidyselect         1.2.0     2022-10-10 [1] CRAN (R 4.3.1)
 tidytree           0.4.5     2023-08-10 [1] CRAN (R 4.3.2)
 treeio             1.24.3    2023-07-24 [1] Bioconductor
 tweenr             2.0.2     2022-09-06 [1] CRAN (R 4.3.2)
 umap               0.2.10.0  2023-02-01 [1] CRAN (R 4.3.2)
 urlchecker         1.0.1     2021-11-30 [1] CRAN (R 4.3.1)
 usethis            2.2.2     2023-07-06 [1] CRAN (R 4.3.1)
 utf8               1.2.3     2023-01-31 [1] CRAN (R 4.3.1)
 uwot               0.1.16    2023-06-29 [1] CRAN (R 4.3.1)
 vctrs              0.6.4     2023-10-12 [1] CRAN (R 4.3.1)
 viridis          * 0.6.4     2023-07-22 [1] CRAN (R 4.3.2)
 viridisLite      * 0.4.2     2023-05-02 [1] CRAN (R 4.3.1)
 withr              2.5.2     2023-10-30 [1] CRAN (R 4.3.2)
 xfun               0.41      2023-11-01 [1] CRAN (R 4.3.2)
 xtable             1.8-4     2019-04-21 [1] CRAN (R 4.3.1)
 XVector            0.40.0    2023-04-25 [1] Bioconductor
 yulab.utils        0.1.0     2023-09-20 [1] CRAN (R 4.3.1)
 zlibbioc           1.46.0    2023-04-25 [1] Bioconductor
 zoo                1.8-12    2023-04-13 [1] CRAN (R 4.3.2)

────────────────────────────────────────────────────────────


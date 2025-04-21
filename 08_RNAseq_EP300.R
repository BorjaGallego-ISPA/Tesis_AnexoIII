################################################################################
######## Transcriptomic Characterization of EP300 Inhibitor Combine with #######
############  DOX Treatment in a DX-R Bone Sarcoma Model: Analysis #############
################################################################################

# Author: Borja Gallego Martínez -----------------------------------------------
  # Sarcomas and Experimental Therapeutics lab (ISPA, Oviedo, Spain)

# Date: 30/10/2024 ----

# Info -------------------------------------------------------------------------

  # This script performs all the steps for the downstream analysis of the RNAseq
    # of 143B DX-R treated with a combination of DOX 1 µM and SGC-CBP30 5 µM, an
    # EP300 inhibitor, (DXSG samples, triplicate) vs. 143B DX-R treated with DOX
    # 1 µM (DOXO samples, triplicate):

      # Proccessing the raw data (from Salmon quant files to DGE limma object).

      # Quality control of the samples (Correlation heatmap and PCA).

      # Differential Expression Analysis (DEGs table, Volcano, Boxplots...).

      # Functional Analysis (GSEA and TF activity inference).

  # IMPORTANT: This script is only to perform the analysis, for more information
    # about the background of this project, information about the samples and
    # visualization of the results please go to the report file (02_Report.Rmd).

################################################################################

# 0. Load Packages and Functions ----------------------------------------------

library(tidyverse)
library(readxl)
library(BiocGenerics)
library(SummarizedExperiment)
library(tximeta)
library(GenomicFeatures)
library(limma)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggrepel)
library(ggVennDiagram)
library(VennDiagram)
library(clusterProfiler)
library(ReactomePA)
library(msigdbr)
library(enrichplot)
library(OmnipathR)
library(dorothea)
library(decoupleR)


select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter

source("./Downstream/00_Functions.R")

set.seed(123)

################################################################################

# 1. Processing Raw Data ------------------------------------------------------

## 1.1. Sample Information ----

Metadata <- data.frame(
  files = list.files("./quants", pattern = "quant.sf", recursive = T,
                     full.names = T),
  names = list.files("./quants/"),
  cell_line = rep("b143", times = length(list.files("./quants/"))),
  condition = rep(c("DOXO", "DXSG"), each = 3),
  replicate = rep(c(1:3), times = 2),
  dox_uM = rep(1, times = 6),
  sgc_uM = rep(c(0, 5), each = 3),
  row.names = list.files("./quants/")
)

## 1.2. Import Transcript Abundance ----

se <- tximeta(Metadata)
gse <- summarizeToGene(se, countsFromAbundance = "lengthScaledTPM")

## 1.3. Create DGE Object and Pre-filter ----

genes <- assays(gse)[["counts"]] %>%
  as.data.frame() %>%
  mutate(
    Symbol = mapIds(org.Hs.eg.db, keys = substr(rownames(.), 1, 15),
                    column = "SYMBOL", keytype = "ENSEMBL",
                    multiVals = "first"),
    Biotype = mapIds(org.Hs.eg.db, keys = substr(rownames(.), 1, 15),
                     column = "GENETYPE", keytype = "ENSEMBL",
                     multiVals = "first"),
    Description = mapIds(org.Hs.eg.db, keys = substr(rownames(.), 1, 15),
                         column = "GENENAME", keytype = "ENSEMBL",
                         multiVals = "first")
  ) %>%
  select(Symbol, Biotype, Description)

y <- DGEList(assays(gse)[["counts"]], samples = Metadata, genes = genes)

design <- model.matrix(~ 0 + condition, data = Metadata)
colnames(design) <- gsub("condition", "", colnames(design))

keep <- filterByExpr(y, design)
y <- y[keep, ]
y <- calcNormFactors(y)

v <- voom(y, design = design)

# 2. Quality Control -----------------------------------------------------------

## 2.1. Correlation Heatmap ----

QC.CorrHeatmap <- QC.MakeCorrHeatmap(
  Counts.norm = v$E,
  heatmap.rows = colnames(v$E),
  main = "Correlation Heatmap",
  output.path = "./Downstream/Results/Figures/QC_CorrHeatmap.tiff"
)

## 2.2. Principal Component Analysis ----

QC.PCA <- QC.MakePCA(
  data = v$E,
  metadata = Metadata,
  main = "PCA",
  output.path = "./Downstream/Results/Figures/QC_PCA.tiff"
)

# 3. Differential Expression Analysis ------------------------------------------

contr.matrix <- makeContrasts(b143.cont = DXSG - DOXO,
                              levels = colnames(design))

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

## 3.1. DEGs Table ----

DEA.Table <- topTable(efit, number = Inf, coef = 1) %>%
  filter(!is.na(Symbol)) %>%
  arrange(adj.P.Val) %>%
  distinct(Symbol, .keep_all = T) %>%
  arrange(Symbol)

write.csv(DEA.Table, "./Downstream/Results/DEGs.csv")
xlsx::write.xlsx2(DEA.Table, "./Downstream/Results/DEGs.xlsx")

## 3.2. Exploratory Plots ----
### 3.2.1. Volcano Plot ----

DEA.Volcano <- DEA.MakeVolcano(
  DEA.table = DEA.Table,
  sym.coln = "Symbol", lfc.coln = "logFC", fdr.coln = "adj.P.Val",
  main = "143B DX-R: DOX+SGC vs. DOX",
  output.path = "./Downstream/Results/Figures/Volcano.tiff"
)

### 3.2.2. Counts Distribution of Top Genes ----

DEA.Boxplot <- DEA.MakeBoxPlot(
  DEA.table = DEA.Table,
  Counts.table = as.data.frame(v$E),
  sorting = "both", sym.coln = "Symbol", lfc.coln = "logFC",
  fdr.coln = "adj.P.Val",
  conditions = c(rep("DOXO", times = 30), rep("DXSG", times = 30)),
  output.file_up = "./Downstream/Results/Figures/Boxplot_up.tiff",
  output.file_down = "./Downstream/Results/Figures/Boxplot_down.tiff"
)

# 4. Functional and Footprint Analysis -----------------------------------------

## 4.1. Gene Set Enrichment Analysis (GSEA) ----

Msigdb.geneset_h <- msigdbr(species = "human", category = "H") %>%
  select(gs_name, entrez_gene)

GSEA.dat <- PathEnrich.MakeGSEA(
  GeneSet = Msigdb.geneset_h,
  Res = DEA.Table, lfc.coln = "logFC",
  table.path = "./Downstream/Results/GSEATable",
)

GSEA.plot <- PathEnrich.MakeBubble(
  GSEA.dat$Results.table,
  main = "GSEA for 143B DX-R: DOX+SGC vs. DOX",
  output.path = "./Downstream/Results/Figures/GSEA_bubble.tiff",
  x.expand = c(.1, 0)
)

RPA.dat <- PathEnrich.MakeGSEA(
  Res = DEA.Table, lfc.coln = "logFC", db = "RPA",
  table.path = "./Downstream/Results/RPATable",
)

RPA.plot <- PathEnrich.MakeBubble(
  RPA.dat$Results.table,
  main = "GSEA-RPA",
  output.path = "./Downstream/Results/Figures/RPA_bubble.tiff",
  x.expand = c(.1, 0)
)

## 4.2. Transcription Factor Activity Inference ----

TFNet <- ChooseTFNet(db = "tri")

TFs <- Get_TFs(
  mat = DEA.Table, condition = "t",
  net = TFNet, method = "ULM", pval.thr = .05, ax.tx.x = 4,
  main = "Significant TFs Regulated in 143B DX-R: DOX+SGC vs. DOX",
  output.path  = "./Downstream/Results/Figures/TFs.tiff"
)

TFs_25 <- Get_TFs(
  mat = DEA.Table, condition = "t",
  net = TFNet, method = "ULM", pval.thr = .05, n_tfs = 25,
  main = "Top 25 TFs Regulated in 143B DX-R: DOX+SGC vs. DOX",
  output.path  = "./Downstream/Results/Figures/TFs_top25.tiff"
)

################################################################################

# 99. Save Results -------------------------------------------------------------

sessInfo <- sessionInfo()

save(
  Metadata, QC.CorrHeatmap, QC.PCA, DEA.Table, DEA.Volcano, DEA.Boxplot,
  GSEA.plot, TFs, TFs_25, sessInfo,
  file = paste0(
    "./Downstream/Results/",
    paste(unlist(strsplit(as.character(Sys.Date()), split = "-")),
          collapse = ""),
    "_Results_RNAseq_EP300.RData")
)

save.image(
  file = paste0(
    "./Downstream/",
    paste(unlist(strsplit(as.character(Sys.Date()), split = "-")),
          collapse = ""),
    "_Workspace_RNAseq_EP300.RData"
  )
)

################################################################################

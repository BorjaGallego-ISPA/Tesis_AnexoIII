################################################################################
########### Transcriptomic Characterization of Doxorubicin Resistance ##########
####################### in Osteosarcoma Organoids Models #######################
################################################################################

# Info #####

# Author: Borja Gallego Martínez
# Sarcomas and Experimental Theraputics Lab (ISPA, Oviedo, Spain)
# Last Modification: 24 Jul 2024

# Setting-Up ####

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
library(msigdbr)
library(enrichplot)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(OmnipathR)
library(dorothea)
library(decoupleR)
library(igraph)
library(STRINGdb)
library(doMC)
library(LEANR)
library(visNetwork)
library(ggraph)
library(tidygraph)
library(Cairo)

select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter
arrange <- dplyr::arrange

source("./Support/HomemadeFunctions_rnaseq.R")

set.seed(123)

# Proccessing Raw Data ####

## Samples Information ####

Metadata <- data.frame(
  files = list.files("./Raw", pattern = "quant.sf", recursive = T, full.names = T),
  names = list.files("./Raw/"),
  cell_line = rep("b143", times = length(list.files("./Raw/"))),
  condition = c(rep("Par", times = 3), rep("DXR", times = 3)),
  replicate = rep(c(1, 2, 3), times = 2),
  row.names = list.files("./Raw/")
) %>%
  mutate(
    group = condition
  )

## Import Transcript Abundance ####

se <- tximeta(Metadata)
gse <- summarizeToGene(se, countsFromAbundance = "lengthScaledTPM")

## Create DGE Object and Prefilter ####

genes <- assays(gse)[["counts"]] %>%
  as.data.frame() %>%
  mutate(
    Symbol = mapIds(org.Hs.eg.db, keys = substr(rownames(.), 1, 15),
                    column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
    Biotype = mapIds(org.Hs.eg.db, keys = substr(rownames(.), 1, 15),
                     column = "GENETYPE", keytype = "ENSEMBL", multiVals = "first"),
    Description = mapIds(org.Hs.eg.db, keys = substr(rownames(.), 1, 15),
                         column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")
  ) %>%
  select(Symbol, Biotype, Description)

y <- DGEList(assays(gse)[["counts"]], samples = Metadata, genes = genes)

keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y)

design <- model.matrix(~ 0 + condition, data = Metadata)
colnames(design) <- gsub("condition", "", colnames(design))

v <- voom(y, design = design)

# Quality Control ####

## Correlation Heatmap ####

QC.CorrHeatmap <- QC.MakeCorrHeatmap(
  Counts.norm = v$E,
  heatmap.rows = colnames(v$E),
  main = "Correlation Heatmap",
  output.path = "./Results/QC_CorrHeatmap.tiff"
)

## Principal Component Analysis ####

QC.PCA <- QC.MakePCA(
  data = v$E,
  metadata = Metadata,
  main = "PCA",
  output.path = "./Results/QC_PCA.tiff",
  fill.color = c("#FF0000", "#0000FF")
)

# Differential Expression Analysis ####

contr.matrix <- makeContrasts(b143.cont = DXR - Par, levels = colnames(design))
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

## DEA Table ####

DEA.Table <- topTable(efit, number = Inf, coef = 1) %>%
  filter(!is.na(Symbol)) %>%
  arrange(adj.P.Val) %>%
  distinct(Symbol, .keep_all = T) %>%
  arrange(Symbol)

write.csv(DEA.Table, "./Results/DEA_Table.csv")
xlsx::write.xlsx2(DEA.Table, "./Results/DEA_Table.xlsx")

## Volcano Plot ####

DEA.Volcano <- DEA.MakeVolcano(
  DEA.table = DEA.Table,
  sym.coln = "Symbol", lfc.coln = "logFC", fdr.coln = "adj.P.Val",
  main = "143B Organoids DX-R vs. Parental",
  output.path = "./Results/DEA_volcano.tiff"
)
DEA.Volcano_big <- DEA.MakeVolcano(
  DEA.table = DEA.Table,
  sym.coln = "Symbol", lfc.coln = "logFC", fdr.coln = "adj.P.Val",
  main = "",
  ax.title.size = 45, ax.tx.size = 30, lg.tx.size = 25, grepel.size = 12,
  pt.size = 4, save.w = 12, save.h = 12, lg.pos = c(.82,.18),
  output.path = "./Results/DEA_volcano_big.tiff"
)

## Counts Distribution of Top Genes ####

DEA.Boxplot <- DEA.MakeBoxPlot(
  DEA.table = filter(DEA.Table, Biotype == "protein-coding"),
  Counts.table = as.data.frame(v$E),
  sorting = "lfc", sym.coln = "Symbol", lfc.coln = "logFC", fdr.coln = "adj.P.Val",
  conditions = c(rep("Par", times = 30), rep("DXR", times = 30)),
  output.file_up = "./Results/DEA_boxplot_up.tiff",
  output.file_down = "./Results/DEA_boxplot_down.tiff"
)

DEA.Boxplot_sig <- DEA.MakeBoxPlot(
  DEA.table = filter(DEA.Table, Biotype == "protein-coding"),
  Counts.table = as.data.frame(v$E),
  sorting = "both", sym.coln = "Symbol", lfc.coln = "logFC", fdr.coln = "adj.P.Val",
  conditions = c(rep("Par", times = 30), rep("DXR", times = 30)),
  output.file_up = "./Results/DEA_boxplot_sig_up.tiff",
  output.file_down = "./Results/DEA_boxplot_sig_down.tiff"
)

# Functional and Footprint Analysis ####

## Gene Set Enrichment Analysis (GSEA) with MSigDB Hallmark Collection ####

Msigdb.geneset_h <- msigdbr(species = "human", category = "H") %>%
  select(gs_name, entrez_gene)

GSEA.dat <- PathEnrich.MakeGSEA(
  GeneSet = Msigdb.geneset_h,
  Res = DEA.Table, lfc.coln = "logFC",
  table.csv = "./Results/GSEATable.csv",
  table.xlsx = "./Results/GSEATable.xlsx"
)

GSEA.plot <- PathEnrich.MakeBubble(
  GSEA.dat$Results.table,
  main = "GSEA for 143B DX-R Organoids",
  output.path = "./Results/GSEA_bubble.tiff"
)
GSEA.plot_mod <- PathEnrich.MakeBubble(
  GSEA.dat$Results.table,
  main = "",
  output.path = "./Results/GSEA_bubble_mod.tiff",
  save.h = 3450, save.w = 2300,
  x.expand = c(.4, 0),
  lg.position = "none"
)

ggsave("./Results/GSEplot_G2M.tiff", plot = GSEA.dat$GSEA.plots[[14]],
       device = "tiff", dpi = 300)

column_ha <- HeatmapAnnotation(
  df = data.frame(
    Condition = Metadata$condition,
    row.names = Metadata$names
  ),
  show_annotation_name = F,
  col = list(Condition = c(DXR = "red", Par = "blue")),
  annotation_legend_param = list(
    Condition = list(
      title = "",
      labels_gp = gpar(fontsize = 15, fontface = "bold")
    )
  )
)

G2M.genelist <- GSEA.dat$Results.table %>%
  filter(Description == "G2M CHECKPOINT") %>%
  pull(core_enrichment) %>%
  strsplit(split = "/") %>%
  unlist()
G2M.genelist <- DEA.Table %>%
  filter(Symbol %in% G2M.genelist) %>%
  select(Symbol) %>%
  arrange(rownames(.))
G2M.counts <- v$E %>%
  as.data.frame() %>%
  filter(rownames(.) %in% rownames(G2M.genelist)) %>%
  cbind(G2M.genelist) %>%
  remove_rownames() %>%
  column_to_rownames(var = "Symbol")
G2M.row_ha <- DEA.Table %>%
  filter(Symbol %in% G2M.genelist$Symbol) %>%
  filter(logFC > 1 & adj.P.Val < .05) %>%
  pull(Symbol)
G2M.heatmap <- Heatmap(
  t(scale(t(G2M.counts))),
  top_annotation = column_ha,
  heatmap_legend_param = list(title = "", labels_gp = gpar(fontsize = 12)),
  show_column_names = F,
  show_row_names = F,
  right_annotation = rowAnnotation(
    siggenes = anno_mark(
      at = unlist(lapply(G2M.row_ha, grep, rownames(G2M.counts))),
      labels = G2M.row_ha,
      labels_gp = gpar(fontsize = 15),
#      extend = unit(10, "mm"),
      padding = 1
    )
  ),
  cluster_rows = T,
  show_row_dend = F,
  cluster_columns = F,
  column_title = "G2M Checkpoint",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  column_title_side = "top",
  column_split = Metadata$condition,
  column_gap = unit(1.5, "mm")
)

tiff("./Results/G2M_heatmap.tiff", width = 1250, height = 2500, units = "px",
     res = 300)
draw(G2M.heatmap)
dev.off()

## Transcription Factor Activity Inference ####

TFNet <- ChooseTFNet(db = "tri")

TFs <- Get_TFs(
  mat = DEA.Table, condition = "t",
  net = TFNet, method = "ULM", pval.thr = .05,
  main = "Significant TFs Regulated in 143B DX-R Organoids",
  output.path  = "./Results/TFs.tiff"
)

TFs_25 <- Get_TFs(
  mat = DEA.Table, condition = "t",
  net = TFNet, method = "ULM", pval.thr = .05, n_tfs = 25,
  main = "Top 25 TFs Regulated in 143B DX-R Organoids",
  output.path  = "./Results/TFs_25.tiff"
)

## Local Enrichment ANalysis (LEAN) ####

string.db <- STRINGdb$new(version = "11.5", species = 9606,
                          network_type = "full",
                          score_threshold = 800,
                          input_directory = "./Support")

lean.metadata <- DEA.Table %>%
  string.db$map("Symbol", removeUnmappedRows = T) %>%
  arrange(adj.P.Val) %>%
  distinct(Symbol, .keep_all = T) %>%
  distinct(STRING_id, .keep_all = T) %>%
  select(Symbol, STRING_id, logFC, FDR = adj.P.Val) %>%
  mutate(
    Sign = sign(logFC),
    Origin = c("Protein")
  )

lean.input <- list(values = lean.metadata %>% pull(FDR, name = STRING_id),
                   net = string.db$get_subnetwork(lean.metadata$STRING_id))

lean.res <- run.lean(
  ranking = lean.input$values, network = lean.input$net,
  n_reps = 10000, ncores = 15, ranked = F
)
lean.res_sig <- lean.res$restab %>%
  as.data.frame() %>%
  rownames_to_column(var = "STRING_id") %>%
  left_join(lean.metadata, by = "STRING_id") %>%
  filter(PLEAN <= .05) %>%
  column_to_rownames(var = "STRING_id") %>%
  arrange(PLEAN)

lean.clusters <- LEAN.MakeClusters(lean.res_sig, Msigdb.geneset_h)
lean.clusters$Network <- visNetwork(
  lean.clusters$Nodes, lean.clusters$Edges,
  main = "Clusterization and Enrichment of Significant LEAN Results",
  submain = "143B DX-R Organoids",
  width = "100%"
) %>%
  visNetwork::visIgraphLayout(randomSeed = 1) %>%
  visNodes(size = 30) %>%
  visLegend(zoom = F) %>%
  visOptions(
    highlightNearest = list(enabled = T, degree = 2, hover = F),
    nodesIdSelection = T
)
visSave(lean.clusters$Network,
        file = "./Results/LEAN_clusters.html")
webshot2::webshot(
  "./Results/LEAN_clusters.html",
  delay = .2, zoom = 1, vwidth = 992, vheight = 744, quiet = T,
  file = "./Results/LEAN_clusters.png"
)

lean.BMP1 <- LEAN.MakeSubNet(
  STRING_id = "9606.ENSP00000305714",
  lean.res = lean.res,
  lean.metadata = lean.metadata,
  main = "BMP1 LEAN Subnetwork", submain = "143B DX-R Organoids"
)

lean.SRC <- LEAN.MakeSubNet(
  STRING_id = "9606.ENSP00000362680",
  lean.res = lean.res,
  lean.metadata = lean.metadata,
  main = "SRC LEAN Subnetwork", submain = "143B DX-R Organoids",
  sel = F
)
lean.SRC.cluster <- LEAN.MakeClusters(select(lean.SRC$Nodes, Symbol = id),
                                      Msigdb.geneset_h)
lean.SRC.cluster$Network <- visNetwork(
  lean.SRC.cluster$Nodes, lean.SRC.cluster$Edges,
  main = "Clusterization and Enrichment of SRC LEAN Subnetwork",
  submain = "143B DX-R Organoids",
  width = "100%"
) %>%
  visNetwork::visIgraphLayout(randomSeed = 1) %>%
  visNodes(size = 30) %>%
  visLegend(zoom = F) %>%
  visOptions(
    highlightNearest = list(enabled = T, degree = 2, hover = F),
    nodesIdSelection = T
  )
visSave(lean.SRC.cluster$Network,
        file = "./Results/LEAN_SRC_clusters.html")
webshot2::webshot(
  "./Results/LEAN_SRC_clusters.html",
  delay = .2, zoom = 1, vwidth = 992, vheight = 744, quiet = T,
  file = "./Results/LEAN_SRC_clusters.png"
)

lean.SRC_sig <- LEAN.MakeSubNet(
  STRING_id = "9606.ENSP00000362680",
  lean.res = lean.res,
  lean.metadata = lean.metadata,
  main = "SRC LEAN Subnetwork", submain = "143B DX-R Organoids",
  sel = T
)
lean.SRC.cluster_sig <- LEAN.MakeClusters(select(lean.SRC_sig$Nodes,
                                                 Symbol = id),
                                          Msigdb.geneset_h)

lean.EGFR <- LEAN.MakeSubNet(
  STRING_id = "9606.ENSP00000275493",
  lean.res = lean.res,
  lean.metadata = lean.metadata,
  main = "EGFR LEAN Subnetwork", submain = "143B DX-R Organoids",
  sel = F
)
lean.EGFR.cluster <- LEAN.MakeClusters(select(lean.EGFR$Nodes, Symbol = id),
                                      Msigdb.geneset_h)
lean.EGFR.cluster$Network <- visNetwork(
  lean.EGFR.cluster$Nodes, lean.EGFR.cluster$Edges,
  main = "Clusterization and Enrichment of EGFR LEAN Subnetwork",
  submain = "143B DX-R Organoids",
  width = "100%"
) %>%
  visNetwork::visIgraphLayout(randomSeed = 1) %>%
  visNodes(size = 30) %>%
  visLegend(zoom = F) %>%
  visOptions(
    highlightNearest = list(enabled = T, degree = 2, hover = F),
    nodesIdSelection = T
  )

lean.EGFR_sel <- LEAN.MakeSubNet(
  STRING_id = "9606.ENSP00000275493",
  lean.res = lean.res,
  lean.metadata = lean.metadata,
  main = "EGFR LEAN Subnetwork", submain = "143B DX-R Organoids",
  sel = T
)

nodes.top <- 20
lean.topnodes <- lean.res_sig %>%
  arrange(desc(k)) %>%
  head(nodes.top) %>%
  ggplot(aes(k, fct_reorder(Symbol, k), fill = PLEAN)) +
  geom_col() +
  labs(title = "Top20 Dense LEAN Networks",
       x = "Nº Significant Nodes",
       y = "Center Node",
       fill = "LEAN\np.value") +
  scale_fill_gradient(low = "darkred", high = "whitesmoke", trans = "log10") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(lean.res_sig$k) + 1),
                     breaks = seq(10, 110, 10)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = .5, face = "bold", size = 22),
    axis.title = element_text(face = "bold", size = 18),
    axis.text.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 14)
  )
ggsave("./Results/LEAN_topnodes.tiff", plot = lean.topnodes, device = "tiff",
       dpi = 300, units = "in", width = 8, height = 8)

### Visualization with ggraph ####

lean.net <- tbl_graph(
  nodes = lean.clusters$Nodes %>%
    mutate(
      group = ifelse(group == "EPITHELIAL MESENCHYMAL TRANSITION", "EMT",
                     group),
      group = ifelse(group == "Cluster 5", "NONE", group),
      logFC = case_when(
        logFC > 1 ~ 1,
        logFC < -1 ~ -1,
        .default = logFC
      )
    ) %>%
    select(Nodes = id, LEAN.pval = PLEAN, logFC, cluster = group),
  edges = lean.clusters$Edges
) %>%
  ggraph(layout = "stress", niter = 500000) +
  geom_edge_fan(edge_alpha = .8, color = "grey50",
                arrow = arrow(length = unit(2, "mm"), type = "open"),
                end_cap = circle(7, "mm")) +
  geom_node_point(aes(fill = logFC, shape = cluster, colour = LEAN.pval),
                  size = 10, stroke = 1) +
  scale_shape_manual(values = c("GLYCOLYSIS" = 24, "EMT" = 23,
                                "TNFA SIGNALING VIA NFKB" = 21,
                                "MYC TARGETS V1" = 22, "NONE" = 25),
                     name = "Cluster") +
  scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4",
                       mid = "whitesmoke", midpoint = 0) +
  scale_color_viridis(name = "LEAN p value") +
  geom_node_text(aes(label = Nodes), size = 2, vjust = .5, hjust = .5,
                 color = "black", fontface = "bold") +
#  ggtitle("LEAN Significant Center Nodes") +
  theme_graph() +
  theme(
    legend.position = c(.2, .6),
#    plot.title = element_text(hjust = .5, face = "bold", size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
ggsave("./Results/LEAN_ggraphnet.tiff", lean.net, device = "tiff", dpi = 300,
       units = "px", width = 6000, height = 2500)

CairoPDF(file = "./Results/LEAN_ggraphnet.pdf", width = 20, height = 9)
print(lean.net)
dev.off()

SRC.net_sig <- tbl_graph(
  nodes = lean.SRC_sig$Nodes %>%
    left_join(lean.SRC.cluster_sig$Nodes, by = "id") %>%
    select(Nodes = id, FDR, logFC, Cluster = group) %>%
    mutate(
      Cluster = ifelse(Cluster == "EPITHELIAL MESENCHYMAL TRANSITION", "EMT",
                       Cluster),
      logFC = case_when(
        logFC > 1 ~ 1,
        logFC < -1 ~ -1,
        .default = logFC
      ),
      FDR = ifelse(FDR > .05, NA, FDR)
    ),
  edges = lean.SRC_sig$Edges
) %>%
  ggraph(layout = "stress", niter = 500000) +
  geom_edge_fan(edge_alpha = .8, color = "grey50",
                arrow = arrow(length = unit(2, "mm"), type = "open"),
                end_cap = circle(7, "mm")) +
  geom_node_point(aes(fill = logFC, colour = FDR), size = 12, shape = 21, stroke = 2) +
  scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4",
                       mid = "whitesmoke", midpoint = 0) +
  scale_color_viridis(na.value = "#FDE725") +
  geom_node_text(aes(label = Nodes), size = 2.5, vjust = .5, hjust = .5,
                 color = "black", fontface = "bold") +
  ggtitle("SRC LEAN Subnetwork (Significant Nodes)") +
  theme_graph() +
  theme(
    plot.title = element_text(hjust = .5, face = "bold", size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
ggsave("./Results/LEAN_SRC_sig.tiff", SRC.net_sig, device = "tiff", dpi = 300,
       units = "px", width = 3500, height = 3000)
CairoPDF(file = "./Results/LEAN_SRC_sig.pdf", width = 12, height = 10)
print(SRC.net_sig)
dev.off()

SRC.net_all <- tbl_graph(
  nodes = lean.SRC$Nodes %>%
    left_join(lean.SRC.cluster$Nodes, by = "id") %>%
    select(Nodes = id, FDR, logFC, Cluster = group) %>%
    mutate(
      Cluster = ifelse(Cluster == "EPITHELIAL MESENCHYMAL TRANSITION", "EMT",
                       Cluster),
      logFC = case_when(
        logFC > 1 ~ 1,
        logFC < -1 ~ -1,
        .default = logFC
      )
    ),
  edges = lean.SRC$Edges
) %>%
  ggraph(layout = "stress", niter = 500000) +
  geom_edge_fan(edge_alpha = .8, color = "grey50",
                arrow = arrow(length = unit(2, "mm"), type = "open"),
                end_cap = circle(7, "mm")) +
  geom_node_point(aes(fill = logFC), size = 10, shape = 21) +
  scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4",
                       mid = "whitesmoke", midpoint = 0) +
  geom_node_text(aes(label = Nodes), size = 2.5, vjust = .5, hjust = .5,
                 color = "black", fontface = "bold") +
  ggtitle("SRC LEAN Subnetwork (All Nodes)") +
  theme_graph() +
  theme(
    plot.title = element_text(hjust = .5, face = "bold", size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )
ggsave("./Results/LEAN_SRC_all.tiff", SRC.net_all, device = "tiff", dpi = 300,
       units = "px", width = 3500, height = 3000)

# Comparison Between 2D and 3D ####

Adh <- read.csv("./Support/b143_DEGsTable.csv", row.names = 1)

OrgvsAdh <- DEA.MakeVenn(
  data.l = list(Organoid = DEA.Table, Adherent = Adh),
  lfc.coln = "logFC", fdr.coln = "adj.P.Val",
  up.main = "Up-Regulated Coincidences",
  up.file = "./Results/OrgvsAdh_up.tiff",
  down.main = "Down-Regulated Coincidences",
  down.file = "./Results/OrgvsAdh_down.tiff"
)

Overlaps <- data.frame(
  Up = OrgvsAdh$Up.data$a3,
  Down = c(
    OrgvsAdh$Down.data$a3,
    rep("", times = length(OrgvsAdh$Up.data$a3) - length(OrgvsAdh$Down.data$a3))
  )
)
xlsx::write.xlsx2(Overlaps, "./Results/OrgvsAdh_overlaps.xlsx")

# Save Results ####

save(
  QC.CorrHeatmap, QC.PCA, DEA.Table, DEA.Volcano, DEA.Boxplot, DEA.Boxplot_sig,
  GSEA.plot, TFs, TFs_25, 
  TFs_25, OrgvsAdh,
  file = paste0("./Results/",
                paste(unlist(strsplit(as.character(Sys.Date()), split = "-")), collapse = ""),
                "_Results_RNAseq_143organoids.RData")
)


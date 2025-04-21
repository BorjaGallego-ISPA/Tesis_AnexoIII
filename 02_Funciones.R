###############################################################################
############### Omic Characterization of Doxorubicin Resistant ################
###################### Bone Sarcoma Models: Functions #########################
###############################################################################

# QC ####

#QC.LibSizeBarplot <- function(dds){
#  data <- as.data.frame(colData(dds)) %>%
#    mutate(LibSize = colSums(assay(dds)) / 1e6) %>%
#    select(cell_line, condition, replicate, LibSize)
#  plot <- ggplot(data, aes(x = forcats::fct_inorder(rownames(data)),
#                           y = LibSize)) +
#    geom_bar(aes(fill = condition), stat = "identity",
#             position = position_dodge(.8), width = .7) +
#    ggtitle("Total Mapped Counts") +
#    labs(y = "Counts (millions)", x = "Sample") +
#    theme_bw() +
#    theme(
#      legend.title = element_blank(),
#      plot.title = element_text(hjust = .6, face = "bold", size = 22),
#      axis.title.y = element_text(size = 18, face = "bold"),
#      axis.title.x = element_blank(),
#      axis.text.x = element_text(size = 14, angle = 45, vjust = .6),
#      axis.text.y = element_text(size = 14)
#    ) +
#    scale_fill_manual(values = c("#0000FF", "#FF0000")) +
#    scale_x_discrete(labels = paste(data$cell_line, data$replicate, sep = "_")) +
#    scale_y_continuous(breaks = seq(0, 30, 5), expand = c(0, 0)) +
#    geom_hline(yintercept = 20, linetype = "dotted", size = .8, color = "darkgrey")
#  return(plot)
#}

QC.MakeCorrHeatmap <- function(Counts.norm, main, output.path, clustering = T,
                            heatmap.rows,
                            colors = colorRampPalette(
                              rev(RColorBrewer::brewer.pal(9, "Blues")))(255)){
  sampleDists <- dist(t(Counts.norm))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- heatmap.rows
  colnames(sampleDistMatrix) <- NULL
  if (clustering == T){
    Heatmap <- pheatmap::pheatmap(sampleDistMatrix,
                                  clustering_distance_rows = sampleDists,
                                  clustering_distance_cols = sampleDists,
                                  col = colors,
                                  main = main)
  } else if (clustering == F) {
    Heatmap <- pheatmap::pheatmap(sampleDistMatrix,
                                  clustering_distance_rows = sampleDists,
                                  clustering_distance_cols = sampleDists,
                                  col = colors,
                                  main = main,
                                  cluster_cols = F, cluster_rows = F)
  } else {
    print("Error: Please set clustering mode in TRUE or FALSE")
  }
  ggsave(output.path, plot = Heatmap, device = "tiff", dpi = 300)
  return(Heatmap)
}

QC.MakePCA <- function(data, metadata, ntop = 500, nPC = 2, main, output.path,
                       fill.color = NULL){
  rv <- matrixStats::rowVars(as.matrix(data))
  select <- order(rv, decreasing = T)[seq_len(min(ntop,length(rv)))]
  pca <- prcomp(t(data[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  intgroup.df <- as.data.frame(metadata[rownames(metadata) %in% colnames(data), ,
                                        drop = F])
  d <- cbind(data.frame(intgroup.df),
             pca$x[,seq_len(min(nPC, ncol(pca$x))), drop = F])
  attr(d, "percentVar") <- percentVar[1:nPC]
  for (i in 1:nPC) {
    percentVar[i] <- round(percentVar[i]*100)
  }
  if (length(unique(d$cell_line)) == 1) {
    if (is.null(fill.color)) {
      p <- ggplot(d, aes(x = PC1, y = PC2, color = condition)) +
        geom_point(size = 5) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle(main) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = .5, face = "bold", size = 24),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(face = "bold", size = 18),
          legend.title = element_blank()
        ) +
        coord_fixed()
    } else {
      p <- ggplot(d, aes(x = PC1, y = PC2, color = condition)) +
        geom_point(size = 5) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle(main) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = .5, face = "bold", size = 24),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(face = "bold", size = 18),
          legend.title = element_blank()
        ) +
        scale_color_manual(values = fill.color)
        coord_fixed()
    }
  } else {
    if (is.null(fill.color)) {
      p <- ggplot(d, aes(x = PC1, y = PC2, color = condition, shape = cell_line)) +
        geom_point(size = 5) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle(main) +
        labs(shape = "Cell Line", colour = "Condition") +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = .5, face = "bold", size = 24),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.background = element_blank(),
          legend.text = element_text(size = 18),
          legend.title = element_text(face = "bold", size = 20)
        ) +
        coord_fixed()
    } else {
      p <- ggplot(d, aes(x = PC1, y = PC2, color = condition, shape = cell_line)) +
        geom_point(size = 5) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle(main) +
        labs(shape = "Cell Line", colour = "Condition") +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = .5, face = "bold", size = 24),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.background = element_blank(),
          legend.text = element_text(size = 18),
          legend.title = element_text(face = "bold", size = 20)
        ) +
        scale_color_manual(values = fill.color)
      coord_fixed()
    }
  }
  PCA.res <- list(Positions = d, PercentVar = percentVar, Plot = p)
  ggsave(output.path, plot = p, device = "tiff", dpi = 300)
  return(PCA.res)
}

# DEA ####

DEA.DESeq2Table <- function(dds, contrast, output.csv, output.xlsx){
  res <- as.data.frame(DESeq2::results(dds, contrast = contrast)) %>%
    mutate(
      Symbol = mapIds(org.Hs.eg.db, keys = rownames(.),
                      column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
      Biotype = mapIds(org.Hs.eg.db, keys = rownames(.),
                       column = "GENETYPE", keytype = "ENSEMBL", multiVals = "first"),
      Description = mapIds(org.Hs.eg.db, keys = rownames(.),
                           column = "GENENAME", keytype = "ENSEMBL", multiVals = "first"),
      pvalue = ifelse(pvalue == 0, .Machine$double.xmin, pvalue),
      padj = ifelse(padj == 0, .Machine$double.xmin, padj)
    ) %>%
    filter(!is.na(padj)) %>%
    filter(!is.na(Symbol)) %>%
    arrange(padj) %>%
    distinct(Symbol, .keep_all = T) %>%
    arrange(Symbol)
  write.csv(res, output.csv)
  xlsx::write.xlsx2(res, output.xlsx)
  return(res)
}

DEA.DEqMS2Table <- function(data, condition, cont, razor, coef_col = 1,
                            output.csv, output.xlsx) {
  design <- model.matrix(~0 + condition)
  colnames(design) <- gsub("condition", "", colnames(design))
  contrast <- makeContrasts(contrasts = cont, levels = design)
  fit1 <- lmFit(data, design)
  fit2 <- contrasts.fit(fit1, contrasts = contrast)
  fit3 <- eBayes(fit2)
  fit3$count <- razor[rownames(fit3$coefficients), "count"]
  fit4 <- spectraCounteBayes(fit3)
  if (length(coef_col) == 1){
    res <- outputResult(fit4, coef_col = coef_col)
  } else{
    res <- list()
    for(i in 1:length(coef_col)){
      res[[i]] <- outputResult(fit4, coef_col = coef_col[i])
    }
  }
  if (length(output.csv) == 1){
    write.csv(res, output.csv)
    xlsx::write.xlsx2(res, output.xlsx)
  } else{
    for(i in 1:length(output.csv)){
      write.csv(res[[i]], output.csv[i])
      xlsx::write.xlsx2(res[[i]], output.xlsx[i])
    }
  }
  return(res)
  }

DEA.limma2Table <- function(data, condition, cont, coef_col = 1, output.csv, output.xlsx) {
  design <- model.matrix(~0 + condition)
  colnames(design) <- gsub("condition", "", colnames(design))
  contrast <- makeContrasts(contrasts = cont, levels = design)
  fit1 <- lmFit(data, design)
  fit2 <- contrasts.fit(fit1, contrasts = contrast)
  fit3 <- eBayes(fit2)
  if (length(coef_col) == 1){
    res <- topTable(fit3, coef = coef_col, number = Inf)
  } else{
    res <- list()
    for(i in 1:length(coef_col)){
      res[[i]] <- topTable(fit3, coef = coef_col[i], number = Inf)
    }
  }
  if (length(output.csv) == 1){
    write.csv(res, output.csv)
    xlsx::write.xlsx2(res, output.xlsx)
  } else{
    for(i in 1:length(output.csv)){
      write.csv(res[[i]], output.csv[i])
      xlsx::write.xlsx2(res[[i]], output.xlsx[i])
    }
  }
  return(res)
}

DEA.MakeVolcano <- function(DEA.table, top = 10, lfc.thr = 1, pval.thr = .05,
                        xlim = Inf, ylim = Inf, main, output.path,
                        top.order = "both", sym.coln, lfc.coln, fdr.coln) {
  if (sym.coln == "row"){
    Data <- DEA.table %>%
      rownames_to_column(var = "Symbol") %>%
      select(Symbol,
             Log2FC = all_of(lfc.coln),
             FDR = all_of(fdr.coln)) %>%
      mutate(
        log10FDR = -log(FDR, 10),
        Expression = case_when(Log2FC >= lfc.thr & FDR <= pval.thr ~
                                 "Up-Regulated",
                               Log2FC <= -lfc.thr & FDR <= pval.thr ~
                                 "Down-Regulated",
                               T ~ "NS")
      )
  } else{
    Data <- DEA.table %>%
      select(Symbol = all_of(sym.coln), Log2FC = all_of(lfc.coln),
             FDR = all_of(fdr.coln)) %>%
      mutate(
        log10FDR = -log(FDR, 10),
        Expression = case_when(Log2FC >= lfc.thr & FDR <= pval.thr ~
                                 "Up-Regulated",
                               Log2FC <= -lfc.thr & FDR <= pval.thr ~
                                 "Down-Regulated",
                               T ~ "NS")
      )
  }
  for (i in 1:nrow(Data)) {
    if (Data$Log2FC[i] > xlim) {
      Data$Log2FC[i] <- xlim
    } else {
      Data$Log2FC[i] <- Data$Log2FC[i]
    }
    if (Data$Log2FC[i] < -xlim) {
      Data$Log2FC[i] <- -xlim
    } else {
      Data$Log2FC[i] <- Data$Log2FC[i]
    }
    if (Data$log10FDR[i] > ylim && is.na(Data$log10FDR[i]) == F) {
      Data$log10FDR[i] <- ylim
    } else {
      Data$log10FDR[i] <- Data$log10FDR[i]
    }
  }
  if (top.order == "both"){
    top_genes <- bind_rows(Data %>%
                             filter(Expression == "Up-Regulated") %>%
                             arrange(FDR, desc(abs(Log2FC))) %>%
                             head(top),
                           Data %>%
                             filter(Expression == "Down-Regulated") %>%
                             arrange(FDR, desc(abs(Log2FC))) %>%
                             head(top))
  } else if (top.order == "lfc") {
    top_genes <- bind_rows(Data %>%
                             filter(Expression == "Up-Regulated") %>%
                             arrange(desc(abs(Log2FC))) %>%
                             head(top),
                           Data %>%
                             filter(Expression == "Down-Regulated") %>%
                             arrange(desc(abs(Log2FC))) %>%
                             head(top))
  } else if (top.order == "pval") {
    top_genes <- bind_rows(Data %>%
                             filter(Expression == "Up-Regulated") %>%
                             arrange(FDR) %>%
                             head(top),
                           Data %>%
                             filter(Expression == "Down-Regulated") %>%
                             arrange(FDR) %>%
                             head(top))
  } else{
    print("Please set an order variable for top genes")
  }
  plot <- ggplot(Data, aes(Log2FC, log10FDR)) +
    geom_point(aes(color = Expression), size = 1) +
    xlab("log2FC") +
    ylab("-log10FDR") +
    scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    geom_label_repel(data = top_genes,
                     mapping = aes(Log2FC, log10FDR, label = Symbol),
                     size = 6, max.overlaps = 2*top) +
    theme(
      legend.title = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 1, colour = "black"),
      plot.title = element_text(hjust = .5, face = "bold", size = 22),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 14, face = "bold"),
      legend.key  = element_blank(),
      legend.background = element_blank()
      ) +
    scale_y_continuous(expand = c(.01,0)) +
    scale_x_continuous(expand = c(.01,0)) +
    geom_vline(aes(xintercept = 0), colour = "black", linewidth = 1) + 
    geom_vline(aes(xintercept = lfc.thr), colour = "gray20", linewidth = 1,
               lty = "dotted") +
    geom_vline(aes(xintercept = -lfc.thr), colour = "gray20", linewidth = 1,
               lty = "dotted") +
    geom_hline(aes(yintercept = -log(pval.thr, 10)), colour = "gray20",
               linewidth = 1, lty = "dotted") +
    ggtitle(main)
  ggsave(output.path, plot = plot, device = "tiff", dpi = 300)
  return(plot)
}

DEA.MakeVenn <- function(data.l, lfc.thr = 1, pval.thr = .05, up.main, up.file,
                         down.main, down.file, lfc.coln, fdr.coln,
                         Symbol = "Symbol"){
  data.l_up <- data.l
  if (Symbol == "row"){
    for (i in 1:length(data.l_up)){
      data.l_up[[i]] <- data.l_up[[i]] %>%
        rownames_to_column(var = "Symbol") %>%
        filter(data.l_up[[i]][lfc.coln] >= lfc.thr & data.l_up[[i]][fdr.coln] <= pval.thr) %>%
        pull(Symbol)
    }
  } else{
    for (i in 1:length(data.l_up)){
      data.l_up[[i]] <- data.l_up[[i]] %>%
        filter(data.l_up[[i]][lfc.coln] >= lfc.thr & data.l_up[[i]][fdr.coln] <= pval.thr) %>%
        pull(Symbol)
    }
  }
  VD.up_plot <- ggVennDiagram(data.l_up, label_alpha = 0, label = "count") +
    scale_fill_gradient(low = "gray80",high = "darkred") +
    #scale_color_manual(values = rep("darkgrey", length(data.l_up))) +
    ggtitle(up.main) +
    scale_x_continuous(expand = expansion(mult = .3)) +
    theme()
  VD.up_data <- calculate.overlap(data.l_up)
  data.l_down <- data.l
  if (Symbol == "row"){
    for (i in 1:length(data.l_down)){
      data.l_down[[i]] <- data.l_down[[i]] %>%
        rownames_to_column(var = "Symbol") %>%
        filter(data.l_down[[i]][lfc.coln] <= -lfc.thr & data.l_down[[i]][fdr.coln] <= pval.thr) %>%
        pull(Symbol)
    } 
  } else{
      for (i in 1:length(data.l_down)){
        data.l_down[[i]] <- data.l_down[[i]] %>%
          filter(data.l_down[[i]][lfc.coln] <= -lfc.thr & data.l_down[[i]][fdr.coln] <= pval.thr) %>%
          pull(Symbol)
      }
  }
  VD.down_plot <- ggVennDiagram(data.l_down, label_alpha = 0, label = "count") +
    scale_fill_gradient(low = "gray80",high = "darkred") +
    #scale_color_manual(values = rep("darkgrey", length(data.l_down))) +
    ggtitle(down.main) +
    scale_x_continuous(expand = expansion(mult = .3)) +
    theme()
  VD.down_data <- calculate.overlap(data.l_down)
  ggsave(up.file, plot = VD.up_plot, device = "tiff", dpi = 300)
  ggsave(down.file, plot = VD.down_plot, device = "tiff", dpi = 300)
  VD.results <- list(Up.plot = VD.up_plot,
                     Down.plot = VD.down_plot,
                     Up.data = VD.up_data,
                     Down.data = VD.down_data)
  return(VD.results)
}

# Functional Analysis ####

## Pathway Enrichment ####

PathEnrich.MakeGSEA <- function(GeneSet, Res, table.csv, table.xlsx, lfc.coln,
                                kegg.ids = NULL, kegg.dir = ".",
                                SymInRow = F){
  set.seed(123)
  if (SymInRow == T) {
    GeneList <- Res %>%
      rownames_to_column(var = "Symbol") %>%
      mutate(EntrezID = mapIds(org.Hs.eg.db,
                               keys = Symbol,
                               column = "ENTREZID",
                               keytype = "SYMBOL", 
                               multiVals = "first")) %>%
      dplyr::filter(!is.na(EntrezID)) %>%
      pull(lfc.coln, name = EntrezID) %>%
      sort(decreasing = T)
  } else {
    GeneList <- Res %>%
      mutate(EntrezID = mapIds(org.Hs.eg.db,
                               keys = Symbol,
                               column = "ENTREZID",
                               keytype = "SYMBOL", 
                               multiVals = "first")) %>%
      dplyr::filter(!is.na(EntrezID)) %>%
      pull(lfc.coln, name = EntrezID) %>%
      sort(decreasing = T)
  }
  gse <- GSEA(GeneList, TERM2GENE = GeneSet, seed = T)
  Count <- NULL
  for (i in 1:nrow(gse@result)){
    Count[i] <- strsplit(gse@result$core_enrichment[i], split = "/") %>%
      unlist() %>%
      length()
  }
  Core <- NULL
  for (i in 1:nrow(gse@result)){
    Core[i] <- mapIds(org.Hs.eg.db,
                      keys = unlist(strsplit(gse@result$core_enrichment[i],
                                             split = "/")),
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first") %>%
      paste(collapse = "/")
  }
  gse.table <- gse@result %>%
    mutate(Description = sapply(strsplit(gse@result$Description,
                                         split = "_"), function(x){
                                           paste(x[2:length(x)], collapse = " ")
                                         }),
           Count = Count,
           Gene_Ratio = Count/setSize,
           core_enrichment = Core) %>%
    arrange(desc(NES)) %>%
    select(Description, NES, p.adjust, Count, setSize, Gene_Ratio, core_enrichment)
  gse.plots <- list()
  for (i in 1:nrow(gse@result)){
    anno <- gse@result[i, c("enrichmentScore", "p.adjust")]
    names(anno) <- c("ES", "FDR")
    lab <- paste0(names(anno), " ", "=", " ", signif(anno, 3), collapse = "\n")
    gse.plots[[i]] <- gseaplot(gse, geneSetID = i, by = "runningScore",
                               title = gse@result$Description[i]) +
      ylab("Running ES") + xlab("Ranked Genes") +
      theme(plot.title = element_text(size = 15, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 8, face = "bold"),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 7)) +
      annotate("text", 2000 + 1500 * sign(gse@result[i, "enrichmentScore"]),
               gse@result[i, "enrichmentScore"] * .75, label = lab,
               hjust = 0, vjust = 0)
  }
  write.csv(gse.table, table.csv)
  xlsx::write.xlsx2(gse.table, table.xlsx)
  if (is.null(kegg.ids)){
    gse.res <- list(Results.table = gse.table,
                    GSEA.plots = gse.plots)
  } else{
    kegg.pathways <- list()
    for(i in 1:length(kegg.ids)){
      kegg.pathways[[i]] <- pathview(gene.data = GeneList,
                                     pathway.id = kegg.ids[i],
                                     species = "hsa",
                                     limit = list(gene = max(abs(GeneList)),
                                                  cpd = 1), kegg.dir = kegg.dir)
      gse.res <-list(Results.table = gse.table,
                     GSEA.plots = gse.plots,
                     KEGG.paths = kegg.pathways)
    }
  }
  return(gse.res)
}

PathEnrich.MakeBubble <- function(Data, main, output.path, x.expand = waiver()){ 
  p <- ggplot(Data, aes(x = NES, y = fct_reorder(Description, NES),
                        size = Gene_Ratio, color = p.adjust)) +
    geom_point() +
    scale_size(range = c(6, 15), name = "Gene Ratio") +
    ggtitle(main) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 20, face = "bold"),
      plot.title = element_text(hjust = .5, face = "bold", size = 22),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18)
    ) +
    scale_x_continuous(expand = x.expand) +
    scale_color_gradient2(trans = "log10")
  ggsave(output.path, p, device = "tiff", dpi = 300,
         units = "px", width = 3883, height = 2399
         )
  return(p)
}

##  TF Activity Inference ####

ChooseTFNet <- function(db = ""){
  if (db == "dorothea"){
    net <- decoupleR::get_dorothea(organism = "human", levels = c("A", "B", "C"))
  } else if (db == "tri"){
    net <- decoupleR::get_collectri(organism = "human", split_complexes = F)
  } else if (db == ""){
    stop("Please choose a database. Databases available: dorothea and tri")
  } else{
    stop("Database not available. Databases available: dorothea and tri")
  }
}

Get_TFs <- function(mat, net, condition, method = "", times = 1000, minsize = 15,
                    cores = 4, pval.thr = NULL, score.thr = NULL, n_tfs = Inf,
                    main, output.path){
  if (any(is.na(mat$Symbol)) == T | any(mat$Symbol == "") == T) {
    mat <- mat %>%
      filter(Symbol != "") %>%
      filter(!is.na(Symbol)) %>% 
      distinct(Symbol, .keep_all = T) %>%
      pull(condition, name = Symbol) %>%
      as.matrix()
  } else {
    mat <- mat %>%
      distinct(Symbol, .keep_all = T) %>%
      pull(condition, name = Symbol) %>%
      as.matrix()
  }
  colnames(mat) <- condition
  if (method == "WMEAN") {
    contrast_acts <- run_wmean(mat = mat, net = net,
                               .source = "source", .target = "target",
                               .mor = "mor", times = times, minsize = minsize)
    if (is.null(pval.thr) & is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(statistic == "norm_wmean") %>%
        mutate(rnk = NA)
    } else if (!is.null(pval.thr) & is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(statistic == "norm_wmean") %>%
        filter(p_value <= pval.thr) %>%
        mutate(rnk = NA)
    } else if (is.null(pval.thr) & !is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(statistic == "norm_wmean") %>%
        filter(abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    } else if(!is.null(pval.thr) & !is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(statistic == "norm_wmean") %>%
        filter(p_value <= pval.thr & abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    }
  } else if (method == "VIPER") {
    contrast_acts <- run_viper(mat = mat, net = net,
                               .source = "source", .target = "target",
                               .mor = "mor", minsize = minsize, pleiotropy = F,
                               eset.filter = F, cores = cores)
    if (is.null(pval.thr) & is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        mutate(rnk = NA)
    } else if (!is.null(pval.thr) & is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(p_value <= pval.thr) %>%
        mutate(rnk = NA)
    } else if (is.null(pval.thr) & !is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    } else if(!is.null(pval.thr) & !is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(p_value <= pval.thr & abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    }
  } else if (method == "ULM"){
    contrast_acts <- run_ulm(mat = mat, net = net,
                             .source = "source", .target = "target",
                             .mor = "mor", minsize = minsize)
    if (is.null(pval.thr) & is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        mutate(rnk = NA)
    } else if (!is.null(pval.thr) & is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(p_value <= pval.thr) %>%
        mutate(rnk = NA)
    } else if (is.null(pval.thr) & !is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    } else if(!is.null(pval.thr) & !is.null(score.thr)){
      f_contrast_acts <- contrast_acts %>%
        filter(p_value <= pval.thr & abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    }
  } else if (method == ""){
    stop("Please choose a method. Methods available: WMEAN, VIPER and ULM")
  } else {
    stop("Method not available. Methods available: WMEAN, VIPER and ULM")
  }
  msk <- f_contrast_acts$score > 0
  f_contrast_acts[msk, "rnk"] <- rank(-f_contrast_acts[msk, "score"])
  f_contrast_acts[!msk, "rnk"] <- rank(-abs(f_contrast_acts[!msk, "score"]))
  tfs <- f_contrast_acts %>%
    arrange(rnk) %>%
    head(n_tfs) %>%
    pull(source)
  f_contrast_acts <- f_contrast_acts %>%
    filter(source %in% tfs)
  Barplot <- ggplot(f_contrast_acts, aes(x = reorder(source, score),
                                         y = score)) +
    geom_bar(aes(fill = score), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "darkred", 
                         mid = "whitesmoke", midpoint = 0) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(face = "bold", size = 18),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 18, face = "bold"),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(hjust = .5, face = "bold", size = 22),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(size = 14)
      ) +
    ylab("NES") +
    ggtitle(main)
  ggsave(output.path, plot = Barplot, device = "tiff", dpi = 300,
         units = "px", width = 3883, height = 2399
         )
  TFs.results <- list(TFs.list = f_contrast_acts, TFs.barplot = Barplot)
  return(TFs.results)
}

## Kinases Activity Inference ####

MakeKinNet <- function(){
  uniprot_kinases <- import_omnipath_annotations(resources = "UniProt_keyword") %>%
    filter(value == "Kinase" & !grepl("COMPLEX", uniprot)) %>%
    distinct() %>%
    pull(genesymbol) %>%
    unique()
  omnipath_ptm <- signed_ptms() %>%
    filter(modification %in% c("dephosphorylation","phosphorylation")) %>%
    filter(!(str_detect(sources, "ProtMapper") & n_resources == 1)) %>%
    mutate(p_site = paste0(substrate_genesymbol, "_", residue_type,
                           residue_offset),
           mor = ifelse(modification == "phosphorylation", 1, -1)) %>%
    transmute(p_site, enzyme_genesymbol, mor) %>%
    filter(enzyme_genesymbol %in% uniprot_kinases)
  omnipath_ptm$id <- paste(omnipath_ptm$p_site, omnipath_ptm$enzyme_genesymbol,
                           sep = "")
  omnipath_ptm <- omnipath_ptm[!duplicated(omnipath_ptm$id), ]
  omnipath_ptm <- omnipath_ptm[, -5]
  names(omnipath_ptm)[c(1, 2)] <- c("target", "source")
  return(omnipath_ptm)
}

Get_Kinase <- function(mat, net, n_kin = Inf, pval.thr = NULL, main, barplot.o,
                       times = 1000, minsize = 3, kin = "", kin_pathway.o = "",
                       score.thr = NULL, method = "", cores = 3, condition){
  mat <- mat %>% mutate(psite_ID = sapply(strsplit(rownames(mat), "_"),
                                          function(x){
                                            gsub("","",
                                                 paste0(x[1],"_",x[2]))
                                          })) %>%
    arrange(adj.P.Val) %>%
    distinct(psite_ID, .keep_all = T) %>%
    remove_rownames() %>%
    column_to_rownames(var = "psite_ID") %>%
    select(t, logFC, adj.P.Val) %>%
    as.matrix()
  if (method == "WMEAN") {
    kin_activity <- run_wmean(mat = mat[, condition, drop = F], net = net,
                              .source = "source", .target = "target",
                              .mor = "mor", times = times, minsize = minsize)
    if (is.null(pval.thr) & is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(statistic == "norm_wmean") %>%
        mutate(rnk = NA)
    } else if (!is.null(pval.thr) & is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(statistic == "norm_wmean") %>%
        filter(p_value <= pval.thr) %>%
        mutate(rnk = NA)
    } else if (is.null(pval.thr) & !is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(statistic == "norm_wmean") %>%
        filter(abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    } else if(!is.null(pval.thr) & !is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(statistic == "norm_wmean") %>%
        filter(p_value <= pval.thr & abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    }
  } else if (method == "VIPER") {
    kin_activity <- run_viper(mat = mat[, condition, drop = F], net = net,
                              .source = "source", .target = "target",
                              .mor = "mor", minsize = minsize, cores = cores,
                              pleiotropy = F, eset.filter = F)
    if (is.null(pval.thr) & is.null(score.thr)){
      kin_activity <- kin_activity %>%
        mutate(rnk = NA)
    } else if (!is.null(pval.thr) & is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(p_value <= pval.thr) %>%
        mutate(rnk = NA)
    } else if (is.null(pval.thr) & !is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    } else if(!is.null(pval.thr) & !is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(p_value <= pval.thr & abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    }
  } else if (method == "ULM"){
    kin_activity <- run_ulm(mat = mat[, condition, drop = F], net = net,
                              .source = "source", .target = "target",
                              .mor = "mor", minsize = minsize)
    if (is.null(pval.thr) & is.null(score.thr)){
      kin_activity <- kin_activity %>%
        mutate(rnk = NA)
    } else if (!is.null(pval.thr) & is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(p_value <= pval.thr) %>%
        mutate(rnk = NA)
    } else if (is.null(pval.thr) & !is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    } else if(!is.null(pval.thr) & !is.null(score.thr)){
      kin_activity <- kin_activity %>%
        filter(p_value <= pval.thr & abs(score) >= score.thr) %>%
        mutate(rnk = NA)
    }
  } else if (method == ""){
    stop("Please choose a method. Methods available: WMEAN, VIPER and ULM")
  } else {
    stop("Method not available. Methods available: WMEAN, VIPER and ULM")
  }
  msk <- kin_activity$score > 0
  kin_activity[msk, "rnk"] <- rank(-kin_activity[msk, "score"])
  kin_activity[!msk, "rnk"] <- rank(-abs(kin_activity[!msk, "score"]))
  kins <- kin_activity %>%
    arrange(rnk) %>%
    head(n_kin) %>%
    pull(source)
  kin_activity <- kin_activity %>%
    filter(source %in% kins)
  Barplot <- ggplot(kin_activity, aes(x = reorder(source, score),
                                      y = score)) +
    geom_bar(aes(fill = score), stat = "identity") +
    scale_fill_gradient2(low = "darkblue", high = "darkred", 
                         mid = "whitesmoke", midpoint = 0) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_text(face = "bold", size = 18),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 18, face = "bold"),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(hjust = .5, face = "bold", size = 22),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    ) +
    ylab("NES") +
    ggtitle(main)
  ggsave(barplot.o, plot = Barplot, device = "tiff", dpi = 300,
         units = "px", width = 3883, height = 2399
         )
  if (kin == "") {
    Kinases.results <- list(Kinases.list = kin_activity,
                            Kinases.barplot = Barplot)
  } else {
    df <- net %>%
      filter(source == kin) %>%
      arrange(target) %>%
      mutate(ID = target, color = "3") %>%
      column_to_rownames("target")
    inter <- sort(intersect(rownames(mat), rownames(df)))
    df <- df[inter, ]
    df[, c("logFC", "t", "adj.P.Val")] <- mat[inter, ]
    df <- df %>%
      mutate(color = if_else(mor > 0 & stat > 0, "1", color)) %>%
      mutate(color = if_else(mor > 0 & stat < 0, "2", color)) %>%
      mutate(color = if_else(mor < 0 & stat > 0, "2", color)) %>%
      mutate(color = if_else(mor < 0 & stat < 0, "1", color))
    Kinase_prot.plot <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val),
                                       color = color, size = 2/5)) +
      geom_point() +
      scale_colour_manual(values = c("red", "royalblue3", "grey")) +
      geom_label_repel(aes(label = ID, size = 1)) + 
      theme_minimal() +
      theme(legend.position = "none") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      ggtitle(kin)
    ggsave(kin_pathway.o, plot = Kinase_prot.plot, device = "tiff", dpi = 300)
    Kinases.results <- list(Kinases.list = kin_activity,
                            Kinases.barplot = Barplot,
                            Kinase.path = Kinase_prot.plot)
  }
  return(Kinases.results)
}

# MultiOmic ####

RNAvsProt <- function(RNA, Prot, pval.thr = 1, lfc.thr_r = 0, lfc.thr_p = 0,
                      top = 10, main = NULL, table.csv = NULL,
                      table.xlsx = NULL, plot.tiff = NULL,
                      makePlot = T, exportData = T, merge = F,
                      pval.thr_p = pval.thr){
  
  RNA <- RNA %>%
    filter(abs(logFC) >= lfc.thr_r & adj.P.Val <= pval.thr) %>%
    remove_rownames() %>%
    select(Symbol, FC_r = logFC, FDR_r = adj.P.Val, t_r = t)
  Prot <- Prot %>%
    filter(abs(logFC) >= lfc.thr_p & adj.P.Val <= pval.thr_p) %>%
    rownames_to_column(var = "Symbol") %>%
    select(Symbol, FC_p = logFC, FDR_p = adj.P.Val, t_p = t)
  RNAvsProt.dat <- RNA %>%
    full_join(Prot, by = "Symbol") %>%
    arrange(Symbol) %>%
#    column_to_rownames(var = "Symbol") %>%
    mutate(
      Agreement = ifelse(sign(FC_r) == sign(FC_p), "Yes", "No"),
      Sign = ifelse(!is.na(FC_r), sign(FC_r), sign(FC_p)),
      CoRegulation = case_when(
        FC_r > lfc.thr_r & FC_p > lfc.thr_p ~ "CoUp",
        FC_r < -lfc.thr_r & FC_p < -lfc.thr_p ~ "CoDown",
        T ~ "No_CoRegulation"
        )
    ) %>%
    mutate(
      Agreement = ifelse(is.na(Agreement), "Unknown", Agreement)
    )
  RNAvsProt.dat <- RNAvsProt.dat %>%
    filter(Symbol != "-") # -> Para filtrar fila rara de prot (solo en saos y tcds) ??
  
  agree <- RNAvsProt.dat %>% filter(Agreement == "Yes") %>% nrow()
  disagree <- RNAvsProt.dat %>% filter(Agreement == "No") %>% nrow()
  unknowkn <- RNAvsProt.dat %>% filter(Agreement == "Unknown") %>% nrow()
  agree_per <- (agree / (agree + disagree)) * 100
  
  print(paste0("Concordance level between RNAseq and proteomic: ",
               round(agree_per, 1), "% (There are ",
               formatC(agree, big.mark = ","), " agreements, ",
               formatC(disagree, big.mark = ","), " disagreements and ",
               formatC(unknowkn, big.mark = ","), " detected in only one of them)"),
        quote = F)
  
  if (exportData == T){
    write.csv(RNAvsProt.dat, table.csv)
    xlsx::write.xlsx2(RNAvsProt.dat, table.xlsx)
  }
  if (makePlot == T){
    TopList <- bind_rows(
      RNAvsProt.dat %>%
        filter(CoRegulation == "CoUp") %>%
        arrange(desc(FC_r), desc(FC_p)) %>%
        head(top),
      RNAvsProt.dat %>%
        filter(CoRegulation == "CoDown") %>%
        arrange(FC_r, FC_p) %>%
        head(top)
    )
    RNAvsProt.plot <- ggplot(RNAvsProt.dat, aes(FC_r, FC_p)) +
      geom_point(aes(color = CoRegulation), size = 1) +
      xlab("RNAseq") +
      ylab("Proteomic") +
      scale_color_manual(values = c("dodgerblue3", "firebrick3", "gray50")) +
      guides(colour = guide_legend(override.aes = list(size = 4))) +
      geom_label_repel(data = TopList,
                       mapping = aes(FC_r, FC_p, label = Symbol),
                       size = 6, max.overlaps = 2*top) +
      theme(
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 1, colour = "black"),
        plot.title = element_text(hjust = .5, face = "bold", size = 22),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14, face = "bold"),
        legend.key  = element_blank(),
        legend.background = element_blank()
      ) +
      scale_y_continuous(expand = c(.01,0)) +
      scale_x_continuous(expand = c(.01,0)) +
      geom_vline(aes(xintercept = lfc.thr_r), colour = "gray20", linewidth = 1,
                 lty = "dotted") +
      geom_vline(aes(xintercept = -lfc.thr_r), colour = "gray20", linewidth = 1,
                 lty = "dotted") +
      geom_hline(aes(yintercept = lfc.thr_p), colour = "gray20",
                 linewidth = 1, lty = "dotted") +
      geom_hline(aes(yintercept = -lfc.thr_p), colour = "gray20",
                 linewidth = 1, lty = "dotted") +
      ggtitle(main)
    ggsave(plot.tiff, plot = RNAvsProt.plot, device = "tiff", dpi = 300)
    RNAvsProt.res <- list(
      RNAvsProt.dat = RNAvsProt.dat,
      RNAvsProt.plot = RNAvsProt.plot
    )
    return(RNAvsProt.res)
  } else if (merge == T) {
    RNAvsProt.dat <- RNAvsProt.dat %>%
      mutate(
        logFC = ifelse(!is.na(FC_r), FC_r, FC_p),
        padj = ifelse(!is.na(FDR_r), FDR_r, FDR_p),
        t = ifelse(!is.na(t_r), t_r, t_p)
      ) %>%
      select(Symbol, logFC, padj, t, Agreement)
    return(RNAvsProt.dat)
  } else {
    return(RNAvsProt.dat)
  }
}

LEAN.MakeSubNet <- function(STRING_id, lean.res, lean.metadata, main = NULL,
                            submain = NULL, sel = F){
  if (sel == T){
    Nodes <- get.ls.info(STRING_id, lean.res) %>%
      left_join(lean.metadata, by = c("ID" = "STRING_id")) %>%
      column_to_rownames(var = "ID") %>%
      filter(selected == T | rownames(.) == STRING_id) %>%
      mutate(
        label = Symbol,
        color = case_when(Sign == 1 ~ "red", Sign == -1 ~ "blue", .default = "grey"),
      ) %>%
      select(id = Symbol, FDR, Sign, Origin, label, color)
    Edges <- string.db$get_interactions(rownames(Nodes)) %>%
      mutate(
        from = Nodes[from, "id"],
        to = Nodes[to, "id"]
      ) %>%
      select(from, to)
    Nodes <- Nodes %>%
      filter(id %in% unique(c(Edges$from, Edges$to)))
    Net <- visNetwork(Nodes, Edges, main = main, submain = submain,
                      width = "100%") %>%
      visNetwork::visIgraphLayout(layout = "layout_with_kk", randomSeed = 1) %>%
      visEdges(
        #arrows = "to",
        color = "dimgrey",
        width = 2
      ) %>%
      visNodes(font = "30px arial black bold") %>%
      #visLegend() %>%
      visOptions(
        highlightNearest = T,
        nodesIdSelection = T,
      )
  } else {
    Nodes <- get.ls.info(STRING_id, lean.res) %>%
      left_join(lean.metadata, by = c("ID" = "STRING_id")) %>%
      column_to_rownames(var = "ID") %>%
      mutate(
        label = Symbol,
        color = case_when(Sign == 1 ~ "red", Sign == -1 ~ "blue", .default = "grey"),
        selected = ifelse(selected == T, "selected", "interacts")
      ) %>%
      select(id = Symbol, FDR, Sign, Origin, selected, label, color)
    Edges <- string.db$get_interactions(rownames(Nodes)) %>%
      mutate(
        from = Nodes[from, "id"],
        to = Nodes[to, "id"]
      ) %>%
      select(from, to)
    Nodes <- Nodes %>%
      filter(id %in% unique(c(Edges$from, Edges$to)))
    Net <- visNetwork(Nodes, Edges, main = main, submain = submain,
                      width = "100%") %>%
      visNetwork::visIgraphLayout(layout = "layout_with_kk", randomSeed = 1) %>%
      visEdges(
        #arrows = "to",
        color = "dimgrey",
        width = 2
      ) %>%
      visNodes(font = "30px arial black bold") %>%
      #visLegend() %>%
      visOptions(
        highlightNearest = T,
        nodesIdSelection = T,
        selectedBy = "selected"
      )
  }
  SubNet <- list(
    Nodes = Nodes,
    Edges = Edges,
    Net = Net
  )
  return(SubNet)
}

LEAN.MakeClusters <- function(data.lean, db){
  Edges <- string.db$get_interactions(rownames(data.lean)) %>%
    select(from, to) %>%
    distinct(from, to, .keep_all = T) %>%
    mutate(
      from = data.lean[from, "Symbol"],
      to = data.lean[to, "Symbol"]
    )
  res_graph <- graph_from_data_frame(Edges)
  cfg <- as.undirected(res_graph, mode = "collapse") %>%
    cluster_fast_greedy()
  clusters <- data.frame(Symbol = cfg$names, Cluster = cfg$membership)
  Nodes <- data.lean %>%
    left_join(clusters, by = "Symbol") %>%
    mutate(
      Cluster = as.factor(Cluster)
      ) %>%
    filter(!is.na(Cluster))
  clusters.l <- list()
  for (i in 1:length(unique(Nodes$Cluster))){
    clusters.l[[i]] <- Nodes %>%
      filter(Cluster == i) %>%
      mutate(
        EntrezID = mapIds(
          org.Hs.eg.db,
          keys = Symbol,
          column = "ENTREZID",
          keytype = "SYMBOL", 
          multiVals = "first"
        )
      ) %>%
      filter(!is.na(EntrezID)) %>%
      pull(EntrezID) %>%
      unname()
  }
  names(clusters.l) <- 1:length(clusters.l)
  clusters.enrich <- compareCluster(
    geneClusters = clusters.l,
    fun = enricher,
    TERM2GENE = db
  ) %>%
    setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  Nodes <- Nodes %>%
    left_join(
      select(
        distinct(clusters.enrich@compareClusterResult, Cluster, .keep_all = T),
        Cluster, ID),
      by = "Cluster") %>%
    mutate(
      ID = ifelse(is.na(ID), paste0("Cluster ", Cluster),
                  sapply(strsplit(ID, split = "_"), function(x){
                    paste(x[2:length(x)], collapse = " ")
                    })
                  ),
      title = Symbol
      ) %>%
    dplyr::rename(id = Symbol, group = ID)
  Cluster.res <- list(
    Nodes = Nodes,
    Edges = Edges,
    Enrichment = clusters.enrich
  )
  return(Cluster.res)
}

Enrich.MakeORA <- function(data, db = c("KEGG", "RPA", "GO", "UNI"), t2g = NULL,
                           ont = "ALL", uni.go = NULL, main, x.expand = waiver(),
                           output.tiff){
  entrez <- data %>%
    mutate(
      EntrezID = mapIds(
        org.Hs.eg.db,
        keys = Symbol,
        column = "ENTREZID",
        keytype = "SYMBOL", 
        multiVals = "first"
      )
    ) %>%
    filter(!is.na(EntrezID)) %>%
    pull(EntrezID) %>%
    unname()
  if (length(db) > 1){
    stop("Error: Analysis method not specified. Choose between \"KEGG\", \"RPA\" (Reactome), \"GO\" or \"UNI\" (Universal)")
  } else if(db == "KEGG"){
    ora.res <- enrichKEGG(
      gene = entrez,
      pvalueCutoff = .05
    ) %>%
      setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  } else if (db == "RPA"){
    ora.res <- enrichPathway(
      gene = entrez,
      pvalueCutoff = .05,
      readable = T
    )
  } else if (db == "GO"){
    ora.res <- enrichGO(
      gene = entrez,
      universe = uni.go,
      ont = ont,
      pvalueCutoff = .05,
      OrgDb = org.Hs.eg.db,
      readable = T
    )
  } else if (db == "UNI"){
    if (is.null(t2g)){
      print("Error: TERM2GENE is empty (t2g = NULL)")
    }
    ora.res <- enricher(
      gene = entrez,
      TERM2GENE = t2g,
      pvalueCutoff = .05
    ) %>%
      setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  } else {
    stop("Error: Analysis method not available Choose between \"KEGG\", \"RPA\" (Reactome), \"GO\" or \"UNI\" (Universal)")
    
  }
  ora.sig <- ora.res@result %>%
    filter(p.adjust < .05) %>%
    mutate(
      ID = sapply(
        strsplit(ID, split = "_"),
        function(x){
          paste(x[2:length(x)], collapse = " ")
        }
      ),
      GeneRatio = sapply(
        strsplit(GeneRatio, split = "/"),
        function(x){
          as.numeric(x[1]) / as.numeric(x[2])
        }
      )
    ) %>%
    select(ID, Count, GeneRatio, p.adjust, geneID)
  ora.plot <- ggplot(ora.sig, aes(x = GeneRatio, y = fct_reorder(ID, GeneRatio),
                               size = Count, color = p.adjust)) +
    geom_point() +
    scale_size(range = c(5, 14), name = "Count") +
    ggtitle(main) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 20, face = "bold"),
      plot.title = element_text(hjust = .5, face = "bold", size = 22),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18)
    ) +
    scale_x_continuous(expand = x.expand) +
    scale_color_gradient2(trans = "log10")
  ggsave(output.tiff, ora.plot, device = "tiff", dpi = 300, units = "px",
         width = 3883, height = 2399)
  ora <- list(
    Raw = ora.res,
    Table = ora.sig,
    Plot = ora.plot
  )
  return(ora)
}

clean_omnipath <- function(){
  PKN_full <- as.data.frame(OmnipathR::import_omnipath_interactions()) %>%
    filter(!is.na(references))
  PKN_clean <- PKN_full %>%
    filter(consensus_stimulation == 1 | consensus_inhibition == 1) %>%
    mutate(
      sign = consensus_stimulation - consensus_inhibition
    ) %>%
    select(
      source = source_genesymbol,
      interaction = sign,
      target = target_genesymbol
    )
  PKN_clean_supp <- PKN_clean %>%
    filter(interaction == 0) %>%
    mutate(
      interaction = -1
    )
  PKN_clean <- PKN_clean %>%
    mutate(
      interaction = ifelse(interaction == 0, 1, interaction)
    ) %>%
    rbind(PKN_clean_supp)
  PKN_filtered <- cosmosR::meta_network_cleanup(PKN_clean)
  return(PKN_filtered)
}

MOON.run <- function(meta_network, downstream_input, upstream_input, RNA_input,
                     n_steps, TFN, downstream_remover = NULL,
                     compressed_mode = T){
  
  if (!is.null(downstream_remover)) {
    
    meta_network <- meta_network[!(meta_network$source %in%
                                     names(downstream_remover)) &
                                   !(meta_network$target %in%
                                       names(downstream_remover)), ]
  }
  
  ## Remove genes that are not expressed from the meta_network
  meta_network <- cosmosR:::filter_pkn_expressed_genes(names(RNA_input),
                                                       meta_pkn = meta_network)
  
  ## Filter inputs and prune the meta_network to only keep nodes that can be
    ## found downstream of the inputs the number of step is quite flexible,
      ## 7 steps already covers most of the network
  upstream_input <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_input,
                                                            meta_network)
  downstream_input <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_input,
                                                              meta_network)
  meta_network <- cosmosR:::keep_controllable_neighbours(meta_network, n_steps,
                                                         names(upstream_input))
  downstream_input <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_input,
                                                              meta_network)
  meta_network <- cosmosR:::keep_observable_neighbours(meta_network, n_steps,
                                                       names(downstream_input))
  upstream_input <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_input,
                                                            meta_network)
  
  if (compressed_mode == T) { ## Compress the network
    
    meta_network_compressed_list <- compress_same_children(
      meta_network, sig_input = upstream_input, metab_input = downstream_input
      )
    meta_network_compressed <- meta_network_compressed_list$compressed_network
    node_signatures <- meta_network_compressed_list$node_signatures
    duplicated_parents <- meta_network_compressed_list$duplicated_signatures
    meta_network_compressed <- meta_network_cleanup(meta_network_compressed)
    
    meta_network_up_to_down <- meta_network_compressed
    
  } else {
    
    meta_network_up_to_down <- meta_network
    
  }
  
  ## MOON iterations
  before <- 1
  after <- 0
  i <- 1
  while (before != after & i < 10) {
    before <- length(meta_network_up_to_down[,1])
    start_time <- Sys.time()
    moon_res <- moon(upstream_input = upstream_input, 
                     downstream_input = downstream_input, 
                     meta_network = meta_network_up_to_down, 
                     n_layers = n_steps, 
                     statistic = "ulm") 
    
    meta_network_up_to_down <- filter_incohrent_TF_target(
      moon_res, TFN, meta_network_up_to_down, RNA_input
    )
    end_time <- Sys.time()
    
    print(end_time - start_time)
    
    after <- length(meta_network_up_to_down[,1])
    i <- i + 1
  }
  
  if (i < 10)
  {
    print(
      paste("Converged after ", paste(i - 1," iterations", sep = ""), sep = "")
      )
  } else
  {
    print(
      paste("Interupted after ",
            paste(i," iterations. Convergence uncertain.", sep = ""),
            sep = "")
      )
  }
  
  if (compressed_mode == T) { ## Decompress MOON Result
    
    source("./Support/support_decompression.R", local = T)
    
    moon_res <- decompress_moon_result(
      moon_res, meta_network_compressed_list, meta_network_up_to_down
      )
    
  }
  
  res <- list(
    moon_res = moon_res,
    meta_network_filtered = meta_network,
    upstream_inputs_filtered = upstream_input
  )
  
  return(res)
}

COSMOS.run <- function(Carni_opt = CARNIVAL_options, PKN_input, expr_input, tf_reg,
                       upstream_input, downstream_input, max.depth = 8, net_cl = T,
                       t_prefor = CARNIVAL_options$timelimit, t_for = 3200,
                       t_prerev = CARNIVAL_options$timelimit, t_rev = 6400,
                       diff_exp_thr = 1, main, submain, ksn_reg){
  # Forward run
  Carni_opt$timelimit <- t_prefor
  Cosmos.prefor <- preprocess_COSMOS_signaling_to_metabolism(
    meta_network = PKN_input,
    signaling_data = upstream_input,
    metabolic_data = downstream_input,
    diff_expression_data = expr_input,
    maximum_network_depth = max.depth,
    remove_unexpressed_nodes = net_cl,
    CARNIVAL_options = Carni_opt,
    tf_regulon = tf_reg,
    diff_exp_threshold = diff_exp_thr
  )
  
  Carni_opt$timelimit <- t_for
  Cosmos.for <- run_COSMOS_signaling_to_metabolism(
    data = Cosmos.prefor,
    CARNIVAL_options = Carni_opt
  )
  Cosmos.for.res <- format_COSMOS_res(Cosmos.for)
  
  # Reverse run
  Carni_opt$timelimit <- t_prerev
  Cosmos.prerev <- preprocess_COSMOS_metabolism_to_signaling(
    meta_network = PKN_input,
    signaling_data = upstream_input,
    metabolic_data = downstream_input,
    diff_expression_data = expr_input,
    maximum_network_depth = max.depth,
    remove_unexpressed_nodes = net_cl,
    CARNIVAL_options = Carni_opt,
    tf_regulon = tf_reg,
    diff_exp_threshold = diff_exp_thr
  )
  
  Carni_opt$timelimit <- t_rev
  Cosmos.rev <- run_COSMOS_metabolism_to_signaling(
    data = Cosmos.prerev,
    CARNIVAL_options = Carni_opt
  )
  Cosmos.rev.res <- format_COSMOS_res(Cosmos.rev)
  
  # Merge and Format Results
  Edges <- as.data.frame(
    rbind(
      Cosmos.for.res[[1]], Cosmos.rev.res[[1]]
    )
  ) %>%
    filter(Weight != 0) %>%
    select(from = Node1, sign = Sign, to = Node2) %>%
    mutate(smooth = T,
           dashes = ifelse(sign > 0, F, T))  %>%
    distinct(from, sign, to, .keep_all = T)
  Nodes <- as.data.frame(
    rbind(
      Cosmos.for.res[[2]], Cosmos.rev.res[[2]]
    )
  ) %>%
    filter(AvgAct != 0) %>%
    unique() %>%
    distinct(Nodes, .keep_all = T) %>%
    mutate(
      label = Nodes,
      measured = ifelse(AvgAct > 0, "active", "inactive")
    ) %>%
    select(id = Nodes, label, measured) %>%
    mutate(
      type = ifelse(id %in% ksn_reg$source, "Kinase", "other"),
      type = ifelse(id %in% tf_reg$tf, "TF", type),
      group = paste(type, measured, sep = "_"),
      title = paste0(label, "_", group)
    ) %>%
    filter(id %in% c(Edges$from, Edges$to))
  
  Network <- visNetwork(
    Nodes, Edges,
    main = main, submain = submain, width = "100%"
    ) %>%
    visIgraphLayout(layout = "layout_with_kk", randomSeed = 1) %>%
    visGroups(groupname = "Kinase_active", color = "mediumseagreen") %>%
    visGroups(groupname = "Kinase_inactive", color = "palegreen") %>%
    visGroups(groupname = "TF_active", color = "#8b0a50") %>%
    visGroups(groupname = "TF_inactive", color = "palevioletred") %>%
    visGroups(groupname = "other_active", color = "burlywood") %>%
    visGroups(groupname = "other_inactive", color = "wheat") %>%
    visEdges(arrows = "to", color = "dimgrey", width = 2) %>%
    visNodes(font = "30px arial black bold") %>%
    visLegend(zoom = F) %>%
    visOptions(highlightNearest = T, nodesIdSelection = T)
  
  Cosmos.res <- list(
    Nodes = Nodes,
    Edges = Edges,
    Net = Network
  )
  return(Cosmos.res)
}

COSMOS.MakeClusters <- function(data.cosmos, db){
  Edges <- data.cosmos$Edges %>%
    select(from, to, sign) %>%
    distinct(from, to, .keep_all = T)
  res_graph <- graph_from_data_frame(Edges)
  cfg <- as.undirected(res_graph, mode = "collapse") %>%
    cluster_fast_greedy()
  clusters <- data.frame(id = cfg$names, Cluster = cfg$membership)
  Nodes <- data.cosmos$Nodes %>%
    select(
      id#, measured, type
      ) %>%
    left_join(clusters, by = "id") %>%
    mutate(
      Cluster = as.factor(Cluster)
    ) %>%
    filter(!is.na(Cluster))
  clusters.l <- list()
  for (i in 1:length(unique(Nodes$Cluster))) {
    clusters.l[[i]] <- Nodes %>%
      filter(Cluster == i) %>%
      mutate(
        EntrezID = mapIds(
          org.Hs.eg.db,
          keys = id,
          column = "ENTREZID",
          keytype = "SYMBOL", 
          multiVals = "first"
        )
      ) %>%
      filter(!is.na(EntrezID)) %>%
      pull(EntrezID) %>%
      unname()
  }
  names(clusters.l) <- 1:length(clusters.l)
  clusters.enrich <- compareCluster(
    geneClusters = clusters.l,
    fun = enricher,
    TERM2GENE = db
  ) %>%
    setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  Nodes <- Nodes %>%
    left_join(
      select(
        distinct(clusters.enrich@compareClusterResult, Cluster, .keep_all = T),
        Cluster, ID),
      by = "Cluster") %>%
    mutate(
      ID = ifelse(is.na(ID), paste0("Cluster ", Cluster),
                  sapply(strsplit(ID, split = "_"), function(x){
                    paste(x[2:length(x)], collapse = " ")
                  })
      ),
      title = id
    ) %>%
    dplyr::rename(group = ID)
  Cluster.res <- list(
    Nodes = Nodes,
    Edges = Edges,
    Enrichment = clusters.enrich
  )
  return(Cluster.res)
}

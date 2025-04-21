################################################################################
############ Synergistic Effect between DOX and Several Inhibitors #############
######################  of Chemoresistance Related Genes #######################
################################################################################

# Author: Borja Gallego Martínez -----------------------------------------------
# Sarcomas and Experimental Therapeutics lab (ISPA, Oviedo, Spain)

# Date: 20/01/2025 ----

# Info -------------------------------------------------------------------------

# This script performs the evaluation of the synergistic effect between DOX and
  # inhibitors of the most interesting genes extracted from MultiOmics analyses
  # of bone sarcoma DX-R models. The objective is to find drug combinations with
  # synergistic effect that could reverse the resistant phenotype.

################################################################################

# 0. Load Packages and Set-Up --------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(scipen = 999)
set.seed(123)

library(tidyverse)
library(synergyfinder)

################################################################################

# 1. Processing Raw Data ------------------------------------------------------

## 1.1. Load Synergy Results Tables ----

b143.data <- readxl::read_xlsx("./Data/Synergy_ResultsTable.xlsx",
                               sheet = 1) %>%
  mutate(
    ConcUnit = "uM",
    PairIndex = case_when(
      Drug1 == "Bortezomib" ~ 1,
      Drug1 == "IKK16" ~ 2,
      Drug1 == "Rapamycin" ~ 3,
      Drug1 == "Torin" ~ 4,
      Drug1 == "Abemaciclib" ~ 5,
      Drug1 == "Palbociclib" ~ 6,
      Drug1 == "CX4945" ~ 7,
      Drug1 == "JX06" ~ 8,
      Drug1 == "SGC" ~ 9,
      .default = NA
    ),
    PairIndex = paste(Sample, PairIndex, sep = "_")
  ) %>%
  select(!(Sample))

Saos2.data <- readxl::read_xlsx("./Data/Synergy_ResultsTable.xlsx",
                                sheet = 2) %>%
  mutate(
    ConcUnit = "uM",
    PairIndex = case_when(
      Drug1 == "Bortezomib" ~ 1,
      Drug1 == "IKK16" ~ 2,
      Drug1 == "Rapamycin" ~ 3,
      Drug1 == "Torin" ~ 4,
      Drug1 == "Abemaciclib" ~ 5,
      Drug1 == "Palbociclib" ~ 6,
      Drug1 == "CX4945" ~ 7,
      Drug1 == "JX06" ~ 8,
      .default = NA
    ),
    PairIndex = paste(Sample, PairIndex, sep = "_")
  ) %>%
  select(!(Sample))

## 1.2. Reshape Tables for Analysis ----

b143.res <- ReshapeData(
  data = b143.data,
  data_type = "viability",
  impute = F,
  noise = T,
  seed = 1
)

Saos2.res <- ReshapeData(
  data = Saos2.data,
  data_type = "viability",
  impute = F,
  noise = T,
  seed = 1
)

## 1.3. Calculate Synergy and Sensitivity Scores ----

b143.res <- CalculateSynergy(
  data = b143.res,
  correct_baseline = "non"
)
b143.res <- CalculateSensitivity(
  data = b143.res,
  correct_baseline = "non"
)

Saos2.res <- CalculateSynergy(
  data = Saos2.res,
  correct_baseline = "non"
)
Saos2.res <- CalculateSensitivity(
  data = Saos2.res,
  correct_baseline = "non"
)

# 2. Visualization of the Results ----------------------------------------------

## 2.1. Heatmaps of Dose-Response Matrices and Synergy Scores ----

Syn.PlotHeatmap <- function(data, plotval, sumstat, explot = F, path = NULL,
                            cellname, w = 10, h = 6, u = "in") {

  plot.list <- list()
  title <- c()
  
  for (i in 1:length(data$drug_pairs$block_id)) {
    if (plotval == "response") {
      title[i] <- paste0(data$drug_pairs$drug1[i], "-", data$drug_pairs$drug2[i],
                         " D/R in ", cellname, " ",
                         unlist(strsplit(data$drug_pairs$block_id[i],
                                         split = "_"))[1],
                         " Model")
    } else {
      title[i] <- paste0(data$drug_pairs$drug1[i], "-", data$drug_pairs$drug2[i],
                         " ", unlist(strsplit(plotval, "_"))[1],
                         " Synergy Score in ", cellname, " ",
                         unlist(strsplit(data$drug_pairs$block_id[i],
                                         split = "_"))[1],
                         " Model")
    }
    
    plot.list[[i]] <- Plot2DrugHeatmap(
      data = data,
      plot_block = data$drug_pairs$block_id[i],
      drugs = c(1, 2),
      plot_value = plotval,
      dynamic = F,
      summary_statistic = sumstat,
      plot_title = title[i],
      text_label_size_scale = 1.5
    ) +
      scale_fill_gradient2(
        low = "royalblue4", mid = "whitesmoke", high = "firebrick", midpoint = 0
      ) +
      labs(
        x = paste0(data$drug_pairs$drug1[i], " (µM)"),
        y = paste0(data$drug_pairs$drug2[i], " (µM)")
      ) +
      theme(
        axis.title.x = element_text(size = 16, face = "bold", hjust = .5,
                                    vjust = 2),
        axis.title.y = element_text(size = 16, face = "bold", hjust = .5,
                                    vjust = -1),
        axis.text.x = element_text(size = 14, vjust = 2, hjust = .5),
        axis.text.y = element_text(size = 14, hjust = .5, angle = 90,
                                   vjust = -1.5),
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 16),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, units = "cm"),
        legend.title = element_text(size = 16, face = "bold", hjust = .5),
        legend.text = element_text(size = 14, hjust = .5)
      )
    
    names(plot.list)[[i]] <- paste(
      unlist(strsplit(data$drug_pairs$block_id[i], split = "_"))[1],
      data$drug_pairs$drug1[i],
      sep = "_"
    )
  }
  
  if (explot == T) {
    for (i in 1: length(plot.list)) {
      ggsave(
        filename = paste0(path, cellname, "_", names(plot.list)[i], "_", plotval,
                          "heatmap.tiff"),
        plot = plot.list[[i]],
        device = "tiff", dpi = 300, width = w, height = h, units = u
      )
    }
  }
  
  return(plot.list)
}

b143.heatmaps_dr <- Syn.PlotHeatmap(
  data = b143.res,
  plotval = "response",
  sumstat = "mean",
  cellname = "143B",
  explot = T,
  path = "./Results/"
)
b143.heatmaps_bliss <- Syn.PlotHeatmap(
  data = b143.res,
  plotval = "Bliss_synergy",
  sumstat = c("mean", "quantile_100"),
  cellname = "143B",
  explot = T,
  path = "./Results/"
)

Saos2.heatmaps_dr <- Syn.PlotHeatmap(
  data = Saos2.res,
  plotval = "response",
  sumstat = "mean",
  cellname = "SaOS2",
  explot = F,
  path = "./Results/"
)
Saos2.heatmaps_bliss <- Syn.PlotHeatmap(
  data = Saos2.res,
  plotval = "Bliss_synergy",
  sumstat = c("mean", "quantile_100"),
  cellname = "SaOS2",
  explot = F,
  path = "./Results/"
)

## 2.2. Contour Plots of Synergy Scores ----

Syn.PlotContour <- function(data, plotval, sumstat, explot = F, path = NULL,
                            cellname, w = 10, h = 5, u = "in") {
  
  plot.list <- list()
  title <- c()
  
  for (i in 1:length(data$drug_pairs$block_id)) {
    if (plotval == "response") {
      title[i] <- paste0(data$drug_pairs$drug1[i], "-", data$drug_pairs$drug2[i],
                         " D/R in ", cellname, " ",
                         unlist(strsplit(data$drug_pairs$block_id[i],
                                         split = "_"))[1],
                         " Model")
    } else {
      title[i] <- paste0(data$drug_pairs$drug1[i], "-", data$drug_pairs$drug2[i],
                         " ", unlist(strsplit(plotval, "_"))[1],
                         " Synergy Score in ", cellname, " ",
                         unlist(strsplit(data$drug_pairs$block_id[i],
                                         split = "_"))[1],
                         " Model")
    }
    
    lim <- max(c(
      abs(data$synergy_scores %>%
            filter(block_id == data$drug_pairs$block_id[i]) %>%
            select(Bliss_synergy) %>%
            min()),
      abs(data$synergy_scores %>%
            filter(block_id == data$drug_pairs$block_id[i]) %>%
            select(Bliss_synergy) %>%
            max())
    ))
    lim <- lim + lim * .2
    
    plot.list[[i]] <- Plot2DrugContour(
      data = data,
      plot_block = data$drug_pairs$block_id[i],
      drugs = c(1, 2),
      plot_value = plotval,
      dynamic = F,
      summary_statistic = sumstat,
      plot_title = title[i],
      interpolate_len = 20
    ) +
      scale_fill_gradientn(
        colours = c("royalblue4", "white", "white", "white", "firebrick"),
        values = scales::rescale(c(-lim, -5, 0, 5, lim)),
        limits = c(-lim, lim)
      ) +
      labs(
        x = paste0(data$drug_pairs$drug1[i], " (µM)"),
        y = paste0(data$drug_pairs$drug2[i], " (µM)")
      ) +
      theme(
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold", hjust = .5),
        legend.text = element_text(size = 14, hjust = .5)
      )
    
    names(plot.list)[[i]] <- paste(
      unlist(strsplit(data$drug_pairs$block_id[i], split = "_"))[1],
      data$drug_pairs$drug1[i],
      sep = "_"
    )
  }
  
  if (explot == T) {
    for (i in 1: length(plot.list)) {
      ggsave(
        filename = paste0(path, cellname, "_", names(plot.list)[i], "_", plotval,
                          "contourplot.tiff"),
        plot = plot.list[[i]],
        device = "tiff", dpi = 300, width = w, height = h, units = u
      )
    }
  }
  
  return(plot.list)
}

b143.contour_bliss <- Syn.PlotContour(
  data = b143.res,
  plotval = "Bliss_synergy",
  sumstat = c("mean", "quantile_100"),
  cellname = "143B",
  explot = T,
  path = "./Results/"
)

Saos2.contour_bliss <- Syn.PlotContour(
  data = Saos2.res,
  plotval = "Bliss_synergy",
  sumstat = c("mean", "quantile_100"),
  cellname = "SaOS2",
  explot = F,
  path = "./Results/"
)

## 2.3. Synergy Scores - CSS Plots ----

b143.SS_bliss <- PlotSensitivitySynergy(
  data = b143.res,
  plot_synergy = "Bliss",
  show_labels = T,
  dynamic = F
)
ggsave("./Results/b143_SSPlotBliss.tiff", plot = b143.SS_bliss, device = "tiff",
       dpi = 300)

Saos2.SS_bliss <- PlotSensitivitySynergy(
  data = Saos2.res,
  plot_synergy = "Bliss",
  show_labels = T,
  dynamic = F
)

################################################################################

# 99. Save Results -------------------------------------------------------------

sessInfo <- sessionInfo()

save.image(
  file = paste0(
    "./Results/",
    paste(unlist(strsplit(as.character(Sys.Date()), split = "-")),
          collapse = ""),
    "_Workspace_SynRanalysis.RData"
  )
)

################################################################################

setwd("/projects/b1038/Pulmonary/cpuritz/PASC/data")
set.seed(0)

library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggsci)
library(reshape2)
source("../code/util/pairwise_wilcox_test.R")

## Read in data
counts <- read.csv("deidentified_data/deidentified_NEP_cell_counts.csv", check.names = FALSE) %>%
    dplyr::select(cell_type, Study_ID, cell_proportion, Status) %>%
    dplyr::mutate(cell_proportion = 100 * cell_proportion) %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = cell_proportion)

# Color palette
palette <- setNames(pal_npg("nrc")(2), c("RPRA", "Healthy"))

# Separate categorical and numerical columns
cat_cols <- c("Study_ID", "Status")
cell_types <- setdiff(colnames(counts), cat_cols)

## Pairwise Wilcoxon rank sum tests with FDR correction
all_comps <- boxplot_signif(counts, cell_types, group.var = "Status",
                            var.name = "cell_type", annot = "value")

df <- melt(counts, id.vars = c("Study_ID", "Status"), variable.name = "cell_type")
boxplots <- ggplot(df, aes(x = value, y = cell_type, fill = Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 0.5, position = position_jitterdodge(jitter.width = 0.1)) +
    scale_fill_manual(name = "", values = palette) +
    xlab("Percent of total cells") +
    ylab(NULL) +
    theme_bw(base_family = "Arial") +
    theme(text = element_text(family = "Arial"),
          axis.title = element_text(size = 11, color = "black"),
          axis.text = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          plot.margin = margin(0, 0, 0, 0),
          legend.position = "bottom")

ggsave(filename = "figures/fig_6/fig_6d.pdf", plot = boxplots,
      width = 8, height = 5, device = cairo_pdf)
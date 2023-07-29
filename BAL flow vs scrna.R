## Here we examine relative abundances of different cell types in BAL fluid,
## comparing flow cytometry with scRNA-seq.

setwd("/projects/b1038/Pulmonary/cpuritz/PASC/data")
set.seed(0)

library(ggplot2)
library(dplyr)

require(tidyr)
require(tibble)
require(patchwork)
require(reshape2)

## Read in scRNA-seq abundance data
scrna <- read.csv("deidentified_data/deidentified_BAL_cell_counts.csv",
                  check.names = FALSE) %>%
    dplyr::select(cell_type, Study_ID, cell_proportion) %>%
    dplyr::mutate(cell_proportion = 100 * cell_proportion)

## Read in flow cytometry data
# RPRA13 (donor) does not have a corresponding scRNA-seq library, omit it
flow <- read.csv("deidentified_data/deidentified_flow_cytometry_data.csv",
                 check.names = FALSE) %>%
    dplyr::filter(Study_ID != "RPRA13 (donor)") %>%
    dplyr::select(-Group) %>%
    tidyr::drop_na()

## Renormalize abundances (rounding errors)
data_cols <- setdiff(colnames(flow), "Study_ID")
flow[data_cols] <- t(apply(flow[data_cols], 1, function(x) { x / sum(x) * 100 }))

## Take the average of the two lungs of RPRA13
lr_ids <- c("RPRA13 (left)", "RPRA13 (right)")
l_ix <- dim(flow)[1] + 1
flow[l_ix, 1] <- "RPRA13"
flow[l_ix, 2:dim(flow)[2]] <- sapply(flow[flow$Study_ID %in% lr_ids, -1], mean)
flow <- dplyr::filter(flow, !Study_ID %in% lr_ids)

## Melt flow data frame
flow <- reshape2::melt(flow, id.vars = "Study_ID", variable.name = "cell_type",
                       value.name = "cell_proportion")

## Only keep samples with scRNA-seq libraries and flow cytometry samples
common_ids <- intersect(unique(flow$Study_ID), unique(scrna$Study_ID))
flow <- flow[flow$Study_ID %in% common_ids,]
scrna <- scrna[scrna$Study_ID %in% common_ids,]

## Map flow and scRNA-seq cell types to coarser cell type labels
cell_type_dict <- list(
    "CD4 T cells" = list(
        flow = c("CD4 T cells"),
        scrna = c("CD4 T cells-1", "CD4 T cells-2")
    ),
    "CD8 T cells" = list(
        flow = c("CD8 T cells"),
        scrna = c("CD8 T cells-1", "CD8 T cells-2", "CD8 T cells-3")
    ),
    "B and plasma cells" = list(
        flow = c("B and plasma cells"),
        scrna = c("B cells", "Plasma cells")
    ),
    "NK cells" = list(
        flow = c("NK cells"),
        scrna = c("gdT cells and NK cells")
    ),
    "Tregs" = list(
        flow = c("Tregs"),
        scrna = c("Tregs")
    ),
    "Macrophages and monocytes" = list(
        flow = c("CD206 low macrophages", "CD206 high macrophages"),
        scrna = c("MoAM-1", "MoAM-2", "MoAM-3", "MoAM-4", "TRAM-1", "TRAM-2",
                  "TRAM-3", "TRAM-4", "TRAM-5", "TRAM-6", "TRAM-7",
                  "Proliferating macrophages", "Perivascular macrophages",
                  "Monocytes-1", "Monocytes-2", "Migratory DC", "DC1", "DC2",
                  "pDC")
    )
)

## Generate plots
plots <- list()
# x,y locations of label
xt <- c(0, 0, 0, 1, 0, 0)
yt <- c(1, 1, 1, 1, 1, 1)
for (ct in names(cell_type_dict)) {
    coarse_flow <- c()
    coarse_scrna <- c()
    
    scrna_sub <- scrna[scrna$cell_type %in% cell_type_dict[[ct]]$scrna,]
    flow_sub <- flow[flow$cell_type %in% cell_type_dict[[ct]]$flow,]
    for (i in 1:length(common_ids)) {
        sid <- common_ids[i]
        coarse_flow[i] <- sum(flow_sub[flow_sub$Study_ID == sid, "cell_proportion"])
        coarse_scrna[i] <- sum(scrna_sub[scrna_sub$Study_ID == sid, "cell_proportion"])
    }
    
    df <- data.frame(x = coarse_flow, y = coarse_scrna)
    xymax <- max(c(max(df$x), max(df$y)))
    cl <- sprintf("%.2f", cor(df$x, df$y, method = "s"))
    plots[[length(plots) + 1]] <- ggplot(data = df, aes(x, y)) +
        geom_point() +
        xlab(NULL) +
        ylab(NULL) +
        xlim(c(0, xymax)) +
        ylim(c(0, xymax)) +
        ggtitle(ct) +
        annotate("text", hjust = 0,
                 x = 0.75 * xymax * xt[length(plots) + 1],
                 y = 0.97 * xymax * yt[length(plots) + 1],
                 label = paste0("rho==", deparse(cl)),
                 parse = TRUE) +
        theme_bw() +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              plot.title = element_text(hjust = 0.5))
}

# Add axis labels
for (i in c(1, 4)) {
    plots[[i]] <- plots[[i]] + ylab("scRNA-seq")
}
for (i in c(4, 5, 6)) {
    plots[[i]] <- plots[[i]] + xlab("Flow cytometry")
}

# Combine and save
pplt <- patchwork::wrap_plots(plots, ncol = 3)
ggsave(filename = "figures/fig_s3/fig_s3_d.pdf", plot = pplt,
       width = 9, height = 6, device = cairo_pdf)
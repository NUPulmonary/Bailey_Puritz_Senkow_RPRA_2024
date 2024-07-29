set.seed(0)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(reshape2)

# Read in scRNA-seq abundance data
scrna <- read.csv("deidentified_data/deidentified_BAL_cell_counts.csv",
                  check.names = FALSE) %>%
    dplyr::select(cell_type, Study_ID, cell_proportion = cell_proportion_scaled) %>%
    dplyr::mutate(cell_proportion = 100 * cell_proportion)

# Read in flow cytometry data
flow <- read.csv("deidentified_data/20230407_PASC_flow_V4_deidentified.csv",
                 check.names = FALSE) %>%
    dplyr::select(-Group) %>%
    drop_na()

# Take the average of the two lungs of RPRA13
lr_ids <- c("RPRA13 (left)", "RPRA13 (right)")
l_ix <- dim(flow)[1] + 1
flow[l_ix, 1] <- "RPRA13"
flow[l_ix, 2:dim(flow)[2]] <- sapply(flow[flow$Study_ID %in% lr_ids, -1], mean)
flow <- dplyr::filter(flow, !Study_ID %in% lr_ids)

flow <- melt(flow, id.vars = "Study_ID",
             variable.name = "cell_type",
             value.name = "cell_proportion")

# Only keep samples with scRNA-seq libraries and flow cytometry samples
common_ids <- intersect(unique(flow$Study_ID), unique(scrna$Study_ID))
flow <- flow[flow$Study_ID %in% common_ids,]
scrna <- scrna[scrna$Study_ID %in% common_ids,]

# Map flow and scRNA-seq cell types to coarser cell type labels
# scrna: exclude Proliferating T cells, pDC, Epithelial cells, Mast cells, SARS-CoV-2
# flow: exclude Neutrophils, C1 FSC_high, C3 HLADR_low
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
    "Dendritic cells" = list(
        flow = c("C4 HLADR_high"),
        scrna = c("DC1", "DC2", "Migratory DC")
    ),
    "Macrophages and monocytes" = list(
        flow = c("CD206 low macrophages", "CD206 high macrophages", "Monocytes"),
        scrna = c("MoAM-1", "MoAM-2", "MoAM-3", "MoAM-4",
                  "TRAM-1", "TRAM-2", "TRAM-3", "TRAM-4", "TRAM-5", "TRAM-6", "TRAM-7",
                  "Proliferating macrophages", "Perivascular macrophages",
                  "Monocytes-1", "Monocytes-2")
    )
)

# Generate plots
plots <- list()
for (ct in names(cell_type_dict)) {
    coarse_flow <- c()
    coarse_scrna <- c()
    
    scrna_sub <- scrna[scrna$cell_type %in% cell_type_dict[[ct]]$scrna, ]
    flow_sub <- flow[flow$cell_type %in% cell_type_dict[[ct]]$flow, ]
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
        xlab("Flow cytometry %") +
        ylab(NULL) +
        xlim(c(0, xymax)) +
        ylim(c(0, xymax)) +
        ggtitle(ct) +
        annotate("text", hjust = 0, x = xymax * 0.75, y = xymax * 0.05,
                 label = paste0("rho==", deparse(cl)), parse = TRUE) +
        theme_bw(base_family = "Arial") +
        theme(text = element_text(family = "Arial"),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              plot.title = element_text(hjust = 0.5))
}
plots[[1]] <- plots[[1]] + ylab("scRNA-seq %")

# Combine and save
ggsave(filename = "figures/fig_s3/fig_s3_c.pdf",
       plot = wrap_plots(plots, nrow = 1),
       width = 21, height = 3.5, device = cairo_pdf)

setwd("/projects/b1038/Pulmonary/cpuritz/PASC/data")
set.seed(1, kind = "L'Ecuyer-CMRG")

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(coin)
library(ggsci)
library(ggsignif)
library(grid)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)

source("../code/util/pairwise_wilcox_test.R")

## Read in data
flow <- read.csv("deidentified_data/deidentified_flow_cytometry_data.csv",
                 check.names = F)
flow <- na.omit(flow)
rownames(flow) <- NULL
colnames(flow) <- gsub("macrophages", "MPs", colnames(flow))

# Group separates RPRA from RPRA (transplant), Status does not
flow["Status"] <- gsub(" (transplant)", "", flow$Group, fixed = TRUE)

## Renormalize abundances (rounding errors)
data_cols <- which(!colnames(flow) %in% c("Study_ID", "Group", "Status"))
flow[data_cols] <- t(apply(flow[data_cols], 1, function(x) { x / sum(x) * 100 }))
cell_types <- colnames(flow[data_cols])

## Perform pairwise Wilcoxon rank sum tests with FDR correction
all_comps <- boxplot_signif(flow, cell_types, group.var = "Status",
                            var.name = "cell_type", annot = "value")

## Create boxplots
palette <- setNames(pal_npg("nrc")(3), c("RPRA", "Healthy", "RPRA (transplant)"))
signif_plots <- list()
nonsignif_plots <- list()
for (i in 1:length(cell_types)) {
    data <- flow %>%
        dplyr::select(dplyr::all_of(c("Status", cell_types[i]))) %>%
        dplyr::rename(cell_type = cell_types[i])
    comparisons <- dplyr::filter(all_comps, cell_type == cell_types[i])
    pplt <- ggplot(data, aes(x = Status, y = cell_type, fill = Status)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(shape = 16, position = position_jitter(width = 0.25)) +
        xlab(NULL) +
        ylab(NULL) +
        ggtitle(gsub("MPs", "macrophages", cell_types[i])) +
        scale_fill_manual(values = palette[1:2]) +
        theme_bw(base_family = "Arial") +
        theme(text = element_text(family = "Arial"),
              legend.position = "none",
              plot.title = element_text(size = 12, hjust = 0.5),
              axis.title = element_text(size = 11, color = "black"),
              axis.text = element_text(size = 10, color = "black"),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank())
    
    if (dim(comparisons)[1] > 0) {
        signif_plots[[length(signif_plots) + 1]] <- pplt +
            ylim(c(0, 1.06 * max(comparisons$y))) +
            geom_signif(xmin = comparisons$xmin,
                        xmax = comparisons$xmax,
                        y_position = 0.97 * comparisons$y,
                        annotation = comparisons$annot,
                        tip_length = 0,
                        vjust = -0.225,
                        textsize = 3.5)
    } else {
        nonsignif_plots[[length(nonsignif_plots) + 1]] <- pplt
    }
}

ylabel <- ggplot(data.frame(l = "Percent of total cells", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 90, size = 4) +
    theme_void() + coord_cartesian(clip = "off")
signif_bps <- ylabel + wrap_plots(signif_plots, nrow = 1) +
    plot_layout(widths = c(1, 30))
nonsignif_bps <- ylabel + wrap_plots(nonsignif_plots, ncol = 2) +
    plot_layout(widths = c(1, 30))


## Hierarchical clustering
flow2 <- column_to_rownames(flow, "Study_ID")
group_annot <- dplyr::select(flow2, Diagnosis = Group)
flow_mat <- t(as.matrix(dplyr::select(flow2, -Group, -Status)))

pheatmap2ggplot <- function (ph) {
    wrap_elements(grid.grabExpr(ComplexHeatmap::draw(ph)))
}
ph_pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(1001)

scale_rows <- function(x) {
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}
fmat_scaled <- scale_rows(flow_mat)

s <- ceiling(max(abs(min(fmat_scaled)), max(fmat_scaled)))
heatmap <- ComplexHeatmap::pheatmap(
    fmat_scaled,
    cluster_rows = T,
    cluster_cols = T,
    show_colnames = T,
    scale = "none",
    angle_col = '45',
    cellheight = 20,
    cellwidth = 20,
    clustering_method = "ward.D2",
    border_color = NA,
    color = ph_pal,
    breaks = seq(-s, s, length.out = length(ph_pal)),
    annotation_col = group_annot,
    annotation_color = list(Diagnosis = palette),
    name = "Percent",
    heatmap_legend_param = list(legend_height = unit(5, "cm"),
                                at = seq(-s, s, s/2))
)
heatmap <- pheatmap2ggplot(heatmap)


## Correlation with CT features
ct_features <- c("Normal", "Fibrotic", "Inflammatory", "Nodularity")
ct_data <- read.csv("deidentified_data/deidentified_ct_scan_1_data.csv",
                    check.names = FALSE)
ct_data <- ct_data %>% dplyr::rename(Study_ID = Subject_ID)
flow_ct <- flow %>%
    dplyr::inner_join(ct_data, by = "Study_ID") %>%
    tibble::column_to_rownames("Study_ID")

cor_mat <- matrix(nrow = length(cell_types),
                  ncol = length(ct_features),
                  dimnames = list(cell_types, ct_features))
p_mat <- matrix(nrow = length(cell_types),
                ncol = length(ct_features),
                dimnames = list(cell_types, ct_features))
for (i in 1:length(cell_types)) {
    for (j in 1:length(ct_features)) {
        x <- flow_ct[,cell_types[i]]
        y <- flow_ct[,ct_features[j]]
        cor_mat[i, j] <- cor(x, y, method = "spearman")
        test <- coin::spearman_test(x ~ y, data = data.frame(x, y),
                                    distribution = coin::approximate(nresample = 1e6,
                                                                     parallel = "multicore",
                                                                     ncpus = 16))
        p_mat[i, j] <- coin::pvalue(test)[1]
    }
}

cor_df <- cor_mat %>% melt() %>% dplyr::rename(Correlation = value)
p_df <- p_mat %>% melt() %>% dplyr::rename(pval = value)
p_adj <- p.adjust(as.vector(p_df$pval, mode = "numeric"), method = "fdr")

# Add significance codes
stars <- sapply(p_adj, function (x) {
    ifelse(x > 0.05, "", ifelse(x > 0.01, "*", ifelse(x > 0.001, "**", "***")))
})
annot <- paste0(sprintf("%.2f", as.vector(cor_mat, mode = "numeric")), stars)
annot_mat <- matrix(annot, nrow = dim(cor_mat)[1], dimnames = dimnames(p_mat))

cor_mat <- t(cor_mat[, rev(ct_features)])
annot_mat <- t(annot_mat[, rev(ct_features)])
ct_corr <- ComplexHeatmap::pheatmap(
    cor_mat,
    cluster_rows = T,
    cluster_cols = T,
    cellwidth = 20,
    cellheight = 20,
    angle_col = '45',
    clustering_method = "ward.D2",
    scale = "none",
    border_color = NA,
    color = ph_pal,
    breaks = seq(-1, 1, length.out = length(ph_pal)),
    name = "Correlation",
    heatmap_legend_param = list(direction = "vertical",
                                legend_height = unit(3, "cm"),
                                title_position = "topcenter",
                                at = seq(-1, 1, 0.5))
)
ct_corr <- pheatmap2ggplot(ct_corr)

## Point and line plots
plplots <- list()
for (i in 1:length(cell_types)) {
    for (j in 1:length(ct_features)) {
        cell <- cell_types[i]
        ct <- ct_features[j]
        fname <- gsub('-|/| ', '_', paste(cell, ct, sep = '_'))
        
        df <- data.frame(x = flow_ct[,ct], y = flow_ct[,cell], Group = flow_ct$Group)
        cr <- deparse(sprintf("%.3f", cor(df$x, df$y, method = "s")))
        pplt <- ggplot(df, aes(x, y)) +
            geom_point() +
            stat_smooth(method = "lm", formula = y ~ x) +
            xlab(paste(ct, "fraction of scan")) +
            ylab(NULL) +
            ggtitle(gsub("MPs", "macrophages", cell)) +
            annotate("text", x = 0.73 * max(df$x), y = 0.98 * max(df$y), hjust = 0,
                     label = paste0("rho == ", cr),
                     parse = TRUE) +
            theme_bw() +
            theme(text = element_text(family = "Arial"),
                  legend.position = "bottom",
                  plot.title = element_text(size = 12, hjust = 0.5),
                  axis.title = element_text(size = 11, color = "black"),
                  axis.text = element_text(size = 10, color = "black"),
                  panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank())
        plplots[[fname]] <- pplt
    }
}

### Build figures
get_annot <- function (x) {
    plot_annotation(title = x, theme = theme(plot.title = element_text(face = 2, size = 20)))
}

## Build Figure 2
cl_bar <- ggplot(data = data.frame(
        x = c(-0.05, -0.05, 2.57, 4.88, 7.2, 10),
        xend = c(10, -0.05, 2.57, 4.88, 7.2, 10),
        y = 0.1,
        yend = c(0.1, 0.14, 0.14, 0.14, 0.14, 0.14))
    ) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend)) +
    annotate("text", x = 1.3, y = 0.07, label = "Cluster 1") +
    annotate("text", x = 3.8, y = 0.07, label = "Cluster 2") +
    annotate("text", x = 6.0, y = 0.07, label = "Cluster 3") +
    annotate("text", x = 8.6, y = 0.07, label = "Cluster 4") +
    xlab(NULL) +
    ylab(NULL) +
    ylim(c(0, 1)) +
    theme_minimal() +
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
fig_2a <- wrap_elements(
    heatmap +
    inset_element(cl_bar, 0.025, -0.16, 0.809, 0.9) +
    get_annot("a")
)
fig_2b <- wrap_elements(signif_bps + get_annot("b"))
fig_2c <- wrap_elements(ct_corr + get_annot("c"))

fig2 <- fig_2a + fig_2b + fig_2c +
    plot_layout(
        design = c(patchwork::area(t = 1, b = 1, l = 1, r = 2),
                   patchwork::area(t = 2, b = 2, l = 1, r = 1),
                   patchwork::area(t = 2, b = 2, l = 2, r = 2)),
        heights = c(1, 0.7),
        widths = c(0.61, 0.39)
    )
ggsave(plot = fig2, filename = "figures/fig_2/fig2.pdf", width = 16, height = 9.5,
       device = cairo_pdf)


## Build Figure S2
fnames <- c("Neutrophils_Normal", "Neutrophils_Inflammatory",
            "Neutrophils_Nodularity", "CD4_T_cells_Fibrotic",
            "B_and_plasma_cells_Fibrotic", "CD206_high_MPs_Nodularity")
ylabel <- ggplot(data.frame(l = "Percent of total cells", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 90, size = 4) +
    theme_void() +
    coord_cartesian(clip = "off")
ct_plots <- ylabel +
    wrap_plots(lapply(fnames, function (x) { plplots[[x]] }), ncol = 2) +
    plot_layout(widths = c(1, 30))

fig_s2a <- wrap_elements(nonsignif_bps + get_annot("a"))
fig_s2b <- wrap_elements(ct_plots + get_annot("b"))
fig_s2 <- fig_s2a + fig_s2b + plot_layout(ncol = 2)

ggsave(plot = fig_s2, filename = "figures/fig_s2/fig_s2.pdf",
       width = 14, height = 10, device = cairo_pdf)
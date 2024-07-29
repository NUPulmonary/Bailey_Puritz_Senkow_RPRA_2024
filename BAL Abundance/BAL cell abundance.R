setwd("/projects/b1038/Pulmonary/cpuritz/PASC")
set.seed(1, kind = "L'Ecuyer-CMRG")

library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(ggdendro)
library(ggsignif)
library(ggsci)
library(patchwork)
library(grid)
library(coin)
library(ComplexHeatmap)
library(reshape2)
source("code/util/pairwise_wilcox_test.R")

# Read in data
counts <- read.csv("data/deidentified_data/deidentified_BAL_cell_counts.csv",
                   check.names = FALSE)
# Used neutrophil-adjusted counts
counts <- dplyr::select(counts, cell_type, Study_ID, Group = is_RPRA,
                        cell_proportion = cell_proportion_scaled, neutrophil_pct)

# Include neutrophil abundances
nphil_counts <- counts %>%
    dplyr::select(Study_ID, Group, cell_proportion = neutrophil_pct) %>%
    unique() %>%
    dplyr::mutate(cell_type = "Neutrophils") %>%
    dplyr::relocate(cell_type, .before = 1)
counts <- counts %>%
    dplyr::select(-neutrophil_pct) %>%
    rbind(nphil_counts) %>%
    dplyr::mutate(cell_proportion = 100 * cell_proportion,
                  Group = dplyr::recode(Group, False = "Healthy", True = "RPRA")) %>%
    pivot_wider(names_from = cell_type, values_from = cell_proportion)

# Separate categorical and numerical columns
cat_cols <- c("Study_ID", "Group")
cell_types_exclude <- c("Epithelial cells", "SARS-CoV-2", "Neutrophils")
cell_types <- setdiff(colnames(counts), c(cat_cols, cell_types_exclude))

# Useful ordering of cells
cell_order <- c(paste0("TRAM-", seq(7)), "Proliferating macrophages",
                paste0("MoAM-", seq(4)), "Perivascular macrophages",
                "Monocytes-1", "Monocytes-2", "CD4 T cells-1", "CD4 T cells-2",
                paste0("CD8 T cells-", seq(3)), "Tregs", "gdT cells and NK cells",
                "Proliferating T cells", "DC1", "DC2", "Migratory DC", "pDC",
                "Mast cells", "B cells", "Plasma cells")

# Pairwise Wilcoxon rank sum tests with FDR correction
all_comps <- boxplot_signif(counts, cell_types, group.var = "Group",
                            var.name = "cell_type", annot = "value")
bp_palette <- setNames(pal_npg("nrc")(2), c("RPRA", "Healthy"))
nonsignif_ct <- c()
signif_plots <- list()
for (i in 1:length(cell_types)) {
    comparisons <- dplyr::filter(all_comps, cell_type == cell_types[i])
    pplt <- ggplot(counts, aes(x = Group, y = !!sym(cell_types[i]), fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = bp_palette) +
        geom_jitter(shape = 16, height = 0) +
        xlab(NULL) +
        ylab(NULL) +
        ggtitle(cell_types[i]) +
        theme_bw(base_family = "Arial") +
        theme(text = element_text(family = "Arial"),
              legend.position = "none",
              axis.title = element_text(size = 14, color = "black"),
              axis.text = element_text(size = 14, color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 15))
    if (dim(comparisons)[1] > 0) {
        signif_plots[[length(signif_plots) + 1]] <- pplt +
            ylim(c(0, 1.07 * max(comparisons$y))) +
            geom_signif(xmin = comparisons$xmin,
                        xmax = comparisons$xmax,
                        y_position = comparisons$y,
                        annotation = comparisons$annot,
                        tip_length = 0,
                        vjust = -0.225,
                        textsize = 4)
    } else {
        nonsignif_ct[length(nonsignif_ct) + 1] <- cell_types[i]
    }
}

# Significant comparisons
ylabel <- ggplot(data.frame(l = "Percent of total cells", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 90, size = 5) +
    theme_void(base_family = "Arial") +
    theme(text = element_text(family = "Arial")) +
    coord_cartesian(clip = "off")
signif_boxplots <- ylabel +
    wrap_plots(signif_plots, nrow = 1) +
    plot_layout(widths = c(0.75, 30))
ggsave(plot = signif_boxplots,
       filename = "data/figures/fig_3/fig_3e.pdf",
       width = 16, height = 4, device = cairo_pdf)

# All comparisons
df <- melt(counts, id.vars = c("Study_ID", "Group"), variable.name = "cell_type")
df$cell_type <- factor(gsub("macrophages", "MPs", df$cell_type),
                       levels = gsub("macrophages", "MPs", c(cell_order, "Epithelial cells", "SARS-CoV-2")),
                       ordered = TRUE)
df <- df[!is.na(df$cell_type), ]

fig_s3h <- ggplot(df, aes(y = value, x = cell_type, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 0.5, position = position_jitterdodge(jitter.height = 0.1)) +
    scale_fill_manual(name = "", values = bp_palette) +
    ylab("Percent of total cells") +
    xlab(NULL) +
    theme_bw(base_family = "Arial") +
    theme(text = element_text(family = "Arial"),
          axis.title = element_text(size = 13, color = "black"),
          axis.text.x = element_text(size = 13, color = "black", angle = 45,
                                     vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          legend.position = "right")
ggsave(plot = fig_s3h,
       filename = "data/figures/fig_s3/fig_s3_h.pdf",
       width = 18, height = 6, device = cairo_pdf)

############################################################

## Load CT data
ct_features <- c("Normal", "Fibrotic", "Inflammatory", "Nodularity")
ct_data <- read.csv("data/deidentified_data/deidentified_ct_scan_1_data.csv",
                    check.names = FALSE) %>%
    dplyr::rename(Study_ID = Subject_ID) %>%
    dplyr::select(-"Emphysema/Cysts")

#### Heatmap ####
zcounts <- counts
zcounts <- dplyr::left_join(zcounts, ct_data, by = "Study_ID")
zcounts$Group[zcounts$Study_ID %in% c("RPRA12", "RPRA13")] <- "RPRA (transplant)"

# Flow cluster assignments
flow_clust <- c(HV01 = "NA", HV02 = 3, HV03 = 1, HV04 = 1, HV05 = 3, HV06 = 1,
                RPRA02 = 4, RPRA03 = 1, RPRA04 = 2, RPRA05 = 3, RPRA06 = 3,
                RPRA07 = 2, RPRA08 = 3, RPRA09 = 1, RPRA10 = 3, RPRA11 = 3,
                RPRA12 = 4, RPRA13 = 4, RPRA18 = 4, RPRA19 = 3, RPRA20 = 3,
                RPRA21 = 2, RPRA22 = 2, RPRA23 = 2, RPRA24 = 2, RPRA26 = 5,
                RPRA27 = 5, RPRA28 = 2, RPRA29 = 2, RPRA30 = 5, RPRA32 = 3,
                RPRA33 = 5, RPRA34 = 5, RPRA35 = 2)
flow_order <- c("RPRA02", "RPRA13", "RPRA12", "RPRA26", "RPRA27", "RPRA30",
                "RPRA33", "RPRA03", "RPRA09", "HV03", "HV04", "HV06", "RPRA23",
                "RPRA29", "RPRA07", "RPRA24", "RPRA28", "RPRA21", "RPRA04",
                "RPRA22", "RPRA08", "HV05", "HV02", "RPRA11", "RPRA32",
                "RPRA05", "RPRA20", "RPRA06", "RPRA19")
flow_clust <- rownames_to_column(as.data.frame(flow_clust), "Study_ID")
zcounts <- dplyr::left_join(zcounts, flow_clust, by = "Study_ID")

ggpheatmap <- function(par, ...) {
    ph <- do.call(ComplexHeatmap::pheatmap, par)
    wrap_elements(grid.grabExpr(ComplexHeatmap::draw(ph, ...)))
}

fcounts <- zcounts[zcounts$flow_clust != "NA", ]
mat <- t(as.matrix(fcounts[, c(cell_types, "Neutrophils")]))
mat <- t(apply(mat, 1, function (x) { (x - mean(x)) / sd(x) }))
colnames(mat) <- fcounts$Study_ID

annot <- fcounts %>%
    dplyr::select(Study_ID, Group, "Flow cluster" = flow_clust) %>%
    column_to_rownames("Study_ID")
grp_palette <- setNames(pal_npg("nrc")(3)[c(1, 3, 2)],
                        c("RPRA", "RPRA (transplant)", "Healthy"))
flow_palette <- setNames(pal_npg("nrc")(8)[4:8], as.character(seq(5)))
annot_color <- list(Group = grp_palette, "Flow cluster" = flow_palette)
ph_pal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(1001)
ph_breaks <- seq(-6, 6, length.out = 1001)
ph_legend <- list(legend_height = unit(5, "cm"), at = seq(-6, 6, 3))

ph_par <- list(cluster_rows = T, show_colnames = T, scale = "none",
               angle_col = '45', cellheight = 20, cellwidth = 20,
               clustering_method = "ward.D2", border_color = NA,
               color = ph_pal, breaks = ph_breaks,
               annotation_color = annot_color, name = "Row\nz-score",
               heatmap_legend_param = ph_legend)

# Cluster subjects
hmap1 <- ggpheatmap(c(list(mat = mat, cluster_cols = T, annotation_col = annot), ph_par),
                    merge_legend = T)
ggsave(plot = hmap1,
       filename = "data/figures/fig_s3/fig_s3_f.pdf",
       width = 12.25, height = 10.5, device = cairo_pdf)

# Use subject clustering from flow
hmap2 <- ggpheatmap(c(list(mat = mat[, flow_order], cluster_cols = F, annotation_col = annot[flow_order, ]), ph_par),
                    merge_legend = T)
ggsave(plot = hmap2,
       filename = "data/figures/fig_s3/fig_s3_g.pdf",
       width = 12.1, height = 9.75, device = cairo_pdf)

############################################################

#### Correlation with CT features ####
counts_ct <- counts %>%
    dplyr::inner_join(ct_data, by = "Study_ID") %>%
    column_to_rownames("Study_ID")

cor_mat <- matrix(nrow = length(cell_types),
                  ncol = length(ct_features),
                  dimnames = list(cell_types, ct_features))
p_mat <- matrix(nrow = length(cell_types),
                ncol = length(ct_features),
                dimnames = list(cell_types, ct_features))
cdst <- coin::approximate(nresample = 1e6,
                          parallel = "multicore",
                          ncpus = 16L)
for (i in 1:length(cell_types)) {
    for (j in 1:length(ct_features)) {
        x <- counts_ct[, cell_types[i]]
        y <- counts_ct[, ct_features[j]]
        cor_mat[i, j] <- cor(x, y, method = "spearman")
        test <- coin::spearman_test(x ~ y, data = data.frame(x, y),
                                    distribution = cdst)
        p_mat[i, j] <- as.numeric(coin::pvalue(test))
    }
}

cor_df <- dplyr::rename(melt(cor_mat), Correlation = value)
p_df <- dplyr::rename(melt(p_mat), pval = value)
p_adj <- p.adjust(as.vector(p_df$pval, mode = "numeric"), method = "fdr")
cor_df_full <- cor_df %>%
    dplyr::full_join(p_df, by = c("Var1", "Var2")) %>%
    dplyr::mutate(padj = p_adj)

# Add annotations
cor_mat_vec <- as.vector(cor_mat, mode = "numeric")
annot <- sapply(seq(length(p_adj)), function (i) {
    if (p_adj[i] < 0.05) sprintf("%.2f", cor_mat_vec[i]) else ""
})
annot_mat <- matrix(annot, nrow = dim(cor_mat)[1], dimnames = dimnames(p_mat))

# Point and line plots for correlations of interest
signif_cells <- as.vector(p_df[which(annot != ""), ][, 1])
signif_ct <- as.vector(p_df[which(annot != ""), ][, 2])
signif_pval <- as.vector(p_adj[which(annot != "")])
cc_plots <- list()
for (i in 1:length(signif_cells)) {
    df <- data.frame(x = counts_ct[, signif_ct[i]],
                     y = counts_ct[, signif_cells[i]])
    cr <- deparse(sprintf("%.3f", cor(df$x, df$y, method = "s")))
    cc_plots[[i]] <- ggplot(df, aes(x, y)) +
        geom_point() +
        stat_smooth(method = "lm", formula = y ~ x) +
        annotate("text", x = min(df$x), y = 0.94 * max(df$y), hjust = 0,
                 label = sprintf("q = %.3f", signif_pval[i])) +
        annotate("text", x = min(df$x), y = max(df$y), hjust = 0,
                 label = paste0("rho == ", cr),
                 parse = TRUE) +
        xlab(paste(signif_ct[i], "fraction of lung")) +
        ylab("Percent of total cells") +
        ggtitle(signif_cells[i]) +
        theme_bw(base_family = "Arial") +
        theme(text = element_text(family = "Arial"),
              legend.position = "none",
              axis.title = element_text(size = 13, color = "black"),
              axis.text = element_text(size = 12, color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5))
}
ggsave(filename = "data/figures/fig_3/fig_3g.pdf",
       plot = wrap_plots(cc_plots, nrow = 2),
       width = 4, height = 7.6, bg = "transparent", device = cairo_pdf)

# Order cells
cor_mat <- t(cor_mat[cell_order, rev(ct_features)])
annot_mat <- t(annot_mat[cell_order, rev(ct_features)])

# Plot CT correlation matrix
ct_corr <- ggpheatmap(list(
    mat = cor_mat, cluster_rows = T, cluster_cols = T, show_colnames = T,
    angle_col = '45', cellheight = 22, cellwidth = 24,
    clustering_method = "ward.D2", scale = "none",
    breaks = seq(-1, 1, length.out = length(ph_pal)), border_color = NA,
    color = ph_pal, display_numbers = annot_mat, number_color = "black",
    name = "Correlation",
    heatmap_legend_param = list(legend_height = unit(3, "cm"))
))
ggsave(filename = "data/figures/fig_3/fig_3f.pdf",
       plot = ct_corr,
       width = 13, height = 3.75, device = cairo_pdf)

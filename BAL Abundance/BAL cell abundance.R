setwd("/projects/b1038/Pulmonary/cpuritz/PASC/data")
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
source("../code/util/pairwise_wilcox_test.R")

## Read in data
counts <- read.csv("deidentified_data/deidentified_BAL_cell_counts.csv",
                   check.names = FALSE) %>%
    dplyr::select(cell_type, Study_ID, Group = is_RPRA,
                  cell_proportion = cell_proportion_scaled) %>%
    dplyr::mutate(cell_proportion = 100 * cell_proportion,
                  Group = dplyr::recode(Group, False = "Healthy", True = "RPRA")) %>%
    dplyr::filter(!cell_type %in% c("Epithelial cells", "SARS-CoV-2")) %>%
    pivot_wider(names_from = cell_type, values_from = cell_proportion)

# Color palette
palette <- setNames(pal_npg("nrc")(2), c("RPRA", "Healthy"))

# Separate categorical and numerical columns
cat_cols <- c("Study_ID", "Group")
cell_types <- setdiff(colnames(counts), cat_cols)

# Useful ordering of cells
cell_order <- c(paste0("TRAM-", seq(7)), "Proliferating macrophages",
                paste0("MoAM-", seq(4)), "Perivascular macrophages",
                "Monocytes-1", "Monocytes-2", "CD4 T cells-1", "CD4 T cells-2",
                paste0("CD8 T cells-", seq(3)), "Tregs", "gdT cells and NK cells",
                "Proliferating T cells", "DC1", "DC2", "Migratory DC", "pDC",
                "Mast cells", "B cells", "Plasma cells")

### Pairwise Wilcoxon rank sum tests with FDR correction ###
all_comps <- boxplot_signif(counts, cell_types, group.var = "Group",
                            var.name = "cell_type", annot = "value")
nonsignif_ct <- c()
signif_plots <- list()
for (i in 1:length(cell_types)) {
    comparisons <- dplyr::filter(all_comps, cell_type == cell_types[i])
    pplt <- ggplot(counts, aes(x = Group, y = !!sym(cell_types[i]), fill = Group)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = palette) +
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
ggsave(plot = signif_boxplots, filename = "figures/fig_3/fig_3e.pdf",
       width = 16, height = 4, device = cairo_pdf)

# Non-significant comparisons
df <- melt(counts, id.vars = c("Study_ID", "Group"), variable.name = "cell_type")
df$cell_type <- factor(df$cell_type, levels = rev(cell_order), ordered = TRUE)
df <- df[df$cell_type %in% nonsignif_ct,]

nonsignif_boxplots <- ggplot(df, aes(x = value, y = cell_type, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size = 0.5, position = position_jitterdodge(jitter.width = 0.1)) +
    scale_fill_manual(name = "", values = palette) +
    xlab("Percent of total cells") +
    ylab(NULL) +
    theme_bw(base_family = "Arial") +
    theme(text = element_text(family = "Arial"),
          axis.title = element_text(size = 13, color = "black"),
          axis.text = element_text(size = 13, color = "black"),
          legend.text = element_text(size = 13, color = "black"),
          plot.margin = margin(0, 0, 0, 0),
          legend.position = "bottom")
ggsave(plot = nonsignif_boxplots, filename = "figures/fig_s4/fig_s4_b.pdf",
       width = 13, height = 12, device = cairo_pdf)

############################################################

## Add in CT data
ct_features <- c("Normal", "Fibrotic", "Inflammatory", "Nodularity")
ct_data <- read.csv("deidentified_data/deidentified_ct_scan_1_data.csv",
                    check.names = FALSE) %>%
    dplyr::rename(Study_ID = Subject_ID) %>%
    dplyr::select(-"Emphysema/Cysts")
counts_ct <- counts %>%
    dplyr::inner_join(ct_data, by = "Study_ID") %>%
    column_to_rownames("Study_ID")

#### Heatmap ####
# z-score abundances
zcounts <- counts
zcounts[cell_types] <- lapply(zcounts[cell_types],
                              function (x) { (x - mean(x)) / sd(x) })
zcounts_mat <- t(as.matrix(dplyr::select(zcounts, -Study_ID, -Group)))
cell_hclust <- hclust(dist(zcounts_mat), "ward.D2")

# Add in CT data
zcounts_ct <- dplyr::left_join(zcounts, ct_data, by = "Study_ID")
zcounts_ct$Group[zcounts_ct$Study_ID %in% c("RPRA12", "RPRA13")] <- "RPRA (transplant)"
df <- melt(zcounts_ct, id.vars = c("Study_ID", "Group", "Normal"))
df$Group <- factor(df$Group)
df <- dplyr::filter(df, variable %in% cell_types)

# Order data
df[is.na(df$Normal), "Normal"] <- -1
df <- dplyr::arrange(df, Normal)
df[df$Normal == -1, "Normal"] <- NA
df$Study_ID <- factor(df$Study_ID, levels = unique(df$Study_ID), ordered = TRUE)
df$variable <- factor(df$variable,
                      levels = rev(rownames(zcounts_mat)[cell_hclust$order]),
                      ordered = TRUE)

# Heatmap
heatmap <-
    ggplot(data = df, aes(x = Study_ID, y = variable)) +
    geom_tile(aes(fill = value)) +
    xlab(NULL) +
    ylab(NULL) +
    scale_fill_gradient2("Abundance", low = "blue", high = "red") +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 0.95, color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.ticks = element_blank(),
          plot.margin = margin(0, 0, 0, 0))

# Diagnosis annotation
palette2 <- setNames(pal_npg("nrc")(3)[c(1, 3, 2)], 
                     c("RPRA", "RPRA (transplant)", "Healthy"))
heatmap_annot1 <-
    ggplot(data = df, aes(x = Study_ID, y = "Diagnosis")) +
    geom_tile(aes(fill = Group)) +
    xlab(NULL) +
    ylab(NULL) +
    scale_fill_manual("Diagnosis", values = palette2) +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold", color = "black"),
          axis.ticks = element_blank(),
          plot.margin = margin(0, 0, 0, 0))

# Normal lung fraction annotation
heatmap_annot2 <-
    ggplot(data = df, aes(x = Study_ID, y = "Normal")) +
    geom_tile(aes(fill = Normal)) +
    scale_fill_gradient("Fraction of lung", na.value = "white", limits = c(0, 1)) +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    xlab(NULL) +
    ylab(NULL) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold", color = "black"),
          axis.ticks = element_blank(),
          plot.margin = margin(0, 0, 0, 0))

# Dendrogram
dend_data <- dendro_data(as.dendrogram(cell_hclust))
segment_data <- dplyr::rename(segment(dend_data), x = y, y = x, xend = yend, yend = xend)
dend_df <- data.frame(y_center = dend_data$labels$x,
                      gene = as.character(dend_data$labels$label),
                      height = 1)
dendro <-
    ggplot(segment_data) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_reverse(expand = c(0, 0.1, 0, 0)) +
    scale_y_reverse(breaks = dend_df$y_center,
                    expand = c(0, 0),
                    limits = c(dim(dend_df)[1], 0) + 0.5) +
    xlab(NULL) +
    ylab(NULL) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(0, 0, 0, 0),
          axis.text = element_blank(),
          axis.ticks = element_blank())

# Put everything together
heatmap_all <-
    heatmap_annot1 + plot_spacer() + heatmap_annot2 +
    plot_spacer() + dendro + plot_spacer() + heatmap +
    plot_layout(
        design = c(patchwork::area(t = 1, b = 1, l = 3, r = 3),
                   patchwork::area(t = 2, b = 2, l = 3, r = 3),
                   patchwork::area(t = 3, b = 3, l = 3, r = 3),
                   patchwork::area(t = 4, b = 4, l = 3, r = 3),
                   patchwork::area(t = 5, b = 5, l = 1, r = 1),
                   patchwork::area(t = 1, b = 5, l = 2, r = 2),
                   patchwork::area(t = 5, b = 5, l = 3, r = 3)),
        widths = c(1, -0.15, 3.5),
        heights = c(1, -1.75, 1, -1.75, 45),
        guides = "collect"
    ) & theme(legend.position = "right")
ggsave(filename = "figures/fig_s4/fig_s4_a.pdf",
       plot = heatmap_all, width = 12, height = 6.5)

############################################################

#### Correlation with CT features ####
cor_mat <- matrix(nrow = length(cell_types),
                  ncol = length(ct_features),
                  dimnames = list(cell_types, ct_features))
p_mat <- matrix(nrow = length(cell_types),
                ncol = length(ct_features),
                dimnames = list(cell_types, ct_features))
for (i in 1:length(cell_types)) {
    for (j in 1:length(ct_features)) {
        x <- counts_ct[,cell_types[i]]
        y <- counts_ct[,ct_features[j]]
        cor_mat[i, j] <- cor(x, y, method = "spearman")
        test <- coin::spearman_test(
            x ~ y, data = data.frame(x, y),
            distribution = coin::approximate(nresample = 1e6,
                                             parallel = "multicore",
                                             ncpus = 16L)
        )
        p_mat[i, j] <- coin::pvalue(test)[1]
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
    ifelse(p_adj[i] < 0.05, sprintf("%.2f", cor_mat_vec[i]), "")
})
annot_mat <- matrix(annot, nrow = dim(cor_mat)[1], dimnames = dimnames(p_mat))

# Point and line plots for correlations of interest
signif_cells <- as.vector(p_df[which(annot != ""),][,1])
signif_ct <- as.vector(p_df[which(annot != ""),][,2])
signif_pval <- as.vector(p_adj[which(annot != "")])
cc_plots <- list()
for (i in 1:length(signif_cells)) {
    cell <- signif_cells[i]
    ct <- signif_ct[i]
    df <- data.frame(x = counts_ct[,ct], y = counts_ct[,cell])
    cr <- deparse(sprintf("%.3f", cor(df$x, df$y, method = "s")))
    cc_plots[[i]] <- ggplot(df, aes(x, y)) +
        geom_point() +
        stat_smooth(method = "lm", formula = y ~ x) +
        annotate("text", x = min(df$x), y = 0.94 * max(df$y), hjust = 0,
                 label = sprintf("q = %.3f", signif_pval[i])) +
        annotate("text", x = min(df$x), y = max(df$y), hjust = 0,
                 label = paste0("rho == ", cr),
                 parse = TRUE) +
        xlab(paste(ct, "fraction of lung")) +
        ylab("Percent of total cells") +
        ggtitle(cell) +
        theme_bw(base_family = "Arial") +
        theme(text = element_text(family = "Arial"),
              legend.position = "none",
              axis.title = element_text(size = 13, color = "black"),
              axis.text = element_text(size = 12, color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5))
}
ggsave(filename = "figures/fig_3/fig_3g.pdf", plot = wrap_plots(cc_plots, nrow = 2),
       width = 4, height = 7.6, device = cairo_pdf, bg = "transparent")

# Order cells
cor_mat <- t(cor_mat[cell_order, rev(ct_features)])
annot_mat <- t(annot_mat[cell_order, rev(ct_features)])

# Plot CT correlation matrix
ct_corr <- ComplexHeatmap::pheatmap(
    cor_mat,
    cluster_rows = T,
    cluster_cols = T,
    show_colnames = T,
    angle_col = '45',
    cellheight = 22,
    cellwidth = 24,
    clustering_method = "ward.D2",
    scale = "none",
    breaks = seq(-1, 1, 2e-3),
    border_color = NA,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(1001),
    display_numbers = annot_mat,
    number_color = "black",
    name = "Correlation",
    heatmap_legend_param = list(legend_height = unit(3, "cm"))
)
ct_corr <- wrap_elements(grid.grabExpr(ComplexHeatmap::draw(ct_corr)))
ggsave(filename = "figures/fig_3/fig_3f.pdf", plot = ct_corr,
       width = 13, height = 3.75, device = cairo_pdf)
setwd("/projects/b1038/Pulmonary/cpuritz/PASC/data")
set.seed(1, kind = "L'Ecuyer-CMRG")

library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

require(coin)
require(ggsci)
require(ggsignif)
require(patchwork)
require(tibble)
require(ComplexHeatmap)

source("../code/util/pairwise_wilcox_test.R")

## Read in data
cyto <- read.csv(
    file = "deidentified_data/deidentified_cytokine_data.csv",
    check.names = FALSE
)
palette <- setNames(ggsci::pal_npg("nrc")(3),
                    c("RPRA", "Healthy", "RPRA (transplant)"))

## Fix analyte names
cyto[cyto$analyte == "Eotaxin", "analyte"] <- "Eotaxin-1"

## Exclude serum samples, rename diagnoses, select columns of interest
cyto <- cyto %>%
    dplyr::filter(sample_origin == "BAL") %>%
    dplyr::select(-sample_origin, -display_name) %>%
    dplyr::mutate(Group = diagnosis) %>%
    tidyr::pivot_wider(names_from = analyte, values_from = mean_concentration)
cyto$Group[cyto$Subject_ID %in% c("RPRA12", "RPRA13")] <- "RPRA (transplant)"

## Remove analytes that weren't measured in any sample
cat_cols <- c("Subject_ID", "diagnosis", "Group")
num_cols <- which(!colnames(cyto) %in% cat_cols)
cyto <- cyto[,c(cat_cols, names(which(colSums(cyto[num_cols]) != 0)))]
analytes <- setdiff(colnames(cyto), cat_cols)

## Perform pairwise Wilcoxon rank sum tests with FDR correction
all_comps <- boxplot_signif(cyto, analytes, group.var = "diagnosis",
                            var.name = "analyte", annot = "value")

## Use Greek letters in cytokine names where necessary
# Used for boxplot titles and CT correlation heatmap
convert_greek <- function(v) {
    l0 <- c("MIP-1a", "MIP-1b", "MIP-1d", "IL-1a", "IL-1b", "TNFa", "TNFb")
    l1 <- expression(paste("MIP-1", alpha, sep = ''),
                     paste("MIP-1", beta, sep = ''),
                     paste("MIP-1", delta, sep = ''),
                     paste("IL-1", alpha, sep = ''),
                     paste("IL-1", beta, sep = ''),
                     paste("TNF", alpha, sep = ''),
                     paste("TNF", beta, sep = ''))
    v[sapply(l0, function (x) { which(v == x) })] <- l1
    return(v)
}
# Used for deconvolution heatmaps
sub_greek <- function(v) {
    sapply(v, function (x) {
        g <- unlist(strsplit(x, "\\(|\\)"))
        gn <- substr(g[2], 1, nchar(g[2]) - 1)
        if (gn %in% c("MIP-1", "IL-1", "TNF")) {
            gl <- substr(g[2], nchar(g[2]), nchar(g[2]))
            g1 <- paste(g[1], "(", gn, sep = '')
            if (gl == 'a') {
                return(as.expression(bquote(.(g1) * alpha * .(")"))))
            } else if (gl == 'b') {
                return(as.expression(bquote(.(g1) * beta * .(")"))))
            } else if (gl == 'd') {
                return(as.expression(bquote(.(g1) * delta * .(")"))))
            }
        }
        return(x)
    })
}

## Create boxplots
signif_plots <- list()
nonsignif_plots <- list()
for (i in 1:length(analytes)) {
    data <- cyto %>%
        dplyr::select(dplyr::all_of(c("diagnosis", analytes[i]))) %>%
        dplyr::rename(analyte = analytes[i])
    comparisons <-  dplyr::filter(all_comps, analyte == analytes[i])
    pplt <- ggplot(data, aes(x = diagnosis, y = analyte, fill = diagnosis)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = palette) +
        geom_jitter(shape = 16, height = 0) +
        xlab(NULL) +
        ylab(NULL) +
        ggtitle(convert_greek(analytes)[i]) +
        theme_bw(base_family = "Arial") +
        theme(text = element_text(family = "Arial"),
              legend.position = "none",
              axis.title = element_text(size = 13, color = "black"),
              axis.text = element_text(size = 12, color = "black"),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              plot.title = element_text(hjust = 0.5))
    if (dim(comparisons)[1] > 0) {
        signif_plots[[length(signif_plots) + 1]] <-
            pplt +
            ylim(c(0, max(comparisons$y))) +
            ggsignif::geom_signif(xmin = comparisons$xmin,
                                  xmax = comparisons$xmax,
                                  y_position = 0.925 * comparisons$y,
                                  annotation = comparisons$annot,
                                  tip_length = 0,
                                  vjust = -0.225)
    } else {
        nonsignif_plots[[length(nonsignif_plots) + 1]] <- pplt
    }
}

ylabel <- ggplot(data.frame(l = "Mean concentration", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 90, size = 4.5) +
    theme_void(base_family = "Arial") +
    theme(text = element_text(family = "Arial")) +
    coord_cartesian(clip = "off")
signif_boxplots <- ylabel +
    patchwork::wrap_plots(signif_plots, ncol = 4) +
    patchwork::plot_layout(widths = c(1, 30))
nonsignif_boxplots <- ylabel +
    patchwork::wrap_plots(nonsignif_plots, ncol = 7) +
    patchwork::plot_layout(widths = c(1, 30))


## Hierarchical clustering
cyto_sub <- tibble::column_to_rownames(cyto, "Subject_ID")
group_annot <- cyto_sub %>% dplyr::select(Group) %>% dplyr::rename(Diagnosis = Group)
cyto_mat <- cyto_sub %>% dplyr::select(-diagnosis, -Group) %>% as.matrix() %>% t()

## z-score matrix rows
scale_rows <- function(x) {
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}
cyto_mat_scaled <- scale_rows(cyto_mat)

legend_param <- function(s, l) {
    list(legend_height = unit(l, "cm"), at = seq(-s, s, s/2))
}

s <- ceiling(max(abs(min(cyto_mat_scaled)), max(cyto_mat_scaled)))
ph_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(1001)

heatmap <- ComplexHeatmap::pheatmap(
    scale_rows(cyto_mat),
    cluster_rows = T,
    cluster_cols = T,
    angle_col = '45',
    cellheight = 16,
    cellwidth = 18,
    clustering_method = "ward.D2",
    scale = "none",
    border_color = NA,
    color = ph_pal,
    breaks = seq(-s, s, length.out = length(ph_pal)),
    labels_row = convert_greek(rownames(cyto_mat)),
    annotation_col = group_annot,
    annotation_color = list(Diagnosis = palette),
    name = "Mean\nconcentration",
    heatmap_legend_param = legend_param(s, 4)
)
heatmap <- patchwork::wrap_elements(
    grid::grid.grabExpr(ComplexHeatmap::draw(heatmap, merge_legend = T))
)

## Correlation with CT features
ct_features <- c("Normal", "Fibrotic", "Inflammatory", "Nodularity")
ct_data <- read.csv("deidentified_data/deidentified_ct_scan_1_data.csv",
                    check.names = FALSE)

cyto_ct <- dplyr::inner_join(cyto, ct_data, by = "Subject_ID") %>%
    dplyr::select(-diagnosis, -Group) %>%
    tibble::column_to_rownames("Subject_ID")
analytes <- setdiff(colnames(cyto_ct), names(ct_data))

cor_mat <- matrix(nrow = length(analytes),
                  ncol = length(ct_features),
                  dimnames = list(analytes, ct_features))
p_mat <- matrix(nrow = length(analytes),
                ncol = length(ct_features),
                dimnames = list(analytes, ct_features))
for (i in 1:length(analytes)) {
    for (j in 1:length(ct_features)) {
        x <- cyto_ct[,analytes[i]]
        y <- cyto_ct[,ct_features[j]]
        cor_mat[i, j] <- cor(x, y, method = "spearman")
        test <- coin::spearman_test(
            x ~ y,
            data = data.frame(x, y),
            distribution = coin::approximate(nresample = 1e6,
                                             parallel = "multicore",
                                             ncpus = 16)
        )
        p_mat[i, j] <- coin::pvalue(test)[1]
    }
}

cor_df <- reshape2::melt(cor_mat) %>% dplyr::rename(Correlation = value)
p_df <- reshape2::melt(p_mat) %>% dplyr::rename(pval = value)
p_adj <- p.adjust(as.vector(p_df$pval, mode = "numeric"), method = "fdr")
cor_df_full <- cor_df %>%
    dplyr::full_join(p_df, by = c("Var1", "Var2")) %>%
    dplyr::mutate(padj = p_adj)

## Add annotations
cor_mat_vec <- as.vector(cor_mat, mode = "numeric")
annot <- sapply(seq(length(p_adj)), function (i) {
    ifelse(p_adj[i] < 0.05, sprintf("%.2f", cor_mat_vec[i]), "")
})
p_adj_mat <- matrix(p_adj, nrow = dim(cor_mat)[1], dimnames = dimnames(p_mat))
annot_mat <- matrix(annot, nrow = dim(cor_mat)[1], dimnames = dimnames(p_mat))

## Plot CT correlation matrix
ct_corr <- ComplexHeatmap::pheatmap(
    t(cor_mat),
    cluster_rows = T,
    cluster_cols = T,
    angle_col = '45',
    cellheight = 20,
    cellwidth = 22,
    clustering_method = "ward.D2",
    scale = "none",
    breaks = seq(-1, 1, length.out = length(ph_pal)),
    border_color = NA,
    color = ph_pal,
    labels_col = convert_greek(rownames(cor_mat)),
    display_numbers = t(annot_mat),
    number_color = "black",
    fontsize_number = 8,
    name = "Correlation",
    heatmap_legend_param = legend_param(1, 3)
)
ct_corr <- patchwork::wrap_elements(
    grid::grid.grabExpr(ComplexHeatmap::draw(ct_corr, merge_legend = T))
)

## Point and line plots
plplots <- list()
for (i in 1:length(analytes)) {
    for (j in 1:length(ct_features)) {
        analyte <- analytes[i]
        ct <- ct_features[j]
        fname <- gsub('-|/| ', '_', paste(analyte, ct, sep = '_'))
        
        df <- data.frame(x = cyto_ct[,ct], y = cyto_ct[,analyte])
        ql <- paste0("q==", deparse(sprintf("%.3f", p_adj_mat[i, j])))
        rl <- paste0("rho==", deparse(sprintf("%.3f", cor(df$x, df$y, method = "s"))))
        pplt <- ggplot(df, aes(x, y)) +
            geom_point() +
            stat_smooth(method = "lm", formula = y ~ x) +
            annotate("text", x = min(df$x), y = c(1.06, 1.13) * max(df$y),
                     hjust = 0, label = c(ql, rl), parse = TRUE) +
            xlab(paste(ct, "fraction of lung")) +
            ylab("Mean concentration") +
            ggtitle(analyte) +
            theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  plot.title = element_text(hjust = 0.5))
        plplots[[fname]] <- pplt
    }
}

## Deconvolution
# Only significant cytokines
mean_expr <- read.csv(
    file = "deidentified_data/deidentified_cytokine_mean_expression.csv",
    check.names = FALSE
)
colnames(mean_expr)[1] <- "Analyte"
mean_expr <- tibble::column_to_rownames(mean_expr, "Analyte")
colnames(mean_expr) <- gsub("macrophages", "MPs", colnames(mean_expr))

mean_expr_scaled <- scale_rows(as.matrix(mean_expr))
s <- ceiling(max(abs(min(mean_expr_scaled)), max(mean_expr_scaled)))

deconv <- ComplexHeatmap::pheatmap(
    mean_expr_scaled,
    cluster_rows = T,
    cluster_cols = F,
    show_colnames = T,
    angle_col = '45',
    cellheight = 18,
    cellwidth = 22,
    clustering_method = "ward.D2",
    scale = "none",
    labels_row = sub_greek(rownames(mean_expr)),
    border_color = NA,
    color = ph_pal,
    breaks = seq(-s, s, length.out = length(ph_pal)),
    name = "Mean\nexpression",
    heatmap_legend_param = legend_param(s, 4)
)
deconv <- patchwork::wrap_elements(
    grid::grid.grabExpr(ComplexHeatmap::draw(deconv, merge_legend = T))
)

# All cytokines
mean_expr_all <- read.csv(
    file = "deidentified_data/deidentified_cytokine_mean_expression_all.csv",
    check.names = FALSE
)
colnames(mean_expr_all)[1] <- "Analyte"
mean_expr_all <- tibble::column_to_rownames(mean_expr_all, "Analyte")

mean_expr_all_scaled <- scale_rows(as.matrix(mean_expr_all))
s <- ceiling(max(abs(min(mean_expr_all_scaled)), max(mean_expr_all_scaled)))



deconv_all <- ComplexHeatmap::pheatmap(
    t(mean_expr_all_scaled),
    cluster_rows = F,
    cluster_cols = T,
    angle_col = '45',
    cellheight = 20,
    cellwidth = 21,
    clustering_method = "ward.D2",
    scale = "none",
    labels_col = sub_greek(rownames(mean_expr_all)),
    border_color = NA,
    color = ph_pal,
    breaks = seq(-s, s, length.out = length(ph_pal)),
    name = "Mean\nexpression",
    heatmap_legend_param = legend_param(s, 5)
)
deconv_all <- patchwork::wrap_elements(
    grid::grid.grabExpr(ComplexHeatmap::draw(deconv_all, merge_legend = T))
)


## Build Figure 5
get_annot <- function (x) {
    patchwork::plot_annotation(
        title = x,
        theme = theme(plot.title = element_text(face = 2, size = 20))
    )
}
thm <- theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

fig_5a <- patchwork::wrap_elements(heatmap + get_annot("a") + thm)
fig_5b <- patchwork::wrap_elements(signif_boxplots + get_annot("b") + thm)
fig_5c <- patchwork::wrap_elements(ct_corr + get_annot("c") + thm)
fig_5d <- patchwork::wrap_elements(plplots[["MCP_1_Fibrotic"]] + get_annot("d") + thm)
fig_5e <- patchwork::wrap_elements(deconv + get_annot("e") + thm)

fig5 <- fig_5a + fig_5b + fig_5c + fig_5d + fig_5e +
    patchwork::plot_layout(
        design = c(patchwork::area(t = 1, b = 1, l = 1, r = 2),
                   patchwork::area(t = 1, b = 1, l = 3, r = 3),
                   patchwork::area(t = 2, b = 2, l = 1, r = 3),
                   patchwork::area(t = 3, b = 3, l = 1, r = 1),
                   patchwork::area(t = 3, b = 4, l = 2, r = 3)),
        heights = c(0.59, 0.15, 0.24, 0.02),
        widths = c(0.5, 0.45, 1)
    ) & theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
ggsave(plot = fig5, filename = "figures/fig_5/fig5.pdf",
       width = 20, height = 20, device = cairo_pdf)

## Build Figure S5
fig_s5_a <- patchwork::wrap_elements(nonsignif_boxplots + get_annot("a") + thm)
fig_s5_b <- patchwork::wrap_elements(deconv_all + get_annot("b") + thm)
fig_s5 <- (fig_s5_a / fig_s5_b) + patchwork::plot_layout(heights = c(0.45, 0.55))
ggsave(plot = fig_s5, filename = "figures/fig_s5/fig_s5.pdf",
       width = 17.5, height = 20, device = cairo_pdf)

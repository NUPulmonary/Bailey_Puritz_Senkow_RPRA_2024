## Generates Figure 4b ##

setwd("/projects/b1038/Pulmonary/cpuritz/PASC/data")
set.seed(2)

library(ggplot2)
library(dplyr)
library(patchwork)

require(ggrepel)

tram1 <- read.csv("01BAL/DEG/output/tram-1/filtered_degs.csv")
tram2 <- read.csv("01BAL/DEG/output/tram-2/filtered_degs.csv")

genes1 <- tram1$X[tram1$sign != ""]
genes2 <- tram2$X[tram2$sign != ""]
tram1_only <- setdiff(genes1, genes2)
tram2_only <- setdiff(genes2, genes1)
common_genes <- intersect(genes1, genes2)

n_up_tram1 <- sum(tram1[tram1$X %in% tram1_only, "sign"] == "Upregulated")
n_down_tram1 <- sum(tram1[tram1$X %in% tram1_only, "sign"] == "Downregulated")
n_up_tram2 <- sum(tram2[tram2$X %in% tram2_only, "sign"] == "Upregulated")
n_down_tram2 <- sum(tram2[tram2$X %in% tram2_only, "sign"] == "Downregulated")
n_up_both <- length(intersect(tram1$X[tram1$sign == "Upregulated"],
                              tram2$X[tram2$sign == "Upregulated"]))
n_down_both <- length(intersect(tram1$X[tram1$sign == "Downregulated"],
                                tram2$X[tram2$sign == "Downregulated"]))

tram1$annot <- "N.S."
tram1[tram1$X %in% tram1_only, "annot"] <- "TRAM-1"
tram1[tram1$X %in% common_genes, "annot"] <- "Both"
tram1$annot <- factor(tram1$annot, levels = c("N.S.", "TRAM-1", "TRAM-2", "Both"))
tram1 <- dplyr::arrange(tram1, annot)
tram1$Name <- tram1$X
tram1[!tram1$X %in% common_genes, "Name"] <- ""
tram1[dim(tram1)[1] + 1,] <- tram1[dim(tram1)[1],]
tram1[dim(tram1)[1], c("X", "annot", "Name")] <- c("", "TRAM-2", "")
tram1$alpha <- c(rep(1, dim(tram1)[1] - 1), 0)

tram2$annot <- "N.S."
tram2[tram2$X %in% tram2_only, "annot"] <- "TRAM-2"
tram2[tram2$X %in% common_genes, "annot"] <- "Both"
tram2$annot <- factor(tram2$annot, levels = c("N.S.", "TRAM-1", "TRAM-2", "Both"))
tram2 <- dplyr::arrange(tram2, annot)
tram2$Name <- tram2$X
tram2[!tram2$X %in% common_genes, "Name"] <- ""
tram2[dim(tram2)[1] + 1,] <- tram2[dim(tram2)[1],]
tram2[dim(tram2)[1], c("X", "annot", "Name")] <- c("", "TRAM-1", "")
tram2$alpha <- c(rep(1, dim(tram2)[1] - 1), 0)

pplt1 <- ggplot(tram1, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(col = annot, alpha = alpha), size = 2) +
    scale_x_log10() +
    ggrepel::geom_label_repel(
        aes(label = Name), show.legend = FALSE, size = 2.5, max.overlaps = Inf,
        force_pull = 1, force = 30, fill = "white", min.segment.length = 0,
        label.padding = 0.15, segment.size = 0.3
    ) +
    labs(x = expression(paste(log["10"], "(base mean)", sep = '')),
         y = expression(paste(log["2"], "(FC)", sep = ''))) +
    ggtitle("TRAM-1") +
    theme(legend.position = "none")

pplt2 <- ggplot(tram2, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(col = annot, alpha = alpha), size = 2) +
    scale_x_log10() +
    ggrepel::geom_label_repel(
        aes(label = Name), show.legend = FALSE, size = 2.5, max.overlaps = Inf,
        force_pull = 1, force = 30, fill = "white", min.segment.length = 0,
        label.padding = 0.15, segment.size = 0.3
    ) +
    labs(x = expression(paste(log["10"], "(base mean)", sep = '')),
         y = expression(paste(log["2"], "(FC)", sep = ''))) +
    ggtitle("TRAM-2") +
    theme(legend.position = "none")

labels <- c(
    "N.S.",
    paste("TRAM-1 (", n_up_tram1, " up, ", n_down_tram1, " down)", sep = ''),
    paste("TRAM-2 (", n_up_tram2, " up, ", n_down_tram2, " down)", sep = ''),
    paste("Shared (", n_up_both, " up, ", n_down_both, " down)", sep = '')
)

pplt <- (pplt1 / pplt2) +
    plot_layout(guides = "collect") &
    scale_color_manual(name = "", labels = labels,
                       values = c("gray", ggsci::pal_npg("nrc")(3))) &
    scale_alpha(guide = 'none') &
    theme_bw(base_family = "Arial") &
    theme(text = element_text(family = "Arial"),
          legend.position = "bottom",
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 12, color = "black"),
          axis.text = element_text(size = 10, color = "black"),
          plot.title = element_text(size = 14, hjust = 0.5))

ggsave(plot = pplt, filename = "figures/fig_4/fig_4b.pdf", width = 13,
       height = 13, device = cairo_pdf, bg = "transparent")

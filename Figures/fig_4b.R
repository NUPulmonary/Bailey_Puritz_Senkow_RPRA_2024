set.seed(0)

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

tram1 <- read.csv("01BAL/DEG/output/tram-1/filtered_degs.csv")
degs <- tram1[tram1$sign != "", "X"]

tram1$annot <- tram1$sign
tram1[tram1$annot == "", "annot"] <- "N.S."
tram1$annot <- factor(tram1$annot,
                      levels = c("N.S.", "Downregulated", "Upregulated"))
tram1 <- dplyr::arrange(tram1, annot)
tram1$Name <- tram1$X
tram1[!tram1$X %in% degs, "Name"] <- ""

pplt <- ggplot(tram1, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(col = annot, alpha = annot), size = 2.2) +
    scale_x_log10() +
    geom_label_repel(
        aes(label = Name), show.legend = FALSE, size = 2.5, max.overlaps = Inf,
        force_pull = 1, force = 30, fill = "white", min.segment.length = 0,
        label.padding = 0.15, segment.size = 0.3
    ) +
    scale_alpha_manual(
        name = "",
        values = c("N.S." = 0.25, "Downregulated" = 1, "Upregulated" = 1)
    ) +
    scale_color_manual(
        name = "",
        values = c("N.S." = "gray", "Downregulated" = "#E64B35FF",
                   "Upregulated" = "#4DBBD5FF")
    ) +
    labs(x = expression(paste(log["10"], "(base mean)", sep = '')),
         y = expression(paste(log["2"], "(FC)", sep = ''))) +
    theme_bw(base_family = "Arial") +
    theme(text = element_text(family = "Arial"),
          legend.position = "bottom",
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12, color = "black"),
          axis.text = element_text(size = 9, color = "black"))

ggsave(plot = pplt,
       filename = "figures/fig_4/fig_4b.pdf",
       width = 8, height = 6,
       device = cairo_pdf, bg = "transparent")

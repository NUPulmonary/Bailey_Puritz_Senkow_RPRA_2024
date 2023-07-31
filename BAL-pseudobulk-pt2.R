library(DESeq2)
library(vsn)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)

set.seed(0)

deg_dir <- "/projects/b1038/Pulmonary/cpuritz/PASC/data/01BAL/DEG"
count_dir <- sprintf("%s/counts", deg_dir)
out_dir <- sprintf("%s/output", deg_dir)
cell_names <- read.csv(sprintf("%s/cell_names.csv", count_dir),
                       check.names = FALSE)

process_cell_type <- function(cell_type, alpha = 0.05) {
    # Read in pseudobulk counts
    counts <- read.csv(sprintf("%s/%s.csv", count_dir, cell_type),
                       row.names = 1, check.names = FALSE, sep = '\t')
    # Read in metadata
    meta <- read.csv(sprintf("%s/%s-meta.csv", count_dir, cell_type),
                     row.names = 1, check.names = FALSE)
    # Create DESeq2 object
    status <- factor(meta$status)
    dds <- DESeqDataSetFromMatrix(counts, DataFrame(status), ~status)
    
    # Create output directory if necessary
    cell_dir <- sprintf("%s/%s", out_dir, cell_type)
    if (!file.exists(cell_dir)) {
        dir.create(cell_dir)
    }
    
    # Plot heatmap of library similarity
    sampleDists <- dist(t(assay(vst(dds))))
    pdf(sprintf("%s/lib-dist.pdf", cell_dir), width = 8, height = 10)
    pheatmap(as.matrix(sampleDists),
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             show_colnames = TRUE,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists)
    dev.off()
    
    # Plot parametric fit dispersion estimates
    # NOTE: in some cases the parametric fit is so bad that DESeq overrides
    # the choice and defaults to the local fit, in which case this plot and
    # the next will be identical
    dds <- DESeq(dds, fitType = "parametric", quiet = TRUE)
    pdf(sprintf("%s/disp-parametric.pdf", cell_dir), width = 6, height = 4)
    plotDispEsts(dds)
    dev.off()
    
    # Plot local fit dispersion estimates
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
    pdf(sprintf("%s/disp-local.pdf", cell_dir), width = 6, height = 4)
    plotDispEsts(dds)
    dev.off()
    
    # Just to ensure a local fit is used
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
    
    # Plot PCA
    pdf(sprintf("%s/pca.pdf", cell_dir), width = 8, height = 6)
    d <- plotPCA(vst(dds, blind = FALSE), ntop = 2000, intgroup = "status",
                 returnData = TRUE)
    percentVar <- round(100 * attr(d, "percentVar"))
    plot(ggplot(d, aes(x = PC1, y = PC2, color = status)) +
             geom_point(size = 3) +
             geom_label_repel(aes(label = colnames(dds)), show.legend = FALSE,
                              fill = "white", size = 1.5) +
             xlab(paste0("PC1: ", percentVar[1], "% variance")) +
             ylab(paste0("PC2: ", percentVar[2], "% variance")) +
             coord_fixed())
    dev.off()
    
    # Extract DEGs
    degs <- results(dds,
                    contrast = c("status", "RPRA", "Healthy"),
                    alpha = alpha)
    degs <- as.data.frame(degs[!is.na(degs$padj), ])
    degs$sign <- ""
    degs_signif <- (degs$padj <= alpha)
    degs$sign[(degs$log2FoldChange < 0) & degs_signif] <- "Downregulated"
    degs$sign[(degs$log2FoldChange > 0) & degs_signif] <- "Upregulated"
    degs <- degs[order(degs$sign),]
    
    # MA plot
    pdf(sprintf("%s/ma.pdf", cell_dir), width = 8, height = 6)
    plot(ggplot(degs) +
             geom_point(aes(x = baseMean, y = log2FoldChange, col = sign)) +
             scale_x_log10() +
             scale_color_manual(values = c("gray", "red", "blue")))
    dev.off()
    
    degs <- degs[order(degs$log2FoldChange, decreasing = TRUE),]
    write.csv(degs, sprintf("%s/degs.csv", cell_dir))
    return(degs)
}

## Process each cell type
degs_all <- NULL
for (i in 1:dim(cell_names)[1]) {
    cell_name <- cell_names[i, 1]
    cell_type <- cell_names[i, 2]
    meta <- read.csv(sprintf("%s/counts/%s-meta.csv", deg_dir, cell_type),
                     row.names = 1, check.names = FALSE)

    # Require at least 3 samples per condition
    num_rpra <- sum(meta$status == "RPRA")
    num_ctrl <- sum(meta$status == "Healthy")
    if (num_rpra < 3) {
        print(sprintf("Only %i RPRA sample(s), skipping %s", num_rpra, cell_name))
    } else if (num_ctrl < 3) {
        print(sprintf("Only %i Healthy sample(s), skipping %s", num_ctrl, cell_name))
    } else {
        print(sprintf("%s/%s: Processing %s", i, dim(cell_names)[1], cell_name))
        degs <- process_cell_type(cell_type) %>%
            dplyr::mutate(cellType = cell_name) %>%
            dplyr::relocate(cellType, .before = dplyr::everything())
        degs_all <- rbind(degs_all, degs)
    }
}

# Save DESeq2 output for all cell types in one file
write.csv(degs_all, sprintf("%s/degs_all.csv", deg_dir))
setwd("/projects/b1038/Pulmonary/cpuritz/PASC/data")
set.seed(0)

library(dplyr)
library(ggplot2)
library(patchwork)

require(ggsignif)
require(ggalluvial)

source("../code/util/pairwise_wilcox_test.R")

##------------ Analysis of CT data ------------##

## Read in data
meta <- read.csv("deidentified_data/deidentified_metadata.csv")
ct_data <- read.csv("deidentified_data/deidentified_ct_data.csv")
ct_data <- ct_data %>%
    dplyr::filter(Type != "Wildcard") %>%
    dplyr::filter(Region == "WholeLung") %>%
    dplyr::select(c("Subject_ID", "Time", "Type", "TypeFrac"))

tf_cols <- unique(ct_data$Type)

## Use new categories
buckets <- list(
    "Normal" = c("NormalNotInflamed"),
    "Fibrotic" = c("Honeycombing", "LinearScar", "Fibronodular", "Reticular",
                   "SubpleuralLine"),
    "Inflammatory" = c("GroundGlass", "NormalInflamed"),
    "Nodularity" = c("Nodule", "Nodular", "CentrilobularNodule"),
    "Emphysema/Cysts" = c("CentrilobularEmphysema", "Emphysematous", "Cyst")
)

group_scans <- function(df, time) {
    df <- df %>%
        dplyr::filter(Time == time) %>%
        tidyr::pivot_wider(names_from = Type, values_from = TypeFrac) %>%
        dplyr::select(-Time)
    for (i in 1:length(buckets)) {
        df <- df %>%
            dplyr::mutate(!!names(buckets)[i] :=
                              rowSums(dplyr::across(dplyr::all_of(buckets[[i]]))))
    }
    dplyr::select(df, -dplyr::all_of(tf_cols))
}

ct_pre <- group_scans(ct_data, "Scan 1")
ct_post <- group_scans(ct_data, "Scan 2")

## Renormalize type fractions
ct_pre[names(buckets)] <- ct_pre[names(buckets)] / rowSums(ct_pre[names(buckets)])
ct_post[names(buckets)] <- ct_post[names(buckets)] / rowSums(ct_post[names(buckets)])

## Statistics of interest
for (f in names(buckets)) {
    print(sprintf("%s: %.1f +/- %.1f", f, 100*mean(ct_pre[[f]]), 100*sd(ct_pre[[f]])))
}

## Boxplots to illustrate compositions of scans
scan_composition <- ggplot(reshape2::melt(ct_pre, id.vars = "Subject_ID"),
                           aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = variable), outlier.shape = NA) +
    geom_jitter(shape = 16, width = 0.15, height = 0) +
    xlab(NULL) +
    ylab("Fraction of scan") +
    scale_fill_manual(values = setNames(ggsci::pal_npg()(5)[c(2, 1, 3, 5, 4)], names(buckets))) +
    theme_bw(base_family = "Arial") +
    theme(text = element_text(family = "Arial"),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"))


## Boxplots of paired pre/post scans
ct_pair <- dplyr::full_join(dplyr::mutate(ct_pre, Time = "Scan 1"),
                            dplyr::mutate(ct_post, Time = "Scan 2"),
                            by = c(names(ct_pre), "Time"))
ct_pair <- dplyr::filter(ct_pair,
                         Subject_ID %in% names(which(table(ct_pair$Subject_ID) == 2)))
ct_pair$Time <- factor(ct_pair$Time, levels = c("Scan 1", "Scan 2"), ordered = TRUE)
all_comps <- boxplot_signif(ct_pair, cols = names(buckets),
                            group.var = "Time", block.var = "Subject_ID",
                            paired = TRUE)
annot_df <- data.frame(name = factor(all_comps$var, levels = names(buckets)),
                       start = paste(all_comps$V1, all_comps$var, sep = '.'),
                       end = paste(all_comps$V2, all_comps$var, sep = '.'),
                       y = 0.97,
                       label = all_comps$annot)

ct_pair_melt <- ct_pair %>%
    tidyr::pivot_longer(cols = -c(Time, Subject_ID)) %>%
    dplyr::mutate(name = factor(name, levels = names(buckets)))

pre_post_boxplots <- 
    ggplot(ct_pair_melt, aes(x = interaction(Time, name), y = value, Time)) +
    geom_boxplot(aes(fill = Time), outlier.shape = NA) +
    geom_line(aes(group = interaction(Subject_ID, name)), linetype = "11",
              color = "gray", linewidth = 0.4) +
    geom_point(shape = 16, size = 1) +
    ggsignif::geom_signif(
        data = annot_df,
        aes(xmin = start, xmax = end, annotations = label, y_position = y),
        manual = TRUE, tip_length = 0, vjust = -0.1
    ) +
    facet_grid(. ~ name, scales = "free_x") + 
    scale_x_discrete(labels = rep(paste("Scan", c(1, 2)), 5)) +
    scale_fill_manual(values = setNames(ggsci::pal_npg("nrc")(9)[c(3, 6)], c("Scan 1", "Scan 2"))) +
    ylim(c(0, 1)) +
    xlab(NULL) +
    ylab("Fraction of scan") +
    theme_bw(base_family = "Arial") +
    theme(text = element_text(family = "Arial"),
          legend.position = "none",
          axis.title = element_text(size = 12, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          strip.text = element_text(size = 12, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())


## Merge pre and post scans into one data frame
ct_all <- dplyr::full_join(ct_pre, ct_post, by = "Subject_ID")
for (x in setdiff(names(ct_pre), "Subject_ID")) {
    ct_all <- ct_all %>%
        mutate(!!paste(x, "(diff)") :=
                   !!sym(paste0(x, ".y")) - !!sym(paste0(x, ".x"))) %>%
        dplyr::rename(!!paste(x, "(pre)") := paste0(x, ".x")) %>%
        dplyr::rename(!!paste(x, "(post)") := paste0(x, ".y"))
}

## Save output
write.csv(ct_pre, file = "deidentified_data/deidentified_ct_scan_1_data.csv", row.names = FALSE)
write.csv(ct_post, file = "deidentified_data/deidentified_ct_scan_2_data.csv", row.names = FALSE)
write.csv(ct_all, file = "deidentified_data/deidentified_ct_data_all.csv", row.names = FALSE)


##---------------------------------------------##


##---------- Build clinical timeline ----------##
## Read in data
clinical <- meta

## Subset to columns of interest
cols <- c("COVID._date", "acute_COVID_admission_date",
          "acute_COVID_intubation_date", "acute_COVID_extubation_date",
          "acute_COVID_discharge_date", "CT.scan.1.date",
          "CT.scan.2.date", "Steroid.start.date", "Steroid.end.date",
          "Bronch.date", "PFT.date")
clinical <- clinical %>%
    tibble::column_to_rownames("Subject_ID") %>%
    dplyr::select(dplyr::all_of(cols)) 

## Fix issues
# Omit this date as it is very old
clinical["RPRA23", "PFT.date"] <- NA
# "x/xx/xx -- other info" -> "x/xx/xx"
for (x in c("RPRA33", "RPRA34")) {
    clinical[x, "Steroid.end.date"] <- strsplit(clinical[x, "Steroid.end.date"], '--')[[1]][1] 
}
# "trach date" -> "date"
clinical$acute_COVID_extubation_date <- gsub("trach ", "",
                                             clinical$acute_COVID_extubation_date)

## Standardize dates and replace remaining non-dates with NA
std_date <- function(x, format = "%m/%d/%y") {
    d <- tryCatch(!is.na(as.Date(x, format)), error = function(err) { FALSE })
    x[!d] <- NA
    x[d] <- as.Date(x[d], format)
    return(x)
}
clinical <- clinical %>% dplyr::mutate(dplyr::across(dplyr::everything(), std_date))

## Center around pivot column
pivot.col <- "COVID._date"
pivot.col.vals <- as.numeric(clinical[[pivot.col]])
clinical <- clinical %>% dplyr::mutate(dplyr::across(dplyr::everything(),
                                                     function(x) { as.numeric(x) - pivot.col.vals }))
## Sort subjects
clinical <- clinical[order(apply(clinical, 1, max, na.rm = TRUE)), ]

## Convert to months
clinical <- clinical / 30

## Record min/max of each row, and add column to specify height
row_mins <- apply(clinical, 1, min, na.rm = TRUE)
row_maxs <- apply(clinical, 1, max, na.rm = TRUE)
clinical$y <- seq(dim(clinical)[1])
minmax_df <- data.frame(xmin = row_mins, xmax = row_maxs, y = clinical$y)

## Compare change in normal fraction with interval between CT scans
ct_interval_df <- clinical %>%
    dplyr::select("CT.scan.1.date", "CT.scan.2.date") %>%
    dplyr::mutate(CT.interval = CT.scan.2.date - CT.scan.1.date) %>%
    dplyr::select(-CT.scan.1.date, -CT.scan.2.date)

print(sprintf("Median time between scans = %.0f +/- %.1f days",
              30 * median(na.omit(ct_interval_df$CT.interval)),
              30 * sd(na.omit(ct_interval_df$CT.interval))))

ct_interval_df <- ct_interval_df %>%
    tibble::rownames_to_column("Subject_ID") %>%
    dplyr::full_join(ct_all, by = "Subject_ID") %>%
    dplyr::select("Normal (diff)", "CT.interval") %>%
    na.omit() %>%
    dplyr::select(Normal.change = "Normal (diff)", CT.interval)

## Correlation of time between scans and normal change
print(sprintf("Correlation between interval and normal change = %.2f",
              cor(ct_interval_df$Normal.change, ct_interval_df$CT.interval, method = "spearman")))
test <- coin::spearman_test(
    Normal.change ~ CT.interval, data = ct_interval_df,
    distribution = coin::approximate(nresample = 1e7,
                                     parallel = "multicore",
                                     ncpus = 16L)
)
print(sprintf("p-value = %.2f", as.numeric(pvalue(test))))

fig_s1 <- ggplot(ct_interval_df, aes(x = CT.interval, y = Normal.change)) +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x) +
    xlab("Months between CT scans") +
    ylab("Change in normal fraction") +
    theme_bw(base_family = "Arial") +
    theme(text = element_text(family = "Arial"),
          legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 13, color = "black"),
          axis.title.y = element_text(size = 13, color = "black"))
ggsave(plot = fig_s1,
       filename = "figures/fig_s1/fig_s1.pdf",
       width = 8, height = 7, device = cairo_pdf)

## Melt data frame
clinical <- clinical %>%
    tibble::rownames_to_column("Study.ID") %>%
    dplyr::relocate("Study.ID", .after = last_col()) %>%
    reshape2::melt(., id.vars = "y", variable.name = "cdata") %>%
    na.omit()

## Set plot aesthetics
cnames <- c("COVID-19 diagnosis", "ICU admission", "Hospital discharge",
            "Intubation", "Extubation/Tracheostomy", "CT Scan 1", "CT Scan 2",
            "Bronchoscopy", "PFT", "Steroid treatment")
colors <- setNames(c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
                     "#17BECFFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF",
                     "#9467BDFF", "#FFAEB9"), cnames)
shapes <- setNames(c(1, 3, 4, 0, 7, 6, 2, 10, 13, 1), cnames)
alphas <- setNames(c(rep(1, length(cnames) - 1), 0), cnames)
linetypes <- setNames(rep(1, length(cnames)), cnames)

## Plot steroid treatment as a line, everything else as a point
ids <- clinical %>% dplyr::filter(cdata == "Study.ID") %>% dplyr::pull(value)
dates <- clinical %>%
    dplyr::filter(!cdata %in% c("Study.ID", "Group")) %>%
    dplyr::mutate(value = as.numeric(value))
points <- dates %>% dplyr::filter(cdata != "Steroid.start.date")
lines <- data.frame(x = dates[dates$cdata == "Steroid.start.date", "value"],
                    xend = dates[dates$cdata == "Steroid.end.date", "value"],
                    y = dates[dates$cdata == "Steroid.start.date", "y"])

## Adjust names in legend
points$cdata <- points$cdata %>%
    dplyr::recode(COVID._date = "COVID-19 diagnosis",
                  acute_COVID_admission_date = "ICU admission",
                  acute_COVID_intubation_date = "Intubation",
                  acute_COVID_extubation_date = "Extubation/Tracheostomy",
                  acute_COVID_discharge_date = "Hospital discharge",
                  Steroid.end.date = "Steroid treatment",
                  CT.scan.1.date = "CT Scan 1",
                  CT.scan.2.date = "CT Scan 2",
                  Bronch.date = "Bronchoscopy",
                  PFT.date = "PFT") %>%
    droplevels() %>%
    factor(levels = cnames)

## Set plot parameters
pt_size <- 2
pt_thickness <- 0.6
line_size <- 1.5
line_alpha <- 0.5
seg_size <- 0.45

timeline <- ggplot(data = points) +
    geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = y),
                 color = colors[["Steroid treatment"]], linewidth = line_size,
                 alpha = line_alpha) +
    geom_segment(data = minmax_df, aes(x = xmin, xend = xmax, y = y, yend = y),
                 linewidth = seg_size, color = "gray") +
    geom_line(aes(x = value, y = y, color = cdata, linetype = cdata),
              linewidth = line_size, alpha = 0) +
    geom_point(aes(x = value, y = y, color = cdata, shape = cdata, alpha = cdata),
               stroke = pt_thickness, size = pt_size) +
    scale_alpha_manual(values = alphas) +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = linetypes) +
    scale_color_manual(values = colors,
                       guide = guide_legend(override.aes = list(
                           alpha = c(rep(1, length(alphas) - 1), line_alpha),
                           size = c(rep(pt_size, length(alphas) - 1), line_size),
                           linetype = c(rep("blank", length(linetypes) - 1), "solid"),
                           shape = c(shapes[-length(shapes)], NA)))) +
    xlab("Months from COVID-19 Diagnosis") +
    ylab(NULL) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 11.5, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 11, color = "black"),
          axis.text.y = element_text(size = 8, color = "black"),
          legend.position = "bottom") +
    scale_x_continuous(breaks = c(0, 6, 12, 18),
                       expand = c(0.0125, 0, 0.0125, 0)) +
    scale_y_continuous(breaks = seq(max(points$y)),
                       labels = ids,
                       expand = c(0.02, 0, 0.02, 0))

##---------------------------------------------##

##--------------- Alluvial plot ---------------##

tbl <- read.csv("deidentified_data/deidentified_metadata.csv", header = TRUE) %>%
    dplyr::select("Flow" = Flow.Cytometry,
                  "Cytokines" = Cytokines,
                  "CT1" = CT.Scan.1,
                  "CT2" = CT.Scan.2,
                  "BAL" = BAL.scRNAseq,
                  "NEP" = NEP.scRNAseq) %>%
    dplyr::count(dplyr::across(dplyr::everything())) %>%
    tidyr::pivot_longer(cols = c("CT1", "CT2", "Flow", "BAL", "Cytokines", "NEP")) %>%
    dplyr::select(Analysis = name, Freq = n, Value = value) %>%
    dplyr::mutate(Analysis = factor(Analysis,
                                    levels = c("CT1", "Flow", "BAL", "Cytokines", "NEP", "CT2"))) %>%
    dplyr::arrange(Analysis) %>%
    dplyr::mutate(Subject = rep(seq(9), 6)) %>%
    dplyr::mutate(Analysis = recode(Analysis, CT1 = "CT Scan 1", CT2 = "CT Scan 2",
                                    Flow = "BAL Flow Cytometry", BAL = "BAL scRNA-seq",
                                    NEP = "Nasal scRNA-seq",
                                    Cytokines = "BAL Cytokines"))
river <- ggplot(tbl, aes(x = Analysis, stratum = Value, alluvium = Subject,
                         y = Freq, fill = Value, label = Value)) +
    scale_x_discrete(expand = c(0, 0.2, 0, 0.3)) +
    scale_color_manual(name = "", values = setNames(ggsci::pal_npg()(2), c("No", "Yes"))) +
    ggalluvial::geom_flow(alpha = 0.6) +
    ggalluvial::geom_stratum(alpha = 0.6) +
    geom_text(stat = "stratum", size = 3) +
    xlab(NULL) +
    ylab("Number of participants") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.margin = unit(c(-0.6, 0.05, 0.05, 0.05), "cm"),
          axis.text.x = element_text(size = 11, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank())

##---------------------------------------------##

##-------------- Build Figure 1 ---------------##

thm <- theme(plot.title = element_text(face = 2, size = 17))
plt_a <- wrap_elements(timeline + plot_annotation(title = "a", theme = thm))
plt_b <- wrap_elements(river + plot_annotation(title = "b", theme = thm))
plt_c <- wrap_elements(scan_composition + plot_annotation(title = "c", theme = thm))
plt_d <- wrap_elements(pre_post_boxplots + plot_annotation(title = "d", theme = thm))

pplt <- plt_a + plt_b + plt_c + plt_d + plot_layout(
    design = c(area(t = 1, b = 1, l = 1, r = 2),
               area(t = 2, b = 2, l = 1, r = 1),
               area(t = 2, b = 2, l = 2, r = 2),
               area(t = 3, b = 3, l = 1, r = 2)),
    widths = c(8.75, 7.25),
    heights = c(7, 4.5, 4.5))

ggsave(plot = pplt,
       filename = "figures/fig_1/fig1.pdf",
       width = 16, height = 16, device = cairo_pdf)
# Perform pairwise Wilcox rank sum tests with FDR correction
#
#' Pairwise Wilcoxon Rank Sum Test
#'
#' Perform pairwise Wilcoxon rank sum tests with FDR correction. This function
#' provides two advantages over the `stats::pairwise.wilcox.test` function:
#'   (1) Exact p-values are computed in the presence of ties using `coin::wilcox_test`
#'   (2) FDR is controlled over all response variables simultaneously
#'
#' @param df data frame
#' @param cols the names of the response variable columns in df
#' @param group.var the name of the column indicating the grouping variable
#' @param block.var the name of the block factor for paired data
#' @param var.name the name of the output variable
#' @param alpha if not NULL, remove comparisons with an adjusted p-value less 
#' than alpha.
#' @param ... further arguments to be passed to `coin::wilcox_test`
pairwise_wilcox_test <- function(df,
                                 cols,
                                 group.var,
                                 block.var = NULL,
                                 var.name = "var",
                                 alpha = NULL,
                                 paired = FALSE,
                                 ...) {
    require(coin)
    require(dplyr)
    
    groups <- unique(as.character(df[[group.var]]))
    relevant_comps <- combn(groups, 2) %>% t() %>% as.data.frame()
    all_comps <- lapply(cols, function(x) {
        p_vals <- vector(mode = "numeric", length = nrow(relevant_comps))
        for (i in 1:nrow(relevant_comps)) {
            sub <- df %>% dplyr::filter(!!sym(group.var) %in% relevant_comps[i,])
            if (!paired) {
                sub <- sub %>%
                    dplyr::select(all_of(eval(c(group.var, x)))) %>%
                    dplyr::rename(value = x)
                sub.df <- data.frame(value = sub$value,
                                     group = factor(sub[[group.var]]))
                test <- coin::wilcox_test(value ~ group,
                                          data = sub.df,
                                          distribution = "exact",
                                          ...)
            } else {
                sub <- sub %>%
                    dplyr::select(all_of(eval(c(group.var, x, block.var)))) %>%
                    dplyr::rename(value = x)
                sub.df <- data.frame(value = sub$value,
                                     group = factor(sub[[group.var]]),
                                     block = factor(sub[[block.var]]))
                test <- coin::wilcoxsign_test(value ~ group | block,
                                              data = sub.df,
                                              distribution = "exact",
                                              ...)
            }
            p_vals[i] <- coin::pvalue(test)
        }
        return(relevant_comps %>% dplyr::mutate(!!var.name := x, pval = p_vals))
    })
    
    # Merge tests and control FDR
    all_comps <- all_comps %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(padj = p.adjust(pval, method = "fdr"))
    
    # Filter comparisons
    if (!is.null(alpha)) {
        all_comps <- all_comps %>% dplyr::filter(padj < alpha)
    }
    
    return(all_comps)
}

boxplot_signif <- function(df,
                           cols,
                           group.var,
                           block.var = NULL,
                           var.name = "var",
                           alpha = 0.05,
                           annot = c("value", "stars"),
                           paired = FALSE) {
    # Perform pairwise Wilcoxon rank sum tests
    all_comps <- pairwise_wilcox_test(df,
                                      cols,
                                      group.var,
                                      block.var = block.var,
                                      var.name = var.name,
                                      alpha = alpha,
                                      paired = paired)
    
    annot <- match.arg(annot)
    if (annot == "value") {
        all_comps <- all_comps %>%
            dplyr::mutate(annot = format(padj, digits = 3, scientific = T))
    } else if (annot == "stars") {
        get_stars <- function(x) {
            if (x <= 0.001) {
                return("***")
            } else if (x <= "0.01") {
                return("**")
            } else if (x <= "0.05") {
                return("*")
            } else {
                return("n.s.")
            }
        }
        all_comps <- all_comps %>%
            dplyr::mutate(annot = sapply(all_comps$padj, get_stars))
    }
    
    # Add in x-coordinates for significance bars
    groups <- unique(df[[group.var]])
    all_comps <- all_comps %>%
        dplyr::left_join(data.frame(V1 = sort(groups), xmin = 1:length(groups)), by = "V1") %>%
        dplyr::left_join(data.frame(V2 = sort(groups), xmax = 1:length(groups)), by = "V2") %>%
        tibble::rownames_to_column(var = "ix")
    
    # Add in y-coordinates for significance bars, 17.5% higher than largest data
    if (dim(all_comps)[1] > 0) {
        all_comps$y <- apply(all_comps, 1, function(x) {
            yvals <- df %>%
                dplyr::filter(!!sym(group.var) %in% !!x[c("V1", "V2")]) %>%
                dplyr::select(x[[var.name]])
            return(1.175 * max(yvals))
        })
    }
    
    # Adjust y-coordinates if there are multiple significance bars
    for (i in 1:length(cols)) {
        comps <- all_comps %>% dplyr::filter(!!sym(var.name) == cols[i])
        if (nrow(comps) > 0) {
            diffs <- abs(comps$xmax - comps$xmin)
            ranks <- rank(diffs, ties.method = "first")
            for (j in 1:nrow(comps)) {
                all_comps[comps[ranks[j],]$ix,]$y <- (1 + 0.375*(j - 1)) * max(comps$y)
            }
        }
    }
    return(all_comps)
}

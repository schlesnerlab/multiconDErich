#' Splits string at -vs- to extract the two contrast groups.
#'
#' @param cont string containing the contrast groups
#' @return A character vector of length two with the contrast split apart
#' @export
split_contrast_groups <- function(cont) {
  stringr::str_split(cont, "-vs-", simplify = TRUE) %>% unlist()
}

#' Joins the DEseq table and the fpkm read table together
#'
#' @param DEseq_tb Deseq result table from DEseq. Saved as TSV table by deseq rule
#' @param fpkm FPKM table from snakemake pipeline. You can use any fpkm table.
#' @return sorted combined tibble between the input files.
#' @importFrom magrittr %>%
#' @export
join_tables <- function(DEseq_tb, fpkm) {
  diff_exp <- DEseq_tb

  filer <- fpkm %>%
    dplyr::filter(gene %in% diff_exp$gene_id) %>%
    dplyr::filter(!duplicated(gene))
  joined_df <-
    dplyr::inner_join(filer, DEseq_tb, by = c("gene" = "gene_id")) %>% dplyr::arrange(padj)

  return(joined_df)
}
#' Calculates mean values from the Four groups, which are still hardcoded. and returns
#' the mean values plus gene names and their adjusted p values
#'
#' @param mat Object that can be safely coerced to tibble with col names from Vascaging
#' @param extra_col Extra column that should be extracted and added
#' @param contrast_groups contrast_groups
#' @param s_map Sample map of samples
#' @param cond_id condition ID to setup by
#' @export
#' @return tibble with mean values, gene names and padj values
#' @importFrom magrittr %>%
mean_tibble_from_mat <-
  function(mat,
           extra_col = NULL,
           contrast_groups = contrast_groups,
           s_map,
           cond_id) {
    mat_tb <- tibble::as_tibble(mat)

    mean_tb <- tibble::tibble(.rows = nrow(mat_tb))
    for (gr in contrast_groups) {
      col_n <- s_map %>%
        dplyr::filter(!!sym(cond_id) == gr) %>%
        dplyr::pull(sample)
      mean_tb[[gr]] <- rowMeans(dplyr::select(mat_tb, col_n))
    }
    if ("gname" %in% colnames(mat_tb)) {
      mean_tb["gname"] <- mat_tb %>% dplyr::select("gname")
    }
    if ("padj" %in% colnames(mat_tb)) {
      mean_tb["padj"] <- mat_tb %>% dplyr::select("padj")
    }
    if (!is.null(extra_col)) {
      mean_tb[extra_col] <- mat_tb %>% dplyr::select(extra_col)
    }

    mean_tb
  }


#' test
#'
#' @param info Info field to check
#' @param index index that is to be checked
#'
#' @return
#' @export
#'
#' @examples
#' NULL
check_against_index <- function(info, index) {
  return(as.numeric(index %in% info))
}


#' Create dynamic table from significant gene names
#'
#' @param sig_gene_names import list of genes to put in matrix
#'
#' @return DT::datatable of genes across comparison groups
#' @export
#'
#' @examples
#' create_overlap_matrix(list("g1" = c("gene1", "gene2", "gene3"), "g2" = c("gene1")))
#'
create_overlap_matrix <- function(sig_gene_names) {
  gene_name_index <- unlist(sig_gene_names) %>%
    unique() %>%
    sort()



  gene_matrix <-
    sapply(sig_gene_names, check_against_index, index = gene_name_index)
  rownames(gene_matrix) <- gene_name_index

  gene_matrix <- cbind(gene_matrix, total = rowSums(gene_matrix))

  ordered_matrix <-
    gene_matrix[order(gene_matrix[, "total"], decreasing = TRUE), ]

  DT::datatable(ordered_matrix)
}
#' filters table so only contrasted samples remain
#'
#' @param rld DESeq2 object that has been normalized with either (rlog2 or vst)
#' @param contrast_groups character vecctor of contrast groups being analyzed.
#' @param j_df joined tibble of deseq + fpkm values to match genes and gnames
#' @param reorder_genes whether the order of genes should be orded by logfolgchange
#' @param cond_id condition variable set in snakeamek
#' @return matrix of epxression values
#' @export
filter_split_table <-
  function(rld,
           contrast_groups,
           j_df,
           reorder_genes = TRUE,
           cond_id) {
    expression_values <-
      as.data.frame(SummarizedExperiment::assay(rld)[j_df$gene,
        rld@colData[, cond_id] %in%
          contrast_groups,
        drop = FALSE
      ])
    rownames(expression_values) <- j_df$gname
    if (reorder_genes) {
      expression_values <-
        expression_values[order(abs(j_df$logFoldChange)), ]
    }
    expression_values
  }

#' Save pheatmapplot to svg  file
#'
#' @param x plot object to save
#' @param filename name where to save file
#' @param width plotwidth
#' @param height plotheight
#' @return nothing
#' @importFrom grDevices dev.off svg
#' @export
save_pheatmap_svg <- function(x,
                              filename,
                              width = 7,
                              height = 7) {
  svglite::svglite(filename, width = width, height = height,fix_text_size = FALSE)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#' Save any plot to svg  file
#'
#' @param x plot object to save
#' @param filename name where to save file
#' @param width plotwidth
#' @param height plotheight
#' @return nothing
#' @importFrom grDevices dev.off svg
#' @export
save_plot_svg <- function(x,
                          filename,
                          width = 7,
                          height = 7) {
  svglite::svglite(filename, width = width, height = height,fix_text_size = FALSE)
  x
  dev.off()
}

#' Save pheatmapplot to svg  file
#'
#' @param x plot object to save
#' @param filename name where to save file
#' @param width plotwidth
#' @param height plotheight
#' @return nothing
#' @importFrom grDevices dev.off svg
#' @import svglite
#' @export
save_cheatmap_svg <- function(x,
                              filename,
                              width = 7,
                              height = 7) {
  svglite::svglite(filename, width = width, height = height,fix_text_size = FALSE)
  ComplexHeatmap::draw(x)
  dev.off()
}

#' Parse deseq2 results
#'
#' @param filepath Path to DESeq2 results object form res saved as tsv
#'
#' @return Tibble of DES2 results.
#' @export
#'
#' @examples
read_deseq_tibble <- function(filepath) {
  readr::read_tsv(
    filepath,
    col_names = c(
      "gene_id",
      "baseMean",
      "logFoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    ),
    skip = 1
  )
}

#' Convert fpkm values to tpm
#'
#' @param fpkm_mat matrix of fpkm values (rows = genes, columns = samples)
#'
#' @return a matrix of tpm values
#' @export
#'
#' @examples
fpkm_to_tpm <- function(fpkm_mat) {
  fpkm_colsum <- colSums(fpkm_mat)
  tpm_mat <- (fpkm_mat/fpkm_colsum) * 10^6 
  tpm_mat
}

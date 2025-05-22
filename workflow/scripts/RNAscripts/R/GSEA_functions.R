#' Get entrezgenenames from ensembl
#'
#' @param gene_names vecotr of gene symbols to be converted
#' @param input_type String denoting the type of input given based on \link[org.Mm.eg.db]{org.Mm.eg.db}
#' @param org_db org db package for the organism analyzed during workflow.
#' @return table of gene names to entrezgene ids
#' @export
#' @examples
#' get_entrezgenes_from_ensembl(c("Aplnr"), org_db = org.Mm.eg.db::org.Mm.eg.db)
get_entrezgenes_from_ensembl <- function(gene_names, input_type = "SYMBOL",
                                         org_db) {
  if (input_type == "ENSEMBL") {
    gene_names <- stringr::str_extract(string = gene_names, "^ENS[A-Z0-9]*")
  }

  ensembl_to_eg <- clusterProfiler::bitr(gene_names,
    fromType = input_type, toType = "ENTREZID",
    OrgDb = org_db
  )
  return(ensembl_to_eg)
}

#' transform table of genes and log fold changes to a vector of LFCs with genes as names
#'
#' @param gene_tb data.frame like object with first col gene names and second col LFC values
#' @param input_type String denoting the type of input given based on \link[org.Mm.eg.db]{org.Mm.eg.db}
#' @param org_db org db package for the organism analyzed during workflow.
#' @return a named vector of column 2 with names from column 1
#'
#' @importFrom magrittr %>%
#' @export
get_entrezgene_vector <- function(gene_tb, input_type = "SYMBOL", org_db) {
  gene_tb <- as.data.frame(gene_tb)
  if (input_type == "ENSEMBL") {
    gene_tb[, 1] <- stringr::str_extract(string = gene_tb[, 1], "^ENS[A-Z0-9]*")
  }
  eg <- get_entrezgenes_from_ensembl(gene_tb[, 1], input_type, org_db = org_db)
  rownames(gene_tb) <- gene_tb[, 1]
  gene_tb <- gene_tb[eg[, 1], ]
  gene_ranks <- stats::setNames(gene_tb[, 2], nm = eg$ENTREZID) %>%
    sort(decreasing = TRUE)
  gene_ranks
}

#' Take two columns from a data frame and move these to a list
#'
#' Useful for dealing with duplicated entries to build a translation table
#' @param tb Tibble containing the two columns to transform
#' @param ref_col Reference column to be used a name
#' @param val_col Value Columns containing values
#' @export
table_to_list <- function(tb, ref_col, val_col) {
  formu <- paste0(val_col, "~", ref_col) %>% as.formula()
  output_list <- tb %>% unstack(formu)

  output_list
}

#' transform gene set using translation table
#'
#' Translation table bulilt by \link{table_to_list}
#' @param gene_set two_column gene_set data_frame to transform
#' @param gene2id Translation table from \link{table_to_list}
#' @export
transform_glist <- function(gene_set, gene2id) {
  index <- gene_set %>% dplyr::pull(2)
  self_hatred <- gene2id[as.character(index)]
  length_info <- sapply(self_hatred, length)
  tibble::tibble(
    gs_name = rep(gene_set$gs_name, length_info),
    ens_id = na.omit(unlist(self_hatred))
  )
}

#' Run over representation
#'
#' @inheritParams gsea_test
#' @param universe Universe of all tested genes from the experiment
#' @param input_type Gene input type used f.e. ENSEMBL or gene symbol
#' @param ... passed to enricher function \link[clusterProfiler]{enricher}
#'
#' @return enrichResult object from \link[DOSE]{DOSE-package}
over_rep_test <- function(DE_genes, T2G, universe, input_type = "ENSEMBL", ...) {
  DE_genes <- as.data.frame(DE_genes)
  if (input_type == "ENSEMBL") {
    DE_genes[, 1] <- stringr::str_extract(string = DE_genes %>% dplyr::pull(1), "^ENS[A-Z0-9]*")
  }
  enrichment_result <- clusterProfiler::enricher(
    gene = DE_genes[, 1],
    pvalueCutoff = 0.05,
    universe = universe,
    TERM2GENE = T2G,
    ...
  )
  return(enrichment_result)
}

#' Run gene set enrichment analysis using \link[fgsea]{fgsea}
#'
#' @param DE_genes Tibble of diferentially expressed genes and their log fold
#' changes
#' @param T2G Term to Gene table of gene set enrichment.
#' @param input_type  Gene ID type used in database
#' @param ... Arguments passed to \link[clusterProfiler]{GSEA}
#' @return enrichment result from \link[fgsea]{fgsea}
#' @importFrom magrittr %>%
gsea_test <- function(DE_genes, T2G, input_type = "gene_symbol", ...) {
  DE_genes <- as.data.frame(DE_genes)
  if (input_type == "ENSEMBL") {
    DE_genes[, 1] <- stringr::str_extract(string = DE_genes %>%
      dplyr::pull(1), "^ENS[A-Z0-9]*")
   
    T2G[, 2] <- stringr::str_extract(string = T2G %>% dplyr::pull(2), "^ENS[A-Z0-9]*")
  }
  glist <- stats::setNames(DE_genes[, 2], nm = DE_genes[, 1]) %>%
    sort(decreasing = TRUE)
  enrichment_result <- clusterProfiler::GSEA(
    geneList = glist,
    TERM2GENE = as.data.frame(T2G),
    eps = 0,
    ...
  )
  enrichment_result
}


#' Run the msig using clusterprofiler on a gene set list.
#'
#' @param gset_list A data.frame/tibble with enriched genes and their log fold changes
#' @param category Msigdb categories to include. Defaults to NULL (all categories)
#' @param subcategory Msigdb subcategory
#' @param species Latin name of species where data comes from. Defaults to Mus musculus.
#' @param GSEA logical argument if GSEA is being used or over representation test
#' @param universe universe of all tested genes to be given when GSEA is FALSE
#' @param translation_table named list translating gene IDs from \link{table_to_list}
#' @param msdb_var Variable to use for gsea from msgdbr table \link[msigdbr]{msigdbr}
#' @param input_type Type of input in gset_list
#' @param custom_geneset Custom gene set to use instead of msigdb. Format is a data.frame with two columns
#' @param ... Parameters passed to enichment functions from clusterProfiler
#'
#' @return \link[DOSE]{enrichResult-class} of given information
#' @export
run_msig_enricher <- function(gset_list, category = NULL, species = "Mus musculus",
                              subcategory = NULL,
                              GSEA = TRUE, universe = NULL, translation_table = NULL,
                              msdb_var = "gene_symbol", input_type = "gene_symbol",
                              custom_geneset = NULL,  # Neuer Parameter für benutzerdefiniertes Gen-Set
                              ...) {
  stopifnot(class(gset_list) == "list")
  col_names <- c("gs_name", msdb_var)
  
  # Überprüfen Sie, ob ein benutzerdefiniertes Gen-Set bereitgestellt wurde
  if (is.null(custom_geneset)) {
    msg_class <- msigdbr::msigdbr(
      species = species,
      category = category,
      subcategory = subcategory
    ) %>%
      dplyr::select(col_names)
  } else {
    # Verwenden Sie das benutzerdefinierte Gen-Set, wenn es bereitgestellt wurde
    msg_class <- custom_geneset
  }
  
  if (!is.null(translation_table)) {
    msg_class <- transform_glist(msg_class, translation_table)
  }
  if (GSEA) {
    enrich_res <- purrr::map(gset_list, gsea_test,
      T2G = msg_class,
      input_type = input_type,
      ...
    )
  } else {
    if (input_type == "ENSEMBL") {
      universe <- stringr::str_extract(string = universe, "^ENS[A-Z0-9]*")
    }
    enrich_res <- purrr::map(gset_list, over_rep_test,
      T2G = msg_class, universe = universe,
      input_type = input_type,
      ...
    )
  }
  return(enrich_res)
}

#' Function to run all gsea queries based on config file
#'
#' @param gsea_genes tibble with gene ids in first column and gsea statistic in second
#' @param de_genes tibble with DE genes in first column
#' @param gset_name Name of gene set in config to analyze
#' @param gset_config gset config object
#' @param species species name
#' @param org_db org.mm.eg.db
#' @param t_table translation_table to convert EnsemblIDs to gene symbols
#' 
#' @return List with results enrichment analysis \link[DOSE]{enrichResult-class}
#' @export
#' 
run_gsea_query <- function(gsea_genes, de_genes, gset_name, 
                          gset_config, species, org_db, t_table = NULL) {
  run_settings <- gset_config[[gset_name]]
  if (run_settings$use_gsea) {
    gene_list <- gsea_genes
  } else {
    gene_list <- de_genes
  }
  if (tolower(run_settings$database) == "msigdb") {
    enrich_obj <- RNAscripts::run_msig_enricher(gset_list = list(gene_list),
                                GSEA = run_settings$use_gsea,
                                category = run_settings$category,
                                subcategory = run_settings$subcategory,
                                species = species,
                                msdb_var = "ensembl_gene", 
                                input_type = "ENSEMBL"
                                )[[1]]
  } else if (tolower(run_settings$database) == "kegg") {
    enrich_obj <- RNAscripts::run_gsea(gene_list,
                          input_type = "ENSEMBL", 
                          p_valcut = 0.05, 
                          species = species)
  } else if (tolower(run_settings$database) == "reactome") {
    g_vec <- RNAscripts::get_entrezgene_vector(gene_list, "ENSEMBL", 
                                               org_db = org_db)
    enrich_obj <- ReactomePA::gsePathway(g_vec,
                                                 tolower(RNAscripts::get_organism_omnipath_name(organism)),
                                                 )
  } else if (tolower(run_settings$database) == "custom_senescence") {
    # TODO: Add Human Version
    if (species == "Mus musculus") {
      senescence_genes <- RNAscripts::senes
    } else {
      stop("Species not supported yet. senescence gene set only available for Mus musculus")
    }

    # Transform from tibble to dataframe for clusterProfiler
    senescence_genes <- as.data.frame(senescence_genes)
    # Convert gene symbols to ENSEMBL IDs using named list t_table
    if (!is.null(t_table)) {
      ensembl_senes <- senescence_genes
      ensembl_senes[,2] <- t_table[as.character(ensembl_senes %>% dplyr::pull(2))] %>% as.character()
      # remove all rows where no ENSEMBL ID was found (NULL in second column)
      ensembl_senes <- ensembl_senes %>% dplyr::filter(!(ensembl_senes[,2])== "NULL")
      senescence_genes <- ensembl_senes
    }
    # Run GSEA

    enrich_obj <- gsea_test(DE_genes =  gene_list,
      T2G = senescence_genes,
      input_type = "ENSEMBL")

  } else if (tolower(run_settings$database) == "mitocarta") {
    if (species == "Mus musculus") {
      mito_genes <- RNAscripts::MitoPathways
    } else {
      stop("Species not supported yet. MitoCarta gene set only available for Mus musculus")
    }

    # Transform from tibble to dataframe for clusterProfiler
    mito_genes <- as.data.frame(mito_genes)
    # Since we are using ENSEMBL IDs we delete gene symbols columns
    mito_genes <- mito_genes %>% dplyr::select(GeneSet, EnsemblGeneID)    

    enrich_obj <- gsea_test(DE_genes =  gene_list,
      T2G = mito_genes,
      input_type = "ENSEMBL")
  } else {
    stop(glue::glue("{run_settings$database} not supported, pleasue use MSigDB, kegg or Reactome"))
  }
  gc()
  enrich_obj
}

#' Run KEGG gene set enrichment on DE_tb table of genes
#'
#' @param DE_tb Data.frame like object of col 1 gene symbols, col2 LFC values
#' @param input_type String denoting the type of input given based on \link[org.Mm.eg.db]{org.Mm.eg.db}
#' @param p_valcut Cutoff passed to gseKEGG
#' @param species Name of species analyzed in workflow
#' @return \link[clusterProfiler]{gseKEGG}
#' @importFrom magrittr %>%
#' @export
run_gsea <- function(DE_tb, input_type, p_valcut = 0.05, species) {
  DE_tb <- as.data.frame(DE_tb)
  if (input_type == "ENSEMBL") {
    DE_tb[, 1] <- stringr::str_extract(string = DE_tb[, 1], "^ENS[A-Z0-9]*")
  }
  org_db <- RNAscripts::get_org_db(species)
  eg <- get_entrezgenes_from_ensembl(DE_tb[, 1],
    input_type = input_type,
    org_db = org_db
  )
  rownames(DE_tb) <- DE_tb[, 1]
  DE_tb <- DE_tb[eg[, 1], ]
  # DE_tb <- DE_tb %>% dplyr::filter(DE_tb[,1] %in% eg[,1])

  gene_ranks <- stats::setNames(DE_tb[, 2], nm = eg$ENTREZID) %>% sort(decreasing = T)

  kegg_species <- get_kegg_name(species)
  fgsea_results <- clusterProfiler::gseKEGG(gene_ranks, organism = kegg_species, pvalueCutoff = p_valcut, eps = 0)

  return(fgsea_results)
}
#' Function to create Reactome DB for gene centric analysis
#'
#' @export
#' @param output_type Type of gene ID to be used as output. Needs to be compatible
#' @param species Name of species being analyzed
#' with org.Mm.eg.db keys
buildReactome <- function(output_type = "ENSEMBL", species) {
  # Get Human gene annotations
  mm <- get_org_db(species)

  # Get Reactome IDs and names
  reac <- as.list(reactome.db::reactomePATHID2EXTID)
  reac.names <- as.list(reactome.db::reactomePATHID2NAME)
  # FIlter for Human Pathways
  index <- grep("R-MMU", names(reac.names), value = TRUE)

  mouse.reac <- reac[index]

  # Get reactome names
  names(mouse.reac) <- unlist(reac.names[index])


  # query unique reactome indices
  query <- unique(unlist(mouse.reac))
  # Create data.frame of ENTREZID to GENE SYMBOLS from
  entrez.to.gene <- AnnotationDbi::select(
    mm,
    keys = query,
    columns = c("ENTREZID", output_type),
    keytype = "ENTREZID"
  )

  # Build Reactome Database by iterating across homo.reac and switch entrez-gene
  # Ids to HGNC gene symbols
  reactomeDB <- sapply(mouse.reac, function(x) {
    index <- which(as.character(entrez.to.gene[, 1]) %in%
      x)

    return(entrez.to.gene[index, 2])
  }, USE.NAMES = TRUE)

  return(reactomeDB)
}

#' Plot Enrichment Barplot
#'
#' @param GSEA_table Table containing results from GSEA analysis
#' @param X Column with Pathways names of GSEA
#' @param Y Column with NES/ES from GSEA
#' @param pval Column containing pvalue
#' @param pval_threshold pvalue filter threshold
#' @param n_max Max number of pathways to plot
#'
#' @return
#' @export
#' @importFrom stats na.omit
#'
#' @examples
#' NULL
plot_enrichment <- function(GSEA_table, X, Y, pval = "pval", pval_threshold = 0.1,
                            n_max = 20) {
  GSEA_table <- GSEA_table %>% dplyr::filter(!duplicated(!!sym(X)))
  GSEA_tb <- tibble::as_tibble(GSEA_table) %>%
    dplyr::arrange(dplyr::desc(!!as.name(Y))) %>%
    dplyr::filter(abs(!!as.name(pval)) < pval_threshold) %>%
    dplyr::mutate(pathway := factor(!!as.name(X), levels = rev(!!as.name(X))))
  if (nrow(GSEA_tb) > n_max) {
    pos_GSEA_tb <- utils::head(GSEA_tb, n = n_max %/% 2) %>% dplyr::filter(!!as.name(Y) > 0)
    neg_GSEA_tb <- utils::tail(GSEA_tb, n = n_max %/% 2) %>% dplyr::filter(!!as.name(Y) < 0)
    GSEA_tb <- rbind(pos_GSEA_tb, neg_GSEA_tb)
  }
  ggplot2::ggplot(GSEA_tb, ggplot2::aes(x = pathway, y = !!rlang::sym(Y))) +
    ggplot2::geom_point(
      size = 3,
      ggplot2::aes(color = ifelse(!!rlang::sym(Y) > 0, "red", "blue"))
    ) +
    ggplot2::geom_segment(ggplot2::aes(x = pathway, xend = pathway, y = 0, yend = !!rlang::sym(Y))) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x =  ggplot2::element_text(size = 11, 
                                                        face = "bold"),
                   axis.text.y =  ggplot2::element_text(size = 11, 
                                                        face = "bold")) 
}

#' Return the correct org_db
#'
#' @param org_name Name of organism used in analysis
#'
#' @return
#' @export
#'
#' @examples NULL
get_org_db <- function(org_name) {
  if (org_name == "Mus musculus") {
    org_db <- org.Mm.eg.db::org.Mm.eg.db
  } else if (org_name == "Homo sapiens") {
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    stop("Species not supported. Pleasue choose Mus musculus or Homo sapiens")
  }
  org_db
}

#' Return the correct kegg name
#'
#' @param org_name Name of organism used in analysis
#'
#' @return
#' @export
#'
#' @examples NULL
get_kegg_name <- function(org_name) {
  if (org_name == "Mus musculus") {
    kegg_name <- "mmu"
  } else if (org_name == "Homo sapiens") {
    kegg_name <- "hsa"
  } else {
    stop("Species not supported. Pleasue choose Mus musculus or Homo sapiens")
  }
  kegg_name
}

#' improved code to perform a dotplot
#' 
#' @param gset GSEA result from Cluster Profiler
#' @param c_groups contrast groups vector of length 2. 
#' 
#' @return list of length two with two dotplot reults
#' @export 
#' 
#' @examples NULL
better_dotplot <- function(gset, c_groups = contrast_groups) {
  pos_gsea <- gset
  pos_gsea@result <- gset@result[gset@result$NES > 0, ]
  if (nrow(pos_gsea) >= 1) {
    num_y_entries <- nrow(pos_gsea)
    # Define the base font size
    base_font_size <- 18
    
    # Define the scaling factor (adjust as needed)
    scaling_factor <- 0.8
    
    # Calculate the new font size
    new_font_size <- base_font_size * scaling_factor / log(num_y_entries + 1)
    
    dp_pos_NES <-
      enrichplot::dotplot(
        pos_gsea,
        size = "NES",
        color = "p.adjust",
        showCategory = 50,
        title = glue::glue("Pathways enriched in {c_groups[1]}")
      ) +
      scale_size(range = c(1, 7), limits = c(1, max(gset@result$NES))) + 
      theme(axis.text.y = element_text(size = new_font_size)) 
  } else {
    dp_pos_NES <- NULL
  }
  neg_gsea <- gset
  neg_gsea@result <- gset[gset@result$NES < 0, ]
  if (nrow(neg_gsea) >= 1) {
    num_y_entries <- nrow(neg_gsea)
    # Define the base font size
    base_font_size <- 18
    
    # Define the scaling factor (adjust as needed)
    scaling_factor <- 0.95
    
    # Calculate the new font size
    new_font_size <- base_font_size * scaling_factor / log(num_y_entries + 1)
    dp_neg_NES <-
      enrichplot::dotplot(
        neg_gsea,
        size = "NES",
        color = "p.adjust",
        showCategory = 50,
        title = glue::glue("Pathways enriched in {c_groups[2]}")
      ) +
      scale_size(range = c(7, 1), limits = c(min(gset@result$NES), -1)) + 
      theme(axis.text.y = element_text(size = new_font_size)) 
  } else {
    dp_neg_NES <- NULL
  }

  return(list(dp_pos_NES, dp_neg_NES))
}
#' Title
#'
#' @param d_plot dotplot object
#' @param p_group plot_group
#' @param p_path plot_path
#' @param gsea_type Type of gsea 
#'
#' @return
#' @export
#'
#' @examples
save_dotplots <- function(d_plot, p_group, p_path = plot_path, gsea_type) {
  ggsave(filename = file.path(p_path, glue::glue("{p_group}_{gsea_type}_dplot.svg")), d_plot, width = 10, height = 10)
}

#' Get BM object via biomart in a stable fashion in case on esource is down
#'
#' @param species 
#'
#' @return
#' @export
#'
#' @examples
stable_get_bm <- function(species) {
  mart <- "www"
  rounds <- 0
  while (class(mart)[[1]] != "Mart") {
    mart <- tryCatch(
      {
        # done here, because error function does not
        # modify outer scope variables, I tried
        if (mart == "www") rounds <- rounds + 1
        # equivalent to useMart, but you can choose
        # the mirror instead of specifying a host
        biomaRt::useEnsembl(
          biomart = "ENSEMBL_MART_ENSEMBL",
          dataset = glue::glue("{species}_gene_ensembl"),
          mirror = mart,
          host = "https://nov2020.archive.ensembl.org"
        )
      },
      error = function(e) {
        # change or make configurable if you want more or
        # less rounds of tries of all the mirrors
        if (rounds >= 3) {
          stop()
        }
        # hop to next mirror
        mart <- switch(mart,
                       useast = "uswest",
                       uswest = "asia",
                       asia = "www",
                       www = {
                         # wait before starting another round through the mirrors,
                         # hoping that intermittent problems disappear
                         Sys.sleep(30)
                         "useast"
                       },
                       host = "https://nov2020.archive.ensembl.org"
        )
      }
    )
  }
  mart
}


#' ENSEMBL IDs to ENSEMBL Gene symbol
#'
#' @param ens_id_vector 
#'
#' @return
#' @export
#'
#' @examples
ensembl_to_symbol <- function(ens_id_vector) {
  ens_short <- stringr::str_extract(ens_id_vector,
                                    pattern = "^ENS[A-Z0-9]*")
  mart <- stable_get_bm("mmusculus")
  g2g <- biomaRt::getBM(
    attributes = c( "ensembl_gene_id",
                    "external_gene_name"),
    filters = "ensembl_gene_id",
    values = ens_short,
    mart = mart,  
  )
  symbol_vec <- setNames(ens_short, nm = ens_short)
  symbol_vec[g2g$ensembl_gene_id] <- g2g$external_gene_name
  
  symbol_vec
}

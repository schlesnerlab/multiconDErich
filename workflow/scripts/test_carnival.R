####

library(progeny)

library(CARNIVAL)
library(OmnipathR)
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(visNetwork)
#library(biomaRt)
library(ggplot2)
library(pheatmap)
library(BiocParallel)
library(DESeq2)
library(decoupleR)
library(RNAscripts)
# library(biomaRt) Basic function to convert human to mouse gene names human = useMart(biomart= 'ensembl', dataset =
# 'hsapiens_gene_ensembl') mouse = useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
#' keep controllable neighbors
#' 
#' keeps the nodes in network that are no more then n_steps away from the starting nodes.
#' @param network network in 3 column data.frame format
#'  (source, interaction, target)
#' @param n_steps largest distance t consider
#' @param input_nodes names of the input nodes in the network
#' @noRd
keep_controllable_neighbours <- function(network, n_steps, input_nodes)
{
  stopifnot(all(c("source","target","interaction") %in% colnames(network)))
  stopifnot(all(input_nodes %in% c(network$source,network$target)))
  
  print(paste("COSMOS: removing nodes that are not reachable from inputs within",n_steps,"steps"))
  meta_g <- igraph::graph_from_data_frame(network[,c("source",'target',"interaction")],directed = TRUE) 
  dn_nbours <- igraph::ego(graph = meta_g, order = n_steps,nodes = input_nodes, mode = "out")
  
  sub_nodes <- c(unique(names(unlist(dn_nbours))), input_nodes)
  
  to_keep = network$source %in% sub_nodes & network$target %in% sub_nodes
  
  print(paste("COSMOS:", sum(!to_keep), "from ", length(to_keep), 
              "interactions are removed from the PKN"))
  
  network <- network[to_keep,]
  
  
  return(network)
}

#' keep observable neighbors
#' 
#' keeps the nodes in network that are no more then n_steps upstreams from the 
#' measured nodes 
#' @param network network in 3 column data.frame format
#'  (source, interaction, target)
#' @param n_steps largest distance t consider
#' @param observed_nodes names of the measured nodes in the network
#' @noRd
keep_observable_neighbours <- function(network, n_steps, observed_nodes)
{
  stopifnot(all(c("source","target","interaction") %in% colnames(network)))
  stopifnot(all(observed_nodes %in% c(network$source,network$target)))
  
  print(paste("COSMOS: removing nodes that are not observable by measurements within",n_steps,"steps"))
  meta_g <- igraph::graph_from_data_frame(network[,c("source",'target',"interaction")],directed = TRUE) 
  up_nbours <- igraph::ego(graph = meta_g, order = n_steps, nodes = observed_nodes, mode = "in")
  
  sub_nodes <- c(unique(names(unlist(up_nbours))), observed_nodes)
  
  to_keep = network$source %in% sub_nodes & network$target %in% sub_nodes
  
  print(paste("COSMOS:", sum(!to_keep), "from ", length(to_keep), 
              "interactions are removed from the PKN"))
  
  network <- network[to_keep,]
  
  
  return(network)
}
mouse_human_homologs <- readr::read_tsv("http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt", col_names = c(
  "hgene",
  "hID", "mgene", "mID", "lcol"
))


if (exists("snakemake")) {
  dds_path <- snakemake@input[["dds_obj"]]
  diffexp_tb_path <- snakemake@input[["table"]]
  fpkm_path <- snakemake@input[["fpkm"]]
  carnival_output <- snakemake@output[["carnival_out"]]
  contrast_groups <- snakemake@wildcards[["contrast"]]
  s_groups <- snakemake@params[["s_groups"]]
  register(MulticoreParam(snakemake@threads))
  contrast_name <- contrast_groups
  plot_path <- snakemake@params[["plot_path"]]
  comp_groups <- snakemake@config[["signif_comp"]]
  color_scheme <- snakemake@config[["group_colors"]]
  cplex_path <- snakemake@config[["cplex_solver"]]
  stopifnot(`Cplexpath doesnt exist please give path` = file.exists(cplex_path))
  temp_path <- snakemake@params[["temp_path"]]
  nr <- snakemake@resources[["nr"]]
  thread_num <- snakemake@threads
  mem_mb <- snakemake@resources[["mem_mb"]]
  time_limit <- (snakemake@resources[["time_min"]] - 20) * 60
  run_vanilla <- snakemake@params[["run_vanilla"]]
  perturbation_gene <- snakemake@params[["perturbation_gene"]] %>% unlist()
  print(perturbation_gene)
  progeny_data <- "../data/progenyMembers.RData"
} else {
  BASE_ANALYSIS_DIR <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/apelin_exp"
  dds_path <- file.path(paste0(BASE_ANALYSIS_DIR), "deseq2/all.rds")
  diffexp_tb_path <- file.path(paste0(BASE_ANALYSIS_DIR), "results/diffexp/condition/YECplusA-vs-YEC.diffexp.tsv")
  fpkm_path <- file.path(BASE_ANALYSIS_DIR, "fpkm/all.tsv")
  contrast_groups <- c("YECplusA", "YEC")
  plot_path <- "."
  register(SerialParam())
  s_groups <- c("YECplusA", "YEC")
  contrast_name <- glue::glue("{contrast_groups[[1]]} vs {contrast_groups[[2]]}")
  the_yaml <- yaml::read_yaml("./configs/VascAge_Apelin_config.yaml")
  comp_groups <- the_yaml$comp_groups
  color_scheme <- the_yaml$group_colors
  carnival_output <- "./test_output.RDS.gz"
  cplex_path <- "/home/heyer/software/external/CPLEX_Studio201/cplex/bin/x86-64_linux/cplex"
  thread_num <- 2
  temp_path <- "./"
  mem_mb <- 8192
  time_limit <- 3600
  run_vanilla <- FALSE
  progeny_data <- "./data/progenyMembers.RData"
  perturbation_gene <- c(Apln = 1)
}

# Read DESeq2 oobject and other tables
dds_obj <- readRDS(dds_path)
diffexp_tb <- read_tsv(diffexp_tb_path,
  col_names = c("gene_id", "baseMean", "logFoldChange", "lfcSE", "stat", "pvalue", "padj"),
  skip = 1
)
# Normalized_counts <- getVarianceStabilizedData(dds_obj)
Normalized_counts <- assay(vst(dds_obj, blind = F))
fpkm <- read_tsv(fpkm_path)


filer <- fpkm %>%
  dplyr::filter(gene %in% rownames(Normalized_counts)) %>%
  dplyr::filter(!duplicated(gname))

joined_df <- join_tables(diffexp_tb, filer) %>%
  dplyr::filter(!duplicated(gname))
Normalized_counts <- Normalized_counts[filer$gene, ]

rownames(Normalized_counts) <- filer$gname

joined_df %>%
  dplyr::select(gname, stat) %>%
  dplyr::filter(!is.na(stat)) %>%
  column_to_rownames(var = "gname") %>%
  as.matrix() -> diffexp_matrix




# regulons <- dorothea_mm %>% dplyr::filter(confidence %in% c('A', 'B'))
organism <- "mouse"
#doro_net <- decoupleR::get_dorothea(organism = organism, levels = c("A", "B", "C"))
doro_net <- decoupleR::get_collectri(organism = organism)
prog_net <- decoupleR::get_progeny(organism = organism, top = 500)
prog_net <- prog_net %>% dplyr::filter(!is.na(target))


PathwayActivity <- PathwayActivity_CARNIVALinput <- decoupleR::run_mlm(
  mat = diffexp_matrix, network = prog_net, .source = "source",
  .target = "target", .mor = "weight",  minsize = 5
) %>%
  dplyr::filter(statistic == "mlm") %>% 
  dplyr::mutate(padj = p.adjust(p_value, method = "BH")) %>%
  dplyr::filter(padj < 0.05) %>%
  as.data.frame()
#if (any(abs(PathwayActivity$score) > 1)) {
#  PathwayActivity$score <- sign(PathwayActivity$score) * (1 - PathwayActivity$p_value)
#  warning("decoupler based enriched failed, falling back on (1-pvalue) * sign(score)")
#}
# PathwayActivity <- PathwayActivity_CARNIVALinput <- progeny(diffexp_matrix, scale = TRUE, organism = 'Mouse', top =
# 100, perm = 10000, z_scores = F ) %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var = 'Pathway')

# colnames(PathwayActivity)[2] <- 'score'
progeny_key <- setNames(object = PathwayActivity$score, nm = PathwayActivity$source)

prog_net %>%
  mutate(progeny_activity = recode(source, !!!progeny_key)) %>%
  mutate(carnival_score = weight * progeny_activity) -> prog_net

if (length(names(progeny_key)) > 1) {
prog_net %>% dplyr::filter(source %in% names(progeny_key)) %>%
    mutate(progeny_activity = recode(source, !!!progeny_key)) %>%
    mutate(carnival_score = weight * progeny_activity) %>%
    group_by(target) %>%
    summarise(carnival_score = mean(carnival_score)) -> prog_net_final
} else {
  prog_net %>% dplyr::filter(source %in% names(progeny_key)) %>%
    group_by(target) ->  prog_net_final
}
tf_activities_stat <- decoupleR::run_ulm(diffexp_matrix, network = doro_net,.source = "source",
                                         .target = "target", .mor = "mor",
                                         minsize = 5) %>%
  filter(statistic == "ulm")
# options = list( minsize = 5, eset.filter = FALSE, cores = 1, verbose = FALSE, nes = TRUE ) )





tf_activities <- tf_activities_CARNIVALinput <- tf_activities_stat %>% dplyr::mutate(padj = p.adjust(p_value, method = "BH")) %>%
  dplyr::filter(padj < 0.1) %>%
  dplyr::select(source, score)
prog_net_final %>%
  filter(!(target %in% tf_activities_stat$source)) -> prog_net_final

## Get Omnipath
if (run_vanilla) {
  omniR <- import_omnipath_interactions(organism = 10090) #signed and directed
  ##omniR <- import_all_interactions(organism = 10090) %>% dplyr::filter(n_references > 0)
} else {
  omniR <- import_omnipath_interactions(organism = 10090) #signed and directed
}
#omniR <- import_omnipath_interactions(organism = 10090) #signed and directed
omnipath_sd <- omniR %>%
  dplyr::filter(consensus_direction == 1 & (consensus_stimulation == 1 | consensus_inhibition == 1))

# changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
omnipath_sd$consensus_stimulation[which(omnipath_sd$consensus_stimulation == 0)] <- -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 1)] <- -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 0)] <- 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omnipath_sd[, c("source_genesymbol", "consensus_stimulation", "consensus_inhibition", "target_genesymbol")] %>%
  dplyr::filter(consensus_stimulation == consensus_inhibition) %>%
  unique.data.frame()

sif$consensus_stimulation <- NULL
colnames(sif) <- c("source", "interaction", "target")

# remove complexes
sif$source <- gsub(":", "_", sif$source)
sif$target <- gsub(":", "_", sif$target)

# dorothea for CARNIVAL
tf_activities_carnival <- data.frame(tf_activities, stringsAsFactors = F)
rownames(tf_activities_carnival) <- tf_activities$source
tf_activities_carnival$source <- NULL
tf_list <- generateTFList(tf_activities_carnival, top = 50, access_idx = 1)
tf_vec <- tf_list$score[, ]

# For mouse change trp53 to tp53 and trp63 to tp63

names(tf_vec) <- gsub("^Trp", "Tp", names(tf_vec))

# progeny for CARNIVAL load(file = progeny_data) progenyMembers$gene$p53 <- 'TP53' progmem_mouse <-
# purrr::map(progenyMembers$gene, convertHumanGeneHomologs, jax_database = mouse_human_homologs) progenyMembers$gene <-
# progmem_mouse progenyMembers$gene$p53 <- 'Tp53'


# PathwayActivity_carnival <- data.frame(PathwayActivity, stringsAsFactors = F) rownames(PathwayActivity_carnival) <-
# PathwayActivity_carnival$Pathway PathwayActivity_carnival$Pathway <- NULL progenylist <- assignPROGENyScores( progeny
# = t(PathwayActivity_carnival), progenyMembers = progenyMembers, id = 'gene', access_idx = 1 ) progeny_vec <-
# progenylist$score
prog_net_final %>% dplyr::filter(target %in% sif$source | target %in% sif$target) -> prog_net_final
progeny_vec <- setNames(prog_net_final$carnival_score, nm = prog_net_final$target)
# scale vec
progeny_vec <- scales::rescale(unlist(progeny_vec), to = c(-1,1))


lp_opts <- CARNIVAL::defaultCplexCarnivalOptions(
  solverPath = cplex_path, cplexMemoryLimit = mem_mb, threads = thread_num,
  timelimit = time_limit, lpFilename = file.path(temp_path, "lptest.lp"), outputFolder = file.path(temp_path, "carnout")
)
cplex_opts <- CARNIVAL::suggestedCplexSpecificOptions()
lp_opts[names(cplex_opts)] <- cplex_opts
# lp_opts$solverPath <- cplex_path
lp_opts$cplexMemoryLimit <- mem_mb
lp_opts$threads <- thread_num
lp_opts$timelimit <- time_limit

# lp_opts$lpFilename <- file.path(temp_path, 'lptest.lp') 
lp_opts$outputFolder <- file.path(temp_path, 'carnout')
dir.create(lp_opts$outputFolder, showWarnings = F, recursive = T)

# run carnival
dir.create(file.path(temp_path, "carnout"), recursive = T, showWarnings = F)
# setwd(file.path(temp_path, 'carnout'))
if (run_vanilla) {
  #lp_opts$mipGap <- 0.1

  lp_opts$poolrelGap <- 0.01
  pert_vec  <- perturbation_gene
  sif<- CARNIVAL::prune_network(sif[c(1,3,2)], upstream_nodes = names(pert_vec), 
  downstream_nodes = names(unlist(tf_vec)))[c(1,3,2)]
  sif <- keep_controllable_neighbours(sif, 5, names(pert_vec) )
  tf_vec <- tf_vec[names(tf_vec) %in% c(sif$source, sif$target)]
  
  
  sif <- keep_observable_neighbours(sif, 5, names(tf_vec))
  tf_vec <- tf_vec[names(tf_vec) %in% c(sif$source, sif$target)]
  progeny_vec <- progeny_vec[names(progeny_vec) %in% c(sif$source, sif$target)]
  
  # check if tf_vec names are in sif source or targets
  
  carnival_result <- runVanillaCarnival(
    perturbations = pert_vec, measurements = unlist(tf_vec), 
    priorKnowledgeNetwork = sif,
    weights = progeny_vec, 
    carnivalOptions = lp_opts
  )
 } else {
   lp_opts$poolrelGap <- 0.01
  tf_vec <- tf_vec[names(tf_vec) %in% c(sif$source, sif$target)]
  sif <- keep_observable_neighbours(sif, 7, names(tf_vec))

  tf_vec <- tf_vec[names(tf_vec) %in% c(sif$source, sif$target)]
  progeny_vec <- progeny_vec[names(progeny_vec) %in% c(sif$source, sif$target)]
  carnival_result <- runInverseCarnival(
    measurements = unlist(tf_vec), priorKnowledgeNetwork = sif, weights = unlist(progeny_vec),
    carnivalOptions = lp_opts
  )
}
if (length(carnival_result$sifAll) < 50) {
  if (nr > 2) {
    write(paste0("Failed to converge after 2 attempts"), file = stderr())
    carnival_result <- list()
  } else {
    stop(glue::glue("Attempt {nr}. Failed to converge restarting"))
  }
}
saveRDS(carnival_result, file = carnival_output)

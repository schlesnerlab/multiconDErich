
library(ComplexHeatmap)
library(circlize)
library(glmmSeq)
library(seriation)
library(DESeq2)
register_DendSer()
if (exists("snakemake")) {
    glmmseq_obj <- snakemake@input[["glmmseq_obj"]]
    corrected_counts <- snakemake@input[["batch_corrected_counts"]]
    png_file <- snakemake@output[["png_file"]]
    var <- snakemake@wildcards[["var"]]
    coef_id <- snakemake@params[["coef"]]
    col_data_to_plot <- snakemake@config[["pca"]][["labels"]] |> unlist()
} else {
    glmmseq_obj <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/glmmseq/glmmseq/glmmseq_obj.rds.gz"
    corrected_counts <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/glmmseq/counts/batch_corrected_counts.rds"
    coef_id <- "EC_statustumor:Aplnr_KOKO"
    var <- "EC_status:Aplnr_KO"
    col_data_to_plot <- c("EC_status", "Aplnr_KO")
}
annotation_colors <- list(
  age = c(
    aged = "#A00101",       # Tomato
    young = "grey"       # SteelBlue
  ),
  EC_status = c(
    healthy = "#98FB98",    # PaleGreen
    tumor = "#FF69B4"       # HotPink
  ),
  experiment = c(
    apelin_2020 = "#FFD700",   # Gold
    APLNR_KO = "#FF8C00",      # DarkOrange
    cre_2022 = "#8A2BE2",      # BlueViolet
    Vasc_age2020 = "#7FFF00",  # Chartreuse
    tumor_vs_ec = "#DC143C"    # Crimson
  ),
  Aplnr_KO = c(
    KO = "blue",
    normal = "grey"
  ),
  Apln_treatment = c(
    up = "#77DD77",
    normal = "grey"
  )
)

plot_heatmap<- function(coef,vst_obj = vst_dds, glmm_obj =  glmmseq_norm_counts,
                        z_score = TRUE, use_vst = FALSE, 
                        ann_col = annotation_colors,
                        coldata_to_plot = c("age", "EC_status","Apln_treatment", "Aplnr_KO", "experiment"),
                        n_genes = 30, qval_filt = 0.01, coef_filt = 0.8, meanExp_cutoff = 7,
                        row_km = 3, col_km = NULL) {
    # Get qvals and coefs for one coeff from glmm_obj
    qvals <- glmm_obj@stats$qvals[,coef[1]]
    coef_values <- glmm_obj@stats$coef[,coef[2]]
    meanexp <- glmm_obj@stats$res[,"meanExp"]
    
    stopifnot(all.equal(names(coef_values) , names(qvals)))
    # get vst normalized expression data 
    vst_data <- assay(vst_obj)

    # retain only genes with qvals < 0.05 and ceef_value > abs(0.5)o
    selected_genes <- which(qvals < 0.05 )
    if (length(selected_genes) > 200) {
    selected_genes <- which(qvals < qval_filt & abs(coef_values) > coef_filt &
                              meanexp > meanExp_cutoff)
    }
    length(selected_genes)
    
    if (use_vst) {
      plot_data <- vst_data#[names(selected_genes),]
      #plot_data <- t(scale(t(plot_data)))
      coldata <- colData(vst_obj)
    } else {
      plot_data <- glmm_obj@countdata[selected_genes,]
      #plot_data <- t(scale(t(plot_data)))

      coldata <- glmm_obj@metadata[,c("age", "EC_status", "experiment")]
    }
    plot_data <- plot_data[names(selected_genes),]
    if (z_score) {
      plot_data <- t(scale(t(plot_data)))
      col_vec <- circlize::colorRamp2(c(-2, 0, 2),  
                                      c("blue", "white", "red")
                                                                     
                                      )
    
    } else {
      col_vec <- circlize::colorRamp2(c(min(plot_data), median(plot_data), quq(plot_data)),
                                      c("blue", "white", "red"))
    }

    # Create Column Heatmapannotation from coldata(vst_dds
    cluster_data <- vst_data[selected_genes, ]
    #row_clust <- hclust(dist(cluster_data),method = "average")
    #col_clust <- hclust(dist(t(cluster_data)), method = "average")
    
    top_anno <- HeatmapAnnotation(df = coldata[, coldata_to_plot],
                                  col = ann_col)
    
    print(dim(plot_data))
    print(dim(cluster_data))
    
    o1 = seriate(dist(plot_data), method = "DendSer")
    o2 = seriate(dist(t(plot_data)), method = "DendSer")
    if (!is.null(col_km)) {
      col_dend <- dendextend::color_branches(as.dendrogram(o2[[1]]), col_km)
    } else {
      col_dend <- as.dendrogram(o2[[1]])
    }
    hmap <- ComplexHeatmap::Heatmap(
      plot_data,
      show_row_names = nrow(plot_data) <= 50,
      top_annotation = top_anno,
      show_column_names = FALSE,
   #   cluster_rows = as.dendrogram(o1[[1]]),
      cluster_columns = col_dend,
      name = "z-scaled expression",
      cluster_rows = dendextend::color_branches(as.dendrogram(o1[[1]]), row_km),
      row_split = row_km,
      column_split = col_km
      # cluster methods for rows and columns
      #     clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
      #      clustering_method_columns = 'ward.D2',
      #      clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
      #     clustering_method_rows = 'ward.D2', col = col_vec,
    )
    # Create row annotation for coefficient values
    anno_df <- data.frame(coef =  coef_values[selected_genes] )
    colnames(anno_df) <- coef[1]
    col_list <- list()
    col_list[[coef[1]]] <-colorRamp2(c(min(coef_values[selected_genes]), 0, 
                                       max(coef_values[selected_genes])), c("blue", "white", "red"))
    coef_values_anno <- rowAnnotation(
      df=anno_df,
      col = col_list,
      annotation_name_rot = 45
    )

    # Create a legend for the coefficient values
    coef_legend <- Legend(
    title = coef[1],
    col_fun = colorRamp2(c(min(coef_values[selected_genes]), 0, max(coef_values[selected_genes])), c("blue", "white", "red")),
     at = c(min(coef_values[selected_genes]), 0, max(coef_values[selected_genes])),
    labels = c(round(min(coef_values[selected_genes]), 2), 0, round(max(coef_values[selected_genes]), 2))
    )

    if (nrow(plot_data) > 1) {
      genes_to_plot <- names(head(sort(abs(coef_values[selected_genes]), decreasing = T),n_genes))
      names_ordered <- rownames(plot_data)[rownames(plot_data) %in% genes_to_plot]
      at_vec <- which(rownames(plot_data) %in% names_ordered)
        genelabels <- rowAnnotation(
        Genes = anno_mark(
        at = at_vec,
        labels = names_ordered,
        labels_gp = gpar(fontsize = 10, fontface = 'bold')),
        width = unit(2.0, 'cm') +
      max_text_width(
        rownames(plot_data)[seq(1, nrow(plot_data), 20)],
        gp = gpar(fontsize = 10,  fontface = 'bold')))

      png(png_file, width = 10, height = 10, units = 'in', res = 300) 
      draw(hmap + genelabels + coef_values_anno,
      heatmap_legend_side = 'left',
      annotation_legend_side = 'right',
      row_sub_title_side = 'left')
      #annotation_legend_list = list(coef_legend))
      dev.off()
    } else {
      png(png_file, width = 10, height = 10, units = 'in', res = 300)
      draw(hmap + coef_values_anno, annotation_legend_list = list(coef_legend))
      dev.off()
    }
}

# Read_glmmseq_obj
glmmseq_obj <- readRDS(glmmseq_obj)
glmmseq_obj$norm_counts <- glmmQvals(glmmseq_obj$norm_counts)
# Read DDs Object 
vst_dds <- readRDS(corrected_counts)
# run_vst
vst_dds <- vst(vst_dds[rownames(glmmseq_obj$norm_counts@countdata),])


plot_heatmap(
  coef = c(var, coef_id),
  vst_obj = vst_dds,
  glmm_obj = glmmseq_obj$norm_counts,
  z_score = TRUE,
  use_vst = TRUE,
  coldata_to_plot = col_data_to_plot,
  n_genes = 30,
  qval_filt = 0.01,
  coef_filt = 0.8,
  meanExp_cutoff = 7,
  row_km = 3,
  col_km = NULL
)

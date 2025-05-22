library(glmmSeq)
library(DESeq2)
library(furrr)
library(snow)
# Insert snakemake boilerplate
if(exists("snakemake")) {
    uncorrected_counts <- snakemake@input[["uncorrected_counts"]]
    the_formula <- snakemake@params[["formula"]] |> as.formula()
    glmmseq_obj <- snakemake@output[["glmmseq_obj"]]
    glmmseq_refit <- snakemake@output[["glmmseq_refit"]]
   # reference_groups <- snakemake@config[["glmmseq"]][["reference_group"]]
    threads <- snakemake@threads
    plan(multicore, workers = snakemake@threads)
} else {
    the_formula <- as.formula("~ age * EC_status + Apln_treatment:age + 
                        Aplnr_KO + EC_status:Aplnr_KO + Apln_treatment +  (1 | condition) + experiment")
    batch_corrected_counts <-  "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/glmmseq/counts/batch_corrected_counts.rds"
    
}
deseq_obj <- readRDS(uncorrected_counts)
#deseq_obj <- DESeq(deseq_obj)
normalized_counts <- counts(deseq_obj, normalized = T)
size_factors <- sizeFactors(deseq_obj)

cpm_filter <- apply(edgeR::cpm(counts(deseq_obj, normalized = T)), 1, function(x) {
                        if(length(which(x > 0.5)) > 0.5 * length(x)) {
                            val <- 1
                        } else {
                            val <- 0 
                        }
                        as.logical(val)
})
dispersions <- setNames(dispersions(deseq_obj), rownames(deseq_obj))
dispersions <- dispersions[cpm_filter]
normalized_counts <- normalized_counts[cpm_filter,]
deseq_obj <- deseq_obj[cpm_filter,]

vst_data <- vst(deseq_obj, blind = F) |> assay()
count_var <- apply(vst_data, 1, var) |> sort()
count_var <- tail(count_var, n = length(count_var)-0.2*length(count_var))
normalized_counts <- normalized_counts[names(count_var),]
deseq_obj <- deseq_obj[names(count_var),]
rm(vst_data)



dispersions <- dispersions[names(count_var)]
met_data <- colData(deseq_obj) |> as.data.frame()

#for (var_name in names(reference_groups)) {
#    sample_mat[[var_name]] <- relevel(sample_mat[[var_name]], ref = reference_groups[[var_name]])
#}
#write(file = stderr(), )
# Run glmmseq
write(file = stderr(), paste0("The formula is: ", the_formula, "\n"))
glmmseq_norm_counts <- glmmSeq(the_formula,
                        countdata =  normalized_counts[,], 
                         metadata = met_data,
                         dispersion = dispersions, 
                         sizeFactors = size_factors,progress = F,
                         cores = threads
)
gc()
saveRDS(list(norm_counts = glmmseq_norm_counts
              ), glmmseq_obj )

refits <- furrr::future_map(rownames(glmmseq_norm_counts@stats$res), 
    glmmRefit, object = glmmseq_norm_counts)
names(refits) <- rownames(glmmseq_norm_counts@stats$res)
saveRDS(refits, glmmseq_refit)






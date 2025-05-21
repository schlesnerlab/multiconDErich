rule batch_correct:
    input:
        counts=get_count_matrix,
    output:
        batch_corrected_counts=join(
            BASE_ANALYSIS_DIR, "counts/batch_corrected_counts.rds"
        ),
        uncorrected_counts=join(BASE_ANALYSIS_DIR, "counts/uncorrected_counts.rds"),
    params:
        samples=config["samples"],
        batch_correct_model=config["glmmseq"]["batch_correct"],
        batch_variable="experiment",
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/glmmseq/batch_correct.log",
    resources:
        mem_mb=8192,
        runtime=59,
    threads: 4
    script:
        "../scripts/glmmseq/batch_correct.R"


rule plot_PCA:
    input:
        batch_corrected_counts=join(
            BASE_ANALYSIS_DIR, "counts/batch_corrected_counts.rds"
        ),
        uncorrected_counts=join(
            BASE_ANALYSIS_DIR, "counts/uncorrected_counts.rds"
        )
    output:
        html_file=join(BASE_ANALYSIS_DIR, "glmmseq/pca.html"),
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/glmmseq/plot_PCA.log",
    resources:
        mem_mb=8192,
        time_min=59,
    threads: 1
    script:
        "../scripts/glmmseq/plot_PCA.Rmd"
rule run_glmmseq:
    input:
        uncorrected_counts=join(
            BASE_ANALYSIS_DIR, "counts/uncorrected_counts.rds"
        ),
    output:
        glmmseq_obj=join(BASE_ANALYSIS_DIR, "glmmseq/glmmseq_obj.rds.gz"),
        glmmseq_refit=join(BASE_ANALYSIS_DIR, "glmmseq/glmmseq_refit.rds.gz"),
    params:
        formula=config["glmmseq"]["formula"],
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/glmmseq/run_glmmseq.log",
    resources:
        mem_mb=16384,
        time_min=59,
    threads: 12
    script:
        "../scripts/glmmseq/glmmseq.R"


rule run_glmmseq_qc:
    input:
        glmmseq_obj=join(BASE_ANALYSIS_DIR, "glmmseq/glmmseq_obj.rds.gz"),
        glmmseq_refit=join(BASE_ANALYSIS_DIR, "glmmseq/glmmseq_refit.rds.gz"),
    output:
        rmd_script=join(BASE_ANALYSIS_DIR, "glmmseq/qc.html"),
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/glmmseq/glmmseq_qc.log",
    resources:
        mem_mb=8192,
        time_min=59,
    threads: 4
    script:
        "../scripts/glmmseq/glmmseq_qc.Rmd"


def get_coef_name(wildcards):
    return config["glmmseq"]["test_group"][wildcards.var]


rule run_glmmseq_enrichments:
    input:
        glmmseq_obj=join(BASE_ANALYSIS_DIR, "glmmseq/glmmseq_obj.rds.gz"),
    output:
        enrichment_obj=join(BASE_ANALYSIS_DIR, "glmmseq/{var}_enrichment_obj.rds.gz"),
    params:
        enrichments=config["glmmseq"]["enrichments"],
        coef_name=get_coef_name,
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/glmmseq/{var}run_glmmseq_enrichments.log",
    resources:
        mem_mb=2 * 8192,
        time_min=59,
    threads: 1
    script:
        "../scripts/glmmseq/glmmseq_enrichments.R"


rule visualize_enrichments:
    input:
        enrichment_obj=expand(
            join(BASE_ANALYSIS_DIR, "glmmseq/{var}_enrichment_obj.rds.gz"),
            var=config["glmmseq"]["test_group"].keys(),
        ),
        glmmseq_obj=join(BASE_ANALYSIS_DIR, "glmmseq/glmmseq_obj.rds.gz"),
    output:
        html_file=join(BASE_ANALYSIS_DIR, "glmmseq/enrichment.html"),
    params:
        coefs=config["glmmseq"]["test_group"],
        enrich_vars=config["glmmseq"]["enrichments"],
    #   coef = get_coef_name
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/glmmseq/visualize_enrichments.log",
    resources:
        mem_mb=8192,
        time_min=59,
    threads: 1
    script:
        "../scripts/glmmseq/visualize_enrichments.Rmd"


rule glmmseq_heatmap:
    input:
        glmmseq_obj=join(BASE_ANALYSIS_DIR, "glmmseq/glmmseq_obj.rds.gz"),
        batch_corrected_counts=join(
            BASE_ANALYSIS_DIR, "counts/batch_corrected_counts.rds"
        ),
    output:
        png_file=join(BASE_ANALYSIS_DIR, "glmmseq/{var}_heatmaps.png"),
    params:
        coef=get_coef_name,
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/glmmseq/{var}_glmmseq_heatmaps.log",
    resources:
        mem_mb=8192,
        time_min=59,
    threads: 1
    script:
        "../scripts/glmmseq/glmmseq_heatmaps.R"

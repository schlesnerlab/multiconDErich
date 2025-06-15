
rule deseq2_init:
    input:
        counts=get_count_matrix,
    output:
        join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    resources:
        mem_mb=8192,
        time_min=59,
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        join(BASE_ANALYSIS_DIR, "deseq2/rlog_transform.RDS.gz"),
    output:
        report(join(BASE_ANALYSIS_DIR, "results/pca.svg"), "../report/pca.rst"),
    params:
        pca_labels=config["pca"]["labels"],
    conda:
        "../envs/deseq2.yaml"
    resources:
        mem_mb=8192,
    log:
        "logs/pca.log",
    script:
        "../scripts/plot-pca.R"


rule deseq2:
    input:
        join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
    output:
        table=report(
            join(
                BASE_ANALYSIS_DIR, "results/diffexp/{condition}/{contrast}.diffexp.tsv"
            ),
            "../report/diffexp.rst",
        ),
        ma_plot=report(
            join(
                BASE_ANALYSIS_DIR, "results/diffexp/{condition}/{contrast}.ma-plot.pdf"
            ),
            "../report/ma.rst",
        ),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    resources:
        mem_mb=8192,
        time_min=29,
    log:
        "logs/deseq2/{condition}/{contrast}.diffexp.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2.R"


rule rlog_transform:
    input:
        dds_obj=join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
    output:
        rld=join(BASE_ANALYSIS_DIR, "deseq2/rlog_transform.RDS.gz"),
        fpkm=join(BASE_ANALYSIS_DIR, "fpkm/all.tsv"),
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/rlog_trans.report.log",
    resources:
        mem_mb=8192,
        time_min=29,
    threads: 4
    script:
        "../scripts/rlog_transform.R"


rule deseq_report:
    input:
        dds_obj=join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
        table=join(
            BASE_ANALYSIS_DIR, "results/diffexp/{condition}/{contrast}.diffexp.tsv"
        ),
        fpkm=get_fpkm,
        featureCounts=get_count_matrix,
        rld=join(BASE_ANALYSIS_DIR, "deseq2/rlog_transform.RDS.gz"),
        gsea_result=join(
            BASE_ANALYSIS_DIR, "results/diffexp/{condition}/{contrast}.gseares.RDS"
        ),
    output:
        join(BASE_ANALYSIS_DIR, "reports/deseq2/{condition}/{contrast}_diffexp.html"),
    params:
        contrast=get_contrast,
        samples=config["samples"],
        plot_path=join(BASE_ANALYSIS_DIR, "reports/deseq2/{condition}/{contrast}_plots"),
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/deseq2/{condition}/{contrast}.report.log",
    threads: 1
    resources:
        time_min=60,
        mem_mb=8192,
        queue="medium-debian"
    script:
        "../scripts/DESeq2_analysis.Rmd"


rule run_gsea:
    input:
        table=join(
            BASE_ANALYSIS_DIR, "results/diffexp/{condition}/{contrast}.diffexp.tsv"
        ),
        fpkm_path=join(BASE_ANALYSIS_DIR, "fpkm/all.tsv"),
    output:
        gsea_result=join(
            BASE_ANALYSIS_DIR, "results/diffexp/{condition}/{contrast}.gseares.RDS"
        ),
    params:
        contrast=get_contrast,
        gsea_use_stat=config["diffexp"]["gsea_use_stat"],
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/run_gsea/{condition}_{contrast}.log",
    threads: 4
    resources:
        time_min=59 * 4,
        mem_mb=8192 * 7,
    script:  #
        "../scripts/run_gsea.R"


rule gsea_report:
    input:
        gsea_result=get_gsea_results,
    output:
        join(BASE_ANALYSIS_DIR, "reports/deseq2/{condition}_joint_gsea_report.html"),
    params:
        contrast_groups=get_all_contrasts,
        plot_path=directory(join(BASE_ANALYSIS_DIR, "reports/deseq2/{condition}_plots")),
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/gsea_report/{condition}.log",
    threads: 1
    resources:
        time_min=59 * 20,
        mem_mb=8192 * 6,
    script:
        "../scripts/RMD_scripts/gsea_report.Rmd"


rule cohort_wide_comparison:
    input:
        tables=get_diffexp_tables,
        fpkm=join(BASE_ANALYSIS_DIR, "fpkm/all.tsv"),
        dds_obj=join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
        rld=join(BASE_ANALYSIS_DIR, "deseq2/rlog_transform.RDS.gz"),
    output:
        join(
            BASE_ANALYSIS_DIR, "reports/deseq2/{condition}_cohort_wide_comparison.html"
        ),
    params:
        contrast=get_all_contrasts,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{condition}_cohort.report.log",
    threads: 1
    resources:
        mem_mb=16384,
        time_min=59
    script:
        "../scripts/DEseq2_cohort.Rmd"


rule run_mitch:
    input:
        tables=get_diffexp_tables,
        fpkm=join(BASE_ANALYSIS_DIR, "fpkm/all.tsv"),
    output:
        mitch_table=join(
            BASE_ANALYSIS_DIR, "results/diffexp/mitch/{condition}_mitch_table.tsv"
        ),
        mitch_rds=join(
            BASE_ANALYSIS_DIR, "results/diffexp/mitch/{condition}_mitch_data.rds"
        ),
        mitch_report=join(
            BASE_ANALYSIS_DIR, "results/diffexp/mitch/{condition}_mitch_report.html"
        ),
    params:
        contrasts=get_all_contrasts,
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/mitch/{condition}_mitch_run.log",
    threads: 4
    resources:
        mem_mb=10240,
        time_min=59 * 5,
    script:
        "../scripts/run_mitch.R"


rule export_diffexp_xlsx:
    input:
        table=get_diffexp_tables,
        fpkm=get_fpkm,
    output:
        outpath=join(
            BASE_ANALYSIS_DIR,
            "results",
            "diffexp",
            "excel_tables",
            "{condition}_diffexp_genes.xlsx",
        ),
    params:
        samp_map=config["samples"],
        contrast_groups=get_all_contrasts,
    conda:
        "../envs/R_4.yaml"
    resources:
        mem_mb=8192,
    log:
        "logs/deseq2/{condition}_export_diffexp__xlsx.log",
    threads: 1
    script:
        "../scripts/export_diffexp_xlsx.R"


rule run_dorothea:
    input:
        dds_path=join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
        fpkm_path=join(BASE_ANALYSIS_DIR, "fpkm/all.tsv"),
    output:
        dorothea_fp=join(BASE_ANALYSIS_DIR, "dorothea/dorothea_results.html"),
    params:
        samp_map=config["samples"],
        plot_path=join(BASE_ANALYSIS_DIR, "dorothea/svg_plots"),
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/run_dorothea/run_dorothea.log",
    resources:
        mem_mb=8192,
    threads: 1
    script:
        "../scripts/dorothea.Rmd"


rule run_dorothea_diffexp:
    input:
        dds_obj=join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
        table=join(
            BASE_ANALYSIS_DIR, "results/diffexp/{condition}/{contrast}.diffexp.tsv"
        ),
        fpkm=join(BASE_ANALYSIS_DIR, "fpkm/all.tsv"),
    output:
        progeny_res=join(
            BASE_ANALYSIS_DIR, "dorothea/{condition}/{contrast}_dorothea.html"
        ),
    params:
        s_groups=get_all_conditions,
        plot_path=join(BASE_ANALYSIS_DIR, "dorothea/svg_{condition}/{contrast}"),
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/run_doroteh_diffexp/{condition}/{contrast}_doroteha.log",
    threads: 1
    resources:
        mem_mb=8192,
        time_min=59,
    script:
        "../scripts/RMD_scripts/dorothea_diffexp.Rmd"


rule run_progeny:
    input:
        dds_obj=join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
        table=join(
            BASE_ANALYSIS_DIR, "results/diffexp/{condition}/{contrast}.diffexp.tsv"
        ),
        fpkm=join(BASE_ANALYSIS_DIR, "fpkm/all.tsv"),
    output:
        progeny_res=join(
            BASE_ANALYSIS_DIR, "progeny/{condition}/{contrast}_progeny.html"
        ),
    params:
        plot_path=join(BASE_ANALYSIS_DIR, "progeny/svg_{condition}/{contrast}"),
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/run_progeny/{condition}/{contrast}_progeny.log",
    threads: 1
    resources:
        mem_mb=8192,
        time_min=59,
        queue="medium-debian"
    script:
        "../scripts/RMD_scripts/progeny_analysis.Rmd"

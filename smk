Using profile /opt/resources/apps/snakemake/slurm for setting default command line arguments.
Building DAG of jobs...
Job stats:
job                             count
----------------------------  -------
all                                 1
create_predictions                  1
filter_predictions                  1
generate_qc_plot_and_summary        1
total                               4


[Tue Mar  4 15:35:12 2025]
rule create_predictions:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Neighborhoods/EnhancerList.txt, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Neighborhoods/GeneList.txt
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz
    jobid: 1
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz
    wildcards: biosample=SPC
    resources: mem_mb=20000, mem_mib=19074, disk_mb=1000, disk_mib=954, tmpdir=<TBD>


		python /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/workflow/scripts/predict.py 			--enhancers /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Neighborhoods/EnhancerList.txt 			--outdir /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions 			--score_column ABC.Score 			--chrom_sizes /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv 			--accessibility_feature ATAC 			--cellType SPC 			--genes /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Neighborhoods/GeneList.txt 			--hic_gamma 1.024238616787792 			--hic_scale 5.9594510043736655 			--hic_pseudocount_distance 5000 			--hic_file /home/jl9324/Sperm_3D_Genome_Project/wo_trimming_fanc/cool/SPC_5kb.hic --hic_type hic --hic_resolution 5000 			--scale_hic_using_powerlaw
		

[Tue Mar  4 15:35:12 2025]
rule filter_predictions:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.bedpe.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictions_threshold0.021_self_promoter.tsv, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/GenePredictionStats_threshold0.021_self_promoter.tsv
    jobid: 8
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv; Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz
    wildcards: biosample=SPC, threshold=0.021, separator=, other_flags=_self_promoter
    resources: mem_mb=<TBD>, disk_mb=<TBD>, tmpdir=<TBD>


		python /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/workflow/scripts/filter_predictions.py 			--output_tsv_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv 			--output_slim_tsv_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictions_threshold0.021_self_promoter.tsv 			--output_bed_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.bedpe.gz 			--output_gene_stats_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/GenePredictionStats_threshold0.021_self_promoter.tsv 			--pred_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutative.tsv.gz 			--pred_nonexpressed_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz 			--score_column ABC.Score 			--threshold 0.021 			--include_self_promoter True 			--only_expressed_genes False
		

[Tue Mar  4 15:35:12 2025]
rule generate_qc_plot_and_summary:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Neighborhoods, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv, /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Metrics/QCSummary_threshold0.021_self_promoter.tsv, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Metrics/QCPlots_threshold0.021_self_promoter.pdf
    jobid: 7
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Metrics/QCSummary_threshold0.021_self_promoter.tsv; Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv
    wildcards: biosample=SPC, threshold=0.021, separator=, other_flags=_self_promoter
    resources: mem_mb=<TBD>, disk_mb=<TBD>, tmpdir=<TBD>


		python /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/workflow/scripts/grabMetrics.py 			--outdir /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Metrics 			--output_qc_summary /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Metrics/QCSummary_threshold0.021_self_promoter.tsv 			--output_qc_plots /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Metrics/QCPlots_threshold0.021_self_promoter.pdf 			--macs_peaks /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed 			--neighborhood_outdir /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Neighborhoods 			--preds_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv 			--chrom_sizes /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv 			--hic_gamma 1.024238616787792 			--hic_scale 5.9594510043736655 
		

[Tue Mar  4 15:35:12 2025]
localrule all:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Metrics/QCSummary_threshold0.021_self_promoter.tsv
    jobid: 0
    reason: Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Metrics/QCSummary_threshold0.021_self_promoter.tsv, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPC/Predictions/EnhancerPredictionsAllPutative.tsv.gz
    resources: mem_mb=<TBD>, disk_mb=<TBD>, tmpdir=/home/jl9324/mnt/scratch/jl9324

Job stats:
job                             count
----------------------------  -------
all                                 1
create_predictions                  1
filter_predictions                  1
generate_qc_plot_and_summary        1
total                               4

Reasons:
    (check individual jobs above for details)
    input files updated by another job:
        all, filter_predictions, generate_qc_plot_and_summary
    missing output files:
        create_predictions, filter_predictions, generate_qc_plot_and_summary

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.

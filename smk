Using profile /opt/resources/apps/snakemake/slurm for setting default command line arguments.
Building DAG of jobs...
Job stats:
job                              count
-----------------------------  -------
all                                  1
call_macs_peaks                      1
create_neighborhoods                 1
create_predictions                   1
filter_predictions                   1
generate_chrom_sizes_bed_file        1
generate_qc_plot_and_summary         1
make_candidate_regions               1
sort_narrowpeaks                     1
total                                9


[Tue Mar  4 15:02:29 2025]
rule call_macs_peaks:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/processed/bigwig/Spermatocytes.tagAlign.gz
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak
    jobid: 5
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak
    wildcards: biosample=SPT
    resources: mem_mb=21073, mem_mib=20097, disk_mb=1317, disk_mib=1256, tmpdir=<TBD>


		if [[ "/home/jl9324/Sperm_3D_Genome_Project/multiome/processed/bigwig/Spermatocytes.tagAlign.gz" == *tagAlign* ]]; then
			FORMAT="BED"
		else
			FORMAT="AUTO"
		fi

		macs2 callpeak 		-f $FORMAT 		-g hs 		-p 0.1 		-n macs2 		--shift -75 		--extsize 150 		--nomodel 		--keep-dup all 		--call-summits 		--outdir /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model//SPT/Peaks 		-t /home/jl9324/Sperm_3D_Genome_Project/multiome/processed/bigwig/Spermatocytes.tagAlign.gz 
		

[Tue Mar  4 15:02:29 2025]
rule generate_chrom_sizes_bed_file:
    input: /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
    jobid: 6
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
    resources: mem_mb=8000, mem_mib=7630, disk_mb=1000, disk_mib=954, tmpdir=<TBD>


		awk 'BEGIN {OFS="	"} {if (NF > 0) print $1,"0",$2 ; else print $0}' /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv > /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
		

[Tue Mar  4 15:02:29 2025]
rule sort_narrowpeaks:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted
    jobid: 4
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted; Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
    wildcards: biosample=SPT
    resources: mem_mb=<TBD>, disk_mb=<TBD>, tmpdir=<TBD>


		# intersect to remove alternate chromosomes
		bedtools intersect -u -a /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak -b /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed | 		bedtools sort -faidx /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv -i stdin > /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted
		

[Tue Mar  4 15:02:29 2025]
rule make_candidate_regions:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted, /home/jl9324/Sperm_3D_Genome_Project/multiome/processed/bigwig/Spermatocytes.tagAlign.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
    jobid: 3
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed; Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
    wildcards: biosample=SPT
    resources: mem_mb=<TBD>, disk_mb=<TBD>, tmpdir=<TBD>


		python /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/workflow/scripts/makeCandidateRegions.py 			--narrowPeak /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted			--accessibility /home/jl9324/Sperm_3D_Genome_Project/multiome/processed/bigwig/Spermatocytes.tagAlign.gz 			--outDir /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks 			--chrom_sizes /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv 			--chrom_sizes_bed /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed 			--regions_blocklist /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_unified_blacklist.bed 			--regions_includelist /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/CollapsedGeneBounds.hg38.TSS500bp.bed 			--peakExtendFromSummit 250 			--nStrongestPeak 150000
		

[Tue Mar  4 15:02:29 2025]
rule create_neighborhoods:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/EnhancerList.txt, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/GeneList.txt, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/processed_genes_file.bed
    jobid: 2
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/EnhancerList.txt, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/GeneList.txt; Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed
    wildcards: biosample=SPT
    resources: mem_mb=32000, mem_mib=30518, disk_mb=<TBD>, tmpdir=<TBD>


		# get sorted & unique gene list
		# intersect first to remove alternate chromosomes
		bedtools intersect -u -a /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/CollapsedGeneBounds.hg38.bed -b /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed | 		bedtools sort -faidx /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv -i stdin | 		uniq > /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/processed_genes_file.bed
						
		python /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/workflow/scripts/run.neighborhoods.py 			--candidate_enhancer_regions /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed 			--DHS  			--ATAC /home/jl9324/Sperm_3D_Genome_Project/multiome/processed/bigwig/Spermatocytes.tagAlign.gz 			--default_accessibility_feature ATAC 			--chrom_sizes /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv 			--chrom_sizes_bed /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/tmp/GRCh38_EBV.no_alt.chrom.sizes.tsv.bed 			--outdir /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods 			--genes /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/processed_genes_file.bed 			--ubiquitously_expressed_genes /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/UbiquitouslyExpressedGenes.txt 			--H3K27ac  			--qnorm /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/EnhancersQNormRef.K562.txt 
		

[Tue Mar  4 15:02:29 2025]
rule create_predictions:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/EnhancerList.txt, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/GeneList.txt
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz
    jobid: 1
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz; Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/EnhancerList.txt, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/GeneList.txt
    wildcards: biosample=SPT
    resources: mem_mb=<TBD>, disk_mb=<TBD>, tmpdir=<TBD>


		python /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/workflow/scripts/predict.py 			--enhancers /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/EnhancerList.txt 			--outdir /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions 			--score_column ABC.Score 			--chrom_sizes /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv 			--accessibility_feature ATAC 			--cellType SPT 			--genes /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods/GeneList.txt 			--hic_gamma 1.024238616787792 			--hic_scale 5.9594510043736655 			--hic_pseudocount_distance 5000 			--hic_file /home/jl9324/Sperm_3D_Genome_Project/wo_trimming_fanc/cool/SPC_5kb.hic --hic_type hic --hic_resolution 5000 			--scale_hic_using_powerlaw
		

[Tue Mar  4 15:02:29 2025]
rule filter_predictions:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.bedpe.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictions_threshold0.021_self_promoter.tsv, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/GenePredictionStats_threshold0.021_self_promoter.tsv
    jobid: 8
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv; Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz
    wildcards: biosample=SPT, threshold=0.021, separator=, other_flags=_self_promoter
    resources: mem_mb=<TBD>, disk_mb=<TBD>, tmpdir=<TBD>


		python /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/workflow/scripts/filter_predictions.py 			--output_tsv_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv 			--output_slim_tsv_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictions_threshold0.021_self_promoter.tsv 			--output_bed_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.bedpe.gz 			--output_gene_stats_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/GenePredictionStats_threshold0.021_self_promoter.tsv 			--pred_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutative.tsv.gz 			--pred_nonexpressed_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz 			--score_column ABC.Score 			--threshold 0.021 			--include_self_promoter True 			--only_expressed_genes False
		

[Tue Mar  4 15:02:29 2025]
rule generate_qc_plot_and_summary:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv, /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv
    output: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Metrics/QCSummary_threshold0.021_self_promoter.tsv, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Metrics/QCPlots_threshold0.021_self_promoter.pdf
    jobid: 7
    reason: Missing output files: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Metrics/QCSummary_threshold0.021_self_promoter.tsv; Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv
    wildcards: biosample=SPT, threshold=0.021, separator=, other_flags=_self_promoter
    resources: mem_mb=<TBD>, disk_mb=<TBD>, tmpdir=<TBD>


		python /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/workflow/scripts/grabMetrics.py 			--outdir /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Metrics 			--output_qc_summary /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Metrics/QCSummary_threshold0.021_self_promoter.tsv 			--output_qc_plots /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Metrics/QCPlots_threshold0.021_self_promoter.pdf 			--macs_peaks /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed 			--neighborhood_outdir /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Neighborhoods 			--preds_file /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsFull_threshold0.021_self_promoter.tsv 			--chrom_sizes /home/jl9324/env/download/ABC-Enhancer-Gene-Prediction/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv 			--hic_gamma 1.024238616787792 			--hic_scale 5.9594510043736655 
		

[Tue Mar  4 15:02:29 2025]
localrule all:
    input: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutative.tsv.gz, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Metrics/QCSummary_threshold0.021_self_promoter.tsv
    jobid: 0
    reason: Input files updated by another job: /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Metrics/QCSummary_threshold0.021_self_promoter.tsv, /home/jl9324/Sperm_3D_Genome_Project/multiome/ABC_model/SPT/Predictions/EnhancerPredictionsAllPutative.tsv.gz
    resources: mem_mb=<TBD>, disk_mb=<TBD>, tmpdir=/home/jl9324/mnt/scratch/jl9324

Job stats:
job                              count
-----------------------------  -------
all                                  1
call_macs_peaks                      1
create_neighborhoods                 1
create_predictions                   1
filter_predictions                   1
generate_chrom_sizes_bed_file        1
generate_qc_plot_and_summary         1
make_candidate_regions               1
sort_narrowpeaks                     1
total                                9

Reasons:
    (check individual jobs above for details)
    input files updated by another job:
        all, create_neighborhoods, create_predictions, filter_predictions, generate_qc_plot_and_summary, make_candidate_regions, sort_narrowpeaks
    missing output files:
        call_macs_peaks, create_neighborhoods, create_predictions, filter_predictions, generate_chrom_sizes_bed_file, generate_qc_plot_and_summary, make_candidate_regions, sort_narrowpeaks

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.

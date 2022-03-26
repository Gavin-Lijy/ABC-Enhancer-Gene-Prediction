# Generates metadata file merging all available DHS/ATAC files with H3K27ac 
python grabDownload.py --dhs ../output/hg38/bam_GRCh38_DHS.tsv --h3k27ac ../output/hg38/bam_GRCh38_H3K27ac.tsv --atac ../output/hg38/bam_GRCh38_ATAC.tsv  --dhs_fastq ../output/hg38/fastq_GRCh38_DHS.tsv --h3k27ac_fastq ../output/hg38/fastq_GRCh38_H3K27ac.tsv --atac_fastq ../output/hg38/fastq_GRCh38_ATAC.tsv --genome_assembly GRCh38 --outdir /users/kmualim/updated_ABC/github/ABC-Enhancer-Gene-Prediction/Snakefiles/workflow/output/hg38/  --expt_file ../output/ABC_ENCODE_hg38_lookup.tsv

# Generates metadata file merging all available DHS/ATAC files with H3K27ac in narrowPeak
python process_dhs_metadata.py metadata/tf-chipseq-metadata-bed.tsv metadata/tf-chipseq-metadata-H3K27ac_bed.tsv
python match_id_to_file.py

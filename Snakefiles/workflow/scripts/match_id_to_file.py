import pandas as pd
import numpy as np

# ABC_ENCODE_hg38_lookup.tsv is the final lookup table generated from grabDownload.py
data = pd.read_csv("ABC_ENCODE_hg38_lookup.tsv", sep="\t", header=None)
# These input files are generated from download.sh, downloads the metadata file for narrowPeak files 
dhs_bam = pd.read_csv("tf-chipseq-metadata-bam.tsv", sep="\t")
h3k27ac_bam = pd.read_csv("tf-chipseq-metadata-h3k27ac-bam.tsv", sep="\t")

dhs_experiment = []
h3k27ac_experiment = []
for idname, dhs, h3k27ac in zip(data[0], data[1], data[2]):
    dhs_1 = dhs_bam.loc[dhs_bam['File accession']==dhs]
    dhs_experiment.append(dhs_1['Experiment accession'].values[0])
    h3k27ac_1 = h3k27ac_bam.loc[h3k27ac_bam['File accession']==h3k27ac]
    h3k27ac_experiment.append(h3k27ac_1['Experiment accession'].values[0])

data['DHS Experiment'] = dhs_experiment
data['H3K27ac Experiment'] = h3k27ac_experiment

data.to_csv("ABC_ENCODE_hg38_lookup.ExperimentAccession.tsv", sep="\t", index=False, header=False)

# File generated from process_dhs_metadata.py
narrowPeak = pd.read_csv("DNase-H3K27ac-lookup.tsv", sep="\t")

indices = []
biosamples = []
for idname, dhs, h3k27ac in zip(data[0], data['DHS Experiment'], data['H3K27ac Experiment']):
    matched = narrowPeak.loc[(narrowPeak['Experiment accession_Accessibility']==dhs) & (narrowPeak['Experiment accession_H3K27ac']==h3k27ac)].index.astype('int')
    if len(matched)>0:
        indices.append(matched)
        biosamples.append([idname]*len(matched))

narrowPeak_subset = narrowPeak.loc[np.concatenate(indices, axis=None), :]
narrowPeak_subset['biosamples'] = np.concatenate(biosamples, axis=None)
to_save = narrowPeak_subset.drop_duplicates(['File accession_Accessibility', 'File accession_H3K27ac', 'Biosample term name'])
to_save.to_csv("ENCODE_narrowPeak_lookup.matched.tsv", sep='\t', index=False)

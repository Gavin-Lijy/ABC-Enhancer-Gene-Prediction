import pandas as pd
import os, sys

# process narrowPeak files
dhs = pd.read_csv(sys.argv[1], sep="\t")
h3k27ac = pd.read_csv(sys.argv[2], sep="\t")

# filter for peaks in narrowPeak files 
dhs1 = dhs.loc[(dhs['File assembly']=="GRCh38") & (dhs['File Status']=='released') & (dhs['File analysis status']=="released") & (dhs['File analysis title'].str.contains("ENCODE4"))]
h3k27ac1 = h3k27ac.loc[(h3k27ac['File assembly']=="GRCh38") & (h3k27ac['File Status']=="released") & (h3k27ac['File analysis status']=="released") & (h3k27ac['File analysis title'].str.contains("ENCODE4"))]
dhs2 = dhs1.loc[dhs1['Output type']=="peaks"]
dhs2 = dhs2.drop_duplicates(['Experiment accession', 'Biological replicate(s)', 'Technical replicate(s)'])

merge_columns = ['Biosample term name'] #'Biosample organism', 'Biosample treatments','Biosample treatments amount', 'Biosample treatments duration','Biosample genetic modifications methods','Biosample genetic modifications categories','Biosample genetic modifications targets', 'Biosample genetic modifications gene targets', 'File assembly', 'File format', 'File type', 'Output type']

h3k27ac2 = h3k27ac1.loc[h3k27ac1['Output type']=="pseudoreplicated peaks"]
h3k27ac2 = h3k27ac1.drop_duplicates(['Experiment accession', 'Biological replicate(s)', 'Technical replicate(s)'])

dhs_lookup_table = (dhs2.groupby('Experiment accession').agg({'File accession' : ','.join})).reset_index()
h3k27ac_lookup_table = (h3k27ac2.groupby('Experiment accession').agg({'File accession' : ','.join})).reset_index()

for exp, filename in zip(dhs_lookup_table['Experiment accession'], dhs_lookup_table['File accession']):
    matched = dhs2.loc[dhs2['Experiment accession']==exp].index.astype('int')
    dhs2.loc[matched, 'File accession'] = filename

for exp, filename in zip(h3k27ac_lookup_table['Experiment accession'], h3k27ac_lookup_table['File accession']):
    matched = h3k27ac2.loc[h3k27ac2['Experiment accession']==exp].index.astype('int')
    h3k27ac2.loc[matched, 'File accession'] = filename

intersected = pd.merge(dhs2.convert_dtypes(), h3k27ac2.convert_dtypes(), how='inner', on=merge_columns, suffixes=('_Accessibility', '_H3K27ac'))
intersected.to_csv("DNase-H3K27ac-lookup.tsv", sep="\t", index=False)

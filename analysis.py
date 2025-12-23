import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
counts_df = pd.read_csv('/home/jiguo/data/data/outdir_rnaseq/counts/counts.txt', 
                        sep='\t', comment='#')
# set gene_id as index
counts_df = counts_df.set_index('Geneid')
# add Genesymbol column
import mygene
mg = mygene.MyGeneInfo()
gene_ids = counts_df.index.str.split(".").str[0]
gene_info = mg.querymany(
    gene_ids,
    scopes="ensembl.gene",
    fields="symbol",
    species="human",
    as_dataframe=True
)

# remove duplicated indices
gene_info = gene_info[~gene_info.index.duplicated(keep="first")]

counts_df["Genesymbol"] = gene_info["symbol"].reindex(gene_ids).values
# fill NaN values in Genesymbol with the original Geneid
counts_df["Genesymbol"] = counts_df["Genesymbol"].fillna(
    pd.Series(gene_ids, index=counts_df.index)
)
# set Genesymbol as index
counts_df = counts_df.set_index("Genesymbol")
# remove annotation columns 
counts_cols = [col for col in counts_df.columns if 'Aligned' in col]
counts = counts_df[counts_cols]
counts.columns = counts.columns.str.replace('.Aligned.sortedByCoord.out.bam', '')
sample_info = pd.DataFrame({
    'sample_id' : ['rep1_0', 'rep1_1', 'rep2_0', 'rep2_1', 'rep3_0', 'rep3_1','control'],
    'cell_line' : ['HeLa', 'HeLa', 'RPE-p53KO', 'RPE-p53KO', 'KE37-WT', 'KE37-WT','RPE-SAS6KO'],
    'condition' : ['with_centrioles'] * 6 + ['without_centrioles'],
    'replicate' : [1, 2] * 3 + [1]
}).set_index('sample_id')
counts.columns = counts.columns.str.replace('.Aligned.sortedByCoord.out.bam', '')
# save the metadata
sample_info.to_csv('/home/jiguo/data/data/outdir_rnaseq/counts/sample_metadata.csv')
counts.to_csv('/home/jiguo/data/data/outdir_rnaseq/counts/annotated_counts.txt', sep='\t')
# keep genes with at least 1 counts in at least 1 samples
keep = (counts>=1).sum(axis=1) >= 1
counts_filtered = counts[keep]
print(f'filtered from {counts.shape[0]} to {counts_filtered.shape[0]} genes')
counts_filtered.to_csv('/home/jiguo/data/data/outdir_rnaseq/counts/filtered_annotated_counts.txt', sep='\t')
library_sizes = counts_filtered.sum(axis=0)
cpm = counts_filtered.div(library_sizes, axis=1) * 1e6 # axis =1 for meaning column-wise division
log_cpm=np.log2(cpm + 1)  # log2 transform with pseudocount of 1
print("library sizes:\n", library_sizes)
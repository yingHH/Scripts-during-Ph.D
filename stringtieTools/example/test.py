#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
from stringtieTools.read_files import read_genes_exp, read_transcripts_exp
from stringtieTools.gfftools import ncbi_gene_id_lst, ncbi_rna_id_lst, ncbi_add_generalID_to_maskID
from stringtieTools.mergeSamples import merge_gene_exp, merge_tx_exp
#%% test 'read_genes_exp'
read_genes_exp('./stringtieTools/example/Chicken_I089-T05_good_1.abund.txt')

#%% test 'read_transcripts_exp'
read_transcripts_exp('stringtieTools\example\Chicken_I089-T05_good_1.abund.gtf')

#%% test 'ncbi_gene_id_lst'
ncbi_gene_id_lst('stringtieTools\example\gal5_ncbi.refseq.head100.gff')

#%% test 'ncbi_rna_id_lst'
ncbi_rna_id_lst('stringtieTools\example\gal5_ncbi.refseq.head100.gff')

#%% test 'ncbi_add_generalID_to_maskID'
ncbi_add_generalID_to_maskID(
    'gene', 
    read_genes_exp(
    './stringtieTools/example/Chicken_I089-T05_good_1.abund.txt'),
    'Gene ID',
    'stringtieTools\example\gal5_ncbi.refseq.head100.gff'
)


#%%
ncbi_add_generalID_to_maskID(
    'tx',
    read_transcripts_exp(
        'stringtieTools\example\Chicken_I089-T05_good_1.abund.gtf'),
    'transcript_id',
    'stringtieTools\example\gal5_ncbi.refseq.head100.gff'
)


#%%
merge_gene_exp(
    'stringtieTools\example', 'stringtieTools\example\gal5_ncbi.refseq.head100.gff',
    exp_type='fpkm',
)


#%%
merge_tx_exp(
    'stringtieTools\example', 'stringtieTools\example\gal5_ncbi.refseq.head100.gff',
    exp_type='fpkm',
)


#%%

# coding: utf-8


import importlib
import sys

import pandas as pd
if not '../utils/' in sys.path:
    sys.path.insert(0,'../utils/')

import dataframes
importlib.reload(dataframes)

# Preformatting of GENCODE V29.
path = "./input_enhancers/v29_regulatory_regions.bed.gz"
df = pd.read_csv(path, header=None, sep="\t", names='chrom,start,end,name'.split(','))
df['name'].head(100).apply(
        lambda v: [dict([(k, v) if k!='gene' else (k,gene_map(v,v)) for k,v in zip(['biotype','gene','strand'], e.split(','))]) for e in v.split('|')])

path = "./gene_tables/gencode_ensemblIDs_to_geneNames.tsv.gz"
gene_map = pd.read_csv(path, header=0, index_col=None, sep="\t")
gene_map = gene_map.set_index('gene_id')['gene_name'].to_dict()

new_name = df['name'].apply(lambda v: '|'.join(
                [dataframes.convert_keyvalue_pairs_to_str([(k, v) if k!='gene' else (k,gene_map.get(v,v))
                    for k,v in zip(['biotype','gene','strand'], e.split(','))],',',':') for e in v.split('|')]))
df['name'] = new_name
df.to_csv("./postprocessed_enhancers/v29_regulatory_regions_NAMES.bed.gz", header=False, index=False, sep="\t", compression="gzip")


# NOW LET'S FORMAT ALL THE DATASETS TO A COMMON FORMAT.
# THIS FORMAT WILL be chrom, start, end, name
# NAME : will combine the annotations of the datasets, with other info, such as the datasetname,
# or a score for the element (that's for GeneHancer mostly)

# GENCODEV29 PROM + UTR
df = pd.read_table("./preprocessed_enhancers/v29_regulatory_regions_NAMES.bed.gz", header=None, names=['chrom','start','end','name'])
df['name'] = 'src=gencodeV2::annots='+df['name']
df.to_csv("./postprocessed_enhancers/FINAL_GENCODE_v29_regulatory_regions_NAMES.bed.gz", header=False, index=False, sep="\t", compression="gzip")

# FANTOM5
df = pd.read_table("./preprocessed_enhancers/fantom5_enh.bed.gz", header=None, names = ['chrom','start','end','name'])
df['name'] = 'src=FANTOM5::annots='+df['name']
df.to_csv("./postprocessed_enhancers/FINAL_fantom5_enh.bed.gz",header=False,index=False,sep="\t",compression="gzip")

# FOCS_F5
df = pd.read_table("./preprocessed_enhancers/FOCS_fantom5_enh.bed.gz", header=None, names = ['chrom','start','end','name'])
df['name'] = 'src=FOCSFANTOM5::annots='+df['name']
df.to_csv("./postprocessed_enhancers/FINAL_FOCSfantom5_enh.bed.gz",header=False,index=False,sep="\t",compression="gzip")

# PEGASUS
df = pd.read_table("./preprocessed_enhancers/pegasus_enh.bed.gz", header=None, names = ['chrom','start','end','name'])
new_name = 'src=PEGASUS::annots='+df['name'].apply(
            lambda v: dataframes.convert_keyvalue_pairs_to_str([(k,v)
                        for k,v in dataframes.transform_to_dict(v,',',':').items() if k!='src']))
df['name'] = new_name
df.to_csv("./postprocessed_enhancers/FINAL_pegasus_enh.bed.gz",header=False,index=False,sep="\t",compression="gzip")


# FOCS_GS
df = pd.read_table("./preprocessed_enhancers/FOCS_groseq_enh.bed.gz", header=None, names = ['chrom','start','end','name'])
df['name'] = 'src=FOCSGROSEQ::annots='+df['name']
df.to_csv("./postprocessed_enhancers/FINAL_FOCSgroseq_enh.bed.gz",header=False,index=False,sep="\t",compression="gzip")

# FOCS_RM
df = pd.read_table("./preprocessed_enhancers/FOCS_roadmap_enh.bed.gz", header=None, names = ['chrom','start','end','name'])
df['name'] = 'src=FOCSROADMAP::annots='+df['name']
df.to_csv("./postprocessed_enhancers/FINAL_FOCSroadmap_enh.bed.gz",header=False,index=False,sep="\t",compression="gzip")


# GENEHANCER
df = pd.read_table("./preprocessed_enhancers/genehancer_enh.bed.gz", header=None, names=['chrom','start','end','name','score','strand','info'])
# First we can try and reduce the redundancy in annotations ;
# for instance : 'src:GH01F000638,targets:MTND2P28,score:0.76|src:GH01F000638,targets:MTCO1P12,score:0.76'
# can be : 'src:GH01F000638,targets:MTND2P28;MTCO1P12,score:0.76'

df['info2'] = df['info'].str.split('|')

def reduce_genehancer_info(info):
    if len(info)==1:
        return info
    else:
        new_info = [dataframes.transform_to_dict(v,',',':') for v in info]
        new_new_info = []
        already_done_idx = []
        for i in range(len(new_info)):
            if i in already_done_idx:
                continue
            el = new_info[i]
            src, score = el['src'], el['score']
            for i2, el2 in enumerate(new_info[i+1:]):
                src2, score2 = el2['src'],el2['score']
                if src==src2 and score==score2:
                    el['targets'] = el['targets']+';'+el2['targets']
                    already_done_idx.append(i2+i+1)
            new_new_info.append(el)
            
        return [dataframes.convert_keyvalue_pairs_to_str(v.items()) for v in new_new_info]
        
df['info2'] = df['info2'].apply(reduce_genehancer_info)
df['name2'] = 'src='+df['name']+'::score='+df['score']+'::annots='+df['info2'].apply(lambda v: '|'.join(v))
df['name2'] = 'src='+df['name']+'::score='+df['score'].astype(float)+'::annots='+df['info2'].apply(lambda v: '|'.join(v))
df['name2'] = 'src='+df['name']+'::score='+df['score'].astype(str)+'::annots='+df['info2'].apply(lambda v: '|'.join(v))

df.loc[:,['chrom','start','end','name2']].to_csv("./postprocessed_enhancers/FINAL_GENEHANCER_enh.bed.gz",header=False,index=False,sep="\t",compression="gzip")

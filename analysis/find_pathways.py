import os
import pandas as pd
import numpy as np
import requests
import pandas as pd
import json
from scipy.stats import hypergeom
import statsmodels.stats.multitest as smm
from tqdm import tqdm
tqdm.pandas()
# import gseapy as gp
import sys
from io import StringIO
from Bio.KEGG import REST
# from Bio.KEGG.KGML import KGML_parser
# from Bio.Graphics.KGML_vis import KGMLCanvas
import io
import urllib
currdir=os.getcwd()
parent = os.path.dirname(currdir)
gparent = os.path.dirname(parent)
import uniprot

# http://togows.org/entry/kegg-pathway/hsa04640/genes.json
## human genes: https://rest.kegg.jp/link/hsa/hsa04640 

# https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/09-KEGG_programming.html
# A bit of code that will help us display the PDF output
def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

# Some code to return a Pandas dataframe, given tabular text
def to_df(result):
    return pd.read_table(io.StringIO(result), header=None)

def all_pathways():
    result = REST.kegg_list("pathway").read()
    result = to_df(result)
    return result

def ontology(pathway): ## get all genes in the pathway
    url = f'https://rest.kegg.jp/link/mmu/mmu{pathway}'
    g = requests.get(url=url)
    if g.text.strip():
        gtext = StringIO(g.text) # convert string to file like format
        gdf = pd.read_csv(gtext, sep= '\t', header=None, names = ['Pathway', 'KEGG_geneIdentifiers'])
        return list(set(gdf['Pathway'].values))[0], list(gdf['KEGG_geneIdentifiers'].values)
    else:
        return pathway, []

def get_KEGGnums(uniprot_id): # conv
    url = f"https://rest.kegg.jp/conv/genes/uniprot:{uniprot_id}"
    response = requests.get(url)
    if response.text.strip():
        vals = StringIO(response.text)
        conv = pd.read_csv(vals, sep='\t', header=None, names=['Uniprot', 'Kegg_Gene'])
        return conv['Kegg_Gene'].values[0]
    else:
        return uniprot_id

def dict2df(d): #dict
    df = pd.DataFrame(columns = ['Genotype', 'adj_pval','Relevant Genes'], index=list(d.keys()))
    for i in df.index:
        genotypes = [x[0] for x in d[i]]
        adjpval = [x[1] for x in d[i]]
        genes = [x[2] for x in d[i]]
        df.loc[i, 'Genotype'] = genotypes
        df.loc[i, 'adj_pval'] = adjpval
        df.loc[i, 'Relevant Genes'] = genes
    return df

def main():
    ALLPATHWAYS = all_pathways()
    PWAYINT = 'Immune'
    immune_pathways  = ['04640', '04610', '04611', '04613', '04620', '04624', 
    '04621', '04622', '04623', '04625', '04650', '04612', 
    '04660', '04658', '04659', '04657', '04662', '04664', 
    '04666', '04670', '04672', '04062']
    ## our data 
    READ_FILE = f'~/Downloads/PPI_list_merged.3_DS (2).xlsx'
    uniprots_json = f'../uniprots.json'
    keggGenes = f'../kegg_gidentifiers.json'
    if os.path.isfile(uniprots_json):
        with open(uniprots_json) as f:
            uniprotd = json.load(f)
    if os.path.isfile(keggGenes):
        with open(keggGenes) as f:
            keggd = json.load(f)
    else:
        keggd = {}
    # READ_FILE = sys.argv [1]
    INDEXCOL = 0 # index = genes?
    if 'xlsx' in READ_FILE:
        data = pd.read_excel(READ_FILE, index_col=INDEXCOL)
    else:
        data = pd.read_csv(READ_FILE)
    data = data.drop(columns=[x for x in data.columns if 'phos' in x.lower() or 'bio' in x.lower()])
    allmygenes = [x.split(';')[0] for x in data.index.values]
    for i in range(0, len(allmygenes), 10):
        subsec = allmygenes[i: i+10]
        if any(gene not in uniprotd for gene in subsec):
            out = uniprot.main(subsec)
            uniprotd.update(out)
            with open(uniprots_json, 'w') as f:
                json.dump(uniprotd, f)
    data.index = data.index.map(uniprotd)
    for x in tqdm(data.index.values):
        if x not in keggd:
            keggd[x] = get_KEGGnums(x)
            with open(keggGenes, 'w') as f:
                json.dump(keggd, f)
    data.index = data.index.map(keggd)
    keggdB = {value: key for key, value in keggd.items()} # backward to map KG 0> uniprot -> gene 
    uniprotdB = {value: key for key, value in uniprotd.items()} # BACKWARDS
    MYGENES = {}
    for gr in data.columns:
        d = data[gr].dropna()
        ggenes = list(d.index)
        MYGENES[gr] = ggenes
    IMMUNEPATHS = ALLPATHWAYS[ALLPATHWAYS[0].apply(lambda x: x.split('map')[-1]).isin(immune_pathways)].reset_index(drop=True)
    ## take intersection between PGENES and immunepaths
    N= 20000 # TOTAL Mouse genes
    PATHWAYS = {f'{IMMUNEPATHS.iloc[i, 1]} ({str(x)})': [] for i, x in enumerate(IMMUNEPATHS[0].values) }
    for i, x in tqdm(enumerate(IMMUNEPATHS[0].values), total=len(IMMUNEPATHS[0].values), desc='Hypergeometric test'):
        pathNum, PGENES = ontology(x.split('map')[-1])
        K = len(PGENES) # NUM OF PATHWAY GENES
        PathwayEnrichment = {}
        for gr,mygenes in MYGENES.items(): #hyp: mygenes not enriched in PGENES, overlap is due to chance.
            n = len(mygenes) # num genes in selec set.
            intersection = list(set(mygenes) & set(PGENES))
            k = len(list(set(mygenes) & set(PGENES)))
            pval = hypergeom.sf(k-1, N, K, n)  # https://calcworkshop.com/discrete-probability-distribution/hypergeometric-distribution/
            PathwayEnrichment[gr] = pval
        TF, newpvals = [True if x < 0.05 else False for x in PathwayEnrichment.values()], list(PathwayEnrichment.values())
        MULTIPLETESTING = False
        if MULTIPLETESTING:
            TF, newpvals, _, _ = smm.multipletests(list(PathwayEnrichment.values()),
                                                    alpha = 0.05, method='holm')
        for ind, (gr, opval) in enumerate(PathwayEnrichment.items()):
            if TF[ind]:
                asuni = [keggdB[mmu] for mmu in intersection]
                asgene = [uniprotdB[mmu] for mmu in asuni]
                PATHWAYS[f'{IMMUNEPATHS.iloc[i, 1]} ({str(x)})'].append((gr, newpvals[ind], asgene))
    df = dict2df(PATHWAYS)
    basename = os.path.basename(READ_FILE).split('.xlsx')[0]
    df.to_csv(f'../results/NoMTadjust{PWAYINT}_GO_{basename}.csv')

if __name__ == "__main__":
    main()

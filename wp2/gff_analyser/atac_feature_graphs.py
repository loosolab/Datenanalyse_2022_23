import os
import numpy as np
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt

def violin_plots(adata, output, group=None, filtered=None, multi_panel=True, log=False):
    
    if filtered == None:
        item_list = []

        for item in adata.obs:
            if adata.obs[item].dtype == 'int64':
                item_list.append(item)
            elif adata.obs[item].dtype == 'float64' and not item.startswith("n_"):
                item_list.append(item)
        
        
    else:
        item_list = filtered
    
    try:
        os.mkdir(f"{output}violins/")
    except:
        print("Folder Exists")
            
    for item in item_list:
        sc.pl.violin(adata, keys=item, groupby=group, multi_panel=multi_panel, rotation=90, log=log)
        
        try:
            tmp_label = group.split(' ')
            label = ''.join(tmp_label)
        except:
            label = group
            
        plt.savefig(f"{output}violins/{item}_{label}_violin.png")
        
def scatter_plots(adata, output, feature_filtered=None, celltype_filtered=None, figsize=[60,3]):
    
    if celltype_filtered != None:
        celltype_of_interest = celltype_filtered
        adata_tmp = adata[adata.obs['cell type'].isin(celltype_of_interest)]
    else:
        adata_tmp = adata
    
    if feature_filtered == None:
        item_list = []
        for item in adata_tmp.obs:
            if adata_tmp.obs[item].dtype == 'float64' and not item.startswith("n_"):
                item_list.append(item)
                
    else:
        item_list = feature_filtered
        
    try:
        os.mkdir(f"{output}images/")
    except:
        print("Skipping creating folder")

    for itemx in item_list:

        fig, ax= plt.subplots(1, len(item_list), figsize=figsize)
        x = adata_tmp.obs[itemx]
        for n, itemy in enumerate(item_list):
            if not itemy == itemx:
                y = adata_tmp.obs[itemy]
                ax = plt.subplot(1, len(item_list), n + 1)
                ax.scatter(x, y)
                plt.title(celltype_of_interest)
                plt.xlabel(itemx, fontsize=10)
                plt.ylabel(itemy, fontsize=15)
                fig.tight_layout()
        plt.savefig(f"{output}images/{itemx}.{itemy}.png")
        plt.show()
        plt.close()
        
def simple_scatter(adata, output, feature=['n_fragments_in_CDS', 'n_total_fragments'], display_all=False):
    item_list = []
    for item in adata.obs:
        if adata.obs[item].dtype == 'int64':
                item_list.append(item)
        elif adata.obs[item].dtype == 'float64' and not item.startswith("n_"):
            item_list.append(item)
    if display_all:
        feature = item_list
        
    try:
        os.mkdir(f"{output}scatter")
    except:
        print("Folder exists")
            
    for itemx in feature:

        x = adata.obs[itemx]
        for n, itemy in enumerate(item_list):
            if not itemy == itemx:
                y = adata.obs[itemy]
                plt.scatter(x, y)
                plt.xlabel(itemx, fontsize=10)
                plt.ylabel(itemy, fontsize=10)
                plt.show()
                plt.savefig(f"{output}scatter/{itemx}_{itemy}.png")
                plt.close()
                
                
def multi_plot(adata, feature1, feature2, group, out, celltype_filtered=None, multi_panel=False):
    
    output= f"{out}{feature1}_{feature2}/"
    try:
        if not out.endswith('/'):
            print('Please specify valid Path ('/')')
            return
        os.mkdir(out)        
    except:
        print('First out exists')
        
    try:
        os.mkdir(output)
    except:
        print('Second too')
    
    if celltype_filtered != None:
        celltype_of_interest = celltype_filtered
        adata_tmp = adata[adata.obs['cell type'].isin(celltype_of_interest)]
    else:
        adata_tmp = adata
        
    sc.settings.figdir = "./"
    sc.pl.violin(adata_tmp, keys=[feature1, feature2], groupby=group, multi_panel=multi_panel, rotation=90)
    plt.savefig(f"{output}violin.png")
    
    for  cell_type_raw in adata_tmp.obs['cell type'].unique():
        
        tmp = cell_type_raw.split(' ')
        cell_type = ''.join(tmp)
        
        # for getting rid of spaces, this object only includes filtered cell types since its filtered before
        tmp_adata = adata_tmp[adata_tmp.obs['cell type'].isin([cell_type_raw])]
        
        
        plt.scatter(tmp_adata.obs[feature1], tmp_adata.obs[feature2])
        plt.title(cell_type_raw)
        plt.xlabel(feature1)
        plt.ylabel(feature2)
        plt.savefig(f"{output}{feature1}_{feature2}_{cell_type}.png")
        plt.show()
        
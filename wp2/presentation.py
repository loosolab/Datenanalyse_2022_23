import os
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
from matplotlib import pyplot as plt

def getList(element) -> list:
    list = []
    if type(element) is list:
        list.extend(element)
    else:
        list.append(element)
    return list

def compareDimensionreductions(adata, key: str or list, comparator: str or list, save=None, legend_loc="right margin"):
    keys = getList(key)
    keys.append(getList(comparator))
    sc.pl.umap(adata, color=keys, save=save, legend_loc=legend_loc)

def compareFeatureToCelltypes(adata, feature: list or str, comparator: str, save=None):
    if type(feature) is list:
        features = feature
        for feature in features:
            print(f'{feature}:')
            renderCompareFeatureToCelltypes(adata, feature, comparator, save=save)
    else:
        renderCompareFeatureToCelltypes(adata, feature, comparator, save=save)

def renderCompareFeatureToCelltypes(adata, feature, comparator, save):
    label = 'Percent' if feature.startswith('pct') else 'Count'
    fig, axes = plt.subplots(nrows=1, ncols=2, gridspec_kw={'wspace':0.4})
    sc.pl.violin(adata, feature, ax = axes[0], show=False, ylabel=label)
    sc.pl.violin(adata, feature, groupby=comparator, ax = axes[1], rotation=90, show=False, ylabel="", save=save)


def createOutputDirectory(out):
    try:
        os.mkdir(f"{out}")
    except:
        print("Folder exists.")

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

    createOutputDirectory(output+"violins/")

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
        adata = adata[adata.obs['cell type'].isin(celltype_of_interest)]
    
    if feature_filtered == None:
        item_list = []
        for item in adata.obs:
            if adata.obs[item].dtype == 'float64' and not item.startswith("n_"):
                item_list.append(item)
                
    else:
        item_list = feature_filtered
        
    try:
        os.mkdir(f"{output}images/")
    except:
        print("Skipping creating folder")

    for itemx in item_list:

        fig, ax= plt.subplots(1, len(item_list), figsize=figsize)
        x = adata.obs[itemx]
        for n, itemy in enumerate(item_list):
            if not itemy == itemx:
                y = adata.obs[itemy]
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
        
    createOutputDirectory(output+"scatter")
            
    for itemx in feature:

        x = adata.obs[itemx]
        for _, itemy in enumerate(item_list):
            if not itemy == itemx:
                y = adata.obs[itemy]
                plt.scatter(x, y)
                plt.xlabel(itemx, fontsize=10)
                plt.ylabel(itemy, fontsize=10)
                plt.show()
                plt.savefig(f"{output}scatter/{itemx}_{itemy}.png")
                plt.close()
                
                
def multi_plot(adata, feature1, feature2, group, out, multi_panel=False):
    # TODO: add pathlib path instead
    output= f"{out}{feature1}_{feature2}/"
    
    if not out.endswith('/'):
        print('Please specify valid Path ('/')')
        return
    createOutputDirectory(out)          
    createOutputDirectory(output)
    
        
    sc.settings.figdir = "./"
    sc.pl.violin(adata, keys=[feature1, feature2], groupby=group, multi_panel=multi_panel, rotation=90)
    plt.savefig(f"{output}violin.png")
    
    for  cell_type_raw in adata.obs['cell type'].unique():
        
        tmp = cell_type_raw.split(' ')
        cell_type = ''.join(tmp)
       
        adata.obs
        tmp_adata = adata[adata.obs['cell type'].isin([cell_type_raw])]
        
        
        plt.scatter(tmp_adata.obs[feature1], tmp_adata.obs[feature2])
        plt.title(cell_type_raw)
        plt.xlabel(feature1)
        plt.ylabel(feature2)
        plt.savefig(f"{output}{feature1}_{feature2}_{cell_type}.png")
        plt.show()
        
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
from matplotlib import pyplot as plt

def compareDimensionreductions(key: str or list, comparator_key: str, save=None, legend_loc="right margin"):
    keys = []
    if type(key) is list:
        keys.extend(key)
    else:
        keys.append(key)
    keys.append(comparator_key)
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

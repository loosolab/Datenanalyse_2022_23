## Data Presentation

After preparing the data for presentation with some normalization and filtering, such as, but not limited to, removing chromosomes X, Y (Optional) and M, removing cells without features or empty features.
These functions are combined into the `presentation.py` script, to avoid unnecessary clutter in any `Jupyter Notebook` you might use them in.

### Different Presentation Styles

- Violinplots
  - display density accross celltypes and featurerate
- Correlation Plots
  - scatterplots that compare two metrics
  - can visualize linear or a parabolical correlation
- Dimension Reduction Maps (umaps)

## Functions

### Dimension Reduction
```py
compareDimesionreductions(adata, key: str or list, comparator: str or list)
```

This function displays the `Dimension Reductions` defined in key and lastly add the `Comperator` or `Comperators` to compare the prior maps to. The comparator will be visible in each row while the dimension reductions labelled by `key` will be one per row.

### Feature Comparison in Violinplots

```py
def compareFeatureToCelltypes(adata, feature: list or str, comparator: str, save=None):
```

The concept is similar to the previous mentioned function, but it focuses on `violin plots`. It expects either a single or a list of strings that represent a feature and a comparator that is used to group these together, in practice we used the predicted celltypes.




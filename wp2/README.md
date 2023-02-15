### Introduction


1. [The data source](#The-data-source)
2. [Reading a General Feature File](#Reading-a-fragment-file)
3. [Calculating Feature Overlap](#Plotting)
   1. [Annotated Features](#Mean-and-Median)
   2. [Other Features](#Distribution-into-groups)
4. [Plotting](#Distribution-per-cell)
   1. [Violin and Correlation](#Calculate-the-fit)
   2. [Dimension Reduction](#Calculate-a-score)
5. [References](#References)
6. [Quick Start](#testing)
7. [Authors](#Authors)

### The data source

### Reading a General Feature File

### Calculating Feature Overlap

#### Annotated Features

#### Other Features

### Plotting

After preparing the data for presentation with some normalization and filtering, such as, but not limited to, removing chromosomes X, Y (Optional) and M, removing cells without features or empty features.
These functions are combined into the `presentation.py` script, to avoid unnecessary clutter in any `Jupyter Notebook` you might use them in.

#### Different Presentation Styles

- Violinplots
  - display density accross celltypes and featurerate
- Correlation Plots
  - scatterplots that compare two metrics
  - can visualize linear or a parabolical correlation
- Dimension Reduction Maps (umaps)

#### Functions

##### Dimension Reduction
```py
compareDimesionreductions(adata, key: str or list, comparator: str or list)
```

This function displays the `Dimension Reductions` defined in key and lastly add the `Comperator` or `Comperators` to compare the prior maps to. The comparator will be visible in each row while the dimension reductions labelled by `key` will be one per row.

##### Feature Comparison in Violinplots

```py
def compareFeatureToCelltypes(adata, feature: list or str, comparator: str, save=None):
```

The concept is similar to the previous mentioned function, but it focuses on `violin plots`. It expects either a single or a list of strings that represent a feature and a comparator that is used to group these together, in practice we used the predicted celltypes.

##### Correlation Plots


### Quick Start


### References

#### <font size=1>- [1] Zhang, K., Hocker, J. D., Miller, M., Hou, X., Chiou, J., Poirion, O. B., Qiu, Y., Li, Y. E., Gaulton, K. J., Wang, A., Preissl, S., &amp; Ren, B. (2021). A single-cell atlas of chromatin accessibility in the human genome. Cell, 184(24). https://doi.org/10.1016/j.cell.2021.10.024 </font>
#### <font size=1>- [2] https://genome.ucsc.edu/FAQ/FAQformat.html </font>
#### <font size=1>- [3] https://scanpy.readthedocs.io/en/stable/ </font> 
#### <font size=1>- [4] https://www.ensembl.org/info/website/upload/gff.html </font>

### Authors

    Daniel Tischler
        daniel.tischler@bioinfsys.uni-giessen.de

    Noah Leon Stürtz

---

### Appendix


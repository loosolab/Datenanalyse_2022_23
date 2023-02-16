# Creating rofile plots.

## Aproach 1: (Using Deeptools)

### Step 1: (filter bed, categorize fragments and create bedgraph)
Filter bed files and create bedgraph using the jupyter notebook: create_filtered_bedgraph.ipynb.
Here the fragments from the bed files can be categorized using various strategies. This may result in multiple bedgraph files.

### Step 2: (sort bedgraph)
    LC_COLLATE=C sort -k1,1 -k2,2n <INPUT_BEDGRAPH_FILE> > <SORTED_BEDGRAPH_FILE>
  
### Step 3: (merge overlapping regions)
    bedtools merge -i <SORTED_BEDGRAPH_FILE> -c 4 -d 0 -o sum > <SORTED_MERGED_OVERLAPS_BEDGRAPH_FILE>

### Step 4: (create bigWig)
    bedGraphToBigWig <SORTED_MERGED_OVERLAPS_BEDGRAPH_FILE> <CHROM_SIZES_FILE> <BIGWIG_FILE>

### Step 5: (calculate fragment count over TSS)
Here depending on the categorization strategies used in Step 1, multiple bigWig files may be passed as input.

    computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -S <BIGWIG_FILE> <BIGWIG_FILE_2> -R <GTF_FILE> -o <MATRIX_FILE_GZ> --outFileSortedRegions <REGIONS_FILE_BED>

### Step 6: (plot using deepTools)
    plotProfile --perGroup -m <MATRIX_FILE_GZ> -out <PLOT_PNG>
  
## Approach 2: (Using HTSeq)
Follow the procedure from the juyter notebook: profile_plots_HTSeq.ipynb.
TODO: Categorization of fragments.

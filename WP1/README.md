# THIS README IS NOT FINISHED AND STILL IN DEVELOPEMENT!!!

### Introduction

This repository contains all the tools and methods developed specifically for the course 
“Applied data analysis in bioinformatics” from the masters program “Bioinformatik und Systembiologie” 
at the Justus-Liebig-University and the Technische Hochschule Mittelhessen in the wintern term 2022/2023.

The goal of this course is to develop a pipeline for the [Max Planck Institute for Heart and Lung Research](https://www.mpg.de/149809/heart-lung-research)
which takes data from [CATLAS](http://catlas.org/humanenhancer/#!/) and performs distinct analyses mainly based on the so-called chromatin accessibility<sup>[1](#--1-zhang-k-hocker-j-d-miller-m-hou-x-chiou-j-poirion-o-b-qiu-y-li-y-e-gaulton-k-j-wang-a-preissl-s-amp-ren-b--2021---a-single-cell-atlas-of-chromatin-accessibility-in-the-human-genome-cell-184--24---httpsdoiorg101016jcell202110024)</sup>.   
Furthermore, this pipeline is organized into two separate packages (WP1/WP2) due to the group distribution in the course. Hence, the main purpose of
these methods is to be used by the second group (WP2), but they can nevertheless be used as a stand-alone tool. A graphical representation of our pipeline is given in the [appendix](#appendix).

In the following we are going to thoroughly illustrate all necessary steps to perform an analysis using this repository.
For this reason we will guide the reader through a basic example reproducing each important step, where we also explain the underlying methods
and their functionality along the way. We also provided a compact [Quick Start](#quick-start) section at the end for anyone who wants to get started right away. The rest of the following content is structured as follows.

1. [The data source](#The-data-source)
2. [Reading a fragment file](#Reading-a-fragment-file)
3. [Plotting](#Plotting)
   1. [Mean and Median](#Mean-and-Median)
   2. [Distribution into groups](#Distribution-into-groups)
4. [Distribution per cell](#Distribution-per-cell)
   1. [Calculate the fit](#Calculate-the-fit)
   2. [Calculate a score](#Calculate-a-score)
5. [Splitting fragments](#Splitting-fragments)
6. [References](#References)
7. [Quick Start](#testing)
8. [Authors](#Authors)

### The data source

As already mentioned in the introduction, the main data source used for developing the pipeline was [CATLAS](http://catlas.org/humanenhancer/#!/). 
This website contains files generated from [Single-Cell ATAC-seq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02075-3) (scATAC-Seq) done on human tissue.
Here, we exclusively focused on the _fragment files_, which consist of row(s) and columns. Each row represents a single fragment, whereas the columns specify the fragment further.
An excerpt of the fragment file [stomach_SM-JF1O3_rep1_fragments.bed](http://yed.ucsd.edu:8787/fragment/) can be seen below.

<p align="center">
   <img src="images/stomach_fragment_excerpt.png" />
</p>

The structure represented above is the same over all provided fragment files. Furthermore, we can specify the
content of each column as follows<sup>[2](#--2-httpsenwikipediaorgwikibedfileformat)</sup>:

| Column Number | Name            | Description                                                                  |
|---------------|-----------------|------------------------------------------------------------------------------|
| 1             | Chromosome      | Abbreviated name of the Chromosome (e.g. chr1 stands for Chromosome 1, etc). |
| 2             | ChromosomeStart | Number of the start location of the fragment (in bp) on the chromosome.      |
| 3             | ChromosomeEnd   | Number of the end location of the fragment (in bp) on the chromosome.        |
| 4             | CellBarcode     | Unique identifier of the cell from which the fragment is extracted.          |
| 5             | Score           | Score of the fragment.                                                       |
| 6             | Strand          | DNA strand orientation ("+" = positive; "-" = negative; "." = no strand)     |

Moreover, in our approach we only use the columns ChromosomeStart, ChromosomeEnd and the CellBarcode
and completely omit the remaining ones. This is further explained in the subsequent chapter.

### Reading a fragment file

In the previous chapter we described the data files which function as input for our pipeline. Now, the first step is
to read these files. For this purpose we implemented the package _fragment_. This package contains the following modules:


- _calc.py_ and 
- _utils.py_

The first module _calc.py_ contains following 4 methods:

    compute_mean(value_list: list, decimal_places=2)
    compute_median(value_list: list)
    is_even(value_list: list)
    calculate_fragment_length(start: int, stop: int)

These methods, as the naming of the module already suggests, are used to _calculate_ different tasks necessary
for the methods in the _utils.py_ module. A short description of the methods above can be seen below.

-  compute_mean(value_list: list, decimal_places=2)
     -  This method computes the **mean** of the provided _value_list_, which should obviusly only contain
   numeric values. Furthermore, the parameter _decimal_places_ defines the decimal places of the resulting **mean**. The standard
   value of this parameter is set to 2.
-  compute_median(value_list: list)
      -  This method computes the **median** of the provided _value_list_. This list should also consist of only numeric
   values. 
-  is_even(value_list: list)
     -  This method simply checks if the length of the list parameter is even and returns a boolean value based on the result.
-  calculate_fragment_length(start: int, stop: int)
     -  This method computes the difference between the integer parameters _start_ and _stop_. The absolute value of this
   difference is then returned.

The second module _utils.py_ holds at least the methods below. Additional methods in this module are not important for our pipeline and will thus
not be discussed here. 

    read_fragment_file(abs_path: str)
    read_fragment_directory(abs_path: str)
    combine_fragment_dictionaries(dictionary_list: list)
    create_dataframe(frag_dictionary: dict, tissue="")

Having all tools necessary for finally reading a fragment file, we can now discuss the first step of the pipeline.
For this purpose we will use the method **read_fragment_file(abs_path: str)**. This method reads a fragment file (.bed) at the
absolute path location provided via parameter and returns a dictionary having the _cellbarcode_ (column 1) as keys and the computed
_fragment lengths_ (column 2, column 3) as the corresponding values. Now, we can use this method with our file _stomach_frag_head_30.bed_, which simply contains
the first 30 lines of [stomach_SM-JF1O3_rep1_fragments.bed](http://yed.ucsd.edu:8787/fragment/).
    
    frag_dictionary = read_fragment_file("~/stomach_frag_head_30.bed")

Formatted printing of the _frag_dictionary_ yields the following output:

    AAACGCAAGCAAACCCGAGATA [32, 47, 311, 45, 153, 393, 161]
    AAACGCAAGCAAACCTAAGTGG [272, 62, 112, 33, 56]
    AAACGCAAGCAAACGGATCAGT [210, 56, 28, 42, 185, 75]
    AAACGCAAGCAAACGTCCCGTT [45, 33, 49, 227, 80, 94, 36, 340, 38]
    AAACGCAAGCAAACTAGCCCTA [48, 479, 85]

Reading a single file is a rather simple case. Thus, we also implemented a convenient wrapper method **read_fragment_directory(abs_path: str)** which
reads all fragment files (.bed) in a directory. Here, the absolute path of the directory is also provided via the parameter of the function. If this directory also
contains other files besides fragment files, they will be ignored by the function. Because the internal datastructure is a dictionary, we strictly prohibit the occurrence of
multiple cellbarcodes and simply extend the corresponding fragment list of the cellbarcode. Hence, enabling an easy-to-use approach for reading fragment files.

We went beyond the steps above and also implemented two additional convenience functions. The first being **combine_fragment_dictionaries(dictionary_list: list)** which takes a list of 
_fragment_dictionaries_ and merges them into one and the second one being **create_dataframe(frag_dictionary: dict, tissue="")**. This method uses a _fragment_dictionary_ and a string of a _tissue_
to create a dataframe<sup>[3](#font-size1---3-httpspandaspydataorgdocsreferenceapipandasdataframehtml-font)</sup>. The _tissue_ string here is used to append the cellbarcode resulting in the cellbarcode of the
form below.

    tissue+cellbarcode

This enables a linking between cellbarcode and tissue the cellbarcode may be taken from and increases distinguishability. Important to mention here is that the resulting dataframe not only consists
of cellbarcodes and fragment length lists, but also the computed _mean_ and _median_ of these lists. Thereby making use of the methods from the _calc.py_ module. To demonstrate this method we can use the earlier 
created _frag_dictionary_ and pass it to our method.

    df = create_dataframe(frag_dictionary)

When printing the dataframe _df_ this results in the following output.

                              Mean Median                               Fragments
    AAACGCAAGCAAACGGATCAGT   99.33   65.5              [210, 56, 28, 42, 185, 75]
    AAACGCAAGCAAACGTCCCGTT  104.67     49  [45, 33, 49, 227, 80, 94, 36, 340, 38]
    AAACGCAAGCAAACCTAAGTGG   107.0     62                  [272, 62, 112, 33, 56]
    AAACGCAAGCAAACCCGAGATA  163.14    153        [32, 47, 311, 45, 153, 393, 161]
    AAACGCAAGCAAACTAGCCCTA   204.0     85                           [48, 479, 85]

These resulting dataframes are used in the next step of our pipeline, where they are processed and analyzed further. 
### Plotting

#### Mean and Median

#### Distribution into groups

### Distribution per cell

#### Calculate the fit

#### Calculate a score

### Splitting fragments

### Quick Start

HIER WENIG ERKLÄRUNG UND GLEICH MIT CODE LOSLEGEN

AUCH AUF TESTDATENSATZ EINGEHEN

### References

#### <font size=1>- [1] Zhang, K., Hocker, J. D., Miller, M., Hou, X., Chiou, J., Poirion, O. B., Qiu, Y., Li, Y. E., Gaulton, K. J., Wang, A., Preissl, S., &amp; Ren, B. (2021). A single-cell atlas of chromatin accessibility in the human genome. Cell, 184(24). https://doi.org/10.1016/j.cell.2021.10.024 </font>
#### <font size=1>- [2] https://en.wikipedia.org/wiki/BED_(file_format) </font>
#### <font size=1>- [3] https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html </font> 

### Authors

    Leon Marvin Geis
      leon.marvin.geis@bioinfsys.uni-giessen.de
      
    Jannik Luebke 
      jannik.luebke@bioinfsys.uni-giessen.de

    Aviral Jain
      aviral.jain@bioinfsys.uni-giessen.de 

---

### Appendix

 INSERT ABSTRACT IMAGE OF PIPELINE 


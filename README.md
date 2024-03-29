# TuBA

New R package of TuBA - Analyzes large graphs to identify biclusters much faster than the [older version](https://github.com/Amartya101/TuBA-Original_Slower). 

## Getting Started

You will need an R version 3.4.0 or more recent in order to use the package.

### Prerequisites

You need the **data.table** and the **Matrix** package. You can install them using the following command:

```
install.packages("data.table",dependencies = T)
install.packages("Matrix",dependencies = T)
```
In order to make the bicluster gene graphs, we also need the following packages - **ggplot2**, **network**, **ggnetwork**. They can be installed using:
```
install.packages("ggplot2",dependencies = T)
install.packages("network",dependencies = T)
install.packages("ggnetwork",dependencies = T)
```

### Installing

Currently, the development version of **TuBA** can be installed using *install_github* from the **devtools** package.

Install **devtools** using:

```
install.packages("devtools",dependencies = T)
```
Once **devtools** is installed, run the following to install **TuBA**:
```
devtools::install_github("Amartya101/TuBA")
```


## Instructions for using TuBA

There are 3 main functions in TuBA, which need to be employed in a sequential manner on the data set of interest. Make sure before running the first function (*DataPrep*) that the data exists in a tabular format (either a .csv or a .txt file), in which the genes (or more generally, features) are along the rows and the samples are along the columns. The first column in the file must contain the IDs or names of the genes (no duplicates allowed). We have included an example file with RP genes in the [ToyExamples](https://github.com/Amartya101/ToyExamples) repository ("*RPGenes.csv*") for reference. 

The descriptions of the functions and instructions on how to use them are provided below:

### DataPrep
The *DataPrep* function requires 2 inputs - *X* and *normalize*. *X* should be the full name of the input file including the specifier for the format of the file (.csv, .txt, or .tsv). It is expected that the dataset in the input file is in the rectangular format, with genes along rows (first column should contain all the gene names) and samples along the columns (the header should contain the sample IDs). It should contain raw counts (for *normalize* = *T*) and should be processed by the user to remove duplicate genes, NAs etc. Alternatively, *X* can be a matrix with Gene IDs as rownames and Sample IDs as colnames. The *DataPrep* function evaluates if there are genes in the dataset that do not have a single non-zero expression value in any sample or have NAs and removes these genes. The *normalize* argument is a logical input. By default, the *normalize* argument is set to *T* and will perform between-samples normalization of the counts using the normalization method used in DESeq [https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106].  The *DataPrep* function returns a matrix and also generates an output file in .csv format with the name of the input file appended by "*Normalized.csv*" (when *normalize* = *T*) or "*_Cleaned.csv*" (when *normalize* = *F*).

Here is an example of a valid function call

```
DataPrep(X = "RPGenes.csv", normalize = F)
```
This will produce an output file with the name "*RPGenes_Cleaned.csv*" in the current working directory.

### GenePairs

The *GenePairs* function has 5 arguments - *X*, *PercSetSize*, *JcdInd*, *highORlow*, *SampleEnrichment*. *X* should be the name of the normalized (or cleaned) file produced by the *DataCleaning* function, or it can be a matrix with gene IDs as its rownames and sample IDs as its colnames. *PercSetSize* requires an input in terms of percentage of the total number of samples that would constitute the extremals for each gene. For example, if you wish to look at the top 10% you should specify this parameter to be 10 (and *highRlow = "h"*, more below). *JcdInd* specifies the Jaccard index cutoff that the overlap between the percentile sets must satisfy for a gene-pair to be a part of the graph. The parameter *highORlow* specifies whether we want the percentile sets for each gene to correspond to samples with the highest expression levels ("h") or lowest expression levels ("l") respectively. The *SampleEnrichment* argument is a logical input that specifies whether we wish for samples that are over-represented in percentile sets to be filtered out (it does this if set to TRUE). By default it is set to TRUE (recommended). 

The *GenePairs* function generates 3 output files. The first file contains all the gene-pairs that have significant overlaps between their percentile sets (Jaccard indices greater than *JcdInd* specified by the user). The name of this file ends with "*_GenePairs.csv*". This file does not contain the full gene IDs or gene names, instead the genes are labelled by their serial number in the input file. The second file contains a binary matrix with genes along the rows and samples along the columns. The first column contains the gene IDs or gene names obtained from the input file. For each gene (row), the presence of a sample in the percentile set is denoted by a 1, while samples not in the percentile set have 0. The third file this function generates contains the names/IDs of the genes in the data set. The name of this file ends with "*_GeneNames.csv*".

Examples of valid function calls are provided below (here we directly used "*RPGenes.csv*", since it was already clean):
```
# For high expression
GenePairs(X = "RPGenes.csv",PercSetSize = 5,JcdInd = 0.2,highORlow = "h")
# For low expression
GenePairs(X = "RPGenes.csv",PercSetSize = 5,JcdInd = 0.2,highORlow = "l")
```
The first one will generate the following 3 files: "*RPGenes_H0.05_JcdInd0.2_GenePairs.csv*", "*RPGenes_H0.05_JcdInd0.2_GenesSamplesBinaryMatrix.csv*", and "*RPGenes_H0.05_JcdInd0.2_GeneNames.csv*". The annotations in the middle of their names indicate the following: *H0.05* indicates that the gene-pairs were obtained for high expression ("H") with the percentile set size of 5% (0.05); *JcdInd0.2* indicates that a Jaccard index threshold of 0.2 was chosen to shortlist the gene-pairs.

The second one will generate the following 3 files: "*RPGenes_L0.05_JcdInd0.2_GenePairs.csv*", "*RPGenes_L0.05_JcdInd0.2_GenesSamplesBinaryMatrix.csv*", and "*RPGenes_L0.05_JcdInd0.2_GeneNames.csv*". The annotations in the middle of their names indicate the following: *L0.05* indicates that the gene-pairs were obtained for low expression ("L") with the percentile set size of 5% (0.05); *JcdInd0.2* indicates that a Jaccard index threshold of 0.2 was chosen to shortlist the gene-pairs.

### Biclustering

The *Biclustering* function has 5 arguments of which 2 inputs are necessary - *GenePairs* and *BinaryMatrix*. *GenePairs* requires the name of the gene-pairs .csv file generated by the *GenePairs* function; *BinaryMatrix* requires the name of the genes-samples binary matrix .csv file also generated by the *GenePairs* function. The 3 optional arguments include - *MinGenes*, *MinSamples*, and *SampleEnrichment*. *MinGenes* specifies the minimum number of genes desired by the user to be present in a bicluster (default is 3); *MinSamples* specifies the minimum number of samples desired by the user to be present in a bicluster (default is 2, but typically the biclusters have much more samples); *SampleEnrichment* is used to specify the level of enrichment that each sample must exhibit within a given bicluster as compared to its presence in the overall graph. The values for *SampleEnrichment* can range between 0 (excluding 0) and 1. Think of *SampleEnrichment* as a *p*-value - smaller values would require stronger sample associations with the bicluster, and would result in fewer samples in the biclusters. The default *SampleEnrichment* is 1.

*Biclustering* function generates 3 files. The file whose name ends with “*_GenesInBiclusters.csv*” contains list of genes contained in each bicluster. It also provides information about the total number of samples present in the bicluster (column 3), as well how many of these samples were present in the percentile sets of each gene within the bicluster (column 4). The file that ends with “*_BiclusterSamplesMatrix.csv*” is a binary matrix with biclusters along the rows and samples along the columns. For each bicluster, only the samples that are present in the bicluster are given a value of 1. The file ending with “*_GenesBiclusterSamplesMatrix.csv*” contains the genes present in each bicluster along with information about which samples are contributed to the bicluster by each gene (1 for samples contributed by a gene, 0 otherwise). The first column of this file contains the gene IDs/names, while the second column contains the bicluster number that these genes belong to; the rest of the file is comprised of the binary matrix.

Example of a valid function call is provided below:
```
Biclustering(GenePairs = "RPGenes_H0.05_JcdInd0.2_GenePairs.csv",BinaryMatrix = "RPGenes_H0.05_JcdInd0.2_GenesSamplesBinaryMatrix.csv")
```
This will generate the following output files: "*RPGenes_H0.05_JcdInd0.2_MinGenesX_MinSamplesY_GenesInBiclusters.csv*", "*RPGenes_H0.05_JcdInd0.2_MinGenesX_MinSamplesY_BiclusterSamplesMatrix.csv*", and "*RPGenes_H0.05_JcdInd0.2_MinGenesX_MinSamplesY_GenesBiclusterSamplesMatrix.csv*". The additional annotations in the middle of their names indicate the following: *MinGenesX* indicates the number of genes in the bicluster with the fewest genes (for the example with default choices above, *X* will be 3), and *MinSamplesY* indicates the number of samples (*Y*) in the bicluster that has the fewest samples.

### Bicluster Genes Graphs
This function can be used to make the graphs/networks showing the genes in the biclusters found by the *Biclustering* function. The *BiclusterGenesGraph* function has 4 input arguments - *BiclusterGenes*, *GenePairs*, *GeneNames*, and *BiclusterNos*. *BiclusterGenes* specifies the name of the "*.._GenesInBicluster.csv*" file produced by the *Biclustering* function for given data set. *GenePairs* specifies the name of the "*.._GenePairs.csv*" produced by the *GenePairs* function for given data set, while *GeneNames* specifies the name of the "*.._GeneNames.csv*" produced by the *GenePairs* function for given data set. The *BiclusterNos* argument can take in a vector of serial numbers for the respective biclusters for which the graphs/networks are desired.

Example of a valid function call is provided below:
```
BiclusterGenesGraph(BiclusterGenes = "RPGenes_H0.05_JcdInd0.2_GenesInBiclusters.csv",GenePairs = "RPGenes_H0.05_JcdInd0.2_GenePairs.csv",GeneNames = "RPGenes_H0.05_JcdInd0.2_GeneNames.csv",BiclusterNos = c(1,2))
```
This will produce following .pdf files - "*RPGenes_H0.05_JcdInd0.2_Bicluster1.pdf*" and "*RPGenes_H0.05_JcdInd0.2_Bicluster2.pdf*" - containing the graphs for the respective biclusters.

## Authors

* **Amartya Singh** - [Amartya101](https://github.com/Amartya101/)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hossein Khiabanian - [Khiabanian Lab](https://khiabanian-lab.org)
* Gyan Bhanot
* Lodovico Terzi di Bergamo

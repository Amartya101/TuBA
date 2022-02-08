# TuBA

New version of TuBA - Analyzes large graphs to identify biclusters much faster than the older version.

## Getting Started

You will need an R version 3.4.0 or more recent in order to use these functions

### Prerequisites

You need the **data.table** package. You can install it using the following command:

```
install.packages("data.table",dependencies = T)
```

### Installing

Currently, the development version of **TuBA** can be installed using *install_github* from the **devtools** package.

Install **devtools** using:

```
install.packages("devtools",dependencies = T)
```
Once **devtools** is installed
```
devtools::install_github("Amartya101/TuBA2")
```


## Instructions for using TuBA

There are 3 functions in TuBA, which need to be employed in a sequential manner on the data set of interest. Make sure before running the first function (*DataCleaning*) that the data exists in a tabular format (either a .csv or a .txt file), in which the genes (or more generally, features) are along the rows and the samples are along the columns. The first column in the file must contain the IDs or names of the genes (no duplicates allowed). We have included an example file with RP genes in the repository (*RPGenes.csv*) for reference. 

The descriptions of the functions and instructions on how to use them are provided below:

### DataCleaning
The *DataCleaning* function requires one input - File. The name of the input file (.txt or .csv) should be specified in full including the specifier for the format of the file. It is expected that the dataset in the input file is in the rectangular format with genes along rows (first column should contain all the gene names) and samples along the columns (the header should contain the sample IDs). It should contain normalized counts and should be processed by the user to remove duplicate genes, NAs etc. The DatasetPreparationFunction evaluates if there are genes in the dataset that do not have a single non-zero expression value in any sample and eliminates these genes. It renames the header of the first column as Gene.ID and generates the output file in .csv format with the name of the input file appended by "_Cleaned.csv". It also generates a file that contains the names of all the genes that have non-zero expression levels in some of the samples. The name of this file ends with “_GeneNames.csv”.

Here is an example of a valid function call

```
DataCleaning(File = "RPGenes.csv")
```

### GenePairs

The SignificantGenePairsFunction requires three inputs - InputFileName, PercentileCutOff, highORlow. InputFileName should be the name of the cleaned file produced by the DatasetPreparationFunction. PercentileCutOff requires an input in terms of percentage of the total sample size that should be considered for comparing the extremals for each gene. For example, if you wish the percentile set size for comparison to be 10% you should specify this parameter to be 10. The input parameter highORlow requires you to specify whether you want the percentile sets to correspond to high expression or low expression. Use this ONLY if the data has been obtained from an RNASeq platform. In order to specify high expression you can assign any one of these character values to highORlow - “h”, “H”, “high”, “High”. Similarly, for low expression you can assign any one of these character values to highORlow - “l”, “L”, “low”, “Low”.

The SignificantGenePairsFunction generates two files. The first file contains a binary matrix with 1 corresponding to the samples in the percentile set for each gene. The name of the output file contains the name of the input file appended by “_H(PercentileCutOff/100)_Genes_Samples_BinaryMatrix.csv” or “_L(PercentileCutOff/100)_Genes_Samples_BinaryMatrix.csv”, depending on whether you specify “h” or “l”. For example, for the input file “GeneExpressionDataset1_Cleaned.csv” and a percentile cutoff of 5% for high expression (highORlow = “h”), the name of the output file will be “GeneExpressionDataset1_H0.05_Genes_Samples_BinaryMatrix.csv”. The second file generated by the SignificantGenePairsFunction contains all the gene pairs that have significant overlaps between their percentile sets (p < 0.05). The name of this file ends with “_SignificantGenePairs.csv”. The genes are labelled by their serial number in the list of genes produced by the DatasetPreparationFunction.

In addition to these two files, this function also generates 6 plots to help guide the choice of the significance level of overlap. The plot with its ending with “_NoOfEdges.pdf” shows the total number of edges in the graph for different levels of the significance of overlap (in terms of -log10(p)). The plot with its name ending with “_NoOfSamples.pdf” represents the total number of samples present in the graph for different levels of significance of overlaps. The plot with its name ending with “_GenesAdded.pdf” shows the the number of new genes added to the graph with each incremental drop in the level of significance of overlap. Similarly, the plot with its name ending with “_SamplesAdded.pdf” shows the number of new samples added to the graph with every incremental drop in the level of significance of overlap. The plot with its name ending with “_GenesPerEdge.pdf” shows the ratio of the number of new genes added to the graph to the number of new edges added to the graph with every incremental drop in the level of significance of overlap. Similarly, the plot with its name ending with “_SamplesPerEdge.pdf” shows the ratio of the number of new samples added to the graph to the number of new edges added to the graph with every incremental drop in the level of significance of overlap.


```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

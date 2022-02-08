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

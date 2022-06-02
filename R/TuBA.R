
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
utils::globalVariables(c("."))
#' @import utils
#' @import stats
#' @import ggplot2
#' @import network
#' @importFrom ggplot2 aes
#' @import ggnetwork

#' @title  DataPrep
#' @description  This function performs basic cleaning of raw counts data sets - removes genes with zero counts across all samples, as well as genes with NAs in some samples. Additionally, it offers the option (default is T) to perform between-sample normalization using the DESeq approach.
#' @export
#' @param X A character variable or a matrix. Specifies the name of the file that contains the raw counts data (should be in .csv, .txt, or .tsv format). Make sure it is in the genes along rows and samples along columns format. Alternatively, you can provide a matrix with raw counts as input. The matrix must contain the Gene IDs as rownames and Sample IDs as colnames.
#' @param normalize A logical variable. Specifies whether the raw counts should be normalized using the DESeq normalization method (default is T).
#' @return A matrix containing the gene expression data after basic cleaning and normalization (if normalize = T). In addition, it generates a .csv file in the working directory that contains the cleaned and normalized (if normalize = T) gene expression data set.
#' @examples
#' \dontrun{
#' DataPrep(X = "RPGenes.csv",normalize = F)
#' }

DataPrep <- function(X,normalize = T)
{
  if (is.character(X)){

    df <- data.table::fread(X)
    Gene.IDs <- as.character(df[,1])
    Expr.Mat <- as.matrix(df[,-1])
    Sample.IDs <- colnames(df)[-1]

    rownames(Expr.Mat) <- Gene.IDs
    colnames(Expr.Mat) <- Sample.IDs

  } else if (sum(dim(as.matrix(X))) > 2) {
    Expr.Mat <- X
    Gene.IDs <- rownames(Expr.Mat)
    Sample.IDs <- colnames(Expr.Mat)
  } else {
    stop("Invalid input. Please provide either the name of the raw counts file (in .csv, .txt, or .tsv format), or a matrix with gene IDs as rownames and sample IDs as colnames.")
  }

  Sums.of.Rows <- rowSums(Expr.Mat)
  if (length(which(Sums.of.Rows == 0 | is.na(Sums.of.Rows) == T)) != 0) {
    Expr.Mat <- Expr.Mat[-which(Sums.of.Rows == 0 | is.na(Sums.of.Rows) == T),]
    Non.Zero.Genes <- Gene.IDs[-which(Sums.of.Rows == 0 | is.na(Sums.of.Rows) == T)]
  }
  else {
    Non.Zero.Genes <- Gene.IDs
  }

  if (normalize == T){
    #Normalize the raw counts using the DESeq method

    DESeqNormalizationFunction <- function(X){

      Rows.With.No.Zeros <- which(apply(X,1,function(x) mean(x == 0)) == 0)

      #Prepare counts matrix with rows with no zeros to calculate the size factors for each sample (column)

      Y <- X[Rows.With.No.Zeros,]

      #Taking geometric mean of gene expression across all samples - DESeq method
      Geo.Mean.Rows <- vector(mode = "numeric", length = nrow(Y))
      for (i in 1:nrow(Y))
      {
        Geo.Mean.Rows[i] <- exp(mean(log(Y[i,])))
      }

      #Compute size factors for each sample
      Size.Factors.Vec <- vector(mode = "numeric",length = ncol(Y))
      for (i in 1:ncol(Y))
      {
        Size.Factors.Vec[i] <- median(Y[,i]/Geo.Mean.Rows)
      }

      #DESeq normalized original counts matrix

      DESeq.Counts.Matrix <- matrix(0,nrow = nrow(X), ncol = ncol(X))
      for (i in 1:ncol(X))
      {
        DESeq.Counts.Matrix[,i] <- X[,i]/Size.Factors.Vec[i]
      }

      colnames(DESeq.Counts.Matrix) <- colnames(X)

      return(DESeq.Counts.Matrix)
    }
    Expr.Mat <- DESeqNormalizationFunction(X = Expr.Mat)
  }

  Expr.df <- data.frame(Non.Zero.Genes,Expr.Mat)
  colnames(Expr.df) <- c("Gene.ID",Sample.IDs)
  if (is.character(X) & normalize == T){
    FileName <- paste0(substr(X,1,nchar(X) -4),"_Normalized.csv")
  } else if (is.character(X) & normalize == F){
    FileName <- paste0(substr(X,1,nchar(X) -4),"_Cleaned.csv")
  } else if(sum(dim(as.matrix(X))) > 2 & normalize == T) {
    FileName <- paste0("NormalizedCounts",".csv")
  } else {
    FileName <- paste0("CleanedCounts",".csv")
  }
  write.table(Expr.df,file = FileName,row.names = F,col.names = T,sep = ",")
  message("Successfully prepared the output .csv file.")
  rownames(Expr.Mat) <- Non.Zero.Genes
  colnames(Expr.Mat) <- Sample.IDs
  return(Expr.Mat)
}

#' @title  GenePairs
#' @description  This function finds gene-pairs that share more than a specified proportion of samples between their percentile sets.
#' @export
#' @param X A character variable or a matrix. Specifies the name of the .csv or .txt file that contains the normalized counts, or a matrix with Gene IDs as rownames and Sample IDs as colnames.
#' @param PercSetSize A numeric variable. Specifies the percentage of samples that should be in the percentile sets (strictly greater than 0 and less than 40).
#' @param JcdInd A numeric variable. Specifies the minimum Jaccard Index for the overlap between the percentile sets of a given gene-pair.
#' @param highORlow A character variable. Specifies whether the percentile sets correspond to the highest expression samples ("h") or the lowest expression samples ("l").
#' @param SampleFilter A logical variable. If TRUE, filters out samples over-represented in percentile sets. Default is TRUE.
#' @return A data frame containing the gene-pairs whose Jaccard indices are greater than the specified threshold (JcdInd). Instead of gene symbols their serial numbers in the input gene expression file are used in order to save space. In addition, this function generates 3 .csv files in the working directory - one contains the gene-pairs and their Jaccard indices, the second contains the binary matrix (genes along rows, samples along columns) in which the presence of a sample in the percentile set of each gene is indicated by a 1. These 2 files are needed as inputs for the \link[TuBA]{Biclustering} function.
#' The third .csv file contains the IDs or names of the genes corresponding to the numbered indices used in the gene-pairs file.
#' @examples
#' \dontrun{
#' # For high expression
#' GenePairs(X = "RPGenes.csv",PercSetSize = 5,JcdInd = 0.2,highORlow = "h")
#' # For low expression
#' GenePairs(X = "RPGenes.csv",PercSetSize = 5,JcdInd = 0.2,highORlow = "l")
#' }
GenePairs <- function(X,PercSetSize,JcdInd,highORlow,SampleFilter = T)
{
  if (is.character(X)){

    df <- data.table::fread(X)
    colnames(df) <- c("Gene.ID",colnames(df)[-1])
    Gene.Names <- df$Gene.ID
    Expr.Mat <- as.matrix(df[,-1])
    Sample.IDs <- colnames(df)[-1]

    rownames(Expr.Mat) <- Gene.Names
    colnames(Expr.Mat) <- Sample.IDs
    rm(df)

  } else if (sum(dim(as.matrix(X))) > 2) {
    Expr.Mat <- X
    Gene.Names <- rownames(Expr.Mat)
    Sample.IDs <- colnames(Expr.Mat)
  } else {
    stop("Invalid input. Please provide either the name of the raw counts file (in .csv, .txt, or .tsv format), or a matrix with gene IDs as rownames and sample IDs as colnames.")
  }

  if (round(PercSetSize) <= 0 | PercSetSize > 40){
    stop("PercSetSize has to be a numeric value strictly greater than 0 and less than 40.")
  }

  CutOffPerc <- round(PercSetSize)/100

  if (JcdInd <= 0 | JcdInd > 1){
    stop("JcdInd has to be a numeric value strictly between 0 and 1.")
  }

  if (is.logical(SampleFilter) == F){
    message("SampleFilter is a logical variable and can only take T or F as valid inputs. Setting it to T by default")
    SampleFilter <- T
  } else if (is.null(SampleFilter)){
    SampleFilter <- T
  }


  if (highORlow == "h" | highORlow == "H" | highORlow == "high" | highORlow == "High" | highORlow == "Hi" | highORlow == "HI" | highORlow == "hi" | highORlow == "HIGH"){

    Threshold <- 1 - CutOffPerc

    ZeroProp.Per.Feature <- apply(Expr.Mat,1,function(x) mean(x==0))

    Relevant.Features <- which(ZeroProp.Per.Feature < Threshold)

    if (length(Relevant.Features) != 0){
      Gene.Names <- Gene.Names[Relevant.Features]
      Expr.Mat <- Expr.Mat[Relevant.Features,]
    }

    Start.Index <- ncol(Expr.Mat) - (ceiling(CutOffPerc*ncol(Expr.Mat))) +1

    message("Preparing genes-samples binary matrix...")

    #Binary matrix with 1 for samples that are in the percentile set for any given gene
    List.Samples.In.PercSet <- vector(mode = "list",length = nrow(Expr.Mat))
    Row.Index.List <- vector(mode = "list",length = nrow(Expr.Mat))
    Values.List <- vector(mode = "list",length = nrow(Expr.Mat))
    for (i in 1:nrow(Expr.Mat))
    {
      Samples.In.Order <-  order(Expr.Mat[i,])

      Samples.In.PercSet <- Samples.In.Order[Start.Index:ncol(Expr.Mat)]

      Samples.With.Zeros <- Samples.In.PercSet[which(Expr.Mat[i,Samples.In.PercSet] == 0)]

      if (length(Samples.With.Zeros) != 0){
        Samples.In.PercSet <- Samples.In.PercSet[!Samples.In.PercSet %in% Samples.With.Zeros]
      }

      if (length(Samples.In.PercSet) != 0){
        List.Samples.In.PercSet[[i]] <- Samples.In.PercSet
        Row.Index.List[[i]] <- rep(i,length(Samples.In.PercSet))
        Values.List[[i]] <- rep(1,length(Samples.In.PercSet))
      }
    }

    Binary.Mat.For.Genes.Outliers <- Matrix::sparseMatrix(i = unlist(Row.Index.List),j = unlist(List.Samples.In.PercSet),x = unlist(Values.List))

    #Find the frequencies for the samples
    Sample.Frequencies <- Matrix::colSums(Binary.Mat.For.Genes.Outliers)

    #Identify samples that only show up in one percentile set
    NonInformativeSamples <- which(Sample.Frequencies == 1)
    if (length(NonInformativeSamples) != 0){
      Binary.Mat.For.Genes.Outliers <- Binary.Mat.For.Genes.Outliers[,-NonInformativeSamples]
      Sample.IDs <- Sample.IDs[-NonInformativeSamples]
    }

    if (SampleFilter == T){
      #Find the threshold for maximum frequency based on percentile set size
      MaxThreshold <- max(Sample.Frequencies)
      n <- choose(n = nrow(Binary.Mat.For.Genes.Outliers),k = 2)
      a <- choose(n = MaxThreshold,k = 2)
      p.val <- a/n
      while (p.val > CutOffPerc) {
        MaxThreshold <- MaxThreshold -1
        a <- choose(n = MaxThreshold,k = 2)
        p.val <- a/n
      }

      NonInformativeSamples <- which(Sample.Frequencies > MaxThreshold)
      if (length(NonInformativeSamples) != 0){
        Binary.Mat.For.Genes.Outliers <- Binary.Mat.For.Genes.Outliers[,-NonInformativeSamples]
        Sample.IDs <- Sample.IDs[-NonInformativeSamples]
      }
    }

    Genes.Samples.Binary.df <- data.frame(Gene.Names,as.matrix(Binary.Mat.For.Genes.Outliers))
    colnames(Genes.Samples.Binary.df) <- c("Gene.ID",Sample.IDs)

    if (is.character(X)){
      FileName <- paste0(substr(X,1,nchar(X) -4),"_H",CutOffPerc,"_JcdInd",JcdInd,"_GeneSamplesBinaryMatrix.csv")
    } else {
      FileName <- paste0("GeneSamplesBinaryMatrix","_H",CutOffPerc,"_JcdInd",JcdInd,".csv")
    }

    write.table(Genes.Samples.Binary.df,file = FileName,row.names = F,col.names = T,sep = ",")

    message("Successfully generated .csv file containing the genes-samples binary matrix.")

    if (is.character(X)){
      FileName <- paste0(substr(X,1,nchar(X) -4),"_H",CutOffPerc,"_JcdInd",JcdInd,"_GeneNames.csv")
    } else {
      FileName <- paste0("GeneNames","_H",CutOffPerc,"_JcdInd",JcdInd,".csv")
    }

    GeneNames.df <- data.frame(Gene.ID = Genes.Samples.Binary.df$Gene.ID)
    write.table(GeneNames.df,file = FileName,row.names = F,col.names = T,sep = ",")

    rm(Genes.Samples.Binary.df)

    message("Finding gene-pairs with significant overlaps in their percentile sets...")

    #Outlier overlaps matrix
    SampleOverlaps.Matrix <- Matrix::tcrossprod(Binary.Mat.For.Genes.Outliers,Binary.Mat.For.Genes.Outliers)

    SampleOverlaps.Matrix <- as.matrix(SampleOverlaps.Matrix)

    SampleUnion <- ceiling(2*CutOffPerc*ncol(Expr.Mat))
    SampleUnion.Matrix <- SampleUnion - SampleOverlaps.Matrix

    Jaccard.Dist.Mat <- SampleOverlaps.Matrix/SampleUnion.Matrix

    #Find relevant gene-pairs

    Relevant.Gene.Pairs <- which(Jaccard.Dist.Mat >= JcdInd,arr.ind = T)

    if (length(Relevant.Gene.Pairs) >= 2){
      Relevant.Gene.Pairs <- Relevant.Gene.Pairs[which(Relevant.Gene.Pairs[,1] < Relevant.Gene.Pairs[,2]),]

      Gene.Pairs.Col1 <- Relevant.Gene.Pairs[,1]

      Gene.Pairs.Col2 <- Relevant.Gene.Pairs[,2]

      Gene.Pairs.df <- data.frame(Gene.Pairs.Col1,Gene.Pairs.Col2,Jaccard.Dist.Mat[Relevant.Gene.Pairs])

      colnames(Gene.Pairs.df) <- c("Gene.1","Gene.2","Jaccard.Index")

      if (is.character(X)){
        FileName <- paste0(substr(X,1,nchar(X) -4),"_H",CutOffPerc,"_JcdInd",JcdInd,"_GenePairs.csv")
      } else {
        FileName <- paste0("GenePairs","_H",CutOffPerc,"_JcdInd",JcdInd,".csv")
      }

      write.table(Gene.Pairs.df,file = FileName, row.names = F,col.names = T,sep = ",")

      if (length(Gene.Pairs.df$Gene.1) != 0)
        message(paste0("Successfully generated .csv file containing ",length(Gene.Pairs.df$Gene.1)," gene-pairs (edges)."))
    } else {
      stop("No gene-pairs found for given choice of parameters.")
    }

  } else if (highORlow == "l" | highORlow == "L" | highORlow == "low" | highORlow == "Low" | highORlow == "Lo" | highORlow == "lo" | highORlow == "LO" | highORlow == "LOW") {

    Threshold <- CutOffPerc

    ZeroProp.Per.Feature <- apply(Expr.Mat,1,function(x) mean(x==0))

    Relevant.Features <- which(ZeroProp.Per.Feature <= Threshold)

    if (length(Relevant.Features) != 0){
      Gene.Names <- Gene.Names[Relevant.Features]
      Expr.Mat <- Expr.Mat[Relevant.Features,]
    }

    Start.Index <- 1

    message("Preparing genes-samples binary matrix...")

    #Binary matrix with 1 for samples that are in the percentile set for any given gene
    List.Samples.In.PercSet <- vector(mode = "list",length = nrow(Expr.Mat))
    Row.Index.List <- vector(mode = "list",length = nrow(Expr.Mat))
    Values.List <- vector(mode = "list",length = nrow(Expr.Mat))
    for (i in 1:nrow(Expr.Mat))
    {
      Samples.In.Order <-  order(Expr.Mat[i,])

      Samples.In.PercSet <- Samples.In.Order[Start.Index:ceiling(CutOffPerc*ncol(Expr.Mat))]

      if (length(Samples.In.PercSet) != 0){
        List.Samples.In.PercSet[[i]] <- Samples.In.PercSet
        Row.Index.List[[i]] <- rep(i,length(Samples.In.PercSet))
        Values.List[[i]] <- rep(1,length(Samples.In.PercSet))
      }
    }

    Binary.Mat.For.Genes.Outliers <- Matrix::sparseMatrix(i = unlist(Row.Index.List),j = unlist(List.Samples.In.PercSet),x = unlist(Values.List))

    #Find the frequencies for the samples
    Sample.Frequencies <- Matrix::colSums(Binary.Mat.For.Genes.Outliers)

    #Identify samples that only show up in one percentile set
    NonInformativeSamples <- which(Sample.Frequencies == 1)
    if (length(NonInformativeSamples) != 0){
      Binary.Mat.For.Genes.Outliers <- Binary.Mat.For.Genes.Outliers[,-NonInformativeSamples]
      Sample.IDs <- Sample.IDs[-NonInformativeSamples]
    }

    if (SampleFilter == T){
      #Find the threshold for maximum frequency based on percentile set size
      MaxThreshold <- max(Sample.Frequencies)
      n <- choose(n = nrow(Binary.Mat.For.Genes.Outliers),k = 2)
      a <- choose(n = MaxThreshold,k = 2)
      p.val <- a/n
      while (p.val > CutOffPerc) {
        MaxThreshold <- MaxThreshold -1
        a <- choose(n = MaxThreshold,k = 2)
        p.val <- a/n
      }

      NonInformativeSamples <- which(Sample.Frequencies > MaxThreshold)
      if (length(NonInformativeSamples) != 0){
        Binary.Mat.For.Genes.Outliers <- Binary.Mat.For.Genes.Outliers[,-NonInformativeSamples]
        Sample.IDs <- Sample.IDs[-NonInformativeSamples]
      }
    }

    Genes.Samples.Binary.df <- data.frame(Gene.Names,as.matrix(Binary.Mat.For.Genes.Outliers))
    colnames(Genes.Samples.Binary.df) <- c("Gene.ID",Sample.IDs)

    if (is.character(X)){
      FileName <- paste0(substr(X,1,nchar(X) -4),"_L",CutOffPerc,"_JcdInd",JcdInd,"_GeneSamplesBinaryMatrix.csv")
    } else {
      FileName <- paste0("GeneSamplesBinaryMatrix","_L",CutOffPerc,"_JcdInd",JcdInd,".csv")
    }

    write.table(Genes.Samples.Binary.df,file = FileName,row.names = F,col.names = T,sep = ",")

    message("Successfully generated .csv file containing the genes-samples binary matrix.")

    if (is.character(X)){
      FileName <- paste0(substr(X,1,nchar(X) -4),"_L",CutOffPerc,"_JcdInd",JcdInd,"_GeneNames.csv")
    } else {
      FileName <- paste0("GeneNames","_L",CutOffPerc,"_JcdInd",JcdInd,".csv")
    }

    GeneNames.df <- data.frame(Gene.ID = Genes.Samples.Binary.df$Gene.ID)
    write.table(GeneNames.df,file = FileName,row.names = F,col.names = T,sep = ",")

    rm(Genes.Samples.Binary.df)

    message("Finding gene-pairs with significant overlaps in their percentile sets...")

    #Outlier overlaps matrix
    SampleOverlaps.Matrix <- Matrix::tcrossprod(Binary.Mat.For.Genes.Outliers,Binary.Mat.For.Genes.Outliers)

    SampleOverlaps.Matrix <- as.matrix(SampleOverlaps.Matrix)

    SampleUnion <- ceiling(2*CutOffPerc*ncol(Expr.Mat))
    SampleUnion.Matrix <- SampleUnion - SampleOverlaps.Matrix

    Jaccard.Dist.Mat <- SampleOverlaps.Matrix/SampleUnion.Matrix

    #Find Jaccard distance between gene-pairs

    Relevant.Gene.Pairs <- which(Jaccard.Dist.Mat >= JcdInd,arr.ind = T)

    if (length(Relevant.Gene.Pairs) >= 2){
      Relevant.Gene.Pairs <- Relevant.Gene.Pairs[which(Relevant.Gene.Pairs[,1] < Relevant.Gene.Pairs[,2]),]
      Gene.Pairs.Col1 <- Relevant.Gene.Pairs[,1]

      Gene.Pairs.Col2 <- Relevant.Gene.Pairs[,2]

      Gene.Pairs.df <- data.frame(Gene.Pairs.Col1,Gene.Pairs.Col2,Jaccard.Dist.Mat[Relevant.Gene.Pairs])

      colnames(Gene.Pairs.df) <- c("Gene.1","Gene.2","Jaccard.Index")

      if (is.character(X)){
        FileName <- paste0(substr(X,1,nchar(X) -4),"_L",CutOffPerc,"_JcdInd",JcdInd,"_GenePairs.csv")
      } else {
        FileName <- paste0("GenePairs","_L",CutOffPerc,"_JcdInd",JcdInd,".csv")
      }

      write.table(Gene.Pairs.df,file = FileName, row.names = F,col.names = T,sep = ",")

      if (length(Gene.Pairs.df$Gene.1) != 0)
        message(paste0("Successfully generated .csv file containing ",length(Gene.Pairs.df$Gene.1)," gene-pairs (edges)."))
    } else {
      stop("No gene-pairs found for given choice of parameters.")
    }
  }
  return(Gene.Pairs.df)
}

#' @title  Biclustering Function
#' @description  This function finds the biclusters using the graph (gene-pairs file) and the genes-samples binary matrix generated by the \link[TuBA]{GenePairs} function.
#' @export
#' @param GenePairs A character variable. Specifies the name of the gene-pairs .csv file generated by the \link[TuBA]{GenePairs} function.
#' @param BinaryMatrix A character variable. Specifies the name of the genes-samples binary matrix .csv file generated by the \link[TuBA]{GenePairs} function.
#' @param MinGenes A numeric (integer) variable. Specifies the minimum number of genes that the biclusters must contain.
#' @param MinSamples A numeric (integer) variable. Specifies the minimum number of samples that the biclusters must contain.
#' @param SampleEnrichment A numeric variable - should be greater than 0 and less than or equal to 1. Specifies to what extent the samples should be enriched in the bicluster. Smaller values indicate higher enrichment. Default is 1.
#' @return A data frame containing the sets of genes in the biclusters along with the information about the number of samples in each bicluster. In addition, it generates 3 .csv files - one contains the list of genes in each bicluster along with information about the total number of samples in the bicluster and how many of them are contributed by each gene within the bicluster; the other contains the bicluster-samples binary matrix (biclusters along rows, samples along columns) in which the presence of a sample within a bicluster is indicated by a 1; the third file contains the genes-biclusters-samples matrix (genes along rows, samples along columns) which contains more detailed information about which sample is present for which gene within a given bicluster. The first column in this genes-biclusters-samples file contains the names of genes in the bicluster, and the second column contains the serial numbers of the biclusters which contains these genes; the rest of the file contains the binary matrix.
#' @examples
#' \dontrun{
#' Biclustering(GenePairs = "RPGenes_H0.05_JcdInd0.2_GenePairs.csv",BinaryMatrix = "RPGenes_H0.05_JcdInd0.2_GenesSamplesBinaryMatrix.csv")
#' }
Biclustering <- function(GenePairs,BinaryMatrix,MinGenes = NULL,MinSamples = NULL,SampleEnrichment = NULL)
{
  if(is.null(MinGenes)){
    MinGenes <- 3
  }

  if(MinGenes < 3){
    MinGenes <- 3
    message("MinGenes cannot be less than 3. Minimum of 3 set as default.")
  }

  if(is.null(MinSamples)){
    MinSamples <- 2
  }

  if (is.null(SampleEnrichment)){
    SampleEnrichment <- 0.05
  }

  #Import data frame that contains the variable pairs found by the significant node pairs function
  Variable.Pairs.df <- data.table::fread(GenePairs)
  colnames(Variable.Pairs.df) <- c("Variable.1","Variable.2","JcdInd")

  JaccardInd <- round(min(Variable.Pairs.df$JcdInd),digits = 2)

  message("Importing files..")

  #Import data frame that contains the binary matrix of variables and their respective percentile set samples
  Variables.Samples.Binary.df <- data.table::fread(BinaryMatrix)
  colnames(Variables.Samples.Binary.df) <- c("Variable.ID",colnames(Variables.Samples.Binary.df[,-1]))

  #Names of variables
  Variable.Names <- as.character(Variables.Samples.Binary.df$Variable.ID)

  #Binary matrix of variables and percentile set samples
  Matrix.For.Variables.Outliers <- as.matrix(Variables.Samples.Binary.df[,-1])
  Matrix.For.Variables.Outliers <- Matrix::Matrix(Matrix.For.Variables.Outliers,sparse = T)

  #IDs of samples (or conditions)
  Sample.IDs <- colnames(Matrix.For.Variables.Outliers)

  rm(Variables.Samples.Binary.df)

  if(length(Variable.Pairs.df$Variable.1) == 0){
    stop("Please check input file. No variable pairs found.")
  } else {
    message(paste0("There are ",length(Variable.Pairs.df$Variable.1)," edges in the graph."))
  }

  #Column 1 of qualified variable pairs data frame
  Qualified.Variable1 <- Variable.Pairs.df$Variable.1

  #Column 2 of qualified variable pairs data frame
  Qualified.Variable2 <- Variable.Pairs.df$Variable.2

  #Summarize graph info in a table
  Variable.Summary.Info <- table(c(Qualified.Variable1,Qualified.Variable2))

  #All variables in graph
  All.Variables.In.Graph <- as.numeric(names(Variable.Summary.Info))

  rm(Variable.Pairs.df)

  message("Preparing graph...")

  #Vector to convert variable serial numbers to the corresponding serial number in the reduced adjacency matrix
  Variable.Annotation.Conversion.Vec <- vector(mode = "numeric", length = length(Variable.Names))
  Variable.Ser.Nos <- 1:length(Variable.Names)
  Matching.Variables.Indices <- match(All.Variables.In.Graph,Variable.Ser.Nos)
  Variable.Annotation.Conversion.Vec[Matching.Variables.Indices] <- 1:length(All.Variables.In.Graph)

  #Prepare adjacency matrix for the graph
  Adjacency.Mat.Variables <- matrix(0,nrow = length(All.Variables.In.Graph),ncol = length(All.Variables.In.Graph))
  Adjacency.Mat.Variables[cbind(Variable.Annotation.Conversion.Vec[Qualified.Variable1],Variable.Annotation.Conversion.Vec[Qualified.Variable2])] <- 1

  #Make the adjacency matrix symmetric
  Adjacency.Mat.Variables <- Adjacency.Mat.Variables + t(Adjacency.Mat.Variables)
  rownames(Adjacency.Mat.Variables) <- All.Variables.In.Graph
  colnames(Adjacency.Mat.Variables) <- All.Variables.In.Graph

  if (length(Qualified.Variable1) > 200000)
    message("This may take several minutes due to the large size of the graph.")

  Total.No.of.Edges.In.Unpruned.Graph <- length(Qualified.Variable1)

  #Find the variables in whose percentile sets each sample shows up
  BinaryMatrix.For.Graph.Variables <- as.matrix(Matrix.For.Variables.Outliers[All.Variables.In.Graph,])
  Variables.Per.Sample <- vector(mode = "list",length = ncol(Matrix.For.Variables.Outliers))
  for (i in 1:ncol(Matrix.For.Variables.Outliers))
  {
    Variables.Per.Sample[[i]] <- All.Variables.In.Graph[which(BinaryMatrix.For.Graph.Variables[,i] == 1)]
  }

  #Find sample background counts along edges in graph
  Samples.Background.Frequencies <- vector(mode = "numeric", length = length(Sample.IDs))
  for (i in 1:length(Samples.Background.Frequencies))
  {
    Nodes.Per.Sample <- Variables.Per.Sample[[i]]
    if (length(Nodes.Per.Sample) >= 2){
      Sub.Adj.Mat <- Adjacency.Mat.Variables[Variable.Annotation.Conversion.Vec[Nodes.Per.Sample],Variable.Annotation.Conversion.Vec[Nodes.Per.Sample]]
      Samples.Background.Frequencies[i] <- sum(Sub.Adj.Mat)/2
    } else{
      Samples.Background.Frequencies[i] <- 0
    }
  }

  message("Finding biclusters...")

  #Find the number of triangular cliques associated with each node-pair in graph
  No.of.Nodes.Associated <- vector(mode = "numeric",length = length(Qualified.Variable1))
  for (i in 1:length(Qualified.Variable1))
  {
    TempI1 <- Qualified.Variable1[i]
    TempI2 <- Qualified.Variable2[i]

    Col.Sums.Vec <- colSums(Adjacency.Mat.Variables[c(Variable.Annotation.Conversion.Vec[TempI1],Variable.Annotation.Conversion.Vec[TempI2]),])
    No.of.Nodes.Associated[i] <- length(which(Col.Sums.Vec == 2))
  }

  #Non-triangular gene-pairs
  Node.Pairs.For.Filtering <- which(No.of.Nodes.Associated == 0)

  if (length(Node.Pairs.For.Filtering) != 0){
    No.of.Nodes.Associated <- No.of.Nodes.Associated[-Node.Pairs.For.Filtering]

    Qualified.Variable1 <- Qualified.Variable1[-Node.Pairs.For.Filtering]
    Qualified.Variable2 <- Qualified.Variable2[-Node.Pairs.For.Filtering]
  }

  if(max(No.of.Nodes.Associated) == 0){
    stop("No clique of size 3 found in input graph.")
  }

  Decreasing.No.of.Nodes <- order(No.of.Nodes.Associated,decreasing = T)

  #Sort column 1 and column2 based on decreasing order of nodes associated
  Qualified.Variable1 <- Qualified.Variable1[Decreasing.No.of.Nodes]
  Qualified.Variable2 <- Qualified.Variable2[Decreasing.No.of.Nodes]

  #Find potential biclusters
  i <- 1
  Temp.Vec1 <- Qualified.Variable1
  Temp.Vec2 <- Qualified.Variable2
  Nodes.In.Subgraph <- vector(mode = "list")
  Temp.Nodes.Vec <- c()
  while (length(Temp.Vec1) != 0){
    Sub.Adj.Mat <- Adjacency.Mat.Variables[c(Variable.Annotation.Conversion.Vec[Temp.Vec1[1]],Variable.Annotation.Conversion.Vec[Temp.Vec2[1]]),]
    Nodes.In.Subgraph[[i]] <- c(Temp.Vec1[1],Temp.Vec2[1],All.Variables.In.Graph[colSums(Sub.Adj.Mat) == 2])
    Nodes.In.Subgraph[[i]] <- Nodes.In.Subgraph[[i]][!Nodes.In.Subgraph[[i]] %in% Temp.Nodes.Vec]
    Temp.Nodes.Vec <- c(Temp.Nodes.Vec,Nodes.In.Subgraph[[i]])

    #Remove those edges that contain the variables in the dense subgraph
    Edges.To.Be.Removed <- which(Temp.Vec1 %in% Nodes.In.Subgraph[[i]] | Temp.Vec2 %in% Nodes.In.Subgraph[[i]])

    Temp.Vec1 <- Temp.Vec1[-Edges.To.Be.Removed]
    Temp.Vec2 <- Temp.Vec2[-Edges.To.Be.Removed]

    i <- i + 1
  }

  if (length(Nodes.In.Subgraph) != 0){
    message("Dense subgraphs identified.")
  } else {
    stop("No subgraph with at least 3 nodes identified.")
  }

  #Reintroduce dense subgraphs back in original graph and add nodes that share edges with at least 2 nodes in dense subgraph
  Nodes.In.Bicluster <- vector(mode = "list",length = length(Nodes.In.Subgraph))
  for (i in 1:length(Nodes.In.Subgraph))
  {
    Col.Ser.Nos.For.Bicluster <- Variable.Annotation.Conversion.Vec[Nodes.In.Subgraph[[i]]]
    Sub.Adj.Mat <- Adjacency.Mat.Variables[,Col.Ser.Nos.For.Bicluster]
    Nodes.In.Bicluster[[i]] <- unique(c(Nodes.In.Subgraph[[i]],as.numeric(rownames(Sub.Adj.Mat)[which(rowSums(Sub.Adj.Mat) >= 2)])))
  }

  #Find biclusters that have nodes that are subsets of other biclusters
  Nested.Biclusters <- vector(mode = "numeric")
  for (i in 1:length(Nodes.In.Bicluster))
  {
    TempI <- length(Nodes.In.Bicluster) - (i-1)
    Temp.Bicluster.I <- Nodes.In.Bicluster[[TempI]]
    j <- 1
    Temp.Intersection <- 1
    while(Temp.Intersection != 0 & TempI != j){
      Temp.Bicluster.J <- Nodes.In.Bicluster[[j]]
      Temp.Intersection <- length(Temp.Bicluster.I) - length(intersect(Temp.Bicluster.I,Temp.Bicluster.J))
      if (Temp.Intersection == 0 & TempI != j & length(Temp.Bicluster.J) >= length(Temp.Bicluster.I)){
        Nested.Biclusters <- c(Nested.Biclusters,TempI)
      }
      j <- j + 1
    }
  }

  #Remove biclusters that are nested within larger biclusters
  if (length(Nested.Biclusters) != 0){
    Nodes.In.Bicluster <- Nodes.In.Bicluster[-Nested.Biclusters]
    Nodes.In.Subgraph <- Nodes.In.Subgraph[-Nested.Biclusters]
  }

  if (length(Nodes.In.Bicluster) > 1){
    Bicluster.Size.Order <- order(unlist(lapply(Nodes.In.Bicluster,length)),decreasing = T)
  } else if (length(Nodes.In.Bicluster) < 2 & length(Nodes.In.Bicluster[[1]]) != 0){
    Bicluster.Size.Order <- 1
  } else {
    stop("No bicluster found with given choice of parameters.")
  }

  Nodes.In.Bicluster <- Nodes.In.Bicluster[Bicluster.Size.Order]
  Nodes.In.Subgraph <- Nodes.In.Subgraph[Bicluster.Size.Order]

  #Find samples preferentially associated with biclusters
  n <- Total.No.of.Edges.In.Unpruned.Graph
  Samples.In.Bicluster <- vector(mode = "list",length = length(Nodes.In.Bicluster))
  for (i in 1:length(Nodes.In.Bicluster))
  {
    Temp.Nodes.In.Bicluster <- Nodes.In.Bicluster[[i]]

    Sub.Adj.Mat <- Adjacency.Mat.Variables[Variable.Annotation.Conversion.Vec[Temp.Nodes.In.Bicluster],Variable.Annotation.Conversion.Vec[Temp.Nodes.In.Bicluster]]
    No.of.Edges.In.Bicluster <- sum(Sub.Adj.Mat)/2

    Bicluster.Samples.Frequencies <- vector(mode = "numeric",length = length(Sample.IDs))
    for (j in 1:length(Bicluster.Samples.Frequencies))
    {
      Nodes.Per.Sample <- Variables.Per.Sample[[j]]
      Nodes.Per.Sample.In.Bicluster <- intersect(Nodes.Per.Sample,Temp.Nodes.In.Bicluster)
      if (length(Nodes.Per.Sample.In.Bicluster) > 0){
        Sub.Adj.Mat <- Adjacency.Mat.Variables[Variable.Annotation.Conversion.Vec[Nodes.Per.Sample.In.Bicluster],Variable.Annotation.Conversion.Vec[Nodes.Per.Sample.In.Bicluster]]
        Bicluster.Samples.Frequencies[j] <- sum(Sub.Adj.Mat)/2
      } else {
        Bicluster.Samples.Frequencies[j] <- 0
      }
    }

    Valid.Sample.Ser.Nos <- which(Bicluster.Samples.Frequencies != 0)
    Bicluster.Samples.Frequencies <- Bicluster.Samples.Frequencies[Valid.Sample.Ser.Nos]
    Order.Bicluster.Samples.Frequencies <- order(Bicluster.Samples.Frequencies,decreasing = T)
    Bicluster.Samples.Frequencies <- Bicluster.Samples.Frequencies[Order.Bicluster.Samples.Frequencies]
    Valid.Sample.Ser.Nos <- Valid.Sample.Ser.Nos[Order.Bicluster.Samples.Frequencies]

    if (SampleEnrichment == 1){
      Sample.Ratios.Per.Bicluster <- vector(mode = "numeric",length = length(Valid.Sample.Ser.Nos))
      for (k in 1:length(Valid.Sample.Ser.Nos))
      {
        t <- Bicluster.Samples.Frequencies[k]
        Sample.Count.In.Graph <- Samples.Background.Frequencies[Valid.Sample.Ser.Nos[k]]

        Sample.Ratios.Per.Bicluster[k] <- (t/Sample.Count.In.Graph)/(No.of.Edges.In.Bicluster/n)
      }

      Temp.ser.nos <- which(Sample.Ratios.Per.Bicluster >= 1)
      Samples.In.Bicluster[[i]] <- Valid.Sample.Ser.Nos[Temp.ser.nos]
    } else {
      p.value.sample <- vector(mode = "numeric",length = length(Valid.Sample.Ser.Nos))
      Sample.Ratios.Per.Bicluster <- vector(mode = "numeric",length = length(Valid.Sample.Ser.Nos))
      for (k in 1:length(Valid.Sample.Ser.Nos))
      {
        Temp.Sample.Frequency <- Bicluster.Samples.Frequencies[k]
        t <- Temp.Sample.Frequency

        Sample.Count.In.Graph <- Samples.Background.Frequencies[Valid.Sample.Ser.Nos[k]]
        if (No.of.Edges.In.Bicluster <= Sample.Count.In.Graph){
          b <- No.of.Edges.In.Bicluster
          a <- Sample.Count.In.Graph
        } else {
          a <- No.of.Edges.In.Bicluster
          b <- Sample.Count.In.Graph
        }
        p.value.sample[k] <- sum(dhyper(t:b,a,n-a,b))
        Sample.Ratios.Per.Bicluster[k] <- (t/Sample.Count.In.Graph)/(No.of.Edges.In.Bicluster/n)
      }
      Temp.ser.nos <- which(p.value.sample <= SampleEnrichment & Sample.Ratios.Per.Bicluster >= 1)
      Samples.In.Bicluster[[i]] <- Valid.Sample.Ser.Nos[Temp.ser.nos]
    }
  }

  #Filter out biclusters that have fewer genes than the specified threshold (MinGenes)
  if (length(which(unlist(lapply(Nodes.In.Bicluster,length)) >= MinGenes)) == 0){
    stop("No biclusters with ",MinGenes," genes found. Please choose different parameters.")
  } else {
    SatisfactoryMinSizeGenes.Biclusters <- which(unlist(lapply(Nodes.In.Bicluster,length)) >= MinGenes)
    Nodes.In.Bicluster <- Nodes.In.Bicluster[SatisfactoryMinSizeGenes.Biclusters]
    Nodes.In.Subgraph <- Nodes.In.Subgraph[SatisfactoryMinSizeGenes.Biclusters]
    Samples.In.Bicluster <- Samples.In.Bicluster[SatisfactoryMinSizeGenes.Biclusters]
  }

  #Find biclusters that contain at least the minimum of samples specified (MinSamples)
  Biclusters.With.Some.Samples <- which(unlist(lapply(Samples.In.Bicluster,length)) >= MinSamples)

  if (length(Biclusters.With.Some.Samples) == 0 & SampleEnrichment != 1){
    stop("No biclusters were found with the given choice of SampleEnrichment. Try with SampleEnrichment = 1")
  } else if (length(Biclusters.With.Some.Samples) == 0 & SampleEnrichment == 1){
    stop("No bicluster found with at least 3 genes")
  } else {
    message("Samples enriched in biclusters identified.")
  }

  #Filter nodes in biclusters based on the samples found enriched in each bicluster - Approach 1 (Remove nodes based on the overlap of their percentile sets with samples found enriched in the bicluster)
  Nodes.In.Bicluster <- Nodes.In.Bicluster[Biclusters.With.Some.Samples]
  Nodes.In.Subgraph <- Nodes.In.Subgraph[Biclusters.With.Some.Samples]
  Samples.In.Bicluster <- Samples.In.Bicluster[Biclusters.With.Some.Samples]

  #Filter nodes based on their association with samples in subgraph
  Matrix.For.Variables.Outliers <- as.matrix(Matrix.For.Variables.Outliers)

  Nodes.In.Final.Biclusters <- vector(mode = "list",length = length(Nodes.In.Bicluster))
  N.PercentileSet <- max(rowSums(Matrix.For.Variables.Outliers))
  for (i in 1:length(Nodes.In.Bicluster))
  {
    Temp.Nodes.In.Bicluster <- Nodes.In.Bicluster[[i]]
    JInd <- vector(mode = "numeric",length = length(Temp.Nodes.In.Bicluster))
    for (j in 1:length(Temp.Nodes.In.Bicluster))
    {
      Samples.In.Percentile.SetJ <- which(Matrix.For.Variables.Outliers[Temp.Nodes.In.Bicluster[j],] == 1)
      Intersecting.Samples <- intersect(Samples.In.Bicluster[[i]],Samples.In.Percentile.SetJ)
      JInd[j] <- length(Intersecting.Samples)/(2*N.PercentileSet -  length(Intersecting.Samples))
    }
    Temp.Filter.Index <- which(JInd < JaccardInd)
    if (length(Temp.Filter.Index != 0)){
      Nodes.In.Final.Biclusters[[i]] <- Nodes.In.Bicluster[[i]][-Temp.Filter.Index]
    } else {
      Nodes.In.Final.Biclusters[[i]] <- Nodes.In.Bicluster[[i]]
    }
  }

  #Find biclusters that have nodes that are subsets of other biclusters
  Nested.Biclusters <- vector(mode = "numeric")
  for (i in 1:length(Nodes.In.Final.Biclusters))
  {
    TempI <- length(Nodes.In.Final.Biclusters) - (i-1)
    Temp.Bicluster.I <- Nodes.In.Final.Biclusters[[TempI]]
    j <- 1
    Temp.Intersection <- 1
    while(Temp.Intersection != 0 & TempI != j){
      Temp.Bicluster.J <- Nodes.In.Final.Biclusters[[j]]
      Temp.Intersection <- length(Temp.Bicluster.I) - length(intersect(Temp.Bicluster.I,Temp.Bicluster.J))
      if (Temp.Intersection == 0 & TempI != j & length(Temp.Bicluster.J) >= length(Temp.Bicluster.I)){
        Nested.Biclusters <- c(Nested.Biclusters,TempI)
      }
      j <- j + 1
    }
  }

  #Remove biclusters that are nested within larger biclusters
  if (length(Nested.Biclusters) != 0){
    Nodes.In.Final.Biclusters <- Nodes.In.Final.Biclusters[-Nested.Biclusters]
    Nodes.In.Subgraph <- Nodes.In.Subgraph[-Nested.Biclusters]
    Samples.In.Bicluster <- Samples.In.Bicluster[-Nested.Biclusters]
  }

  if (length(Nodes.In.Final.Biclusters) > 1){
    Bicluster.Size.Order <- order(unlist(lapply(Nodes.In.Final.Biclusters,length)),decreasing = T)
  } else if (length(Nodes.In.Final.Biclusters) < 2 & length(Nodes.In.Final.Biclusters[[1]]) != 0){
    Bicluster.Size.Order <- 1
  } else {
    stop("No bicluster found with given choice of parameters.")
  }

  Nodes.In.Final.Biclusters <- Nodes.In.Final.Biclusters[Bicluster.Size.Order]
  Nodes.In.Subgraph <- Nodes.In.Subgraph[Bicluster.Size.Order]
  Samples.In.Bicluster <- Samples.In.Bicluster[Bicluster.Size.Order]

  #Filter out biclusters that have fewer genes than the specified threshold (MinGenes)
  if (length(which(unlist(lapply(Nodes.In.Final.Biclusters,length)) >= MinGenes)) == 0){
    stop("No biclusters with ",MinGenes," genes found. Please choose different parameters.")
  } else {
    SatisfactoryMinSizeGenes.Biclusters <- which(unlist(lapply(Nodes.In.Final.Biclusters,length)) >= MinGenes)
    Nodes.In.Final.Biclusters <- Nodes.In.Final.Biclusters[SatisfactoryMinSizeGenes.Biclusters]
    Nodes.In.Subgraph <- Nodes.In.Subgraph[SatisfactoryMinSizeGenes.Biclusters]
    Samples.In.Bicluster <- Samples.In.Bicluster[SatisfactoryMinSizeGenes.Biclusters]
  }

  message("Nodes in final biclusters identified.")

  #Make bicluster-samples matrix
  Bicluster.Samples.Matrix <- matrix(0,nrow = length(Samples.In.Bicluster),ncol = ncol(Matrix.For.Variables.Outliers))
  for (i in 1:length(Samples.In.Bicluster))
  {
    Bicluster.Samples.Matrix[i,Samples.In.Bicluster[[i]]] <- 1
  }

  if (length(dim(Bicluster.Samples.Matrix)) == 0){   #For the case where only one bicluster is found
    MinBiclusterSamples <- sum(Bicluster.Samples.Matrix)
    No.of.Samples <- length(Bicluster.Samples.Matrix)
  } else if (length(rowSums(Bicluster.Samples.Matrix)) != 0){
    MinBiclusterSamples <- min(rowSums(Bicluster.Samples.Matrix))
    No.of.Samples <- ncol(Bicluster.Samples.Matrix)
  } else {
    MinBiclusterSamples <- 0
  }

  Temp.No.of.Nodes.In.Bicluster <- unlist(lapply(Nodes.In.Final.Biclusters,length))
  if (length(Temp.No.of.Nodes.In.Bicluster) != 0){
    MinBiclusterGenes <- min(Temp.No.of.Nodes.In.Bicluster)
  } else {
    MinBiclusterGenes <- 0
  }

  OriginalMinSamples <- MinSamples
  if (MinBiclusterSamples > MinSamples)
    MinSamples <- MinBiclusterSamples

  if (OriginalMinSamples != 2){
    message(paste0("Found ",length(Nodes.In.Final.Biclusters)," biclusters with at least ",MinGenes," genes and ",OriginalMinSamples," samples."))
  } else {
    message(paste0("Found ",length(Nodes.In.Final.Biclusters)," biclusters with at least ",MinGenes," genes and ",MinSamples," samples."))
  }

  #Prepare the output files
  if (length(Nodes.In.Final.Biclusters) != 0){

    Nodes.Biclusters.Info.df <- data.frame(Variable.Names[unlist(Nodes.In.Final.Biclusters)],rep(1:length(Nodes.In.Final.Biclusters),unlist(lapply(Nodes.In.Final.Biclusters,length))))
    colnames(Nodes.Biclusters.Info.df) <- c("Gene.ID","Bicluster.No")

    No.of.Samples.In.Bicluster <- vector(mode = "numeric",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Samples.Per.Gene.In.Bicluster <- vector(mode = "list",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    No.of.Samples.Per.Gene.In.Bicluster <- vector(mode = "numeric",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Proportion.of.Samples.In.Bicluster <- vector(mode = "numeric",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Genes.Bicluster.Samples.Matrix <- matrix(0,nrow = length(Nodes.Biclusters.Info.df$Gene.ID),ncol = No.of.Samples)
    for (i in 1:length(Nodes.Biclusters.Info.df$Gene.ID))
    {
      Gene.Ser.No <- which(Variable.Names == Nodes.Biclusters.Info.df$Gene.ID[i])
      Bicluster.No <- Nodes.Biclusters.Info.df$Bicluster.No[i]
      if (length(dim(Bicluster.Samples.Matrix)) == 0){
        No.of.Samples.In.Bicluster[i] <- length(which(Bicluster.Samples.Matrix == 1))
        Samples.Per.Gene.In.Bicluster[[i]] <- intersect(which(Matrix.For.Variables.Outliers[Gene.Ser.No,] == 1),which(Bicluster.Samples.Matrix == 1))
      } else {
        No.of.Samples.In.Bicluster[i] <- length(which(Bicluster.Samples.Matrix[Bicluster.No,] == 1))
        Samples.Per.Gene.In.Bicluster[[i]] <- intersect(which(Matrix.For.Variables.Outliers[Gene.Ser.No,] == 1),which(Bicluster.Samples.Matrix[Bicluster.No,] == 1))
      }

      No.of.Samples.Per.Gene.In.Bicluster[i] <- length(Samples.Per.Gene.In.Bicluster[[i]])
      Proportion.of.Samples.In.Bicluster[i] <- No.of.Samples.Per.Gene.In.Bicluster[i]/No.of.Samples.In.Bicluster[i]
      Genes.Bicluster.Samples.Matrix[i,Samples.Per.Gene.In.Bicluster[[i]]] <- 1
    }

    Nodes.Biclusters.Info.df$Samples.In.Bicluster <- No.of.Samples.In.Bicluster
    Nodes.Biclusters.Info.df$Samples.Per.Gene.In.Bicluster <- No.of.Samples.Per.Gene.In.Bicluster
    Nodes.Biclusters.Info.df$Proportion.of.Samples <- Proportion.of.Samples.In.Bicluster

    File.Name <- paste0(substr(GenePairs,1,nchar(GenePairs)-13),"MinGenes",MinBiclusterGenes,"_MinSamples",MinBiclusterSamples,"_GenesInBiclusters.csv")

    data.table::fwrite(Nodes.Biclusters.Info.df,file = File.Name,row.names = F,col.names = T,sep = ",")

    message("Successfully generated .csv file containing the list of genes in biclusters.")

    Bicluster.Samples.df <- data.frame(paste0("Bicluster.",1:length(Nodes.In.Final.Biclusters)),Bicluster.Samples.Matrix)
    colnames(Bicluster.Samples.df) <- c("Bicluster.No",Sample.IDs)

    File.Name <- paste0(substr(GenePairs,1,nchar(GenePairs)-13),"MinGenes",MinBiclusterGenes,"_MinSamples",MinBiclusterSamples,"_BiclusterSamplesMatrix.csv")

    data.table::fwrite(Bicluster.Samples.df,file = File.Name,row.names = F,col.names = T,sep = ",")

    if (nrow(Bicluster.Samples.Matrix) != 0)
      message("Successfully generated .csv file containing the bicluster-samples binary matrix.")

    Genes.Bicluster.Samples.df <- data.frame(Nodes.Biclusters.Info.df$Gene.ID,Nodes.Biclusters.Info.df$Bicluster.No,Genes.Bicluster.Samples.Matrix)
    colnames(Genes.Bicluster.Samples.df) <- c("Gene.ID","Bicluster.No",Sample.IDs)

    File.Name <- paste0(substr(GenePairs,1,nchar(GenePairs)-13),"MinGenes",MinBiclusterGenes,"_MinSamples",MinBiclusterSamples,"_GenesBiclusterSamplesMatrix.csv")

    data.table::fwrite(Genes.Bicluster.Samples.df,file = File.Name,row.names = F,col.names = T,sep = ",")

    if (nrow(Bicluster.Samples.Matrix) != 0)
      message("Successfully generated .csv file containing the genes-bicluster-samples binary matrix.")
  }
  return(Nodes.Biclusters.Info.df)
}

#' @title  Bicluster Genes Graph Function
#' @description  This function can be used to make graphs/networks containing the genes in the biclusters found by the \link[TuBA]{Biclustering} function.
#' @export
#' @param BiclusterGenes A character variable. Specifies the name of the .csv file containing the lists of genes in the biclusters generated by the \link[TuBA]{Biclustering} function.
#' @param GenePairs A character variable. Specifies the name of the genes-pairs .csv file generated by the \link[TuBA]{GenePairs} function.
#' @param GeneNames A character variable. Specifies the name of the .csv file containing the gene names/IDs generated by the \link[TuBA]{GenePairs} function.
#' @param BiclusterNos A numeric vector. Specifies the serial numbers of the biclusters for which the bicluster graphs/networks are desired.
#' @return Generates .pdf files containing the graphs/networks of genes corresponding to the biclusters specified by the user.
#' @examples
#' \dontrun{
#' BiclusterGenesGraph(BiclusterGenes = "RPGenes_H0.05_JcdInd0.2_GenesInBiclusters.csv",GenePairs = "RPGenes_H0.05_JcdInd0.2_GenePairs.csv",GeneNames = "RPGenes_H0.05_JcdInd0.2_GeneNames.csv",BiclusterNos = c(1,2))
#' }

BiclusterGenesGraph <- function(BiclusterGenes,GenePairs,GeneNames,BiclusterNos)
{
  if(BiclusterGenes == ""){
    stop("Please specify the name of the .csv file containing the lists of genes in biclusters.")
  }

  if(GenePairs == ""){
    stop("Please specify the name of the .csv file containing the gene-pairs.")
  }

  if(GeneNames == ""){
    stop("Please specify the name of the .csv file containing the names/IDs of the genes.")
  }

  if(length(BiclusterNos) == 0 | length(which(BiclusterNos <= 0) != 0)){
    stop("Bicluster serial number(s) must be specified, and should be present in the bicluster results.")
  }

  #Import bicluster genes result
  BiclusterGenes.df <- data.table::fread(BiclusterGenes)

  #Import gene-pairs file
  GenePairs.df <- data.table::fread(GenePairs)

  #Import gene names
  GeneNames.Vec <- data.table::fread(GeneNames)$Gene.ID

  #Replace serial numbers with corresponding gene names
  GenePairs.df$Gene.1 <- GeneNames.Vec[GenePairs.df$Gene.1]
  GenePairs.df$Gene.2 <- GeneNames.Vec[GenePairs.df$Gene.2]

  #Make graphs for the bicluster serial numbers provided by user

  for (i in 1:length(BiclusterNos))
  {
    TempI <- BiclusterNos[i]
    BiclusterGenes <- BiclusterGenes.df$Gene.ID[BiclusterGenes.df$Bicluster.No == TempI]

    Col1.Edgelist <- as.character(GenePairs.df$Gene.1[GenePairs.df$Gene.1 %in% BiclusterGenes & GenePairs.df$Gene.2 %in% BiclusterGenes])
    Col2.Edgelist <- as.character(GenePairs.df$Gene.2[GenePairs.df$Gene.1 %in% BiclusterGenes & GenePairs.df$Gene.2 %in% BiclusterGenes])

    Bicluster.Edgelist <- cbind(Col1.Edgelist,Col2.Edgelist)

    n <- network::network(Bicluster.Edgelist,directed = F)

    n <- ggnetwork::ggnetwork(n, layout = "fruchtermanreingold", cell.jitter = 0.75)

    FileName <- paste0(substr(GeneNames,1,nchar(GeneNames)-13),"Bicluster",TempI,".pdf")

    x <- y <- xend <- yend <- vertex.names <- NULL

    if (length(BiclusterGenes) < 50){
      p <- ggplot2::ggplot(n, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggnetwork::geom_edges(color = "grey",curvature = 0.1) +
        ggnetwork::geom_nodes(ggplot2::aes(color = "peachpuff"),size = 15,alpha = 0.5) +
        ggnetwork::geom_nodes(ggplot2::aes(x, y, color = "peachpuff"), size = 10) +
        ggnetwork::geom_nodetext(ggplot2::aes(label = vertex.names),
                                 fontface = "bold") +
        ggnetwork::theme_blank() +
        ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(filename = FileName,p,device = "pdf",width = 10,height = 10,units = "in")
    } else if (length(BiclusterGenes) >= 50 & length(BiclusterGenes) < 100){
      p <- ggplot2::ggplot(n, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggnetwork::geom_edges(color = "grey",curvature = 0.1) +
        ggnetwork::geom_nodes(ggplot2::aes(color = "peachpuff"),size = 15,alpha = 0.5) +
        ggnetwork::geom_nodes(ggplot2::aes(x, y, color = "peachpuff"), size = 10) +
        ggnetwork::geom_nodetext(ggplot2::aes(label = vertex.names),
                                 fontface = "bold") +
        ggnetwork::theme_blank() +
        ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(filename = FileName,p,device = "pdf",width = 12,height = 12,units = "in")
    } else if (length(BiclusterGenes) >= 100 & length(BiclusterGenes) < 150){
      p <- ggplot2::ggplot(n, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggnetwork::geom_edges(color = "grey",curvature = 0.1) +
        ggnetwork::geom_nodes(ggplot2::aes(color = "peachpuff"),size = 15,alpha = 0.5) +
        ggnetwork::geom_nodes(ggplot2::aes(x, y, color = "peachpuff"), size = 10) +
        ggnetwork::geom_nodetext(ggplot2::aes(label = vertex.names),
                                 fontface = "bold") +
        ggnetwork::theme_blank() +
        ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(filename = FileName,p,device = "pdf",width = 14,height = 14,units = "in")
    } else if (length(BiclusterGenes) >= 150 & length(BiclusterGenes) < 200){
      p <- ggplot2::ggplot(n, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggnetwork::geom_edges(color = "grey",curvature = 0.1) +
        ggnetwork::geom_nodes(ggplot2::aes(color = "peachpuff"),size = 15,alpha = 0.5) +
        ggnetwork::geom_nodes(ggplot2::aes(x, y, color = "peachpuff"), size = 10) +
        ggnetwork::geom_nodetext(ggplot2::aes(label = vertex.names),
                                 fontface = "bold") +
        ggnetwork::theme_blank() +
        ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(filename = FileName,p,device = "pdf",width = 16,height = 16,units = "in")
    } else if (length(BiclusterGenes) >= 200 & length(BiclusterGenes) < 250){
      p <- ggplot2::ggplot(n, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggnetwork::geom_edges(color = "grey",curvature = 0.1) +
        ggnetwork::geom_nodes(ggplot2::aes(color = "peachpuff"),size = 15,alpha = 0.5) +
        ggnetwork::geom_nodes(ggplot2::aes(x, y, color = "peachpuff"), size = 10) +
        ggnetwork::geom_nodetext(ggplot2::aes(label = vertex.names),
                                 fontface = "bold") +
        ggnetwork::theme_blank() +
        ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(filename = FileName,p,device = "pdf",width = 18,height = 18,units = "in")
    } else if(length(BiclusterGenes) >= 250){
      p <- ggplot2::ggplot(n, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggnetwork::geom_edges(color = "grey",curvature = 0.1) +
        ggnetwork::geom_nodes(ggplot2::aes(color = "peachpuff"),size = 15,alpha = 0.5) +
        ggnetwork::geom_nodes(ggplot2::aes(x, y, color = "peachpuff"), size = 10) +
        ggnetwork::geom_nodetext(ggplot2::aes(label = vertex.names),
                                 fontface = "bold") +
        ggnetwork::theme_blank() +
        ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(filename = FileName,p,device = "pdf",width = 20,height = 20,units = "in")
    }
  }
}

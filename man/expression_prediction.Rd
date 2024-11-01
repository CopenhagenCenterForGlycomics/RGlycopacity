% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{expression_prediction}
\alias{expression_prediction}
\title{single cell gene expression prediction}
\usage{
expression_prediction(ncells, pseudobulk)
}
\arguments{
\item{ncells}{number of cells in the cluster}

\item{pseudobulk}{mean expression of the gene in the cluster}
}
\value{
a pvalue corresponding to the probability that the gene is expressed in the cluster
}
\description{
predicted expression of a gene in a single cell cluster
}

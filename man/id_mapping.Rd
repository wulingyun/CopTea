% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapping.R
\name{id_mapping}
\alias{id_mapping}
\title{Mapping between different gene ids}
\usage{
id_mapping(ids, db, from = 1, to = 2)
}
\arguments{
\item{ids}{The vector of input gene ids.}

\item{db}{The gene id mapping database.}

\item{from}{The index of input gene id in the mapping database.}

\item{to}{The index of output gene id in the mapping database.}
}
\value{
This function will return a vector of output gene ids mapped from given gene ids.
}
\description{
Map gene ids to other gene ids by using organism annotation package in Bioconductor.
}
\seealso{
\code{\link{id_mapping_from_symbol}}, \code{\link{id_mapping_to_symbol}}, 
\code{\link{id_mapping_species}}
}

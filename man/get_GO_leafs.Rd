% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{get_GO_leafs}
\alias{get_GO_leafs}
\title{Get GO leaf terms}
\usage{
get_GO_leafs(domains = c("BP", "CC", "MF"))
}
\arguments{
\item{domains}{A vector of strings indicated the specific GO domains, possible values include
\code{"BP"}, \code{"CC"}, and \code{"MF"}.}
}
\value{
This function will return a vector of strings, which is the GO IDs of extracted terms.
}
\description{
Extract the leaf terms from the specific GO domains.
}

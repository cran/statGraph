% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.kmeans.R
\name{graph.kmeans}
\alias{graph.kmeans}
\title{K-means for Graphs}
\usage{
graph.kmeans(Graphs, k, nstart = 2, dist = "JS", ...)
}
\arguments{
\item{Graphs}{a list of undirected graphs.
If each graph has the  attribute \code{eigenvalues} containing its
eigenvalues , such values will be used to
compute their spectral density.}

\item{k}{an integer specifying the number of clusters.}

\item{nstart}{the number of trials of k-means clusterizations. The algorithm
returns the clusterization with the best silhouette.}

\item{dist}{string indicating if you want to use the 'JS' (default), 'L1' or 'L2'
distances. 'JS' means Jensen-Shannon divergence.}

\item{...}{Other relevant parameters for \code{\link{graph.spectral.density}}.}
}
\value{
A list with class 'statGraph' containing the following components:
\item{\code{method:}}{ a string indicating the used method.}
\item{\code{info:}}{ a string showing details about the method.}
\item{\code{data.name:}}{ a string with the data's name(s).}
\item{\code{cluster:}}{ a vector of the same length of \code{Graphs} containing the clusterization
labels.}
\item{\code{centers:}}{a list containing the centroids of each cluster.}
}
\description{
\code{graph.kmeans} clusters graphs following a k-means algorithm based on the
Jensen-Shannon divergence between the spectral densities of the graphs.
}
\examples{
set.seed(1)
g <- list()
for(i in 1:5){
  g[[i]] <- igraph::sample_gnp(30, p=0.2)
}
for(i in 6:10){
  g[[i]] <- igraph::sample_gnp(30, p=0.5)
}
res <- graph.kmeans(g, k=2, nstart=2)
res

}
\references{
MacQueen, James. 'Some methods for classification and analysis of
multivariate observations.' Proceedings of the fifth Berkeley symposium on
mathematical statistics and probability. Vol. 1. No. 14. 1967.

Lloyd, Stuart. 'Least squares quantization in PCM.' IEEE transactions on
information theory 28.2 (1982): 129-137.
}
\keyword{k-means}

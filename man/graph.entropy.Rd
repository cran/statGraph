% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.entropy.R
\name{graph.entropy}
\alias{graph.entropy}
\title{Graph Spectral Entropy}
\usage{
graph.entropy(Graph, ...)
}
\arguments{
\item{Graph}{the undirected graph (igraph object).
If \code{Graph} has the  attribute \code{eigenvalues} containing
the eigenvalues of \code{Graph}, such values will be used to
compute its spectral density.}

\item{...}{Other relevant parameters for \code{\link{graph.spectral.density}}.}
}
\value{
A list with class 'statGraph' containing the following components:
\item{\code{method:}}{ a string indicating the used method.}
\item{\code{info:}}{ a string showing details about the method.}
\item{\code{data.name:}}{ a string with the data's name(s).}
\item{\code{entropy:}}{ a real number corresponding to the graph spectral entropy.}
}
\description{
\code{graph.entropy} returns the spectral entropy of an undirected graph.
}
\examples{
set.seed(1)
G <- igraph::sample_gnp(n=100, p=0.5)
entropy <- graph.entropy(Graph = G)
entropy

}
\references{
Takahashi, D. Y., Sato, J. R., Ferreira, C. E. and Fujita A. (2012)
Discriminating Different Classes of Biological Networks by Analyzing the
Graph Spectra  Distribution. _PLoS ONE_, *7*, e49949.
doi:10.1371/journal.pone.0049949.
}
\keyword{spectral_entropy}

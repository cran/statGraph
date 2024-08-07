% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp.anogva.R
\name{sp.anogva}
\alias{sp.anogva}
\title{Semi-parametric Analysis of Graph Variability (SP-ANOGVA)}
\usage{
sp.anogva(Graphs, model, maxBoot = 100, ...)
}
\arguments{
\item{Graphs}{a list of undirected graphs.
If each graph has the  attribute \code{eigenvalues} containing its
eigenvalues , such values will be used to
compute their spectral density.}

\item{model}{A string that indicates one of the following models: 'ER'
(Erdos-Renyi random graph model), 'GRG' (geometric random graph model), 'WS'
(Watts-Strogatz random graph model), and 'BA' (Barabási-Albert random graph
model).}

\item{maxBoot}{integer indicating the number of bootstrap resamples (default
is \code{500}).}

\item{...}{Other relevant parameters for \code{\link{graph.param.estimator}}.}
}
\value{
A list with class 'htest' containing the following components:
\item{\code{statistic:}}{ the F statistic of the test.}
\item{\code{p.value:}}{ the p-value of the test.}
\item{\code{method:}}{ a string indicating the used method.}
\item{\code{data.name:}}{ a string with the data's name(s).}
\item{\code{estimates:}}{ a vector containing the estimated parameters for each graph.}
}
\description{
\code{sp.anogva} statistically tests whether two or more graphs are
generated by the same model and set of parameters.
}
\examples{
\donttest{
set.seed(1)
model <- 'ER'
G <- list()

# Under H0
G[[1]] <- igraph::sample_gnp(50, 0.5)
G[[2]] <- igraph::sample_gnp(50, 0.5)
G[[3]] <- igraph::sample_gnp(50, 0.5)
result1 <- sp.anogva(G, model, maxBoot = 10,eps=0.1)
result1

# Under H1
G[[1]] <- igraph::sample_gnp(50, 0.5)
G[[2]] <- igraph::sample_gnp(50, 0.75)
G[[3]] <- igraph::sample_gnp(50, 0.5)
result2 <- sp.anogva(G, model, maxBoot = 10,eps=0.1)
result2
}


}
\references{
Andre Fujita, Eduardo Silva Lira, Suzana de Siqueira Santos, Silvia Yumi
Bando, Gabriela Eleuterio Soares, Daniel Yasumasa Takahashi. A
semi-parametric statistical test to compare complex networks, Journal of
Complex Networks, cnz028, https://doi.org/10.1093/comnet/cnz028

Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
selection method for kernel density estimation.
_Journal of the Royal Statistical Society series B_, 53, 683-690.
http://www.jstor.org/stable/2345597.
}
\keyword{semi_parametric_analysis_of_graph_variability}

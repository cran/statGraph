% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.spectral.density.R
\name{graph.spectral.density}
\alias{graph.spectral.density}
\title{Graph Spectral Density}
\usage{
graph.spectral.density(Graph, method = "diag", ...)
}
\arguments{
\item{Graph}{the undirected graph (igraph object).
If \code{Graph} has the  attribute \code{eigenvalues} containing
the eigenvalues of \code{Graph}, such values will be used to
compute its spectral density.}

\item{method}{String that specifies the method to obtain the spectral density. It can
take two possible values 'diag' (Default) and 'fast'. If 'diag' is used then
the exact spectral density is obtained, otherwise the degree-based spectral
density is obtained.}

\item{...}{Other relevant parameters to obtain the spectral density such as \code{from}, \code{to},
an \code{npoints}.  \code{from}, \code{to} specify  the lower and upper bound of the eigenvalues' support
(automatically computed if not given); and \code{npoints} is the number of discretization points (default \code{1024}) of
the interval <\code{from},\code{to}>.
There are other parameters that depend on the value of the parameter \code{method}:
If \code{method='diag'}, then the parameter \code{bandwidth} can be used.
This parameter is a string that specifies the criterion to choose the
bandwidth during the spectral density estimation. Choose between the
following criteria: 'Silverman' (default), 'Sturges', 'bcv', 'ucv' and 'SJ'.
'bcv' is an abbreviation of biased cross-validation, while 'ucv' means
unbiased cross-validation. 'SJ'  implements the methods of Sheather & Jones
(1991) to select the bandwidth using pilot estimation of derivatives.
Otherwise, if \code{method='fast'}, then the parameter \code{numCores} can be used. This parameter
specifies the number of cores (default \code{1}) to use for parallelization.}
}
\value{
A list with class 'statGraph' containing the following components:
\item{\code{method:}}{ a string indicating the used method.}
\item{\code{info:}}{ a string showing details about the method.}
\item{\code{data.name:}}{ a string with the data's name(s).}
\item{\code{x:}}{ a vector corresponding to the x axis coordinates of the density function.}
\item{\code{y:}}{ a vector corresponding to the y axis coordinates of the density function.}
\item{\code{from:}}{ a real number corresponding to the smallest value of the x axis.}
\item{\code{to:}}{ a real number corresponding to the largest value of the x axis.}
}
\description{
\code{graph.spectral.density} returns the exact or degree-based spectral density
in the interval <\code{from},\code{to}> by using \code{npoints} discretization points.
}
\examples{
set.seed(1)
G <- igraph::sample_smallworld(dim = 1, size = 50, nei = 2, p = 0.2)

# Obtain the spectral density
density <- graph.spectral.density(Graph = G)
density

}
\references{
#' Takahashi, D. Y., Sato, J. R., Ferreira, C. E. and Fujita A. (2012)
Discriminating Different Classes of Biological Networks by Analyzing the
Graph Spectra  Distribution. _PLoS ONE_, *7*, e49949.
doi:10.1371/journal.pone.0049949.

Silverman, B. W. (1986) _Density Estimation_.  London: Chapman and Hall.

Sturges, H. A. The Choice of a Class Interval. _J. Am. Statist. Assoc._,
*21*, 65-66.

Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
selection method for kernel density estimation.
_Journal of the Royal Statistical Society series B_, 53, 683-690.
http://www.jstor.org/stable/2345597.

Newman, M. E. J., Zhang, X., & Nadakuditi, R. R. (2019).
Spectra of random networks with arbitrary degrees.
Physical Review E, 99(4), 042309.
}
\keyword{eigenvalue_density}

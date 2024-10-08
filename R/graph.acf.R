#' Autocorrelation Function Estimation for Graphs
#'
#' The function \code{graph.acf} computes estimates of the autocorrelation
#' function for graphs.
#'
#' @param Graphs a list of undirected graphs.
#' If each graph has the  attribute \code{eigenvalues} containing its
#' eigenvalues , such values will be used to
#' compute their spectral density.
#'
#' @param plot logical. If \code{TRUE} (default) the graph.acf is plotted.
#'
#' @return An object of class acf.
#'
#' @keywords autocorrelation
#'
#' @references
#' Fujita, A., Takahashi, D. Y., Balardin, J. B., Vidal, M. C. and Sato, J. R.
#' (2017) Correlation between graphs with an application to brain network
#' analysis. _Computational Statistics & Data Analysis_ *109*, 76-92.
#'
#' @examples
#' set.seed(1)
#' G <- list()
#' p <- array(0, 100)
#' p[1:3] <- rnorm(3)
#' for (t in 4:100) {
#'   p[t] <- 0.5*p[t-3] + rnorm(1)
#' }
#' ma <- max(p)
#' mi <- min(p)
#' p <- (p - mi)/(ma-mi)
#' for (t in 1:100) {
#'   G[[t]] <- igraph::sample_gnp(100, p[t])
#' }
#' graph.acf(G, plot=TRUE)
#'
#' @export
graph.acf <- function(Graphs, plot = TRUE) {
    if (!valid.input(Graphs, level = 1)) {
        stop("The input should be a list of igraph objects!")
    }

    G.radius <- unlist(Map(f = function(G) {
        get.largest.eigenvalue(G)
    }, Graphs))
    res <- acf(G.radius, plot = plot)
    return(res)
}

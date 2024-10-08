#' Graph Model Selection
#'
#' \code{graph.model.selection} selects the graph model that best approximates the
#' observed graph according to the Graph Information Criterion (GIC).
#'
#' @param Graph the undirected graph (igraph object).
#' If \code{Graph} has the  attribute \code{eigenvalues} containing
#' the eigenvalues of \code{Graph}, such values will be used to
#' compute its spectral density.
#'
#' @param models either a vector of strings, or a list of functions:
#'
#' A vector of strings containing some of the following models: 'ER' (Erdos-Renyi
#' random graph), 'GRG' (geometric random graph), 'KR' (k regular random graph),
#' 'WS' (Watts-Strogatz model), and 'BA' (Barabási-Albert model).
#'
#' A list of functions. Each function returns a graph (igraph object)
#' generated by a graph model and has two arguments: the graph
#' size and the model parameter, in this order.
#'
#' If the argument \code{models} is \code{NULL}, then the 'ER', 'WS', and 'BA' models will
#' be considered for the model selection.
#'
#' @param parameters list of numeric vectors or list of lists. If a list of numeric vectors is given,
#' then each vector contains the values that will be considered for the parameter estimation of the
#' corresponding model. If a list of lists is given, then each list contains \code{lo} and \code{hi} elements
#' that indicate the model's parameter search interval <\code{lo},\code{hi}>.
#' If the user does not provide the argument \code{parameters}, then default values
#' are used for the predefined models ('ER', 'GRG', 'KR', 'WS', and 'BA') as done in \code{graph.param.estimator}.
#'
#' @param ... Other relevant parameters for \code{\link{graph.param.estimator}}.
#'
#' @return A list with class 'statGraph' containing the following components:
#' \item{\code{method:}}{ a string indicating the used method.}
#' \item{\code{info:}}{ a string showing details about the method.}
#' \item{\code{model:}}{ the indice(s) or name(s) of the selected model(s), i. e. the
#' model(s) that minimize(s) the Graph Information Criterion (GIC).}
#' \item{\code{estimates:}}{ a matrix in which each row corresponds to a model, the
#' column 'param' corresponds to the parameter estimate, and the column 'GIC'
#' corresponds to the Graph Information Criterion (GIC), i. e. the
#' distance measure (Kullback-Leibler divergence by default) between the observed graph and the model.}
#'
#' @keywords model_selection
#'
#' @references
#' Takahashi, D. Y., Sato, J. R., Ferreira, C. E. and Fujita A. (2012)
#' Discriminating Different Classes of Biological Networks by Analyzing the
#' Graph Spectra  Distribution. _PLoS ONE_, *7*, e49949.
#' doi:10.1371/journal.pone.0049949.
#'
#' Silverman, B. W. (1986) _Density Estimation_.  London: Chapman and Hall.
#'
#' Sturges, H. A. The Choice of a Class Interval. _J. Am. Statist. Assoc._,
#' *21*, 65-66.
#'
#' Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
#' selection method for kernel density estimation.
#' _Journal of the Royal Statistical Society series B_, 53, 683-690.
#' http://www.jstor.org/stable/2345597.
#'
#' @examples
#'
#' ## Example using an igraph object as input data
#' set.seed(1)
#' G <- igraph::sample_gnp(n=30, p=0.5)
#'
#' # Using strings to indicate the graph models
#' result1 <- graph.model.selection(G, models=c('ER', 'WS'), eps = 0.5)
#' result1
#'
#' \donttest{
#' ## Using functions to describe the graph models
#' # Erdos-Renyi graph
#' model1 <- function(n, p) {
#'   return(igraph::sample_gnp(n, p))
#' }
#' # Watts-Strogatz small-world graph
#' model2 <- function(n, pr, K=8) {
#'   return(igraph::sample_smallworld(1, n, K, pr))
#' }
#' parameters <- list(seq(0.01, 0.99, 0.49), seq(0.01, 0.99, 0.49))
#' result2 <- graph.model.selection(G, list(model1, model2), parameters)
#' result2
#' }
#'
#' @export
graph.model.selection <- function(Graph, models = NULL, parameters = NULL, ...) {

    if (!valid.input(Graph)) {
        stop("The input should be an igraph object!")
    }

    data.name <- deparse(substitute(Graph))
    if (is.null(models)) {
        models <- list("ER", "WS", "BA")
    }

    Graph$density <- graph.spectral.density(Graph = Graph, ...)

    results <- matrix(NA, length(models), 2)
    # if models is a list of functions return the position/ otherwise return the model name
    model_names <- c()
    for (idx in 1:length(models)) {
        if (methods::is(models[idx], "function"))
            model_names <- c(model_names, as.character(idx)) else model_names <- c(model_names, models[idx])
    }
    rownames(results) <- unlist(model_names)
    colnames(results) <- c("param", "GIC")
    for (i in 1:length(models)) {
        param <- NULL
        if (!is.null(parameters)) {
            param <- parameters[[i]]
        }
        r <- graph.param.estimator(Graph = Graph, model = models[[i]], interval = param, ...)
        results[i, "param"] <- r$param
        results[i, "GIC"] <- r$dist
    }
    m <- which(results[, "GIC"] == min(results[, "GIC"]))
    if (!is.null(rownames(results))) {
        m <- rownames(results)[m]
    }
    ###
    method_info <- "Graph Model Selection"
    info <- "Selects the graph model that best approximates the observed graph"
    output <- list(method = method_info, info = info, data.name = data.name, model = m, estimates = results)
    class(output) <- "statGraph"
    return(output)
}

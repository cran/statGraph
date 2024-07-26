#' @exportS3Method plot statGraph
plot.statGraph <- function(x, ...) {
  # only works when x is the output of 'spectral.density'
  if (any(names(x) == "x") && any(names(x) == "y")) {
    plot(x$x,x$y,xlab = expression("Eigenvalues ("~lambda~")" ),ylab = expression(rho~"("~lambda~")"),type = "l",col = "red",main = x$method,lwd = 2.0)
  }
}

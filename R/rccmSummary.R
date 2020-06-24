#' BIC for RCM
#'
#' This function calculates the BIC for the random covariance model (RCM)
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}.
#' @param Omegas \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level precision matrices.
#' @param Gk_est \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level networks.
#' @return Numeric BIC value
#'
#' @author
#' Lin Zhang
#'
#' @examples
#' # Generate data
#' set.seed(1994)
#' myData <- rccSim(G = 1, clustSize = 10, p = 10, n = 100, overlap = 0.50, rho = 0.10)
#'
#' # Analyze with RCM
#' rcmRes <- randCov(myData$simDat, lambda1 = 0.01, lambda2 = 0.01, lambda3 = 0, delta = 0.0001)
#'
#' # Calculate BIC for the RCM
#' bic_cal(x = myData$simDat, Omegas = rcmRes$Omegas)
#'
#' @export
bic_cal <- function(x, Omegas, Gk_est = NULL) {
  nk <- sapply(x, FUN = nrow)
  S <- lapply(x, FUN = cov)

  p <- dim(Omegas)[1]
  K <- length(x)

  if(is.null(Gk_est)) {
    Gk_est <- (round(Omegas,3) != 0) - array(diag(p), c(p, p, K))
  }

  nedges <- apply(Gk_est, MARGIN = 3, FUN = sum) / 2

  bic <- mapply(FUN = function(x1, x2, x3, x4) {(x1 - 1) * sum(diag(x2 %*% x3)) - x1 * log(det(x3)) + x4 * log(x1)},
                nk, S, lapply(1:K, FUN = function(k){Omegas[, , k]}), nedges + p)

  return(sum(bic))
}

#' Modified BIC for RCM
#'
#' This function calculates the modified BIC for the random covariance model (RCM)
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}.
#' @param Omega0 \eqn{p} x \eqn{p} group-level precision matrix estimate.
#' @param Omegas \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level precision matrices.
#' @param lambda2 Non-negative scalar. Induces similarity between subject-level matrices and group-level matrix.
#' @param G0_est \eqn{p} x \eqn{p} group-level network estimate.
#' @param Gk_est \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level networks.
#' @return Numeric Modified BIC value
#'
#' @author
#' Lin Zhang
#'
#' @examples
#' # Generate data
#' set.seed(1994)
#' myData <- rccSim(G = 1, clustSize = 10, p = 10, n = 100, overlap = 0.50, rho = 0.10)
#'
#' # Analyze with RCM
#' rcmRes <- randCov(myData$simDat, lambda1 = 0.01, lambda2 = 0.01, lambda3 = 0, delta = 0.0001)
#'
#' # Calculate modified BIC for the RCM
#' mbic_cal(x = myData$simDat, Omega0 = rcmRes$Omega0, Omegas = rcmRes$Omegas,
#'          lambda2 = 0.01)
#'
#' @export
mbic_cal <- function(x, Omega0, Omegas, lambda2, G0_est = NULL, Gk_est = NULL) {
  nk <- sapply(x, FUN = nrow)
  S <- lapply(x, FUN = cov)

  p <- dim(Omegas)[1]
  K <- length(x)

  if(is.null(Gk_est)) {Gk_est = (round(Omegas,3) != 0) - array(diag(p),c(p,p,K))}
  if(is.null(G0_est)) { G0_est = (round(Omega0,3) != 0) - diag(p)}
  nedges <- apply(Gk_est, MARGIN = 3, FUN = sum) / 2

  df.r <- (nedges + p) / (1 + lambda2)
  df.f <- (sum(G0_est) / 2 + p) * lambda2 / (1 + lambda2)

  mbic <- mapply(FUN = function(x1,x2,x3) {(x1-1)*sum(diag(x2%*%x3)) - x1*log(det(x3))},
                 nk, S, lapply(1:K, FUN = function(k){Omegas[, , k]}))

  return(sum(mbic) + (sum(df.r) + df.f) * log(sum(nk)))
}

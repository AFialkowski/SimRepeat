#' @title Calculate Expected Correlation Matrix of Outcomes (Y) for Correlated Systems of Continuous Variables
#'
#' @description This function calculates the expected correlation matrix for outcomes (Y) in a correlated system of continuous variables.
#'     This system is generated with \code{\link[SimRepeat]{nonnormsys}} using the techniques of Headrick and Beasley (\doi{10.1081/SAC-120028431}).
#'     These correlations are determined based on the beta (slope) coefficients calculated with \code{\link[SimRepeat]{calc_betas}}, the correlations
#'     between independent variables \eqn{X_{(pj)}} for a given outcome \eqn{Y_p}, for \code{p = 1, ..., M}, the
#'     correlations between error terms, and the variances.  The result can be used to compare the simulated correlation matrix to the theoretical correlation matrix.
#'     If there are continuous mixture variables and the betas are specified in terms of non-mixture and mixture variables
#'     and/or \code{error_type = "mix"}, then the correlations in \code{corr.x} and/or \code{corr.e} will be calculated in terms of non-mixture and mixture variables using
#'     \code{\link[SimCorrMix]{rho_M1M2}} and \code{\link[SimCorrMix]{rho_M1Y}}.  In this case, the dimensions of the matrices in \code{corr.x} should not
#'     match the number of columns of \code{betas}.  The vignette \strong{Theory and Equations for Correlated Systems of Continuous Variables}
#'     gives the equations, and the vignette \strong{Correlated Systems of Statistical Equations with Non-Mixture and Mixture Continuous
#'     Variables} gives examples.  There are also vignettes in \code{\link[SimCorrMix]{SimCorrMix}} which provide more details on continuous
#'     non-mixture and mixture variables.
#'
#' @param betas a matrix of the slope coefficients calculated with \code{\link[SimRepeat]{calc_betas}}, rows represent the outcomes
#' @param corr.x list of length \code{M}, each component a list of length \code{M}; \code{corr.x[[p]][[q]]} is matrix of correlations
#'     for independent variables in equations p (\eqn{X_{(pj)}} for outcome \eqn{Y_p}) and q (\eqn{X_{(qj)}} for outcome \eqn{Y_q});
#'     if p = q, \code{corr.x[[p]][[q]]} is a correlation matrix with \code{nrow(corr.x[[p]][[q]])} = # \eqn{X_{(pj)}} for outcome \eqn{Y_p};
#'     if p != q, \code{corr.x[[p]][[q]]} is a non-symmetric matrix of correlations where rows correspond to covariates for \eqn{Y_p}
#'     so that \code{nrow(corr.x[[p]][[q]])} = # \eqn{X_{(pj)}} for outcome \eqn{Y_p} and
#'     columns correspond to covariates for \eqn{Y_q} so that \code{ncol(corr.x[[p]][[q]])} = # \eqn{X_{(qj)}} for outcome \eqn{Y_q};
#'     order is 1st continuous non-mixture and 2nd components of continuous mixture variables
#' @param corr.e correlation matrix for continuous non-mixture or components of mixture error terms
#' @param vars a list of same length as \code{corr.x} of vectors of variances for \eqn{X_{(pj)}, E}; E term should be last;
#'     order should be the same as in \code{corr.x}
#' @param mix_pis a list of same length as \code{corr.x}, where \code{mix_pis[[p]][[j]]} is a vector of mixing probabilities for \eqn{X_{mix(pj)}} that sum to 1,
#'     the j-th mixture covariate for outcome \eqn{Y_p}; the last element of \code{mix_pis[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_pis[[p]] = NULL}
#' @param mix_mus a list of same length as \code{corr.x}, where \code{mix_mus[[p]][[j]]} is a vector of means for \eqn{X_{mix(pj)}},
#'     the j-th mixture covariate for outcome \eqn{Y_p}; the last element of \code{mix_mus[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_mus[[p]] = NULL}
#' @param mix_sigmas a list of same length as \code{corr.x}, where \code{mix_sigmas[[p]][[j]]} is a vector of standard deviations for \eqn{X_{mix(pj)}},
#'     the j-th mixture covariate for outcome \eqn{Y_p}; the last element of \code{mix_sigmas[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_sigmas[[p]] = NULL}
#' @param error_type "non_mix" if all error terms have continuous non-mixture distributions, "mix" if all error terms have continuous mixture distributions,
#'     defaults to "non_mix"
#'
#' @return \code{corr.y} the correlation matrix for the outcomes \eqn{Y}
#'
#' @import SimCorrMix
#' @export
#' @keywords continuous mixture Headrick Beasley
#' @seealso \code{\link[SimRepeat]{nonnormsys}}, \code{\link[SimRepeat]{calc_betas}}, \code{\link[SimCorrMix]{rho_M1M2}}, \code{\link[SimCorrMix]{rho_M1Y}}
#' @references
#' Headrick TC, Beasley TM (2004).  A Method for Simulating Correlated Non-Normal Systems of Linear Statistical Equations.
#'     Communications in Statistics - Simulation and Computation, 33(1).  \doi{10.1081/SAC-120028431}
#'
#' @examples \dontrun{
#' # Example: system of three equations for 2 independent variables, where each
#' # error term has unit variance, from Headrick & Beasley (2002)
#' corr.yx <- list(matrix(c(0.4, 0.4), 1), matrix(c(0.5, 0.5), 1),
#'   matrix(c(0.6, 0.6), 1))
#' corr.x <- list()
#' corr.x[[1]] <- corr.x[[2]] <- corr.x[[3]] <- list()
#' corr.x[[1]][[1]] <- matrix(c(1, 0.1, 0.1, 1), 2, 2)
#' corr.x[[1]][[2]] <- matrix(c(0.1974318, 0.1859656, 0.1879483, 0.1858601),
#'   2, 2, byrow = TRUE)
#' corr.x[[1]][[3]] <- matrix(c(0.2873190, 0.2589830, 0.2682057, 0.2589542),
#'   2, 2, byrow = TRUE)
#' corr.x[[2]][[1]] <- t(corr.x[[1]][[2]])
#' corr.x[[2]][[2]] <- matrix(c(1, 0.35, 0.35, 1), 2, 2)
#' corr.x[[2]][[3]] <- matrix(c(0.5723303, 0.4883054, 0.5004441, 0.4841808),
#'   2, 2, byrow = TRUE)
#' corr.x[[3]][[1]] <- t(corr.x[[1]][[3]])
#' corr.x[[3]][[2]] <- t(corr.x[[2]][[3]])
#' corr.x[[3]][[3]] <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
#' corr.e <- matrix(0.4, nrow = 3, ncol = 3)
#' diag(corr.e) <- 1
#' vars <- list(rep(1, 3), rep(1, 3), rep(1, 3))
#' betas <- calc_betas(corr.yx, corr.x, vars)
#' calc_corr_y(betas, corr.x, corr.e, vars)
#' }
calc_corr_y <- function(betas = NULL, corr.x = list(), corr.e = NULL,
                        vars = list(), mix_pis = list(), mix_mus = list(),
                        mix_sigmas = list(),
                        error_type = c("non_mix", "mix")) {
  M <- length(corr.x)
  K.x <- numeric(M)
  for (p in 1:M) {
    if (!is.null(corr.x[[p]])) K.x[p] <- ncol(corr.x[[p]][[p]])
  }
  if (!isTRUE(all.equal(K.x, apply(betas, 1,
    function(x) length(x) - sum(round(x, 10) == 0)),
    check.attributes = FALSE)) & length(mix_pis) == 0)
    stop("The dimensions of betas should match the dimensions of corr.x if
         there are no mixture variables.")
  if (length(error_type) == 2) error_type <- "non_mix"
  K.mix <- rep(0, M)
  if (length(mix_pis) > 0) K.mix <- lengths(mix_pis)
  K.cont <- lengths(vars) - K.mix
  if (error_type == "mix") K.mix2 <- K.mix - 1 else K.mix2 <- K.mix
  vars0 <- vars
  if (isTRUE(all.equal(K.x, apply(betas, 1,
    function(x) length(x) - sum(round(x, 10) == 0)),
    check.attributes = FALSE)) & length(mix_pis) > 0) {
    vars <- list()
    for (i in 1:M) {
      vars <- append(vars, list(NULL))
      if (error_type == "non_mix") {
        if ((K.cont[i] - 1) > 0) {
          vars[[i]] <- append(vars[[i]], vars0[[i]][1:(K.cont[i] - 1)])
        }
        if (K.mix2[i] > 0) {
          vars[[i]] <- append(vars[[i]], unlist(mix_sigmas[[i]])^2)
        }
      } else {
        if (K.cont[i] > 0) {
          vars[[i]] <- append(vars[[i]], vars0[[i]][1:K.cont[i]])
        }
        if (K.mix2[i] > 0) {
          vars[[i]] <- append(vars[[i]],
            unlist(mix_sigmas[[i]][-length(mix_sigmas[[i]])])^2)
        }
      }
      vars[[i]] <- append(vars[[i]], vars0[[i]][length(vars0[[i]])])
    }
  }
  if (error_type == "mix") {
    epis <- mapply('[', mix_pis, lengths(mix_pis))
    emus <- mapply('[', mix_mus, lengths(mix_mus))
    esigmas <- mapply('[', mix_sigmas, lengths(mix_sigmas))
    k.error <- c(0, cumsum(lengths(epis)))
    corr.e0 <- corr.e
    corr.e <- diag(1, M)
    for (i in 1:(M - 1)) {
      for (j in (i + 1):M) {
        corr.e[i, j] <- rho_M1M2(list(epis[[i]], epis[[j]]),
          list(emus[[i]], emus[[j]]), list(esigmas[[i]], esigmas[[j]]),
          corr.e0[(k.error[i] + 1):(k.error[i + 1]),
          (k.error[j] + 1):(k.error[j + 1])])
        corr.e[j, i] <- corr.e[i, j]
      }
    }
    mix_pis <- lapply(mix_pis, function(x) if (length(x) %in% c(0, 1))
      list(NULL) else x[-length(x)])
    mix_mus <- lapply(mix_mus, function(x) if (length(x) %in% c(0, 1))
      list(NULL) else x[-length(x)])
    mix_sigmas <- lapply(mix_sigmas, function(x) if (length(x) %in% c(0, 1))
      list(NULL) else x[-length(x)])
  }
  if (!isTRUE(all.equal(K.x, apply(betas, 1,
    function(x) length(x) - sum(round(x, 10) == 0)),
    check.attributes = FALSE)) & length(mix_pis) > 0) {
    corr.x0 <- corr.x
    K.cont <- numeric(M)
    K.cmix <- lapply(mix_pis, function(x) c(0, cumsum(sapply(x, length))))
    for (p in 1:M) {
      if (is.null(corr.x[[p]])) next
      K.cont[p] <- K.x[p] - length(unlist(mix_pis[[p]]))
      for (q in 1:M) {
        if (is.null(corr.x[[p]][[q]])) next
        K.cont[q] <- K.x[q] - length(unlist(mix_pis[[q]]))
        if (q >= p) {
          corr.x[[p]][[q]] <- matrix(1, length(vars0[[p]]) - 1,
                                     length(vars0[[q]]) - 1)
          for (i in 1:nrow(corr.x[[p]][[q]])) {
            for (j in 1:ncol(corr.x[[p]][[q]])) {
              if (i <= K.cont[p] & j <= K.cont[q])
                corr.x[[p]][[q]][i, j] <- corr.x0[[p]][[q]][i, j]
              if (i <= K.cont[p] & j > K.cont[q])
                corr.x[[p]][[q]][i, j] <-
                  rho_M1Y(mix_pis[[q]][[j - K.cont[q]]],
                    mix_mus[[q]][[j - K.cont[q]]],
                    mix_sigmas[[q]][[j - K.cont[q]]],
                    corr.x0[[p]][[q]][i, (K.cont[q] + K.cmix[[q]][j -
                      K.cont[q]] + 1):(K.cont[q] + K.cmix[[q]][j -
                      K.cont[q] + 1])])
              if (i > K.cont[p] & j <= K.cont[q])
                corr.x[[p]][[q]][i, j] <-
                  rho_M1Y(mix_pis[[p]][[i - K.cont[p]]],
                    mix_mus[[p]][[i - K.cont[p]]],
                    mix_sigmas[[p]][[i - K.cont[p]]],
                    corr.x0[[p]][[q]][(K.cont[p] + K.cmix[[p]][i -
                      K.cont[p]] + 1):(K.cont[p] + K.cmix[[p]][i -
                      K.cont[p] + 1]), j])
              if (i > K.cont[p] & j > K.cont[q])
                corr.x[[p]][[q]][i, j] <-
                  rho_M1M2(list(mix_pis[[p]][[i - K.cont[p]]],
                    mix_pis[[q]][[j - K.cont[q]]]),
                    list(mix_mus[[p]][[i - K.cont[p]]],
                         mix_mus[[q]][[j - K.cont[q]]]),
                    list(mix_sigmas[[p]][[i - K.cont[p]]],
                         mix_sigmas[[q]][[j - K.cont[q]]]),
                    corr.x0[[p]][[q]][(K.cont[p] + K.cmix[[p]][i -
                      K.cont[p]] + 1):(K.cont[p] + K.cmix[[p]][i -
                      K.cont[p] + 1]),
                      (K.cont[q] + K.cmix[[q]][j -
                      K.cont[q]] + 1):(K.cont[q] + K.cmix[[q]][j -
                      K.cont[q] + 1])])
            }
          }
        } else {
          corr.x[[p]][[q]] <- t(corr.x[[q]][[p]])
        }
      }
    }
  }
  corr.y <- matrix(0, nrow = M, ncol = M)
  for (p in 1:M) {
    if (is.null(corr.x[[p]])) {
      for (q in 1:M) {
        if (is.null(corr.x[[q]])) {
          corr.y[p, q] <- corr.e[p, q]
          next
        }
        K2 <- ncol(corr.x[[q]][[q]])
        sum2 <- 0
        betas_q <- betas[q, which(round(betas[q, ], 10) != 0), drop = FALSE]
        if (K2 > 1) {
          for (i in 1:(K2 - 1)) {
            for (j in (i + 1):K2) {
              sum2 <- sum2 + betas_q[1, i] * betas_q[1, j] *
                sqrt(vars[[q]][i]) * sqrt(vars[[q]][j]) *
                corr.x[[q]][[q]][i, j]
            }
          }
        }
        denom <- sqrt(vars[[q]][length(vars[[q]])] +
          sum(betas_q[1, ]^2 * vars[[q]][-length(vars[[q]])]) + 2 * sum2)
        corr.y[p, q] <- sqrt(vars[[q]][length(vars[[q]])]) *
          corr.e[p, q]/denom
      }
    } else {
      betas_p <- betas[p, which(round(betas[p, ], 10) != 0), drop = FALSE]
      for (q in 1:M) {
        K1 <- ncol(corr.x[[p]][[p]])
        if (is.null(corr.x[[q]])) {
          sum1 <- 0
          if (K1 > 1) {
            for (i in 1:(K1 - 1)) {
              for (j in (i + 1):K1) {
                sum1 <- sum1 + betas_p[1, i] * betas_p[1, j] *
                  sqrt(vars[[p]][i]) * sqrt(vars[[p]][j]) *
                  corr.x[[p]][[p]][i, j]
              }
            }
          }
          denom <- sqrt(vars[[p]][length(vars[[p]])] +
            sum(betas_p[1, ]^2 * vars[[p]][-length(vars[[p]])]) + 2 * sum1)
          corr.y[p, q] <- sqrt(vars[[p]][length(vars[[p]])]) *
            corr.e[p, q]/denom
          next
        }
        if (q > p) {
          betas_q <- betas[q, which(round(betas[q, ], 10) != 0), drop = FALSE]
          K2 <- ncol(corr.x[[q]][[q]])
          sum1 <- 0
          sum2 <- 0
          if (K1 > 1) {
            for (i in 1:(K1 - 1)) {
              for (j in (i + 1):K1) {
                sum1 <- sum1 + betas_p[1, i] * betas_p[1, j] *
                  sqrt(vars[[p]][i]) * sqrt(vars[[p]][j]) *
                  corr.x[[p]][[p]][i, j]
              }
            }
          }
          if (K2 > 1) {
            for (i in 1:(K2 - 1)) {
              for (j in (i + 1):K2) {
                sum2 <- sum2 + betas_q[1, i] * betas_q[1, j] *
                  sqrt(vars[[q]][i]) * sqrt(vars[[q]][j]) *
                  corr.x[[q]][[q]][i, j]
              }
            }
          }
          denom <- sqrt(vars[[p]][length(vars[[p]])] +
            sum(betas_p[1, ]^2 * vars[[p]][-length(vars[[p]])]) +
            2 * sum1) * sqrt(vars[[q]][length(vars[[q]])] +
            sum(betas_q[1, ]^2 * vars[[q]][-length(vars[[q]])]) +
            2 * sum2)
          corr.y[p, q] <- (sqrt(vars[[p]][length(vars[[p]])] *
            vars[[q]][length(vars[[q]])]) * corr.e[p, q] +
              (betas_p[1, ] * sqrt(vars[[p]][-length(vars[[p]])])) %*%
              (corr.x[[p]][[q]] %*% matrix((betas_q[1, ] *
                sqrt(vars[[q]][-length(vars[[q]])])), ncol = 1)))/denom
          corr.y[q, p] <- corr.y[p, q]
        }
      }
    }
  }
  diag(corr.y) <- 1
  rownames(corr.y) <- paste("Y", 1:M, sep = "")
  colnames(corr.y) <- paste("Y", 1:M, sep = "")
  return(corr.y)
}

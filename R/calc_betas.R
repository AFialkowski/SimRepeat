#' @title Calculate Beta Coefficients for Correlated Systems of Continuous Variables
#'
#' @description This function calculates the beta (slope) coefficients used in \code{\link[SimRepeat]{nonnormsys}}
#'     by the techniques of Headrick and Beasley (\doi{10.1081/SAC-120028431}).  These coefficients are determined based on the
#'     correlations between independent variables \eqn{X_{(pj)}} for a given outcome \eqn{Y_p}, for \code{p = 1, ..., M}, the
#'     correlations between that outcome \eqn{Y_p} and the \eqn{X_{(pj)}} terms, and the variances.  If there are continuous mixture variables and the
#'     matrices in \code{corr.yx} are specified in terms of correlations between outcomes and non-mixture and mixture variables, then the
#'     solutions are the slope coefficients for the non-mixture and mixture variables.  In this case, the number of columns of the matrices of \code{corr.yx}
#'     should not match the dimensions of the matrices in \code{corr.x}.  The correlations in \code{corr.x} will be calculated in terms of non-mixture and
#'     mixture variables using \code{\link[SimCorrMix]{rho_M1M2}} and \code{\link[SimCorrMix]{rho_M1Y}}.  If there are continuous mixture variables and the
#'     matrices in \code{corr.yx} are specified in terms of correlations between outcomes and non-mixture and components of mixture variables,
#'     then the solutions are the slope coefficients for the non-mixture and components of mixture variables.  In this case, the number of columns of the matrices of \code{corr.yx}
#'     should match the dimensions of the matrices in \code{corr.x}.  The vignette \strong{Theory and Equations for Correlated Systems of
#'     Continuous Variables} gives the equations, and the vignette \strong{Correlated Systems of Statistical Equations with Non-Mixture and
#'     Mixture Continuous Variables} gives examples.  There are also vignettes in \code{\link[SimCorrMix]{SimCorrMix}} which provide more details on continuous
#'     non-mixture and mixture variables.
#'
#' @param corr.yx a list of length \code{M} = # of equations, where the p-th component is a 1 row matrix of correlations between \eqn{Y_p} and \eqn{X_{(pj)}};
#'     if there are mixture variables and the \code{betas} are desired in terms of these (and not the components), then \code{corr.yx}
#'     should be specified in terms of correlations between outcomes and non-mixture or mixture variables, and the number of columns of the matrices
#'     of \code{corr.yx} should not match the dimensions of the matrices in \code{corr.x}; if the \code{betas} are desired in terms of
#'     the components, then \code{corr.yx} should be specified in terms of correlations between outcomes and non-mixture or components of
#'     mixture variables, and the number of columns of the matrices of \code{corr.yx} should match the dimensions of the matrices in \code{corr.x}
#' @param corr.x list of length \code{M}, each component a list of length \code{M}; \code{corr.x[[p]][[q]]} is matrix of correlations
#'     for independent variables in equations p (\eqn{X_{(pj)}} for outcome \eqn{Y_p}) and q (\eqn{X_{(qj)}} for outcome \eqn{Y_q});
#'     if p = q, \code{corr.x[[p]][[q]]} is a correlation matrix with \code{nrow(corr.x[[p]][[q]])} = # \eqn{X_{(pj)}} for outcome \eqn{Y_p};
#'     if p != q, \code{corr.x[[p]][[q]]} is a non-symmetric matrix of correlations where rows correspond to covariates for \eqn{Y_p}
#'     so that \code{nrow(corr.x[[p]][[q]])} = # \eqn{X_{(pj)}} for outcome \eqn{Y_p} and
#'     columns correspond to covariates for \eqn{Y_q} so that \code{ncol(corr.x[[p]][[q]])} = # \eqn{X_{(qj)}} for outcome \eqn{Y_q};
#'     order is 1st continuous non-mixture and 2nd components of continuous mixture variables
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
#' @param n the number of sets of random uniform(0, 1) numbers used as starting values in \code{\link[nleqslv]{nleqslv}} to find the betas
#' @param seed the seed for random number generation
#'
#' @return \code{betas} a matrix of slope coefficients where rows represent the outcomes; extra zeros are appended at the end of a row
#'     if that outcome has fewer \eqn{X_{(pj)}} terms
#'
#' @import nleqslv
#' @import SimCorrMix
#' @export
#' @keywords continuous mixture Headrick Beasley
#' @seealso \code{\link[SimRepeat]{nonnormsys}}, \code{\link[SimCorrMix]{rho_M1M2}}, \code{\link[SimCorrMix]{rho_M1Y}}
#' @references
#' Headrick TC, Beasley TM (2004).  A Method for Simulating Correlated Non-Normal Systems of Linear Statistical Equations.
#'     Communications in Statistics - Simulation and Computation, 33(1).  \doi{10.1081/SAC-120028431}
#'
#' @examples
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
#' vars <- list(rep(1, 3), rep(1, 3), rep(1, 3))
#' calc_betas(corr.yx, corr.x, vars)
#'
calc_betas <- function(corr.yx = list(), corr.x = list(), vars = list(),
                       mix_pis = list(), mix_mus = list(), mix_sigmas = list(),
                       error_type = c("non_mix", "mix"), n = 25, seed = 1234) {
  M <- length(corr.x)
  K.x <- numeric(M)
  for (p in 1:M) {
    if (!is.null(corr.x[[p]])) K.x[p] <- ncol(corr.x[[p]][[p]])
  }
  if (!isTRUE(all.equal(K.x, sapply(corr.yx, function(x) if (is.null(x)) 0
    else ncol(x)), check.attributes = FALSE)) & length(mix_pis) == 0)
    stop("The dimensions of corr.yx should match the dimensions of corr.x if
         there are no mixture variables.")
  if (length(error_type) == 2) error_type <- "non_mix"
  K.mix <- rep(0, M)
  if (length(mix_pis) > 0) K.mix <- lengths(mix_pis)
  K.cont <- lengths(vars) - K.mix
  if (error_type == "mix") K.mix2 <- K.mix - 1 else K.mix2 <- K.mix
  vars0 <- vars
  if (isTRUE(all.equal(K.x, sapply(corr.yx, function(x) if (is.null(x)) 0
    else ncol(x)), check.attributes = FALSE)) & length(mix_pis) > 0) {
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
    mix_pis <- lapply(mix_pis, function(x) if (length(x) %in% c(0, 1))
      list(NULL) else x[-length(x)])
    mix_mus <- lapply(mix_mus, function(x) if (length(x) %in% c(0, 1))
      list(NULL) else x[-length(x)])
    mix_sigmas <- lapply(mix_sigmas, function(x) if (length(x) %in% c(0, 1))
      list(NULL) else x[-length(x)])
  }
  if (!isTRUE(all.equal(K.x, sapply(corr.yx, function(x) if (is.null(x)) 0
    else ncol(x)), check.attributes = FALSE)) & length(mix_pis) > 0) {
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
          corr.x[[p]][[q]] <- matrix(1, ncol(corr.yx[[p]]),
                                     ncol(corr.yx[[q]]))
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
  K.max <- max(sapply(Filter(Negate(is.null), corr.yx), ncol))
  betas <- NULL
  for (p in 1:M) {
    if (is.null(corr.x[[p]])) {
      betas <- rbind(betas, rep(0, K.max))
      colnames(betas) <- 1:K.max
      next
    }
    K <- ncol(corr.x[[p]][[p]])
    f.corr.yx <- list()
    for (i in 1:K) {
      f.corr.yx[[i]] <- "("
      for (j in 1:K) {
        if (j == i) {
          f.corr.yx[[i]] <- paste(f.corr.yx[[i]], "b[", j, "] * ",
            sqrt(vars[[p]][j]), " + ", sep = "")
        } else {
          f.corr.yx[[i]] <- paste(f.corr.yx[[i]], "b[", j, "] * ",
            sqrt(vars[[p]][j]) * corr.x[[p]][[p]][i, j], " + ", sep = "")
        }
      }
      f.corr.yx[[i]] <- paste(f.corr.yx[[i]], "0)/(sqrt(",
        vars[[p]][length(vars[[p]])], " + ", sep = "")
      for (j in 1:K) {
        f.corr.yx[[i]] <- paste(f.corr.yx[[i]], "b[", j, "]^2 * ",
                                vars[[p]][j], " + ", sep = "")
      }
      f.corr.yx[[i]] <- paste(f.corr.yx[[i]], "2 * (", sep = "")
      if (K > 1) {
        for (j in 1:(K - 1)) {
          for (k in (j + 1):K) {
            f.corr.yx[[i]] <- paste(f.corr.yx[[i]], "b[", j, "] * b[", k,
              "] * ", sqrt(vars[[p]][j]) * sqrt(vars[[p]][k]) *
                corr.x[[p]][[p]][j, k], " + ", sep = "")
          }
        }
      }
      f.corr.yx[[i]] <- paste(f.corr.yx[[i]], "0))) - ", corr.yx[[p]][1, i],
                              sep = "")
    }
    f.corr.yx2 <- function(b) {
      f <- numeric(K)
      for (i in 1:K) {
        f[i] <- eval(parse(text = f.corr.yx[[i]]))
      }
      return(f)
    }
    set.seed(seed)
    converged <- NULL
    cstart <- matrix(runif(n * K), nrow = n, ncol = K)
    for (i in 1:nrow(cstart)) {
      nl_solution <- nleqslv(x = cstart[i, ], fn = f.corr.yx2,
                             method = "Broyden",
                             control = list(ftol = 1e-05))
      if (nl_solution$termcd == 1) {
        if (K.max - K > 0) {
          converged <- as.data.frame(rbind(converged, c(nl_solution$x,
            rep(0, K.max - K), sum(nl_solution$fvec^2), p)))
        } else {
          converged <- as.data.frame(rbind(converged, c(nl_solution$x,
            sum(nl_solution$fvec^2), p)))
        }
      }
    }
    if (!is.null(converged)) {
      colnames(converged)[(ncol(converged) - 1):ncol(converged)] <-
        c("fnorm", "Y")
      converged <- subset(converged, fnorm == min(fnorm))
      colnames(converged) <- 1:K.max
      betas <- rbind(betas, converged[, 1:K.max])
    } else {
      warning(paste("Solutions could not be found for outcome ", p, sep = ""))
    }
  }
  if (is.null(betas)) stop("Solutions could not be found.")
  rownames(betas) <- paste("Y", 1:M, sep = "")
  colnames(betas) <- paste("B", 1:ncol(betas), sep = "")
  return(as.matrix(betas))
}

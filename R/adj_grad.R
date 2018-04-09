#' @title Convert Non-Positive-Definite Correlation Matrix to Positive-Definite Matrix Using the Adjusted Gradient Updating Method
#'
#' @description This function converts a non-positive-definite correlation matrix to a positive-definite matrix using the
#'     adjusted gradient updating method with initial matrix \code{B1}.
#'
#' @param Sigma the non-PD correlation matrix
#' @param B1 the initial matrix for algorithm; if NULL, uses a scaled initial matrix with diagonal elements \code{sqrt(nrow(Sigma))/2}
#' @param tau parameter used to calculate theta
#' @param tol maximum error for Frobenius norm distance between new matrix and original matrix
#' @param steps maximum number of steps for k (default = 100)
#' @param msteps maximum number of steps for m (default = 10)
#'
#' @return list with \code{Sigma2} the new correlation matrix, \code{dist} the Frobenius norm distance between \code{Sigma2} and \code{Sigma},
#'     \code{eig0} original eigenvalues of \code{Sigma}, \code{eig2} eigenvalues of \code{Sigma2}
#'
#' @importFrom Matrix norm
#' @export
#' @references
#' S Maree (2012). Correcting Non Positive Definite Correlation Matrices. BSc Thesis Applied Mathematics, TU Delft.
#'   \url{http://resolver.tudelft.nl/uuid:2175c274-ab03-4fd5-85a9-228fe421cdbf}.
#'
#' JF Yin and Y Zhang (2013). Alternative gradient algorithms for computing the nearest correlation matrix.
#'     Applied Mathematics and Computation, 219(14): 7591-7599. \url{https://doi.org/10.1016/j.amc.2013.01.045}.
#'
#' Y Zhang and JF Yin. Modified alternative gradients algorithm for computing the nearest correlation matrix.
#'     Internal paper of the Tongji University, Shanghai.
#'
#' @examples
#' Sigma <- matrix(c(1, 0, 0.8, 0, 1, 0.8, 0.8, 0.8, 1), 3, 3, byrow = TRUE)
#' adj_grad(Sigma)
#'
adj_grad <- function(Sigma = NULL, B1 = NULL, tau = 0.5, tol = 0.1,
                     steps = 100, msteps = 10) {
  eig0 <- eigen(Sigma, symmetric = TRUE)$values
  if (min(eig0) >= 0) stop("Sigma is positive-definite.")
  if (is.null(B1)) {
    B <- Sigma
    diag(B) <- sqrt(nrow(B))/2
    B <- diag(1/sqrt(diag(B %*% t(B))), nrow(B)) %*% B
  } else B <- B1
  BB <- B %*% t(B)
  r <- numeric(steps)
  r[1] <- norm((BB - Sigma), type = "F")
  mstart <- 1
  for (k in 2:steps) {
    gradB <- 2 * (BB - Sigma) %*% B
    for (m in mstart:msteps) {
      theta <- 0.5 * (tau^m) * (norm(gradB, type = "F")/norm(gradB %*% t(B),
                                                             type = "F"))^2
      B2 <- B - theta * gradB
      B2 <- diag(1/sqrt(diag(B2 %*% t(B2))), nrow(B2)) %*% B2
      BB <- B2 %*% t(B2)
      r[k] <- norm((BB - Sigma), type = "F")
      if ((r[k] - r[k - 1]) < 0) {
        mstart <- m
        B <- B2
        break
      }
    }
    if (abs(r[k] - r[k - 1])/r[k] <= tol) break
  }
  Sigma2 <- B %*% t(B)
  eig2 <- eigen(Sigma2, symmetric = TRUE)$values
  return(list(Sigma2 = Sigma2, dist = norm((Sigma2 - Sigma), type = "F"),
              eig0 = eig0, eig2 = eig2))
}

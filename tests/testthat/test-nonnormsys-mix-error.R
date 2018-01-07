skip_on_cran()
library("SimRepeat")
context("Simulate correlated system of continuous variables with mixture error terms")

options(scipen = 999)
tol <- 1e-5

seed <- 276
n <- 10000
M <- 3

L <- calc_theory("Logistic", c(0, 1))
C <- calc_theory("Chisq", 4)
B <- calc_theory("Beta", c(4, 1.5))

skews <- lapply(seq_len(M), function(x) c(L[3], B[3]))
skurts <- lapply(seq_len(M), function(x) c(L[4], B[4]))
fifths <- lapply(seq_len(M), function(x) c(L[5], B[5]))
sixths <- lapply(seq_len(M), function(x) c(L[6], B[6]))
Six <- lapply(seq_len(M), function(x) list(1.75, 0.03))

mix_pis <- lapply(seq_len(M), function(x) list(c(0.4, 0.6), c(0.3, 0.2, 0.5)))
mix_mus <- lapply(seq_len(M), function(x) list(c(-2, 2), c(L[1], C[1], B[1])))
mix_sigmas <- lapply(seq_len(M),
  function(x) list(c(1, 1), c(L[2], C[2], B[2])))
mix_skews <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[3], C[3], B[3])))
mix_skurts <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[4], C[4], B[4])))
mix_fifths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[5], C[5], B[5])))
mix_sixths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[6], C[6], B[6])))
mix_Six <- lapply(seq_len(M), function(x) list(NULL, NULL, 1.75, NULL, 0.03))
Nstcum <- calc_mixmoments(mix_pis[[1]][[1]], mix_mus[[1]][[1]],
  mix_sigmas[[1]][[1]], mix_skews[[1]][[1]], mix_skurts[[1]][[1]],
  mix_fifths[[1]][[1]], mix_sixths[[1]][[1]])
Mstcum <- calc_mixmoments(mix_pis[[2]][[2]], mix_mus[[2]][[2]],
  mix_sigmas[[2]][[2]], mix_skews[[2]][[2]], mix_skurts[[2]][[2]],
  mix_fifths[[2]][[2]], mix_sixths[[2]][[2]])

same.var <- 1
betas.0 <- rep(0, M)

means <- lapply(seq_len(M),
  function(x) c(L[1], B[1], Nstcum[1], Mstcum[1]))
vars <- lapply(seq_len(M),
  function(x) c(L[2]^2, B[2]^2, Nstcum[2]^2, Mstcum[2]^2))

corr.x <- list()
corr.x[[1]] <- corr.x[[2]] <- corr.x[[3]] <- list()
corr.x[[1]][[1]] <- matrix(0.1, 4, 4)
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[2]] <- matrix(0.2, 4, 4)
corr.x[[1]][[3]] <- matrix(0.3, 4, 4)
colnames(corr.x[[1]][[1]]) <- rownames(corr.x[[1]][[1]]) <-
  rownames(corr.x[[1]][[2]]) <- rownames(corr.x[[1]][[3]]) <-
  c("C1", "C12", "M11_1", "M11_2")
colnames(corr.x[[1]][[2]]) <- c("C1", "C22", "M21_1", "M21_2")
colnames(corr.x[[1]][[3]]) <- c("C1", "C32", "M31_1", "M31_2")

# set correlations between comp of same mixture variable to 0
corr.x[[1]][[1]]["M11_1", "M11_2"] <- corr.x[[1]][[1]]["M11_2", "M11_1"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]][[1]] <- t(corr.x[[1]][[2]])
corr.x[[2]][[2]] <- matrix(0.35, 4, 4)
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[3]] <- matrix(0.4, 4, 4)
rownames(corr.x[[2]][[2]]) <- rownames(corr.x[[2]][[3]]) <-
  colnames(corr.x[[2]][[2]]) <- colnames(corr.x[[1]][[2]])
colnames(corr.x[[2]][[3]]) <- colnames(corr.x[[1]][[3]])
# set correlations between comp of same mixture variable to 0
corr.x[[2]][[2]]["M21_1", "M21_2"] <- corr.x[[2]][[2]]["M21_2", "M21_1"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[2]][[2]][, same.var] <- corr.x[[2]][[3]][, same.var] <-
  corr.x[[1]][[2]][same.var, ]
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[2]][[2]][same.var, ] <- corr.x[[2]][[2]][, same.var]

corr.x[[3]][[1]] <- t(corr.x[[1]][[3]])
corr.x[[3]][[2]] <- t(corr.x[[2]][[3]])
corr.x[[3]][[3]] <- matrix(0.5, 4, 4); diag(corr.x[[3]][[3]]) <- 1
rownames(corr.x[[3]][[3]]) <- colnames(corr.x[[3]][[3]]) <-
  colnames(corr.x[[1]][[3]])
# set correlations between comp of same mixture variable to 0
corr.x[[3]][[3]]["M31_1", "M31_2"] <- corr.x[[3]][[3]]["M31_2", "M31_1"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[3]][[3]][, same.var] <- corr.x[[1]][[3]][same.var, ]
corr.x[[3]][[3]][same.var, ] <- corr.x[[3]][[3]][, same.var]

corr.yx <- list(matrix(0.3, 1, 4), matrix(0.4, 1, 4), matrix(0.5, 1, 4))
corr.yx2 <- lapply(seq_len(M), function(x) matrix(1, 1, 3))

for (i in 1:M) {
  if (is.null(corr.yx[[i]])) next
  corr.yx2[[i]][1, 1:length(skews[[i]])] <-
    corr.yx[[i]][1, 1:length(skews[[i]])]
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]][-length(mix_pis[[i]])])))
  for (j in 1:length(mix_pis[[i]][-length(mix_pis[[i]])])) {
    corr.yx2[[i]][1, length(skews[[i]]) + j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (length(skews[[i]]) + k.cmix[j] + 1):(length(skews[[i]])
                                                        + k.cmix[j + 1])])
  }
}

corr.e <- matrix(0.4, 9, 9)
colnames(corr.e) <- rownames(corr.e) <- c("E1_1", "E1_2", "E1_3", "E2_1",
  "E2_2", "E2_3", "E3_1", "E3_2", "E3_3")
corr.e[c("E1_1", "E1_2", "E1_3"), c("E1_1", "E1_2", "E1_3")] <- diag(3)
corr.e[c("E2_1", "E2_2", "E2_3"), c("E2_1", "E2_2", "E2_3")] <- diag(3)
corr.e[c("E3_1", "E3_2", "E3_3"), c("E3_1", "E3_2", "E3_3")] <- diag(3)
error_type = "mix"

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y <- calc_corr_y(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                 error_type)
YE <- calc_corr_ye(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
YX <- calc_corr_yx(N$betas, corr.x, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 2 continuous, 1 mixture,
          error mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.2978052, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.10259813, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[1]][1, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.2348418, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.11041695, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[1]][1, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 2 continuous, 1 mixture,
          error mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
})

fifths <- lapply(seq_len(M), function(x) c(L[5], B[5]))
sixths <- lapply(seq_len(M), function(x) c(L[6], B[6]))
Six <- lapply(seq_len(M), function(x) list(1.75, 0.03))
mix_fifths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[5], C[5], B[5])))
mix_sixths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[6], C[6], B[6])))
mix_Six <- lapply(seq_len(M), function(x) list(NULL, NULL, 1.75, NULL, 0.03))

skews <- append(list(NULL), skews[2:3])
skurts <- append(list(NULL), skurts[2:3])
fifths <- append(list(NULL), fifths[2:3])
sixths <- append(list(NULL), sixths[2:3])
Six <- append(list(NULL), Six[2:3])

mix_pis[[1]] <- mix_pis[[1]][-1]
mix_mus[[1]] <- mix_mus[[1]][-1]
mix_sigmas[[1]] <- mix_sigmas[[1]][-1]
mix_skews[[1]] <- mix_skews[[1]][-1]
mix_skurts[[1]] <- mix_skurts[[1]][-1]
mix_fifths[[1]] <- mix_fifths[[1]][-1]
mix_sixths[[1]] <- mix_sixths[[1]][-1]
mix_Six[[1]] <- mix_Six[[1]][-c(1, 2)]

same.var <- matrix(c(2, 1, 3, 1), 1, 4)
betas.0 <- rep(0, M)
same <- "C1"

means <- append(list(Mstcum[1]), lapply(seq_len(M - 1),
  function(x) c(L[1], B[1], Nstcum[1], Mstcum[1])))
vars <- append(list(Mstcum[2]^2), lapply(seq_len(M - 1),
  function(x) c(L[2]^2, B[2]^2, Nstcum[2]^2, Mstcum[2]^2)))

corr.x <- list(NULL, list(), list())
corr.x[[2]] <- list(NULL, matrix(0.25, 4, 4), matrix(0.3, 4, 4))
diag(corr.x[[2]][[2]]) <- 1
rownames(corr.x[[2]][[2]]) <- rownames(corr.x[[2]][[3]]) <-
  colnames(corr.x[[2]][[2]]) <- c("C1", "C22", "M21_1", "M21_2")
colnames(corr.x[[2]][[3]]) <- c("C1", "C32", "M31_1", "M31_2")
# set correlations between comp of same mixture variable to 0
corr.x[[2]][[2]]["M21_1", "M21_2"] <- corr.x[[2]][[2]]["M21_2", "M21_1"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[2]][[3]][, same] <- corr.x[[2]][[2]][, same]

corr.x[[3]] <- list(NULL, t(corr.x[[2]][[3]]), matrix(0.35, 4, 4))
diag(corr.x[[3]][[3]]) <- 1
rownames(corr.x[[3]][[3]]) <- colnames(corr.x[[3]][[3]]) <-
  colnames(corr.x[[2]][[3]])
# set correlations between comp of same mixture variable to 0
corr.x[[3]][[3]]["M31_1", "M31_2"] <- corr.x[[3]][[3]]["M31_2", "M31_1"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[3]][[3]][, same] <- corr.x[[2]][[3]][same, ]
corr.x[[3]][[3]][same, ] <- corr.x[[3]][[3]][, same]

corr.yx <- list(NULL, matrix(0.3, 1, 4), matrix(0.4, 1, 4))
corr.yx2 <- list(NULL, matrix(1, 1, 3), matrix(1, 1, 3))

for (i in 1:M) {
  if (is.null(corr.yx[[i]])) next
  corr.yx2[[i]][1, 1:length(skews[[i]])] <-
    corr.yx[[i]][1, 1:length(skews[[i]])]
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]][-length(mix_pis[[i]])])))
  for (j in 1:length(mix_pis[[i]][-length(mix_pis[[i]])])) {
    corr.yx2[[i]][1, length(skews[[i]]) + j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (length(skews[[i]]) + k.cmix[j] + 1):(length(skews[[i]])
                                                        + k.cmix[j + 1])])
  }
}

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y <- calc_corr_y(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                 error_type)
YE <- calc_corr_ye(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
YX <- calc_corr_yx(N$betas, corr.x, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 2 continuous, 1 mixture,
          error mix, Y1 has no X terms", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[2]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[2]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[2]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.1070921, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.12164848, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[2]][2, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.1120992, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.1216485, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[2]][2, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 2 continuous, 1 mixture,
          error mix, Y1 has no X terms", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[2]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[2]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[2]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
})

fifths <- lapply(seq_len(M), function(x) c(L[5], B[5]))
sixths <- lapply(seq_len(M), function(x) c(L[6], B[6]))
Six <- lapply(seq_len(M), function(x) list(1.75, 0.03))
mix_fifths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[5], C[5], B[5])))
mix_sixths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[6], C[6], B[6])))
mix_Six <- lapply(seq_len(M), function(x) list(NULL, NULL, 1.75, NULL, 0.03))

skews <- list(skews[[2]], NULL, skews[[3]])
skurts <- list(skurts[[2]], NULL, skurts[[3]])
fifths <- list(fifths[[2]], NULL, fifths[[3]])
sixths <- list(sixths[[2]], NULL, sixths[[3]])
Six <- list(Six[[2]], NULL, Six[[3]])

mix_pis <- list(mix_pis[[2]], mix_pis[[1]], mix_pis[[3]])
mix_mus <- list(mix_mus[[2]], mix_mus[[1]], mix_mus[[3]])
mix_sigmas <- list(mix_sigmas[[2]], mix_sigmas[[1]], mix_sigmas[[3]])
mix_skews <- list(mix_skews[[2]], mix_skews[[1]], mix_skews[[3]])
mix_skurts <- list(mix_skurts[[2]], mix_skurts[[1]], mix_skurts[[3]])
mix_fifths <- list(mix_fifths[[2]], list(c(L[5], C[5], B[5])), mix_fifths[[3]])
mix_sixths <- list(mix_sixths[[2]], list(c(L[6], C[6], B[6])), mix_sixths[[3]])
mix_Six <- list(mix_Six[[2]], list(1.75, NULL, 0.03), mix_Six[[3]])

same.var <- matrix(c(1, 1, 3, 1), 1, 4)
betas.0 <- rep(0, M)
same <- "C1"

means <- list(means[[2]], Mstcum[1], means[[3]])
vars <- list(vars[[2]], Mstcum[2]^2, vars[[3]])

corr.x <- list(list(), NULL, list())
corr.x[[1]] <- list(matrix(0.25, 4, 4), NULL, matrix(0.3, 4, 4))
diag(corr.x[[1]][[1]]) <- 1
rownames(corr.x[[1]][[1]]) <- rownames(corr.x[[1]][[3]]) <-
  colnames(corr.x[[1]][[1]]) <- c("C1", "C12", "M11_1", "M11_2")
colnames(corr.x[[1]][[3]]) <- c("C1", "C32", "M31_1", "M31_2")
# set correlations between comp of same mixture variable to 0
corr.x[[1]][[1]]["M11_1", "M11_2"] <- corr.x[[1]][[1]]["M11_2", "M11_1"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[1]][[3]][, same] <- corr.x[[1]][[1]][, same]

corr.x[[3]] <- list(t(corr.x[[1]][[3]]), NULL, matrix(0.35, 4, 4))
diag(corr.x[[3]][[3]]) <- 1
rownames(corr.x[[3]][[3]]) <- colnames(corr.x[[3]][[3]]) <-
  colnames(corr.x[[1]][[3]])
# set correlations between comp of same mixture variable to 0
corr.x[[3]][[3]]["M31_1", "M31_2"] <- corr.x[[3]][[3]]["M31_2", "M31_1"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[3]][[3]][, same] <- corr.x[[1]][[3]][same, ]
corr.x[[3]][[3]][same, ] <- corr.x[[3]][[3]][, same]

corr.yx <- list(matrix(0.3, 1, 4), NULL, matrix(0.4, 1, 4))
corr.yx2 <- list(matrix(1, 1, 3), NULL, matrix(1, 1, 3))

for (i in 1:M) {
  if (is.null(corr.yx[[i]])) next
  corr.yx2[[i]][1, 1:length(skews[[i]])] <-
    corr.yx[[i]][1, 1:length(skews[[i]])]
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]][-length(mix_pis[[i]])])))
  for (j in 1:length(mix_pis[[i]][-length(mix_pis[[i]])])) {
    corr.yx2[[i]][1, length(skews[[i]]) + j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (length(skews[[i]]) + k.cmix[j] + 1):(length(skews[[i]])
                                                        + k.cmix[j + 1])])
  }
}

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y <- calc_corr_y(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                 error_type)
YE <- calc_corr_ye(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
YX <- calc_corr_yx(N$betas, corr.x, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 2 continuous, 1 mixture,
          error mix, Y2 has no X terms", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.10709213, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.10709213, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[1]][1, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.1120992, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.1120992, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[1]][1, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 2 continuous, 1 mixture,
          error mix, Y2 has no X terms", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
})

fifths <- lapply(seq_len(M), function(x) c(L[5], B[5]))
sixths <- lapply(seq_len(M), function(x) c(L[6], B[6]))
Six <- lapply(seq_len(M), function(x) list(1.75, 0.03))
mix_fifths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[5], C[5], B[5])))
mix_sixths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[6], C[6], B[6])))
mix_Six <- lapply(seq_len(M), function(x) list(NULL, NULL, 1.75, NULL, 0.03))

skews <- list(skews[[1]], skews[[3]], NULL)
skurts <- list(skurts[[1]], skurts[[3]], NULL)
fifths <- list(fifths[[1]], fifths[[3]], NULL)
sixths <- list(sixths[[1]], sixths[[3]], NULL)
Six <- list(Six[[1]], Six[[3]], NULL)

mix_pis <- list(mix_pis[[1]], mix_pis[[3]], mix_pis[[2]])
mix_mus <- list(mix_mus[[1]], mix_mus[[3]], mix_mus[[2]])
mix_sigmas <- list(mix_sigmas[[1]], mix_sigmas[[3]], mix_sigmas[[2]])
mix_skews <- list(mix_skews[[1]], mix_skews[[3]], mix_skews[[2]])
mix_skurts <- list(mix_skurts[[1]], mix_skurts[[3]], mix_skurts[[2]])
mix_fifths <- list(mix_fifths[[1]], mix_fifths[[3]], list(c(L[5], C[5], B[5])))
mix_sixths <- list(mix_sixths[[1]], mix_sixths[[3]], list(c(L[6], C[6], B[6])))
mix_Six <- list(mix_Six[[1]], mix_Six[[3]], list(1.75, NULL, 0.03))

same.var <- matrix(c(1, 1, 2, 1), 1, 4)
betas.0 <- rep(0, M)
same <- "C1"

means <- list(means[[1]], means[[3]], means[[2]])
vars <- list(vars[[1]], vars[[3]], vars[[2]])

corr.x <- list(list(), list(), NULL)
corr.x[[1]] <- list(matrix(0.25, 4, 4), matrix(0.3, 4, 4), NULL)
diag(corr.x[[1]][[1]]) <- 1
rownames(corr.x[[1]][[1]]) <- rownames(corr.x[[1]][[2]]) <-
  colnames(corr.x[[1]][[1]]) <- c("C1", "C12", "M11_1", "M11_2")
colnames(corr.x[[1]][[2]]) <- c("C1", "C22", "M21_1", "M21_2")
# set correlations between comp of same mixture variable to 0
corr.x[[1]][[1]]["M11_1", "M11_2"] <- corr.x[[1]][[1]]["M11_2", "M11_1"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[1]][[2]][, same] <- corr.x[[1]][[1]][, same]

corr.x[[2]] <- list(t(corr.x[[1]][[2]]), matrix(0.35, 4, 4), NULL)
diag(corr.x[[2]][[2]]) <- 1
rownames(corr.x[[2]][[2]]) <- colnames(corr.x[[2]][[2]]) <-
  colnames(corr.x[[1]][[2]])
# set correlations between comp of same mixture variable to 0
corr.x[[2]][[2]]["M21_1", "M21_2"] <- corr.x[[2]][[2]]["M21_2", "M21_1"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[2]][[2]][, same] <- corr.x[[1]][[2]][same, ]
corr.x[[2]][[2]][same, ] <- corr.x[[2]][[2]][, same]

corr.yx <- list(matrix(0.3, 1, 4), matrix(0.4, 1, 4), NULL)
corr.yx2 <- list(matrix(1, 1, 3), matrix(1, 1, 3), NULL)

for (i in 1:M) {
  if (is.null(corr.yx[[i]])) next
  corr.yx2[[i]][1, 1:length(skews[[i]])] <-
    corr.yx[[i]][1, 1:length(skews[[i]])]
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]][-length(mix_pis[[i]])])))
  for (j in 1:length(mix_pis[[i]][-length(mix_pis[[i]])])) {
    corr.yx2[[i]][1, length(skews[[i]]) + j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (length(skews[[i]]) + k.cmix[j] + 1):(length(skews[[i]])
                                                        + k.cmix[j + 1])])
  }
}

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y <- calc_corr_y(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                 error_type)
YE <- calc_corr_ye(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
YX <- calc_corr_yx(N$betas, corr.x, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 2 continuous, 1 mixture,
          error mix, Y3 has no X terms", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.30430330, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.1070921, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[1]][1, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.2373529, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.1120992, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[1]][1, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 2 continuous, 1 mixture,
          error mix, Y3 has no X terms", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- skurts <- fifths <- sixths <- Six <- list()
mix_pis[[3]] <- mix_pis[[2]]
mix_mus[[3]] <- mix_mus[[2]]
mix_sigmas[[3]] <- mix_sigmas[[2]]
mix_skews[[3]] <- mix_skews[[2]]
mix_skurts[[3]] <- mix_skurts[[2]]
mix_fifths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[5], C[5], B[5])))
mix_sixths <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[6], C[6], B[6])))
mix_Six <- lapply(seq_len(M), function(x) list(NULL, NULL, 1.75, NULL, 0.03))

same.var <- c(1, 2)
betas.0 <- rep(0, M)

means <- lapply(seq_len(M), function(x) c(Nstcum[1], Mstcum[1]))
vars <- lapply(seq_len(M), function(x) c(Nstcum[2]^2, Mstcum[2]^2))

corr.x <- lapply(seq_len(M), function(x) list(diag(2), diag(2), diag(2)))

corr.yx <- list(matrix(0.3, 1, 2), matrix(0.35, 1, 2), matrix(0.4, 1, 2))
corr.yx2 <- list(matrix(1, 1, 1), matrix(1, 1, 1), matrix(1, 1, 1))

for (i in 1:M) {
  if (is.null(corr.yx[[i]])) next
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]][-length(mix_pis[[i]])])))
  for (j in 1:length(mix_pis[[i]][-length(mix_pis[[i]])])) {
    corr.yx2[[i]][1, j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (k.cmix[j] + 1):(k.cmix[j + 1])])
  }
}

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y <- calc_corr_y(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                 error_type)
YE <- calc_corr_ye(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
YX <- calc_corr_yx(N$betas, corr.x, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 1 mixture for all M,
          error mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][3, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][3, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][3, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.3057165, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.1101574, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[1]][1, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.1406715, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.1205121, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[1]][1, 1],
    0.1363636, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 1 mixture for all M,
          error mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][3, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][3, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][3, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
})

M <- 2
skews <- skurts <- fifths <- sixths <- Six <- list()
mix_pis <- list(list(c(0.4, 0.6), c(0.3, 0.2, 0.5)),
                list(c(0.3, 0.2, 0.5)))
mix_mus <- list(list(c(-2, 2), c(L[1], C[1], B[1])),
                list(c(L[1], C[1], B[1])))
mix_sigmas <- list(list(c(1, 1), c(L[2], C[2], B[2])),
                   list(c(L[2], C[2], B[2])))
mix_skews <- list(list(c(0, 0), c(L[3], C[3], B[3])),
                  list(c(L[3], C[3], B[3])))
mix_skurts <- list(list(c(0, 0), c(L[4], C[4], B[4])),
                   list(c(L[4], C[4], B[4])))
mix_fifths <- list(list(c(0, 0), c(L[5], C[5], B[5])),
                   list(c(L[5], C[5], B[5])))
mix_sixths <- list(list(c(0, 0), c(L[6], C[6], B[6])),
                   list(c(L[6], C[6], B[6])))
mix_Six <- list(list(NULL, NULL, 1.75, NULL, 0.03), list(1.75, NULL, 0.03))

betas.0 <- rep(0, M)

means <- list(c(Nstcum[1], Mstcum[1]), Mstcum[1])
vars <- list(c(Nstcum[2]^2, Mstcum[2]^2), Mstcum[2]^2)

corr.x <- list(list(diag(2), NULL), NULL)

corr.yx <- list(matrix(0.3, 1, 2), NULL)
corr.yx2 <- list(matrix(1, 1, 1), NULL)

for (i in 1:M) {
  if (is.null(corr.yx[[i]])) next
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]][-length(mix_pis[[i]])])))
  for (j in 1:length(mix_pis[[i]][-length(mix_pis[[i]])])) {
    corr.yx2[[i]][1, j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (k.cmix[j] + 1):(k.cmix[j + 1])])
  }
}

corr.e <- corr.e[1:6, 1:6]
same.var <- NULL

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y <- calc_corr_y(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                 error_type)
YE <- calc_corr_ye(N$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
YX <- calc_corr_yx(N$betas, corr.x, vars, mix_pis, mix_mus, mix_sigmas,
                   error_type)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 1 mixture for Y1, no X for Y2,
          error mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][3, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][3, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][3, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.1101574, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.1101574, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[1]][1, 2],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.1205121, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.1205121, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[1]][1, 1],
    0.1363636, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, N3$E_mix, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 1 mixture for Y1, no X for Y2,
          error mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx2, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][3, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][3, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N3$constants[[1]][3, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
})

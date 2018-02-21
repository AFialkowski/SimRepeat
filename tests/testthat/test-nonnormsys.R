skip_on_cran()
library("SimRepeat")
context("Simulate correlated system of continuous variables")

options(scipen = 999)
tol <- 1e-4

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
  function(x) c(L[1], Nstcum[1], Mstcum[1], B[1]))
vars <- lapply(seq_len(M),
  function(x) c(L[2]^2, Nstcum[2]^2, Mstcum[2]^2, B[2]^2))

corr.x <- list()
corr.x[[1]] <- corr.x[[2]] <- corr.x[[3]] <- list()
corr.x[[1]][[1]] <- matrix(0.1, 6, 6)
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[2]] <- matrix(0.2, 6, 6)
corr.x[[1]][[3]] <- matrix(0.3, 6, 6)
colnames(corr.x[[1]][[1]]) <- rownames(corr.x[[1]][[1]]) <-
  rownames(corr.x[[1]][[2]]) <- rownames(corr.x[[1]][[3]]) <-
  c("C1", "M11_1", "M11_2", "M12_1", "M12_2", "M12_3")
colnames(corr.x[[1]][[2]]) <- c("C1", "M21_1", "M21_2", "M22_1",
  "M22_2", "M22_3")
colnames(corr.x[[1]][[3]]) <- c("C1", "M31_1", "M31_2", "M32_1",
  "M32_2", "M32_3")

# set correlations between comp of same mixture variable to 0
corr.x[[1]][[1]]["M11_1", "M11_2"] <- corr.x[[1]][[1]]["M11_2", "M11_1"] <-
  corr.x[[1]][[1]]["M12_1", "M12_2"] <- corr.x[[1]][[1]]["M12_2", "M12_1"] <-
  corr.x[[1]][[1]]["M12_1", "M12_3"] <- corr.x[[1]][[1]]["M12_3", "M12_1"] <-
  corr.x[[1]][[1]]["M12_2", "M12_3"] <- corr.x[[1]][[1]]["M12_3", "M12_2"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]][[1]] <- t(corr.x[[1]][[2]])
corr.x[[2]][[2]] <- matrix(0.35, 6, 6)
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[3]] <- matrix(0.4, 6, 6)
rownames(corr.x[[2]][[2]]) <- rownames(corr.x[[2]][[3]]) <-
  colnames(corr.x[[2]][[2]]) <- colnames(corr.x[[1]][[2]])
colnames(corr.x[[2]][[3]]) <- colnames(corr.x[[1]][[3]])
# set correlations between comp of same mixture variable to 0
corr.x[[2]][[2]]["M21_1", "M21_2"] <- corr.x[[2]][[2]]["M21_2", "M21_1"] <-
  corr.x[[2]][[2]]["M22_1", "M22_2"] <- corr.x[[2]][[2]]["M22_2", "M22_1"] <-
  corr.x[[2]][[2]]["M22_1", "M22_3"] <- corr.x[[2]][[2]]["M22_3", "M22_1"] <-
  corr.x[[2]][[2]]["M22_2", "M22_3"] <- corr.x[[2]][[2]]["M22_3", "M22_2"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[2]][[2]][, same.var] <- corr.x[[2]][[3]][, same.var] <-
  corr.x[[1]][[2]][same.var, ]
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[2]][[2]][same.var, ] <- corr.x[[2]][[2]][, same.var]

corr.x[[3]][[1]] <- t(corr.x[[1]][[3]])
corr.x[[3]][[2]] <- t(corr.x[[2]][[3]])
corr.x[[3]][[3]] <- matrix(0.5, 6, 6); diag(corr.x[[3]][[3]]) <- 1
rownames(corr.x[[3]][[3]]) <- colnames(corr.x[[3]][[3]]) <-
  colnames(corr.x[[1]][[3]])
# set correlations between comp of same mixture variable to 0
corr.x[[3]][[3]]["M31_1", "M31_2"] <- corr.x[[3]][[3]]["M31_2", "M31_1"] <-
  corr.x[[3]][[3]]["M32_1", "M32_2"] <- corr.x[[3]][[3]]["M32_2", "M32_1"] <-
  corr.x[[3]][[3]]["M32_1", "M32_3"] <- corr.x[[3]][[3]]["M32_3", "M32_1"] <-
  corr.x[[3]][[3]]["M32_2", "M32_3"] <- corr.x[[3]][[3]]["M32_3", "M32_2"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[3]][[3]][, same.var] <- corr.x[[1]][[3]][same.var, ]
corr.x[[3]][[3]][same.var, ] <- corr.x[[3]][[3]][, same.var]

corr.yx <- list(matrix(0.3, 1, 6), matrix(0.4, 1, 6), matrix(0.5, 1, 6))
corr.yx2 <- lapply(seq_len(M), function(x) matrix(1, 1, 3))

for (i in 1:M) {
  corr.yx2[[i]][1, 1:(length(skews[[i]]) - 1)] <-
    corr.yx[[i]][1, 1:(length(skews[[i]]) - 1)]
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]])))
  for (j in 1:length(mix_pis[[i]])) {
    corr.yx2[[i]][1, (length(skews[[i]]) - 1) + j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (length(skews[[i]]) + k.cmix[j]):(length(skews[[i]]) -
        1 + k.cmix[j + 1])])
  }
}

corr.e <- matrix(0.4, nrow = M, ncol = M)
diag(corr.e) <- 1
error_type = "non_mix"

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
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
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, E_mix = NULL, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 1 continuous, 2 mixture,
          error non_mix", {
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
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.5620107, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.310582, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[1]][1, 1],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.456203, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.3738404, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[1]][1, 1],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, E_mix = NULL, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 1 continuous, 2 mixture,
          error non_mix", {
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
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- lapply(seq_len(M), function(x) c(L[3], B[3]))
skurts <- lapply(seq_len(M), function(x) c(L[4], B[4]))
fifths <- lapply(seq_len(M), function(x) c(L[5], B[5]))
sixths <- lapply(seq_len(M), function(x) c(L[6], B[6]))
Six <- lapply(seq_len(M), function(x) list(1.75, 0.03))
mix_pis <- mix_mus <- mix_sigmas <- mix_skews <- mix_skurts <-
  mix_fifths <- mix_sixths <- mix_Six <- list()
same.var <- 1
betas.0 <- rep(0, M)

means <- lapply(seq_len(M), function(x) c(L[1], B[1]))
vars <- lapply(seq_len(M), function(x) c(L[2]^2, B[2]^2))

corr.x <- lapply(seq_len(M), function(x) list(matrix(1, 1, 1), matrix(1, 1, 1),
                                              matrix(1, 1, 1)))

corr.yx <- list(matrix(0.3, 1, 1), matrix(0.4, 1, 1), matrix(0.5, 1, 1))

corr.e <- matrix(0.4, nrow = M, ncol = M)
diag(corr.e) <- 1
error_type = "non_mix"

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y <- calc_corr_y(N$betas, corr.x, corr.e, vars)
YE <- calc_corr_ye(N$betas, corr.x, corr.e, vars)
YX <- calc_corr_yx(N$betas, corr.x, vars)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Polynomial method: 1 continuous all M,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.4697199, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.3815757, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[1]][1, 1],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 1 continuous all M,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

M <- 2
skews <- list(B[3], c(L[3], B[3]))
skurts <- list(B[4], c(L[4], B[4]))
fifths <- list(B[5], c(L[5], B[5]))
sixths <- list(B[6], c(L[6], B[6]))
Six <- list(0.03, list(1.75, 0.03))
mix_pis <- mix_mus <- mix_sigmas <- mix_skews <- mix_skurts <-
  mix_fifths <- mix_sixths <- mix_Six <- list()
same.var <- NULL
betas.0 <- rep(0, M)

means <- list(B[1], c(L[1], B[1]))
vars <- list(B[2]^2, c(L[2]^2, B[2]^2))

corr.x <- list(NULL, list(NULL, matrix(1, 1, 1)))

corr.yx <- list(NULL, matrix(0.3, 1, 1))

corr.e <- matrix(0.4, nrow = M, ncol = M)
diag(corr.e) <- 1
error_type = "non_mix"

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y <- calc_corr_y(N$betas, corr.x, corr.e, vars)
YE <- calc_corr_ye(N$betas, corr.x, corr.e, vars)
YX <- calc_corr_yx(N$betas, corr.x, vars)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Polynomial method: 1 continuous for Y2, no X for Y1,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[2]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[2]][1, "c3"],
    0.03605225, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.3815757, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.4, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[2]][2, 1],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 1 continuous for Y2, no X for Y1,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means,
    vars, skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six,
    same.var = same.var, betas.0 = betas.0, corr.x = corr.x,
    corr.yx = corr.yx, corr.e = corr.e), TRUE)
  expect_equal(all.equal(N$constants[[2]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[2]][1, "c3"],
    0.03605955, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

M <- 3
skews <- append(list(B[3]), lapply(seq_len(M - 1), function(x) c(L[3], B[3])))
skurts <- append(list(B[4]), lapply(seq_len(M - 1), function(x) c(L[4], B[4])))
fifths <- append(list(B[5]), lapply(seq_len(M - 1), function(x) c(L[5], B[5])))
sixths <- append(list(B[6]), lapply(seq_len(M - 1), function(x) c(L[6], B[6])))
Six <- append(list(0.03), lapply(seq_len(M - 1), function(x) list(1.75, 0.03)))

mix_pis <- append(list(NULL), lapply(seq_len(M - 1),
  function(x) list(c(0.4, 0.6), c(0.3, 0.2, 0.5))))
mix_mus <- append(list(NULL), lapply(seq_len(M - 1),
  function(x) list(c(-2, 2), c(L[1], C[1], B[1]))))
mix_sigmas <- append(list(NULL), lapply(seq_len(M - 1),
  function(x) list(c(1, 1), c(L[2], C[2], B[2]))))
mix_skews <- append(list(NULL), lapply(seq_len(M - 1),
  function(x) list(rep(0, 2), c(L[3], C[3], B[3]))))
mix_skurts <- append(list(NULL), lapply(seq_len(M - 1),
  function(x) list(rep(0, 2), c(L[4], C[4], B[4]))))
mix_fifths <- append(list(NULL), lapply(seq_len(M - 1),
  function(x) list(rep(0, 2), c(L[5], C[5], B[5]))))
mix_sixths <- append(list(NULL), lapply(seq_len(M - 1),
  function(x) list(rep(0, 2), c(L[6], C[6], B[6]))))
mix_Six <- append(list(NULL), lapply(seq_len(M - 1),
  function(x) list(NULL, NULL, 1.75, NULL, 0.03)))
Nstcum <- calc_mixmoments(mix_pis[[2]][[1]], mix_mus[[2]][[1]],
  mix_sigmas[[2]][[1]], mix_skews[[2]][[1]], mix_skurts[[2]][[1]],
  mix_fifths[[2]][[1]], mix_sixths[[2]][[1]])
Mstcum <- calc_mixmoments(mix_pis[[2]][[2]], mix_mus[[2]][[2]],
  mix_sigmas[[2]][[2]], mix_skews[[2]][[2]], mix_skurts[[2]][[2]],
  mix_fifths[[2]][[2]], mix_sixths[[2]][[2]])

same.var <- matrix(c(2, 1, 3, 1), 1, 4)
same <- "C1"
betas.0 <- rep(0, M)

means <- append(list(B[1]), lapply(seq_len(M - 1),
  function(x) c(L[1], Nstcum[1], Mstcum[1], B[1])))
vars <- append(list(B[2]^2), lapply(seq_len(M - 1),
  function(x) c(L[2]^2, Nstcum[2]^2, Mstcum[2]^2, B[2]^2)))

corr.x <- list(NULL, list(), list())
corr.x[[2]] <- list(NULL, matrix(0.35, 6, 6), matrix(0.4, 6, 6))
diag(corr.x[[2]][[2]]) <- 1
rownames(corr.x[[2]][[2]]) <- rownames(corr.x[[2]][[3]]) <-
  colnames(corr.x[[2]][[2]]) <- c("C1", "M21_1", "M21_2", "M22_1",
  "M22_2", "M22_3")
colnames(corr.x[[2]][[3]]) <- c("C1", "M31_1", "M31_2", "M32_1",
  "M32_2", "M32_3")
# set correlations between comp of same mixture variable to 0
corr.x[[2]][[2]]["M21_1", "M21_2"] <- corr.x[[2]][[2]]["M21_2", "M21_1"] <-
  corr.x[[2]][[2]]["M22_1", "M22_2"] <- corr.x[[2]][[2]]["M22_2", "M22_1"] <-
  corr.x[[2]][[2]]["M22_1", "M22_3"] <- corr.x[[2]][[2]]["M22_3", "M22_1"] <-
  corr.x[[2]][[2]]["M22_2", "M22_3"] <- corr.x[[2]][[2]]["M22_3", "M22_2"] <- 0
# since C1 same across outcomes, set correlation the same
corr.x[[2]][[3]][, same] <- corr.x[[2]][[2]][, same]

corr.x[[3]] <- list(NULL, t(corr.x[[2]][[3]]), matrix(0.5, 6, 6))
diag(corr.x[[3]][[3]]) <- 1
rownames(corr.x[[3]][[3]]) <- colnames(corr.x[[3]][[3]]) <-
  colnames(corr.x[[2]][[3]])
# set correlations between comp of same mixture variable to 0
corr.x[[3]][[3]]["M31_1", "M31_2"] <- corr.x[[3]][[3]]["M31_2", "M31_1"] <-
  corr.x[[3]][[3]]["M32_1", "M32_2"] <- corr.x[[3]][[3]]["M32_2", "M32_1"] <-
  corr.x[[3]][[3]]["M32_1", "M32_3"] <- corr.x[[3]][[3]]["M32_3", "M32_1"] <-
  corr.x[[3]][[3]]["M32_2", "M32_3"] <- corr.x[[3]][[3]]["M32_3", "M32_2"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[3]][[3]][, same] <- corr.x[[2]][[3]][same, ]
corr.x[[3]][[3]][same, ] <- corr.x[[3]][[3]][, same]

corr.yx <- list(NULL, matrix(0.4, 1, 6), matrix(0.5, 1, 6))
corr.yx2 <- append(list(NULL), lapply(seq_len(M - 1),
  function(x) matrix(1, 1, 3)))

for (i in 1:M) {
  if (is.null(corr.yx[[i]])) next
  corr.yx2[[i]][1, 1:(length(skews[[i]]) - 1)] <-
    corr.yx[[i]][1, 1:(length(skews[[i]]) - 1)]
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]])))
  for (j in 1:length(mix_pis[[i]])) {
    corr.yx2[[i]][1, (length(skews[[i]]) - 1) + j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (length(skews[[i]]) + k.cmix[j]):(length(skews[[i]]) -
        1 + k.cmix[j + 1])])
  }
}

corr.e <- matrix(0.4, nrow = M, ncol = M)
diag(corr.e) <- 1
error_type = "non_mix"

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
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
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, E_mix = NULL, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 1 continuous, 2 mixture,
          error non_mix, Y1 has no X terms", {
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
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.2873793, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.4, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[2]][2, 1],
    0.4, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.3591677, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.4, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[2]][2, 1],
    0.4, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, E_mix = NULL, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 1 continuous, 2 mixture,
          error non_mix, Y1 has no X terms", {
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
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- list(c(L[3], B[3]), B[3], c(L[3], B[3]))
skurts <- list(c(L[4], B[4]), B[4], c(L[4], B[4]))
fifths <- list(c(L[5], B[5]), B[5], c(L[5], B[5]))
sixths <- list(c(L[6], B[6]), B[6], c(L[6], B[6]))
Six <- list(list(1.75, 0.03), list(0.03), list(1.75, 0.03))

mix_pis <- list(list(c(0.4, 0.6), c(0.3, 0.2, 0.5)), NULL,
  list(c(0.4, 0.6), c(0.3, 0.2, 0.5)))
mix_mus <- list(list(c(-2, 2), c(L[1], C[1], B[1])), NULL,
  list(c(-2, 2), c(L[1], C[1], B[1])))
mix_sigmas <- list(list(c(1, 1), c(L[2], C[2], B[2])), NULL,
  list(c(1, 1), c(L[2], C[2], B[2])))
mix_skews <- list(list(rep(0, 2), c(L[3], C[3], B[3])), NULL,
  list(rep(0, 2), c(L[3], C[3], B[3])))
mix_skurts <- list(list(rep(0, 2), c(L[4], C[4], B[4])), NULL,
  list(rep(0, 2), c(L[4], C[4], B[4])))
mix_fifths <- list(list(rep(0, 2), c(L[5], C[5], B[5])), NULL,
  list(rep(0, 2), c(L[5], C[5], B[5])))
mix_sixths <- list(list(rep(0, 2), c(L[6], C[6], B[6])), NULL,
  list(rep(0, 2), c(L[6], C[6], B[6])))
mix_Six <- list(list(NULL, NULL, 1.75, NULL, 0.03), NULL,
  list(NULL, NULL, 1.75, NULL, 0.03))
Nstcum <- calc_mixmoments(mix_pis[[1]][[1]], mix_mus[[1]][[1]],
  mix_sigmas[[1]][[1]], mix_skews[[1]][[1]], mix_skurts[[1]][[1]],
  mix_fifths[[1]][[1]], mix_sixths[[1]][[1]])
Mstcum <- calc_mixmoments(mix_pis[[1]][[2]], mix_mus[[1]][[2]],
  mix_sigmas[[1]][[2]], mix_skews[[1]][[2]], mix_skurts[[1]][[2]],
  mix_fifths[[1]][[2]], mix_sixths[[1]][[2]])

same.var <- matrix(c(1, 1, 3, 1), 1, 4)
same <- "C1"
betas.0 <- rep(0, M)

means <- list(c(L[1], Nstcum[1], Mstcum[1], B[1]), B[1],
  c(L[1], Nstcum[1], Mstcum[1], B[1]))
vars <- list(c(L[2]^2, Nstcum[2]^2, Mstcum[2]^2, B[2]^2), B[2]^2,
  c(L[2]^2, Nstcum[2]^2, Mstcum[2]^2, B[2]^2))

corr.x <- list(list(), NULL, list())
corr.x[[1]] <- list(matrix(0.25, 6, 6), NULL, matrix(0.4, 6, 6))
diag(corr.x[[1]][[1]]) <- 1
rownames(corr.x[[1]][[1]]) <- rownames(corr.x[[1]][[3]]) <-
  colnames(corr.x[[1]][[1]]) <- c("C1", "M11_1", "M11_2", "M12_1",
  "M12_2", "M12_3")
colnames(corr.x[[1]][[3]]) <- c("C1", "M31_1", "M31_2", "M32_1",
  "M32_2", "M32_3")
# set correlations between comp of same mixture variable to 0
corr.x[[1]][[1]]["M11_1", "M11_2"] <- corr.x[[1]][[1]]["M11_2", "M11_1"] <-
  corr.x[[1]][[1]]["M12_1", "M12_2"] <- corr.x[[1]][[1]]["M12_2", "M12_1"] <-
  corr.x[[1]][[1]]["M12_1", "M12_3"] <- corr.x[[1]][[1]]["M12_3", "M12_1"] <-
  corr.x[[1]][[1]]["M12_2", "M12_3"] <- corr.x[[1]][[1]]["M12_3", "M12_2"] <- 0
# since C1 same across outcomes, set correlation the same
corr.x[[1]][[3]][, same] <- corr.x[[1]][[1]][, same]

corr.x[[3]] <- list(t(corr.x[[1]][[3]]), NULL, matrix(0.5, 6, 6))
diag(corr.x[[3]][[3]]) <- 1
rownames(corr.x[[3]][[3]]) <- colnames(corr.x[[3]][[3]]) <-
  colnames(corr.x[[1]][[3]])
# set correlations between comp of same mixture variable to 0
corr.x[[3]][[3]]["M31_1", "M31_2"] <- corr.x[[3]][[3]]["M31_2", "M31_1"] <-
  corr.x[[3]][[3]]["M32_1", "M32_2"] <- corr.x[[3]][[3]]["M32_2", "M32_1"] <-
  corr.x[[3]][[3]]["M32_1", "M32_3"] <- corr.x[[3]][[3]]["M32_3", "M32_1"] <-
  corr.x[[3]][[3]]["M32_2", "M32_3"] <- corr.x[[3]][[3]]["M32_3", "M32_2"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[3]][[3]][, same] <- corr.x[[1]][[3]][same, ]
corr.x[[3]][[3]][same, ] <- corr.x[[3]][[3]][, same]

corr.yx <- list(matrix(0.3, 1, 6), NULL, matrix(0.5, 1, 6))
corr.yx2 <- list(matrix(1, 1, 3), NULL, matrix(1, 1, 3))

for (i in 1:M) {
  if (is.null(corr.yx[[i]])) next
  corr.yx2[[i]][1, 1:(length(skews[[i]]) - 1)] <-
    corr.yx[[i]][1, 1:(length(skews[[i]]) - 1)]
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]])))
  for (j in 1:length(mix_pis[[i]])) {
    corr.yx2[[i]][1, (length(skews[[i]]) - 1) + j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (length(skews[[i]]) + k.cmix[j]):(length(skews[[i]]) -
        1 + k.cmix[j + 1])])
  }
}

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
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
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, E_mix = NULL, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 1 continuous, 2 mixture,
          error non_mix, Y2 has no X terms", {
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
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.3364521, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.3364521, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[1]][1, 1],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.3762730, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.3762730, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[1]][1, 1],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, E_mix = NULL, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 1 continuous, 2 mixture,
          error non_mix, Y2 has no X terms", {
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
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- list(c(L[3], B[3]), c(L[3], B[3]), B[3])
skurts <- list(c(L[4], B[4]), c(L[4], B[4]), B[4])
fifths <- list(c(L[5], B[5]), c(L[5], B[5]), B[5])
sixths <- list(c(L[6], B[6]), c(L[6], B[6]), B[6])
Six <- list(list(1.75, 0.03), list(1.75, 0.03), list(0.03))

mix_pis <- list(list(c(0.4, 0.6), c(0.3, 0.2, 0.5)),
  list(c(0.4, 0.6), c(0.3, 0.2, 0.5)), NULL)
mix_mus <- list(list(c(-2, 2), c(L[1], C[1], B[1])),
  list(c(-2, 2), c(L[1], C[1], B[1])), NULL)
mix_sigmas <- list(list(c(1, 1), c(L[2], C[2], B[2])),
  list(c(1, 1), c(L[2], C[2], B[2])), NULL)
mix_skews <- list(list(rep(0, 2), c(L[3], C[3], B[3])),
  list(rep(0, 2), c(L[3], C[3], B[3])), NULL)
mix_skurts <- list(list(rep(0, 2), c(L[4], C[4], B[4])),
  list(rep(0, 2), c(L[4], C[4], B[4])), NULL)
mix_fifths <- list(list(rep(0, 2), c(L[5], C[5], B[5])),
  list(rep(0, 2), c(L[5], C[5], B[5])), NULL)
mix_sixths <- list(list(rep(0, 2), c(L[6], C[6], B[6])),
  list(rep(0, 2), c(L[6], C[6], B[6])), NULL)
mix_Six <- list(list(NULL, NULL, 1.75, NULL, 0.03),
  list(NULL, NULL, 1.75, NULL, 0.03), NULL)
Nstcum <- calc_mixmoments(mix_pis[[1]][[1]], mix_mus[[1]][[1]],
  mix_sigmas[[1]][[1]], mix_skews[[1]][[1]], mix_skurts[[1]][[1]],
  mix_fifths[[1]][[1]], mix_sixths[[1]][[1]])
Mstcum <- calc_mixmoments(mix_pis[[1]][[2]], mix_mus[[1]][[2]],
  mix_sigmas[[1]][[2]], mix_skews[[1]][[2]], mix_skurts[[1]][[2]],
  mix_fifths[[1]][[2]], mix_sixths[[1]][[2]])

same.var <- matrix(c(1, 1, 2, 1), 1, 4)
same <- "C1"
betas.0 <- rep(0, M)

means <- list(c(L[1], Nstcum[1], Mstcum[1], B[1]),
  c(L[1], Nstcum[1], Mstcum[1], B[1]), B[1])
vars <- list(c(L[2]^2, Nstcum[2]^2, Mstcum[2]^2, B[2]^2),
  c(L[2]^2, Nstcum[2]^2, Mstcum[2]^2, B[2]^2), B[2]^2)

corr.x <- list(list(), list(), NULL)
corr.x[[1]] <- list(matrix(0.2, 6, 6), matrix(0.3, 6, 6), NULL)
diag(corr.x[[1]][[1]]) <- 1
rownames(corr.x[[1]][[1]]) <- rownames(corr.x[[1]][[2]]) <-
  colnames(corr.x[[1]][[1]]) <- c("C1", "M11_1", "M11_2", "M12_1",
  "M12_2", "M12_3")
colnames(corr.x[[1]][[2]]) <- c("C1", "M21_1", "M21_2", "M22_1",
  "M22_2", "M22_3")
# set correlations between comp of same mixture variable to 0
corr.x[[1]][[1]]["M11_1", "M11_2"] <- corr.x[[1]][[1]]["M11_2", "M11_1"] <-
  corr.x[[1]][[1]]["M12_1", "M12_2"] <- corr.x[[1]][[1]]["M12_2", "M12_1"] <-
  corr.x[[1]][[1]]["M12_1", "M12_3"] <- corr.x[[1]][[1]]["M12_3", "M12_1"] <-
  corr.x[[1]][[1]]["M12_2", "M12_3"] <- corr.x[[1]][[1]]["M12_3", "M12_2"] <- 0
# since C1 same across outcomes, set correlation the same
corr.x[[1]][[2]][, same] <- corr.x[[1]][[1]][, same]

corr.x[[2]] <- list(t(corr.x[[1]][[2]]), matrix(0.35, 6, 6), NULL)
diag(corr.x[[2]][[2]]) <- 1
rownames(corr.x[[2]][[2]]) <- colnames(corr.x[[2]][[2]]) <-
  colnames(corr.x[[1]][[2]])
# set correlations between comp of same mixture variable to 0
corr.x[[2]][[2]]["M21_1", "M21_2"] <- corr.x[[2]][[2]]["M21_2", "M21_1"] <-
  corr.x[[2]][[2]]["M22_1", "M22_2"] <- corr.x[[2]][[2]]["M22_2", "M22_1"] <-
  corr.x[[2]][[2]]["M22_1", "M22_3"] <- corr.x[[2]][[2]]["M22_3", "M22_1"] <-
  corr.x[[2]][[2]]["M22_2", "M22_3"] <- corr.x[[2]][[2]]["M22_3", "M22_2"] <- 0

# since C1 same across outcomes, set correlation the same
corr.x[[2]][[2]][, same] <- corr.x[[1]][[2]][same, ]
corr.x[[2]][[2]][same, ] <- corr.x[[2]][[2]][, same]

corr.yx <- list(matrix(0.3, 1, 6), matrix(0.4, 1, 6), NULL)
corr.yx2 <- list(matrix(1, 1, 3), matrix(1, 1, 3), NULL)

for (i in 1:M) {
  if (is.null(corr.yx[[i]])) next
  corr.yx2[[i]][1, 1:(length(skews[[i]]) - 1)] <-
    corr.yx[[i]][1, 1:(length(skews[[i]]) - 1)]
  k.cmix <- c(0, cumsum(lengths(mix_pis[[i]])))
  for (j in 1:length(mix_pis[[i]])) {
    corr.yx2[[i]][1, (length(skews[[i]]) - 1) + j] <-
      rho_M1Y(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
      mix_sigmas[[i]][[j]],
      corr.yx[[i]][1, (length(skews[[i]]) + k.cmix[j]):(length(skews[[i]]) -
        1 + k.cmix[j + 1])])
  }
}

method <- "Polynomial"
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
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
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, E_mix = NULL, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
Y3 <- calc_corr_y(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus, mix_sigmas,
  error_type)
YE3 <- calc_corr_ye(N3$betas, corr.x, corr.e, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)
YX3 <- calc_corr_yx(N3$betas, corr.x, vars, mix_pis, mix_mus,
  mix_sigmas, error_type)

test_that("works for Polynomial method: 1 continuous, 2 mixture,
          error non_mix, Y3 has no X terms", {
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
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y[1, 2],
    0.6219752, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE[1, 2],
    0.3302623, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX[[1]][1, 1],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(Y3[1, 2],
    0.4601323, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YE3[1, 2],
    0.3755321, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(YX3[[1]][1, 1],
    0.3, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
fifths <- sixths <- Six <- mix_fifths <- mix_sixths <- mix_Six <- list()
N <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N2 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx, corr.e,
  seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)
N3 <- nonnormsys(n, M, method, error_type, means, vars, skews, skurts,
  fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
  mix_fifths, mix_sixths, mix_Six, same.var, betas.0, corr.x, corr.yx2,
  corr.e, seed, use.nearPD = FALSE)
S3 <- summary_sys(N3$Y, N3$E, E_mix = NULL, N3$X, N3$X_all, M, method,
  means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
  mix_skews, mix_skurts, mix_fifths, mix_sixths, corr.x = corr.x,
  corr.e = corr.e)

test_that("works for Fleishman method: 1 continuous, 2 mixture,
          error non_mix, Y3 has no X terms", {
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
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S3$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

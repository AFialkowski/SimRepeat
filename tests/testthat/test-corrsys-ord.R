skip_on_cran()
library("SimRepeat")
context("Simulate correlated system using method 1: ordinal only")

options(scipen = 999)
tol <- 1e-5

seed <- 276
n <- 10000
M <- 3
Time <- 1:M

# Error terms have a beta(4, 1.5) distribution with an AR(1, p = 0.4)
# correlation structure
error_type <- "non_mix"
B <- calc_theory("Beta", c(4, 1.5))
corr.e <- matrix(c(1, 0.4, 0.4^2,
                   0.4, 1, 0.4,
                   0.4^2, 0.4, 1), M, M, byrow = TRUE)
means <- lapply(seq_len(M), function(x) B[1])
vars <- lapply(seq_len(M), function(x) B[2]^2)
skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

mix_pis <- mix_mus <- mix_sigmas <- mix_skews <- mix_skurts <- mix_fifths <-
  mix_sixths <- mix_Six <- list()

# Ordinal variables
marginal <- lapply(seq_len(M),
  function(x) list(0.4, c(0.2, 0.5), c(0.1, 0.4, 0.7)))
support <- lapply(seq_len(M), function(x) list(c(0, 1), NULL, NULL))

lam <- p_zip <- size <- prob <- mu <- p_zinb <- list()
pois_eps <- nb_eps <- list()

same.var <- 1

K <- 3
corr.x <- list()
corr.x[[1]] <- corr.x[[2]] <- corr.x[[3]] <- list()
corr.x[[1]][[1]] <- matrix(0.3, K, K)
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[2]] <- matrix(0.35, K, K)
corr.x[[1]][[3]] <- matrix(0.4, K, K)
colnames(corr.x[[1]][[1]]) <- rownames(corr.x[[1]][[1]]) <-
  rownames(corr.x[[1]][[2]]) <- rownames(corr.x[[1]][[3]]) <-
  paste("X1_", 1:ncol(corr.x[[1]][[1]]), sep = "")
colnames(corr.x[[1]][[2]]) <- paste("X2_", 1:ncol(corr.x[[1]][[2]]), sep = "")
colnames(corr.x[[1]][[3]]) <- paste("X3_", 1:ncol(corr.x[[1]][[3]]), sep = "")
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]][[1]] <- t(corr.x[[1]][[2]])
corr.x[[2]][[2]] <- matrix(0.4, K, K)
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[3]] <- matrix(0.45, K, K)
colnames(corr.x[[2]][[2]]) <- rownames(corr.x[[2]][[2]]) <-
  rownames(corr.x[[2]][[3]]) <- colnames(corr.x[[1]][[2]])
colnames(corr.x[[2]][[3]]) <- colnames(corr.x[[1]][[3]])
corr.x[[2]][[2]][same.var, ] <- corr.x[[1]][[2]][same.var, ]
corr.x[[2]][[2]][, same.var] <- t(corr.x[[2]][[2]][same.var, ])
corr.x[[2]][[3]][, same.var] <- t(corr.x[[1]][[2]][same.var, ])
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]

corr.x[[3]][[1]] <- t(corr.x[[1]][[3]])
corr.x[[3]][[2]] <- t(corr.x[[2]][[3]])
corr.x[[3]][[3]] <- matrix(0.5, K, K)
diag(corr.x[[3]][[3]]) <- 1
colnames(corr.x[[3]][[3]]) <- rownames(corr.x[[3]][[3]]) <-
  colnames(corr.x[[1]][[3]])
corr.x[[3]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[3]][[3]][, same.var] <- t(corr.x[[3]][[3]][same.var, ])

subj.var <- matrix(c(1, 2, 1, 3, 2, 2, 2, 3, 3, 2, 3, 3), 6, 2, byrow = TRUE)
int.var <- tint.var <- NULL
betas.0 <- 0
betas <- list(seq(0.5, 0.5 + (K - 1) * 0.25, 0.25))
betas.subj <- list(seq(0.5, 0.5 + (K - 2) * 0.1, 0.1))
betas.int <- list()
betas.t <- 1
betas.tint <- list(0.25)

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 3 ordinal, same.var, subj.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 3 ordinal, same.var, subj.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

int.var <- matrix(c(1, 2, 3, 2, 2, 3, 3, 2, 3), 3, 3, byrow = TRUE)
betas.int <- list(1)
betas.subj <- list(seq(0.5, 0.5 + (K - 1) * 0.1, 0.1))

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 3 ordinal, same.var, subj.var, int.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 3 ordinal, same.var, subj.var, int.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

tint.var = matrix(c(1, 2, 1, 3, 1, 4, 2, 2, 2, 3, 2, 4, 3, 2, 3, 3, 3, 4), 9,
  2, byrow = TRUE)
betas.tint <- list(c(0.25, 0.5, 0.75))

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 3 ordinal, same.var, subj.var, int.var,
          tint.var, error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 3 ordinal, same.var, subj.var, int.var,
          tint.var, error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- skurts <- fifths <- sixths <- Six <- list()
L <- calc_theory("Logistic", c(0, 1))
C <- calc_theory("Chisq", 4)
B <- calc_theory("Beta", c(4, 1.5))

mix_pis <- lapply(seq_len(M), function(x) list(c(0.3, 0.2, 0.5)))
mix_mus <- lapply(seq_len(M), function(x) list(c(L[1], C[1], B[1])))
mix_sigmas <- lapply(seq_len(M), function(x) list(c(L[2], C[2], B[2])))
mix_skews <- lapply(seq_len(M), function(x) list(c(L[3], C[3], B[3])))
mix_skurts <- lapply(seq_len(M), function(x) list(c(L[4], C[4], B[4])))
mix_fifths <- lapply(seq_len(M), function(x) list(c(L[5], C[5], B[5])))
mix_sixths <- lapply(seq_len(M), function(x) list(c(L[6], C[6], B[6])))
mix_Six <- lapply(seq_len(M), function(x) list(1.75, NULL, 0.03))
Mstcum <- calc_mixmoments(mix_pis[[1]][[1]], mix_mus[[1]][[1]],
  mix_sigmas[[1]][[1]], mix_skews[[1]][[1]], mix_skurts[[1]][[1]],
  mix_fifths[[1]][[1]], mix_sixths[[1]][[1]])

means <- lapply(seq_len(M), function(x) Mstcum[1])
vars <- lapply(seq_len(M), function(x) Mstcum[2]^2)

corr.e <- matrix(0.4, 9, 9)
colnames(corr.e) <- rownames(corr.e) <- c("E1_1", "E1_2", "E1_3", "E2_1",
  "E2_2", "E2_3", "E3_1", "E3_2", "E3_3")
corr.e[c("E1_1", "E1_2", "E1_3"), c("E1_1", "E1_2", "E1_3")] <- diag(3)
corr.e[c("E2_1", "E2_2", "E2_3"), c("E2_1", "E2_2", "E2_3")] <- diag(3)
corr.e[c("E3_1", "E3_2", "E3_3"), c("E3_1", "E3_2", "E3_3")] <- diag(3)
error_type = "mix"

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 3 ordinal, same.var, subj.var, int.var,
          tint.var, error mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][3, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][3, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 3 ordinal, same.var, subj.var, int.var,
          tint.var, error mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][3, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][3, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

# Error terms have a beta(4, 1.5) distribution with an AR(1, p = 0.4)
# correlation structure
error_type <- "non_mix"
B <- calc_theory("Beta", c(4, 1.5))
corr.e <- matrix(c(1, 0.4, 0.4^2,
                   0.4, 1, 0.4,
                   0.4^2, 0.4, 1), M, M, byrow = TRUE)
means <- lapply(seq_len(M), function(x) B[1])
vars <- lapply(seq_len(M), function(x) B[2]^2)
skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

mix_pis <- mix_mus <- mix_sigmas <- mix_skews <- mix_skurts <- mix_fifths <-
  mix_sixths <- mix_Six <- list()

# Ordinal variables
marginal <- lapply(seq_len(M), function(x) list(0.4))
support <- lapply(seq_len(M), function(x) list(c(0, 1)))

corr.x <- lapply(seq_len(M), function(x) list(matrix(1, 1, 1), matrix(1, 1, 1),
  matrix(1, 1, 1)))

same.var <- 1
subj.var <- int.var <- tint.var <- NULL
betas <- list(0.5)
betas.subj <- betas.int <- list()
betas.tint <- list(0.25)

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 1 ordinal all M, same.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 1 ordinal all M, same.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

marginal <- list(NULL, list(0.4), list(0.4))
support <- list(NULL, list(c(0, 1)), list(c(0, 1)))

corr.x <- list(NULL, list(NULL, matrix(1, 1, 1), matrix(1, 1, 1)),
  list(NULL, matrix(1, 1, 1), matrix(1, 1, 1)))

same.var <- matrix(c(2, 1, 3, 1), 1, 4, byrow = TRUE)
subj.var <- int.var <- tint.var <- NULL
betas <- list(NULL, 0.5, 0.5)
betas.subj <- betas.int <- list()
betas.tint <- list(NULL, 0.25, 0.25)

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: Y1 no X, same.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: Y1 no X, same.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

marginal <- list(list(0.4), NULL, list(0.4))
support <- list(list(c(0, 1)), NULL, list(c(0, 1)))

corr.x <- list(list(matrix(1, 1, 1), NULL, matrix(1, 1, 1)), NULL,
  list(matrix(1, 1, 1), NULL, matrix(1, 1, 1)))

same.var <- matrix(c(1, 1, 3, 1), 1, 4, byrow = TRUE)
subj.var <- int.var <- tint.var <- NULL
betas <- list(0.5, NULL, 0.5)
betas.subj <- betas.int <- list()
betas.tint <- list(0.25, NULL, 0.25)

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: Y2 no X, same.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: Y2 no X, same.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

marginal <- list(list(0.4), list(0.4), NULL)
support <- list(list(c(0, 1)), list(c(0, 1)), NULL)

corr.x <- list(list(matrix(1, 1, 1), matrix(1, 1, 1), NULL),
  list(matrix(1, 1, 1), matrix(1, 1, 1), NULL), NULL)

same.var <- matrix(c(1, 1, 2, 1), 1, 4, byrow = TRUE)
subj.var <- int.var <- tint.var <- NULL
betas <- list(0.5, 0.5, NULL)
betas.subj <- betas.int <- list()
betas.tint <- list(0.25, 0.25, NULL)

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: Y3 no X, same.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: Y3 no X, same.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

M <- 2
Time <- 1:M
means <- lapply(seq_len(M), function(x) B[1])
vars <- lapply(seq_len(M), function(x) B[2]^2)
skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

marginal <- list(list(0.4), NULL)
support <- list(list(c(0, 1)), NULL)

corr.x <- list(list(matrix(1, 1, 1), NULL), NULL)

same.var <- NULL
subj.var <- int.var <- tint.var <- NULL
betas <- list(0.5, NULL)
betas.subj <- betas.int <- list()
betas.tint <- list(0.25, NULL)
corr.e <- corr.e[1:2, 1:2]

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 1 X for Y1 only,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 1 X for Y1 only,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

###############################
## Random effects
###############################
seed <- 276
n <- 10000
M <- 3
Time <- 1:M

error_type <- "non_mix"
B <- calc_theory("Beta", c(4, 1.5))
corr.e <- matrix(c(1, 0.4, 0.4^2,
                   0.4, 1, 0.4,
                   0.4^2, 0.4, 1), M, M, byrow = TRUE)
means <- lapply(seq_len(M), function(x) B[1])
vars <- lapply(seq_len(M), function(x) B[2]^2)
skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

mix_pis <- mix_mus <- mix_sigmas <- mix_skews <- mix_skurts <- mix_fifths <-
  mix_sixths <- mix_Six <- list()

# Ordinal variables
marginal <- lapply(seq_len(M),
  function(x) list(0.4, c(0.2, 0.5), c(0.1, 0.4, 0.7)))
support <- lapply(seq_len(M), function(x) list(c(0, 1), NULL, NULL))

lam <- p_zip <- size <- prob <- mu <- p_zinb <- list()
pois_eps <- nb_eps <- list()

same.var <- 1

K <- 3
corr.x <- list()
corr.x[[1]] <- corr.x[[2]] <- corr.x[[3]] <- list()
corr.x[[1]][[1]] <- matrix(0.3, K, K)
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[2]] <- matrix(0.35, K, K)
corr.x[[1]][[3]] <- matrix(0.4, K, K)
colnames(corr.x[[1]][[1]]) <- rownames(corr.x[[1]][[1]]) <-
  rownames(corr.x[[1]][[2]]) <- rownames(corr.x[[1]][[3]]) <-
  paste("X1_", 1:ncol(corr.x[[1]][[1]]), sep = "")
colnames(corr.x[[1]][[2]]) <- paste("X2_", 1:ncol(corr.x[[1]][[2]]), sep = "")
colnames(corr.x[[1]][[3]]) <- paste("X3_", 1:ncol(corr.x[[1]][[3]]), sep = "")
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]][[1]] <- t(corr.x[[1]][[2]])
corr.x[[2]][[2]] <- matrix(0.4, K, K)
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[3]] <- matrix(0.45, K, K)
colnames(corr.x[[2]][[2]]) <- rownames(corr.x[[2]][[2]]) <-
  rownames(corr.x[[2]][[3]]) <- colnames(corr.x[[1]][[2]])
colnames(corr.x[[2]][[3]]) <- colnames(corr.x[[1]][[3]])
corr.x[[2]][[2]][same.var, ] <- corr.x[[1]][[2]][same.var, ]
corr.x[[2]][[2]][, same.var] <- t(corr.x[[2]][[2]][same.var, ])
corr.x[[2]][[3]][, same.var] <- t(corr.x[[1]][[2]][same.var, ])
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]

corr.x[[3]][[1]] <- t(corr.x[[1]][[3]])
corr.x[[3]][[2]] <- t(corr.x[[2]][[3]])
corr.x[[3]][[3]] <- matrix(0.5, K, K)
diag(corr.x[[3]][[3]]) <- 1
colnames(corr.x[[3]][[3]]) <- rownames(corr.x[[3]][[3]]) <-
  colnames(corr.x[[1]][[3]])
corr.x[[3]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[3]][[3]][, same.var] <- t(corr.x[[3]][[3]][same.var, ])

subj.var <- matrix(c(1, 2, 1, 3, 2, 2, 2, 3, 3, 2, 3, 3), 6, 2, byrow = TRUE)
int.var <- tint.var <- NULL
betas.0 <- 0
betas <- list(seq(0.5, 0.5 + (K - 1) * 0.25, 0.25))
betas.subj <- list(seq(0.5, 0.5 + (K - 2) * 0.1, 0.1))
betas.int <- list()
betas.t <- 1
betas.tint <- list(0.25)

# 1) rand.int = none, rand.tsl = none, non_mix rand.var 2 per Y
rand.int <- "none"
rand.tsl <- "none"
rand.var <- matrix(c(1, 1, 1, 2, 2, 1, 2, 2, 3, 1, 3, 2), 6, 2, byrow = TRUE)

Stcum5 <- calc_theory("Logistic", c(0, 1))
Stcum6 <- calc_theory("t", 10)

rmeans <- lapply(seq_len(M), function(x) c(Stcum5[1], Stcum6[1]))
rvars <- lapply(seq_len(M), function(x) c(Stcum5[2]^2, Stcum6[2]^2))
rskews <- lapply(seq_len(M), function(x) c(Stcum5[3], Stcum6[3]))
rskurts <- lapply(seq_len(M), function(x) c(Stcum5[4], Stcum6[4]))
rfifths <- lapply(seq_len(M), function(x) c(Stcum5[5], Stcum6[5]))
rsixths <- lapply(seq_len(M), function(x) c(Stcum5[6], Stcum6[6]))
rSix <- list(list(1.75, NULL), list(1.75, NULL), list(1.75, NULL))

means <- append(means, rmeans)
vars <- append(vars, rvars)
skews <- append(skews, rskews)
skurts <- append(skurts, rskurts)
fifths <- append(fifths, rfifths)
sixths <- append(sixths, rsixths)
Six <- append(Six, rSix)

corr.u <- list(list(matrix(c(1, 0.4, 0.4, 1), 2, 2), matrix(0.4, 2, 2),
                    matrix(0.4, 2, 2)),
               list(matrix(0.4, 2, 2), matrix(c(1, 0.4, 0.4, 1), 2, 2),
                    matrix(0.4, 2, 2)),
               list(matrix(0.4, 2, 2), matrix(0.4, 2, 2),
                    matrix(c(1, 0.4, 0.4, 1), 2, 2)))

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N$U, N$U_all, rand.int, rand.tsl,
  corr.u, N$rmeans2, N$rvars2)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE,
  errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Polynomial method: 3 ordinal, same.var, subj.var,
          error non_mix; rand.int = none, rand.tsl = none,
          non_mix rand.var 2 per Y", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N$U, N$U_all, rand.int, rand.tsl,
  corr.u, N$rmeans2, N$rvars2)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE,
  errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Fleishman method: 3 ordinal, same.var, subj.var,
          error non_mix; rand.int = none, rand.tsl = none,
          non_mix rand.var 2 per Y", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

# 2) rand.int = none, non_mix, mix; rand.tsl = none, non_mix, mix;
# non_mix rand.var 1 per Y

rand.int <- c("none", "non_mix", "mix")
rand.tsl <- c("none", "non_mix", "mix")
rand.var <- matrix(c(1, 1, 2, 1, 3, 1), 3, 2, byrow = TRUE)

mix_pis <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, NULL, list(c(0.4, 0.6), c(0.3, 0.2, 0.5))))
mix_mus <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, NULL, list(c(-2, 2), c(L[1], C[1], B[1]))))
mix_sigmas <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, NULL, list(c(1, 1), c(L[2], C[2], B[2]))))
mix_skews <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, NULL, list(c(0, 0), c(L[3], C[3], B[3]))))
mix_skurts <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, NULL, list(c(0, 0), c(L[4], C[4], B[4]))))
mix_fifths <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, NULL, list(c(0, 0), c(L[5], C[5], B[5]))))
mix_sixths <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, NULL, list(c(0, 0), c(L[6], C[6], B[6]))))
mix_Six <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, NULL, list(NULL, NULL, 1.75, NULL, 0.03)))
Nstcum <- calc_mixmoments(mix_pis[[6]][[1]], mix_mus[[6]][[1]],
  mix_sigmas[[6]][[1]], mix_skews[[6]][[1]], mix_skurts[[6]][[1]],
  mix_fifths[[6]][[1]], mix_sixths[[6]][[1]])
Mstcum <- calc_mixmoments(mix_pis[[6]][[2]], mix_mus[[6]][[2]],
  mix_sigmas[[6]][[2]], mix_skews[[6]][[2]], mix_skurts[[6]][[2]],
  mix_fifths[[6]][[2]], mix_sixths[[6]][[2]])

means <- lapply(seq_len(M), function(x) B[1])
vars <- lapply(seq_len(M), function(x) B[2]^2)
skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

skews <- append(skews, list(Stcum5[3], c(Stcum5[3], Stcum6[3], Stcum5[3]),
                            Stcum5[3]))
skurts <- append(skurts, list(Stcum5[4], c(Stcum5[4], Stcum6[4], Stcum5[4]),
                              Stcum5[4]))
fifths <- append(fifths, list(Stcum5[5], c(Stcum5[5], Stcum6[5], Stcum5[5]),
                              Stcum5[5]))
sixths <- append(sixths, list(Stcum5[6], c(Stcum5[6], Stcum6[6], Stcum5[6]),
                              Stcum5[6]))
Six <- append(Six, list(list(1.75), list(1.75, NULL, 1.75), list(1.75)))

means <- append(means, list(Stcum5[1], c(Stcum5[1], Stcum6[1], Stcum5[1]),
  c(Nstcum[1], Mstcum[1], Stcum5[1])))
vars <- append(vars, list(Stcum5[2]^2,
  c(Stcum5[2]^2, Stcum6[2]^2, Stcum5[2]^2),
  c(Nstcum[2]^2, Mstcum[2]^2, Stcum5[2]^2)))

corr.u <- list(list(matrix(0.4, 1, 1), matrix(0.4, 1, 3), matrix(0.4, 1, 6)),
  list(matrix(0.4, 3, 1), matrix(0.4, 3, 3), matrix(0.4, 3, 6)),
  list(matrix(0.4, 6, 1), matrix(0.4, 6, 3), matrix(0.4, 6, 6)))

diag(corr.u[[1]][[1]]) <- diag(corr.u[[2]][[2]]) <-
  diag(corr.u[[3]][[3]]) <- 1

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N$U, N$U_all, rand.int, rand.tsl,
  corr.u, N$rmeans2, N$rvars2)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE,
  errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Polynomial method: 3 ordinal, same.var, subj.var,
          error non_mix; rand.int = none, non_mix, mix;
          rand.tsl = none, non_mix, mix; non_mix rand.var 1 per Y", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N$U, N$U_all, rand.int, rand.tsl,
  corr.u, N$rmeans2, N$rvars2)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE,
  errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Fleishman method: 3 ordinal, same.var, subj.var,
          error non_mix; rand.int = non_mix, rand.tsl = non_mix,
          non_mix rand.var 2 per Y", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

# 3) rand.int = none, mix, non_mix; rand.tsl = none, non_mix, mix;
# mix, non_mix rand.var 2 per Y2, Y3

rand.int <- c("none", "mix", "non_mix")
rand.tsl <- c("none", "non_mix", "mix")
rand.var <- matrix(c(2, 2, 2, 1, 3, 2, 3, 1), 4, 2, byrow = TRUE)

mix_pis <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, list(c(0.4, 0.6), c(0.3, 0.2, 0.5)),
       list(c(0.4, 0.6), c(0.3, 0.2, 0.5))))
mix_mus <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, list(c(-2, 2), c(L[1], C[1], B[1])),
       list(c(-2, 2), c(L[1], C[1], B[1]))))
mix_sigmas <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, list(c(1, 1), c(L[2], C[2], B[2])),
       list(c(1, 1), c(L[2], C[2], B[2]))))
mix_skews <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, list(c(0, 0), c(L[3], C[3], B[3])),
       list(c(0, 0), c(L[3], C[3], B[3]))))
mix_skurts <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, list(c(0, 0), c(L[4], C[4], B[4])),
       list(c(0, 0), c(L[4], C[4], B[4]))))
mix_fifths <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, list(c(0, 0), c(L[5], C[5], B[5])),
       list(c(0, 0), c(L[5], C[5], B[5]))))
mix_sixths <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, list(c(0, 0), c(L[6], C[6], B[6])),
       list(c(0, 0), c(L[6], C[6], B[6]))))
mix_Six <- append(lapply(seq_len(M), function(x) NULL),
  list(NULL, list(NULL, NULL, 1.75, NULL, 0.03),
  list(NULL, NULL, 1.75, NULL, 0.03)))
Nstcum <- calc_mixmoments(mix_pis[[6]][[1]], mix_mus[[6]][[1]],
  mix_sigmas[[6]][[1]], mix_skews[[6]][[1]], mix_skurts[[6]][[1]],
  mix_fifths[[6]][[1]], mix_sixths[[6]][[1]])
Mstcum <- calc_mixmoments(mix_pis[[6]][[2]], mix_mus[[6]][[2]],
  mix_sigmas[[6]][[2]], mix_skews[[6]][[2]], mix_skurts[[6]][[2]],
  mix_fifths[[6]][[2]], mix_sixths[[6]][[2]])

means <- lapply(seq_len(M), function(x) B[1])
vars <- lapply(seq_len(M), function(x) B[2]^2)
skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

skews <- append(skews, list(NULL, c(Stcum5[3], Stcum6[3]),
                            c(Stcum5[3], Stcum6[3])))
skurts <- append(skurts, list(NULL, c(Stcum5[4], Stcum6[4]),
                              c(Stcum5[4], Stcum6[4])))
fifths <- append(fifths, list(NULL, c(Stcum5[5], Stcum6[5]),
                              c(Stcum5[5], Stcum6[5])))
sixths <- append(sixths, list(NULL, c(Stcum5[6], Stcum6[6]),
                              c(Stcum5[6], Stcum6[6])))
Six <- append(Six, list(NULL, list(1.75, NULL), list(1.75, NULL)))

means <- append(means, list(NULL,
  c(Nstcum[1], Stcum5[1], Stcum6[1], Mstcum[1]),
  c(Stcum5[1], Nstcum[1], Stcum6[1], Mstcum[1])))
vars <- append(vars, list(NULL,
  c(Nstcum[2]^2, Stcum5[2]^2, Stcum6[2]^2, Mstcum[2]^2),
  c(Stcum5[2]^2, Nstcum[2]^2, Stcum6[2]^2, Mstcum[2]^2)))

corr.u <- list(NULL,
  list(NULL, matrix(0.4, 7, 7), matrix(0.4, 7, 7)),
  list(NULL, matrix(0.4, 7, 7), matrix(0.4, 7, 7)))
diag(corr.u[[2]][[2]]) <- diag(corr.u[[3]][[3]]) <- 1

method <- "Polynomial"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N$U, N$U_all, rand.int, rand.tsl,
  corr.u, N$rmeans2, N$rvars2)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE,
  errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Polynomial method: 3 ordinal, same.var, subj.var,
          error non_mix; rand.int = none, mix, non_mix;
          rand.tsl = none, non_mix, mix;
          mix, non_mix rand.var 2 per Y2, Y3", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N$U, N$U_all, rand.int, rand.tsl,
  corr.u, N$rmeans2, N$rvars2)
N2 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE,
  errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Fleishman method: 3 ordinal, same.var, subj.var,
          error non_mix; rand.int = none, mix, non_mix;
          rand.tsl = none, non_mix, mix;
          mix, non_mix rand.var 2 per Y2, Y3", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][1, "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

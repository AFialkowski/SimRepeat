skip_on_cran()
library("SimRepeat")
context("Simulate correlated system using method 1: continuous only")

options(scipen = 999)
tol <- 1e-5

seed <- 276
n <- 10000
M <- 3
Time <- 1:M

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

means <- lapply(seq_len(M),
  function(x) c(L[1], Nstcum[1], Mstcum[1], B[1]))
vars <- lapply(seq_len(M),
  function(x) c(L[2]^2, Nstcum[2]^2, Mstcum[2]^2, B[2]^2))

same.var <- c(2, 3)

corr.x <- list()
corr.x[[1]] <- list(matrix(0.1, 6, 6), matrix(0.2, 6, 6), matrix(0.3, 6, 6))
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[1]][2:3, 2:3] <- diag(2)
corr.x[[1]][[1]][4:6, 4:6] <- diag(3)
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]] <- list(t(corr.x[[1]][[2]]), matrix(0.35, 6, 6),
                    matrix(0.4, 6, 6))
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[2]][4:6, 4:6] <- diag(3)
corr.x[[2]][[2]][, same.var] <- corr.x[[2]][[3]][, same.var] <-
  t(corr.x[[1]][[2]][same.var, ])
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[2]][[2]][same.var, ] <- t(corr.x[[2]][[2]][, same.var])

corr.x[[3]] <- list(t(corr.x[[1]][[3]]), t(corr.x[[2]][[3]]),
                    matrix(0.5, 6, 6))
diag(corr.x[[3]][[3]]) <- 1
corr.x[[3]][[3]][4:6, 4:6] <- diag(3)
corr.x[[3]][[3]][, same.var] <- t(corr.x[[1]][[3]][same.var, ])
corr.x[[3]][[3]][same.var, ] <- t(corr.x[[3]][[3]][, same.var])

# Error terms have a beta(4, 1.5) distribution with an AR(1, p = 0.4)
# correlation structure
error_type <- "non_mix"
corr.e <- matrix(c(1, 0.4, 0.4^2,
                   0.4, 1, 0.4,
                   0.4^2, 0.4, 1), M, M, byrow = TRUE)

marginal <- support <- list()
lam <- p_zip <- size <- prob <- mu <- p_zinb <- list()
pois_eps <- nb_eps <- list()

K <- 3
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

test_that("works for Polynomial method: 1 cont, 2 mix, same.var, subj.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Fleishman method: 1 cont, 2 mix, same.var, subj.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Polynomial method: 1 cont, 2 mix, same.var, subj.var, int.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Fleishman method: 1 cont, 2 mix, same.var, subj.var, int.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Polynomial method: 1 cont, 2 mix, same.var, subj.var, int.var,
          tint.var, error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Fleishman method: 1 cont, 2 mix, same.var, subj.var, int.var,
          tint.var, error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

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

means <- lapply(seq_len(M),
                function(x) c(L[1], B[1], Nstcum[1], Mstcum[1]))
vars <- lapply(seq_len(M),
               function(x) c(L[2]^2, B[2]^2, Nstcum[2]^2, Mstcum[2]^2))

same.var <- 1

corr.x <- list()
corr.x[[1]] <- list(matrix(0.1, 4, 4), matrix(0.2, 4, 4), matrix(0.3, 4, 4))
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[1]][3:4, 3:4] <- diag(2)
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]] <- list(t(corr.x[[1]][[2]]), matrix(0.35, 4, 4),
                    matrix(0.4, 4, 4))
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[2]][3:4, 3:4] <- diag(2)
corr.x[[2]][[2]][, same.var] <- corr.x[[2]][[3]][, same.var] <-
  t(corr.x[[1]][[2]][same.var, ])
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[2]][[2]][same.var, ] <- t(corr.x[[2]][[2]][, same.var])

corr.x[[3]] <- list(t(corr.x[[1]][[3]]), t(corr.x[[2]][[3]]),
                    matrix(0.5, 4, 4))
diag(corr.x[[3]][[3]]) <- 1
corr.x[[3]][[3]][3:4, 3:4] <- diag(2)
corr.x[[3]][[3]][, same.var] <- t(corr.x[[1]][[3]][same.var, ])
corr.x[[3]][[3]][same.var, ] <- t(corr.x[[3]][[3]][, same.var])

K <- 3
subj.var <- matrix(c(1, 2, 1, 3, 2, 2, 2, 3, 3, 2, 3, 3), 6, 2, byrow = TRUE)
int.var <- tint.var <- NULL
betas.0 <- 0
betas <- list(seq(0.5, 0.5 + (K - 1) * 0.25, 0.25))
betas.subj <- list(seq(0.5, 0.5 + (K - 2) * 0.1, 0.1))
betas.int <- list()
betas.t <- 1
betas.tint <- list(0.25)

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

test_that("works for Polynomial method: 1 cont, 2 mix, same.var, subj.var, int.var,
          tint.var, error mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Fleishman method: 1 cont, 2 mix, same.var, subj.var, int.var,
          tint.var, error mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

M <- 3
Time <- 1:M
skews <- lapply(seq_len(M), function(x) c(L[3], B[3]))
skurts <- lapply(seq_len(M), function(x) c(L[4], B[4]))
fifths <- lapply(seq_len(M), function(x) c(L[5], B[5]))
sixths <- lapply(seq_len(M), function(x) c(L[6], B[6]))
Six <- lapply(seq_len(M), function(x) list(1.75, 0.03))

mix_pis <- mix_mus <- mix_sigmas <- mix_skews<- mix_skurts <- mix_fifths <-
  mix_sixths <- mix_six <- list()

means <- lapply(seq_len(M), function(x) c(L[1], B[1]))
vars <- lapply(seq_len(M), function(x) c(L[2]^2, B[2]^2))

corr.x <- list(list(matrix(1, 1, 1), matrix(0.4, 1, 1), matrix(0.4, 1, 1)),
               list(matrix(0.4, 1, 1), matrix(1, 1, 1), matrix(0.4, 1, 1)),
               list(matrix(0.4, 1, 1), matrix(0.4, 1, 1), matrix(1, 1, 1)))

error_type <- "non_mix"
corr.e <- matrix(c(1, 0.4, 0.4^2,
                   0.4, 1, 0.4,
                   0.4^2, 0.4, 1), M, M, byrow = TRUE)
same.var <- NULL
subj.var <- matrix(c(1, 1, 3, 1), 2, 2, byrow = TRUE)
int.var <- tint.var <- NULL
betas.0 <- 0
betas <- list(0.5)
betas.subj <- list()
betas.int <- list()
betas.t <- 1
betas.tint <- list(NULL, 0.25, NULL)

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
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.0001)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 1 cont all M, subj.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.0001)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 1 cont all M, subj.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- list(B[3], c(L[3], B[3]), c(L[3], B[3]))
skurts <- list(B[4], c(L[4], B[4]), c(L[4], B[4]))
fifths <- list(B[5], c(L[5], B[5]), c(L[5], B[5]))
sixths <- list(B[6], c(L[6], B[6]), c(L[6], B[6]))
Six <- list(list(0.03), list(1.75, 0.03), list(1.75, 0.03))

means <- list(B[1], c(L[1], B[1]), c(L[1], B[1]))
vars <- list(B[2]^2, c(L[2]^2, B[2]^2), c(L[2]^2, B[2]^2))

corr.x <- list(NULL,
               list(NULL, matrix(1, 1, 1), matrix(0.4, 1, 1)),
               list(NULL, matrix(0.4, 1, 1), matrix(1, 1, 1)))

same.var <- NULL
subj.var <- matrix(c(2, 1), 1, 2, byrow = TRUE)
int.var <- NULL
tint.var <- matrix(c(2, 1), 1, 2, byrow = TRUE)
betas <- list(NULL, 0.5, 0.5)
betas.subj <- list()
betas.int <- list()
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
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.0001)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: Y1 no X, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Fleishman method: Y1 no X, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- list(c(L[3], B[3]), B[3], c(L[3], B[3]))
skurts <- list(c(L[4], B[4]), B[4], c(L[4], B[4]))
fifths <- list(c(L[5], B[5]), B[5], c(L[5], B[5]))
sixths <- list(c(L[6], B[6]), B[6], c(L[6], B[6]))
Six <- list(list(1.75, 0.03), list(0.03), list(1.75, 0.03))

means <- list(c(L[1], B[1]), B[1], c(L[1], B[1]))
vars <- list(c(L[2]^2, B[2]^2), B[2]^2, c(L[2]^2, B[2]^2))

corr.x <- list(list(matrix(1, 1, 1), NULL, matrix(1, 1, 1)), NULL,
               list(matrix(1, 1, 1), NULL, matrix(1, 1, 1)))

same.var <- matrix(c(1, 1, 3, 1), 1, 4)
subj.var <- matrix(c(3, 1), 1, 2, byrow = TRUE)
int.var <- NULL
tint.var <- matrix(c(1, 1), 1, 2, byrow = TRUE)
betas <- list(0.5, NULL, 0.5)
betas.subj <- list()
betas.int <- list()
betas.tint <- list(0.25, NULL, NULL)

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

test_that("works for Polynomial method: Y2 no X, same.var, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
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

test_that("works for Fleishman method: Y2 no X, same.var, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- list(c(L[3], B[3]), c(L[3], B[3]), B[3])
skurts <- list(c(L[4], B[4]), c(L[4], B[4]), B[4])
fifths <- list(c(L[5], B[5]), c(L[5], B[5]), B[5])
sixths <- list(c(L[6], B[6]), c(L[6], B[6]), B[6])
Six <- list(list(1.75, 0.03), list(1.75, 0.03), list(0.03))

means <- list(c(L[1], B[1]), c(L[1], B[1]), B[1])
vars <- list(c(L[2]^2, B[2]^2), c(L[2]^2, B[2]^2), B[2]^2)

corr.x <- list(list(matrix(1, 1, 1), matrix(0.4, 1, 1), NULL),
               list(matrix(0.4, 1, 1), matrix(1, 1, 1), NULL), NULL)

same.var <- NULL
subj.var <- matrix(c(2, 1), 1, 2, byrow = TRUE)
int.var <- NULL
tint.var <- matrix(c(2, 1), 1, 2, byrow = TRUE)
betas <- list(0.5, 0.5, NULL)
betas.subj <- list()
betas.int <- list()
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
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.0001)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: Y3 no X, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Fleishman method: Y3 no X, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

M <- 2
Time <- 1:M
means <- list(c(L[1], B[1]), B[1])
vars <- list(c(L[2]^2, B[2]^2), B[2]^2)
skews <- list(c(L[3], B[3]), B[3])
skurts <- list(c(L[4], B[4]), B[4])
fifths <- list(c(L[5], B[5]), B[5])
sixths <- list(c(L[6], B[6]), B[6])
Six <- list(list(1.75, 0.03), NULL)

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

test_that("works for Polynomial method: 1 X for Y1 only,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
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

test_that("works for Fleishman method: 1 X for Y1 only,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

M <- 3
Time <- 1:M
skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

mix_pis <- lapply(seq_len(M), function(x) list(c(0.3, 0.2, 0.5)))
mix_mus <- lapply(seq_len(M), function(x) list(c(L[1], C[1], B[1])))
mix_sigmas <- lapply(seq_len(M), function(x) list(c(L[2], C[2], B[2])))
mix_skews<- lapply(seq_len(M), function(x) list(c(L[3], C[3], B[3])))
mix_skurts <- lapply(seq_len(M), function(x) list(c(L[4], C[4], B[4])))
mix_fifths <- lapply(seq_len(M), function(x) list(c(L[5], C[5], B[5])))
mix_sixths <- lapply(seq_len(M), function(x) list(c(L[6], C[6], B[6])))
mix_Six <- lapply(seq_len(M), function(x) list(1.75, NULL, 0.03))

means <- lapply(seq_len(M), function(x) c(Mstcum[1], B[1]))
vars <- lapply(seq_len(M), function(x) c(Mstcum[2]^2, B[2]^2))

corr.x <- list(list(diag(3), matrix(0.4, 3, 3), matrix(0.4, 3, 3)),
               list(matrix(0.4, 3, 3), diag(3), matrix(0.4, 3, 3)),
               list(matrix(0.4, 3, 3), matrix(0.4, 3, 3), diag(3)))

error_type <- "non_mix"
corr.e <- matrix(c(1, 0.4, 0.4^2,
                   0.4, 1, 0.4,
                   0.4^2, 0.4, 1), M, M, byrow = TRUE)
same.var <- NULL
subj.var <- matrix(c(1, 1, 3, 1), 2, 2, byrow = TRUE)
int.var <- NULL
tint.var <- matrix(c(1, 1, 3, 1), 2, 2, byrow = TRUE)
betas <- list(0.5)
betas.subj <- list()
betas.int <- list()
betas.tint <- list(0.25, 0.25, 0.25)

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

test_that("works for Polynomial method: 1 mix all M, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Fleishman method: 1 cont all M, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

mix_pis <- list(NULL, list(c(0.3, 0.2, 0.5)), list(c(0.3, 0.2, 0.5)))
mix_mus <- list(NULL, list(c(L[1], C[1], B[1])), list(c(L[1], C[1], B[1])))
mix_sigmas <- list(NULL, list(c(L[2], C[2], B[2])), list(c(L[2], C[2], B[2])))
mix_skews<- list(NULL, list(c(L[3], C[3], B[3])), list(c(L[3], C[3], B[3])))
mix_skurts <- list(NULL, list(c(L[4], C[4], B[4])), list(c(L[4], C[4], B[4])))
mix_fifths <- list(NULL, list(c(L[5], C[5], B[5])), list(c(L[5], C[5], B[5])))
mix_sixths <- list(NULL, list(c(L[6], C[6], B[6])), list(c(L[6], C[6], B[6])))
mix_Six <- list(NULL, list(1.75, NULL, 0.03), list(1.75, NULL, 0.03))

means <- list(B[1], c(Mstcum[1], B[1]), c(Mstcum[1], B[1]))
vars <- list(B[2]^2, c(Mstcum[2]^2, B[2]^2), c(Mstcum[2]^2, B[2]^2))

corr.x <- list(NULL,
               list(NULL, diag(3), diag(3)),
               list(NULL, diag(3), diag(3)))

same.var <- matrix(c(2, 1, 3, 1, 2, 2, 3, 2, 2, 3, 3, 3), 3, 4, byrow = TRUE)
subj.var <- matrix(c(3, 1), 1, 2)
int.var <- NULL
tint.var <- matrix(c(3, 1), 1, 2)
betas <- list(NULL, 0.5, 0.5)
betas.subj <- list()
betas.int <- list()
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

test_that("works for Polynomial method: Y1 no X, same.var, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
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

test_that("works for Fleishman method: Y1 no X, same.var, subj.var, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

mix_pis <- list(list(c(0.3, 0.2, 0.5)), NULL, list(c(0.3, 0.2, 0.5)))
mix_mus <- list(list(c(L[1], C[1], B[1])), NULL, list(c(L[1], C[1], B[1])))
mix_sigmas <- list(list(c(L[2], C[2], B[2])), NULL, list(c(L[2], C[2], B[2])))
mix_skews<- list(list(c(L[3], C[3], B[3])), NULL, list(c(L[3], C[3], B[3])))
mix_skurts <- list(list(c(L[4], C[4], B[4])), NULL, list(c(L[4], C[4], B[4])))
mix_fifths <- list(list(c(L[5], C[5], B[5])), NULL, list(c(L[5], C[5], B[5])))
mix_sixths <- list(list(c(L[6], C[6], B[6])), NULL, list(c(L[6], C[6], B[6])))
mix_Six <- list(list(1.75, NULL, 0.03), NULL, list(1.75, NULL, 0.03))

means <- list(c(Mstcum[1], B[1]), B[1], c(Mstcum[1], B[1]))
vars <- list(c(Mstcum[2]^2, B[2]^2), B[2]^2, c(Mstcum[2]^2, B[2]^2))

corr.x <- list(list(diag(3), NULL, matrix(0.4, 3, 3)), NULL,
               list(matrix(0.4, 3, 3), NULL, diag(3)))

same.var <- NULL
subj.var <- matrix(c(1, 1), 1, 2, byrow = TRUE)
int.var <- NULL
tint.var <- matrix(c(1, 1), 1, 2, byrow = TRUE)
betas <- list(0.5, NULL, 0.5)
betas.subj <- list()
betas.int <- list()
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

test_that("works for Polynomial method: Y2 no X, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
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

test_that("works for Fleishman method: Y2 no X, tint.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

mix_pis <- list(list(c(0.3, 0.2, 0.5)), list(c(0.3, 0.2, 0.5)), NULL)
mix_mus <- list(list(c(L[1], C[1], B[1])), list(c(L[1], C[1], B[1])), NULL)
mix_sigmas <- list(list(c(L[2], C[2], B[2])), list(c(L[2], C[2], B[2])), NULL)
mix_skews<- list(list(c(L[3], C[3], B[3])), list(c(L[3], C[3], B[3])), NULL)
mix_skurts <- list(list(c(L[4], C[4], B[4])), list(c(L[4], C[4], B[4])), NULL)
mix_fifths <- list(list(c(L[5], C[5], B[5])), list(c(L[5], C[5], B[5])), NULL)
mix_sixths <- list(list(c(L[6], C[6], B[6])), list(c(L[6], C[6], B[6])), NULL)
mix_Six <- list(list(1.75, NULL, 0.03), list(1.75, NULL, 0.03), NULL)

means <- list(c(Mstcum[1], B[1]), c(Mstcum[1], B[1]), B[1])
vars <- list(c(Mstcum[2]^2, B[2]^2), c(Mstcum[2]^2, B[2]^2), B[2]^2)

corr.x <- list(list(diag(3), matrix(0.4, 3, 3), NULL),
               list(matrix(0.4, 3, 3), diag(3), NULL), NULL)

same.var <- NULL
subj.var <- matrix(c(2, 1), 1, 2, byrow = TRUE)
int.var <- NULL
tint.var <- NULL
betas <- list(0.5, 0.5, NULL)
betas.subj <- list()
betas.int <- list()
betas.tint <- list(0.25, NULL, NULL)

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

test_that("works for Polynomial method: Y3 no X, subj.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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

test_that("works for Fleishman method: Y3 no X, subj.var,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

M <- 2
Time <- 1:M
mix_pis <- list(NULL, list(c(0.3, 0.2, 0.5)))
mix_mus <- list(NULL, list(c(L[1], C[1], B[1])))
mix_sigmas <- list(NULL, list(c(L[2], C[2], B[2])))
mix_skews<- list(NULL, list(c(L[3], C[3], B[3])))
mix_skurts <- list(NULL, list(c(L[4], C[4], B[4])))
mix_fifths <- list(NULL, list(c(L[5], C[5], B[5])))
mix_sixths <- list(NULL, list(c(L[6], C[6], B[6])))
mix_Six <- list(NULL, list(1.75, NULL, 0.03))

means <- list(B[1], c(Mstcum[1], B[1]))
vars <- list(B[2]^2, c(Mstcum[2]^2, B[2]^2))
skews <- list(B[3], B[3])
skurts <- list(B[4], B[4])
fifths <- list(B[5], B[5])
sixths <- list(B[6], B[6])
Six <- list(list(0.03), list(0.03))
corr.x <- list(NULL, list(NULL, diag(3)))

same.var <- NULL
subj.var <- int.var <- tint.var <- NULL
betas <- list(NULL, 0.5)
betas.subj <- betas.int <- list()
betas.tint <- list()
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

test_that("works for Polynomial method: 1 X for Y1 only,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
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

test_that("works for Fleishman method: 1 X for Y1 only,
          error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

###############################
## Random effects
###############################
seed <- 276
n <- 10000
M <- 3
Time <- 1:M

skews0 <- lapply(seq_len(M), function(x) c(L[3], B[3]))
skurts0 <- lapply(seq_len(M), function(x) c(L[4], B[4]))
fifths0 <- lapply(seq_len(M), function(x) c(L[5], B[5]))
sixths0 <- lapply(seq_len(M), function(x) c(L[6], B[6]))
Six0 <- lapply(seq_len(M), function(x) list(1.75, 0.03))

mix_pis0 <- lapply(seq_len(M), function(x) list(c(0.4, 0.6), c(0.3, 0.2, 0.5)))
mix_mus0 <- lapply(seq_len(M), function(x) list(c(-2, 2), c(L[1], C[1], B[1])))
mix_sigmas0 <- lapply(seq_len(M),
  function(x) list(c(1, 1), c(L[2], C[2], B[2])))
mix_skews0 <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[3], C[3], B[3])))
mix_skurts0 <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[4], C[4], B[4])))
mix_fifths0 <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[5], C[5], B[5])))
mix_sixths0 <- lapply(seq_len(M),
  function(x) list(rep(0, 2), c(L[6], C[6], B[6])))
mix_Six0 <- lapply(seq_len(M), function(x) list(NULL, NULL, 1.75, NULL, 0.03))

means0 <- lapply(seq_len(M),
                function(x) c(L[1], B[1], Nstcum[1], Mstcum[1]))
vars0 <- lapply(seq_len(M),
               function(x) c(L[2]^2, B[2]^2, Nstcum[2]^2, Mstcum[2]^2))

same.var <- 1

corr.x <- list()
corr.x[[1]] <- list(matrix(0.1, 4, 4), matrix(0.2, 4, 4), matrix(0.3, 4, 4))
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[1]][3:4, 3:4] <- diag(2)
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]] <- list(t(corr.x[[1]][[2]]), matrix(0.35, 4, 4),
                    matrix(0.4, 4, 4))
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[2]][3:4, 3:4] <- diag(2)
corr.x[[2]][[2]][, same.var] <- corr.x[[2]][[3]][, same.var] <-
  t(corr.x[[1]][[2]][same.var, ])
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[2]][[2]][same.var, ] <- t(corr.x[[2]][[2]][, same.var])

corr.x[[3]] <- list(t(corr.x[[1]][[3]]), t(corr.x[[2]][[3]]),
                    matrix(0.5, 4, 4))
diag(corr.x[[3]][[3]]) <- 1
corr.x[[3]][[3]][3:4, 3:4] <- diag(2)
corr.x[[3]][[3]][, same.var] <- t(corr.x[[1]][[3]][same.var, ])
corr.x[[3]][[3]][same.var, ] <- t(corr.x[[3]][[3]][, same.var])

K <- 3
subj.var <- matrix(c(1, 2, 1, 3, 2, 2, 2, 3, 3, 2, 3, 3), 6, 2, byrow = TRUE)
int.var <- tint.var <- NULL
betas.0 <- 0
betas <- list(seq(0.5, 0.5 + (K - 1) * 0.25, 0.25))
betas.subj <- list(seq(0.5, 0.5 + (K - 2) * 0.1, 0.1))
betas.int <- list()
betas.t <- 1
betas.tint <- list(0.25)

corr.e <- matrix(0.4, 9, 9)
colnames(corr.e) <- rownames(corr.e) <- c("E1_1", "E1_2", "E1_3", "E2_1",
  "E2_2", "E2_3", "E3_1", "E3_2", "E3_3")
corr.e[c("E1_1", "E1_2", "E1_3"), c("E1_1", "E1_2", "E1_3")] <- diag(3)
corr.e[c("E2_1", "E2_2", "E2_3"), c("E2_1", "E2_2", "E2_3")] <- diag(3)
corr.e[c("E3_1", "E3_2", "E3_3"), c("E3_1", "E3_2", "E3_3")] <- diag(3)
error_type = "mix"

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

means <- append(means0, rmeans)
vars <- append(vars0, rvars)
skews <- append(skews0, rskews)
skurts <- append(skurts0, rskurts)
fifths <- append(fifths0, rfifths)
sixths <- append(sixths0, rsixths)
Six <- append(Six0, rSix)
mix_pis <- mix_pis0
mix_mus <- mix_mus0
mix_sigmas <- mix_sigmas0
mix_skews <- mix_skews0
mix_skurts <- mix_skurts0
mix_fifths <- mix_fifths0
mix_sixths <- mix_sixths0
mix_Six <- mix_Six0
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
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
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
  errorloop = TRUE, epsilon = 0.001)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Polynomial method: 2 cont, 1 mix, same.var, subj.var,
          error mix; rand.int = none, rand.tsl = none,
          non_mix rand.var 2 per Y", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
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
  errorloop = TRUE, epsilon = 0.001)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Fleishman method: 2 cont, 1 mix, same.var, subj.var,
          error mix; rand.int = none, rand.tsl = none,
          non_mix rand.var 2 per Y", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

# 2) rand.int = none, non_mix, mix; rand.tsl = none, non_mix, mix;
# non_mix rand.var 1 per Y

rand.int <- c("none", "non_mix", "mix")
rand.tsl <- c("none", "non_mix", "mix")
rand.var <- matrix(c(1, 1, 2, 1, 3, 1), 3, 2, byrow = TRUE)

rmix_pis <- list(NULL, NULL, list(c(0.4, 0.6), c(0.3, 0.2, 0.5)))
rmix_mus <- list(NULL, NULL, list(c(-2, 2), c(L[1], C[1], B[1])))
rmix_sigmas <- list(NULL, NULL, list(c(1, 1), c(L[2], C[2], B[2])))
rmix_skews <- list(NULL, NULL, list(c(0, 0), c(L[3], C[3], B[3])))
rmix_skurts <- list(NULL, NULL, list(c(0, 0), c(L[4], C[4], B[4])))
rmix_fifths <- list(NULL, NULL, list(c(0, 0), c(L[5], C[5], B[5])))
rmix_sixths <- list(NULL, NULL, list(c(0, 0), c(L[6], C[6], B[6])))
rmix_Six <- list(NULL, NULL, list(NULL, NULL, 1.75, NULL, 0.03))

rskews <- list(Stcum5[3], c(Stcum5[3], Stcum6[3], Stcum5[3]),
                            Stcum5[3])
rskurts <- list(Stcum5[4], c(Stcum5[4], Stcum6[4], Stcum5[4]),
                              Stcum5[4])
rfifths <- list(Stcum5[5], c(Stcum5[5], Stcum6[5], Stcum5[5]),
                              Stcum5[5])
rsixths <- list(Stcum5[6], c(Stcum5[6], Stcum6[6], Stcum5[6]),
                              Stcum5[6])
rSix <- list(list(1.75), list(1.75, NULL, 1.75), list(1.75))

rmeans <-list(Stcum5[1], c(Stcum5[1], Stcum6[1], Stcum5[1]),
  c(Nstcum[1], Mstcum[1], Stcum5[1]))
rvars <- list(Stcum5[2]^2,
  c(Stcum5[2]^2, Stcum6[2]^2, Stcum5[2]^2),
  c(Nstcum[2]^2, Mstcum[2]^2, Stcum5[2]^2))

means <- append(means0, rmeans)
vars <- append(vars0, rvars)
skews <- append(skews0, rskews)
skurts <- append(skurts0, rskurts)
fifths <- append(fifths0, rfifths)
sixths <- append(sixths0, rsixths)
Six <- append(Six0, rSix)
mix_pis <- append(mix_pis0, rmix_pis)
mix_mus <- append(mix_mus0, rmix_mus)
mix_sigmas <- append(mix_sigmas0, rmix_sigmas)
mix_skews <- append(mix_skews0, rmix_skews)
mix_skurts <- append(mix_skurts0, rmix_skurts)
mix_fifths <- append(mix_fifths0, rmix_fifths)
mix_sixths <- append(mix_sixths0, rmix_sixths)
mix_Six <- append(mix_Six0, rmix_Six)

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
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
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
  errorloop = TRUE, epsilon = 0.001)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Polynomial method: 2 cont, 1 mix, same.var, subj.var,
          error mix; rand.int = none, non_mix, mix;
          rand.tsl = none, non_mix, mix; non_mix rand.var 1 per Y", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
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
  errorloop = TRUE, epsilon = 0.001)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Fleishman method: 2 cont, 1 mix, same.var, subj.var,
          error mix; rand.int = non_mix, rand.tsl = non_mix,
          non_mix rand.var 2 per Y", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

# 3) rand.int = none, mix, non_mix; rand.tsl = none, non_mix, mix;
# mix, non_mix rand.var 2 per Y2, Y3

rand.int <- c("none", "mix", "non_mix")
rand.tsl <- c("none", "non_mix", "mix")
rand.var <- matrix(c(2, 2, 2, 1, 3, 2, 3, 1), 4, 2, byrow = TRUE)

rmix_pis <- list(NULL, list(c(0.4, 0.6), c(0.3, 0.2, 0.5)),
       list(c(0.4, 0.6), c(0.3, 0.2, 0.5)))
rmix_mus <- list(NULL, list(c(-2, 2), c(L[1], C[1], B[1])),
       list(c(-2, 2), c(L[1], C[1], B[1])))
rmix_sigmas <- list(NULL, list(c(1, 1), c(L[2], C[2], B[2])),
       list(c(1, 1), c(L[2], C[2], B[2])))
rmix_skews <- list(NULL, list(c(0, 0), c(L[3], C[3], B[3])),
       list(c(0, 0), c(L[3], C[3], B[3])))
rmix_skurts <- list(NULL, list(c(0, 0), c(L[4], C[4], B[4])),
       list(c(0, 0), c(L[4], C[4], B[4])))
rmix_fifths <- list(NULL, list(c(0, 0), c(L[5], C[5], B[5])),
       list(c(0, 0), c(L[5], C[5], B[5])))
rmix_sixths <- list(NULL, list(c(0, 0), c(L[6], C[6], B[6])),
       list(c(0, 0), c(L[6], C[6], B[6])))
rmix_Six <- list(NULL, list(NULL, NULL, 1.75, NULL, 0.03),
  list(NULL, NULL, 1.75, NULL, 0.03))

rskews <- list(NULL, c(Stcum5[3], Stcum6[3]),
                            c(Stcum5[3], Stcum6[3]))
rskurts <- list(NULL, c(Stcum5[4], Stcum6[4]),
                              c(Stcum5[4], Stcum6[4]))
rfifths <- list(NULL, c(Stcum5[5], Stcum6[5]),
                              c(Stcum5[5], Stcum6[5]))
rsixths <- list(NULL, c(Stcum5[6], Stcum6[6]),
                              c(Stcum5[6], Stcum6[6]))
rSix <- list(NULL, list(1.75, NULL), list(1.75, NULL))

rmeans <- list(NULL,
  c(Nstcum[1], Stcum5[1], Stcum6[1], Mstcum[1]),
  c(Stcum5[1], Nstcum[1], Stcum6[1], Mstcum[1]))
rvars <- list(NULL,
  c(Nstcum[2]^2, Stcum5[2]^2, Stcum6[2]^2, Mstcum[2]^2),
  c(Stcum5[2]^2, Nstcum[2]^2, Stcum6[2]^2, Mstcum[2]^2))
means <- append(means0, rmeans)
vars <- append(vars0, rvars)
skews <- append(skews0, rskews)
skurts <- append(skurts0, rskurts)
fifths <- append(fifths0, rfifths)
sixths <- append(sixths0, rsixths)
Six <- append(Six0, rSix)
mix_pis <- append(mix_pis0, rmix_pis)
mix_mus <- append(mix_mus0, rmix_mus)
mix_sigmas <- append(mix_sigmas0, rmix_sigmas)
mix_skews <- append(mix_skews0, rmix_skews)
mix_skurts <- append(mix_skurts0, rmix_skurts)
mix_fifths <- append(mix_fifths0, rmix_fifths)
mix_sixths <- append(mix_sixths0, rmix_sixths)
mix_Six <- append(mix_Six0, rmix_Six)

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
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
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
  errorloop = TRUE, epsilon = 0.001)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Polynomial method: 2 cont, 1 mix, same.var, subj.var,
          error mix; rand.int = none, mix, non_mix;
          rand.tsl = none, non_mix, mix;
          mix, non_mix rand.var 2 per Y2, Y3", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04495033, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
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
  rand.int, rand.tsl, rand.var, corr.u, seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
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
  errorloop = TRUE, epsilon = 0.001)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e, N2$U, N2$U_all, rand.int, rand.tsl,
  corr.u, N2$rmeans2, N2$rvars2)

test_that("works for Fleishman method: 2 cont, 1 mix, same.var, subj.var,
          error mix; rand.int = none, mix, non_mix;
          rand.tsl = none, non_mix, mix;
          mix, non_mix rand.var 2 per Y2, Y3", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
    rand.int, rand.tsl, rand.var, corr.u), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565185, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[3, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

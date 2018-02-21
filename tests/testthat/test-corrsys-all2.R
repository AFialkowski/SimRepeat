skip_on_cran()
library("SimRepeat")
context("Simulate correlated system using method 2: all types")

options(scipen = 999)
tol <- 1e-4

seed <- 276
n <- 10000
M <- 3
Time <- 1:M

L <- calc_theory("Logistic", c(0, 1))
C <- calc_theory("Chisq", 4)
B <- calc_theory("Beta", c(4, 1.5))

skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

mix_pis <- lapply(seq_len(M), function(x) list(c(0.4, 0.6)))
mix_mus <- lapply(seq_len(M), function(x) list(c(-2, 2)))
mix_sigmas <- lapply(seq_len(M), function(x) list(c(1, 1)))
mix_skews <- lapply(seq_len(M), function(x) list(c(0, 0)))
mix_skurts <- lapply(seq_len(M), function(x) list(c(0, 0)))
mix_fifths <- lapply(seq_len(M), function(x) list(c(0, 0)))
mix_sixths <- lapply(seq_len(M), function(x) list(c(0, 0)))
mix_Six <- list()
Nstcum <- calc_mixmoments(mix_pis[[1]][[1]], mix_mus[[1]][[1]],
  mix_sigmas[[1]][[1]], mix_skews[[1]][[1]], mix_skurts[[1]][[1]],
  mix_fifths[[1]][[1]], mix_sixths[[1]][[1]])

means <- lapply(seq_len(M), function(x) c(Nstcum[1], B[1]))
vars <- lapply(seq_len(M), function(x) c(Nstcum[2]^2, B[2]^2))

marginal <- lapply(seq_len(M), function(x) list(0.4))
support <- list(NULL, list(c(0, 1)), NULL)
lam <- list(1, 5, 10)
p_zip <- list(NULL, 0.05, 0.1)
size <- list(10, 15, 20)
prob <- list(0.3, 0.4, 0.5)
mu <- mapply(function(x, y) x * (1 - y)/y, size, prob, SIMPLIFY = FALSE)
p_zinb <- list(NULL, 0.05, 0.1)
pois_eps <- nb_eps <- list()

same.var <- c(2, 3)

K <- 5
corr.x <- list()
corr.x[[1]] <- list(matrix(0.1, K, K), matrix(0.2, K, K), matrix(0.3, K, K))
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[1]][2:3, 2:3] <- diag(2)
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]] <- list(t(corr.x[[1]][[2]]), matrix(0.35, K, K),
                    matrix(0.4, K, K))
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[2]][2:3, 2:3] <- diag(2)
corr.x[[2]][[2]][, same.var] <- corr.x[[2]][[3]][, same.var] <-
  t(corr.x[[1]][[2]][same.var, ])
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[2]][[2]][same.var, ] <- t(corr.x[[2]][[2]][, same.var])

corr.x[[3]] <- list(t(corr.x[[1]][[3]]), t(corr.x[[2]][[3]]),
                    matrix(0.5, K, K))
diag(corr.x[[3]][[3]]) <- 1
corr.x[[3]][[3]][2:3, 2:3] <- diag(2)
corr.x[[3]][[3]][, same.var] <- t(corr.x[[1]][[3]][same.var, ])
corr.x[[3]][[3]][same.var, ] <- t(corr.x[[3]][[3]][, same.var])

# Error terms have a beta(4, 1.5) distribution with an AR(1, p = 0.4)
# correlation structure
error_type <- "non_mix"
corr.e <- matrix(c(1, 0.4, 0.4^2,
                   0.4, 1, 0.4,
                   0.4^2, 0.4, 1), M, M, byrow = TRUE)

subj.var <- matrix(c(1, 2, 1, 3, 2, 2, 2, 3, 3, 2, 3, 3), 6, 2, byrow = TRUE)
int.var <- tint.var <- NULL
betas.0 <- 0
betas <- list(seq(0.5, 0.5 + (K - 2) * 0.25, 0.25))
betas.subj <- list(seq(0.5, 0.5 + (K - 2) * 0.1, 0.1))
betas.int <- list()
betas.t <- 1
betas.tint <- list(c(0.25, 0.5))

method <- "Polynomial"
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 1 ord, 1 mix, 1 Poisson, 1 NB,
          same.var, subj.var, error non_mix", {
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
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 1 ord, 1 mix, 1 Poisson, 1 NB,
          same.var, subj.var, error non_mix", {
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

int.var <- matrix(c(1, 1, 4, 2, 2, 3), 2, 3, byrow = TRUE)
betas.int <- list(1, 1, NULL)
betas.subj <- list(c(0.25, 0.5, 0.75, 0.25, 0.5, 0.75),
                   c(0.25, 0.5, 0.75, 0.25, 0.5, 0.75),
                   c(0.25, 0.5, 0.75, 0.25))
betas.tint <- list(c(0.25, 0.5, 0.75), c(0.25, 0.5), c(0.25, 0.5))

method <- "Polynomial"
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 1 ord, 1 mix, 1 Poisson, 1 NB,
          same.var, subj.var, int.var, error non_mix", {
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
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 1 ord, 1 mix, 1 Poisson, 1 NB,
          same.var, subj.var, int.var, error non_mix", {
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

betas <- list(seq(0.5, 0.5 + (K - 2) * 0.25, 0.25))
int.var <- matrix(c(1, 1, 4, 1, 2, 3, 3, 2, 3), 3, 3, byrow = TRUE)
betas.int <- list(c(1, 1), NULL, 1)
subj.var <- matrix(c(1, 2, 1, 3), 2, 2, byrow = TRUE)
betas.subj <- list(c(0.25, 0.5, 0.75, 0.25, 0.5, 0.75, 0.25, 0.5, 0.75),
                   NULL, NULL)
tint.var = matrix(c(1, 5, 1, 6, 3, 1, 3, 5), 4, 2, byrow = TRUE)
betas.tint <- list(c(0.25, 0.5), c(0.25, 0.5, 0.75, 1), c(0.25, 0.5))

method <- "Polynomial"
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 1 ord, 1 mix, 1 Poisson, 1 NB,
          same.var, subj.var, int.var, tint.var, error non_mix", {
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
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, E_mix = NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 1 ord, 1 mix, 1 Poisson, 1 NB,
          same.var, subj.var, int.var, tint.var, error non_mix", {
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

means <- lapply(seq_len(M), function(x) c(B[1], Nstcum[1]))
vars <- lapply(seq_len(M), function(x) c(B[2]^2, Nstcum[2]^2))

same.var <- 1

K <- 4
corr.x <- list()
corr.x[[1]] <- list(matrix(0.1, K, K), matrix(0.2, K, K), matrix(0.3, K, K))
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]] <- list(t(corr.x[[1]][[2]]), matrix(0.35, K, K),
                    matrix(0.4, K, K))
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[2]][, same.var] <- corr.x[[2]][[3]][, same.var] <-
  t(corr.x[[1]][[2]][same.var, ])
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[2]][[2]][same.var, ] <- t(corr.x[[2]][[2]][, same.var])

corr.x[[3]] <- list(t(corr.x[[1]][[3]]), t(corr.x[[2]][[3]]),
                    matrix(0.5, K, K))
diag(corr.x[[3]][[3]]) <- 1
corr.x[[3]][[3]][, same.var] <- t(corr.x[[1]][[3]][same.var, ])
corr.x[[3]][[3]][same.var, ] <- t(corr.x[[3]][[3]][, same.var])

subj.var <- matrix(c(2, 4), 1, 2)
int.var <- tint.var <- NULL
betas <- list(seq(0.5, 0.5 + (K - 1) * 0.25, 0.25))
betas.subj <- list(NULL, c(0.5, 0.75, 1), NULL)
betas.int <- list()
betas.tint <- list(rep(0.25, K), rep(0.25, K - 1), rep(0.25, K))

corr.e <- matrix(0.4, 6, 6)
corr.e[1:2, 1:2] <- corr.e[3:4, 3:4] <- corr.e[5:6, 5:6] <- diag(2)
error_type = "mix"

method <- "Polynomial"
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: 1 ord, 1 cont, 1 Poisson, 1 NB,
          same.var, subj.var, error mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[3, "Mean"],
    -2, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[3, "Mean"],
    -2, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: 1 ord, 1 cont, 1 Poisson, 1 NB,
          same.var, subj.var, error mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[3, "Mean"],
    -2, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[3, "Mean"],
    -2, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- list(NULL, B[3], B[3])
skurts <- list(NULL, B[4], B[4])
fifths <- list(NULL, B[5], B[5])
sixths <- list(NULL, B[6], B[6])
Six <- list(NULL, list(c(0.01, 0.02, 0.03)), list(0.03))

means <- list(Nstcum[1], c(B[1], Nstcum[1]), c(B[1], Nstcum[1]))
vars <- list(Nstcum[2]^2, c(B[2]^2, Nstcum[2]^2), c(B[2]^2, Nstcum[2]^2))

marginal <- list(NULL, list(0.4), list(0.4))
support <- list(NULL, list(c(0, 1)), list(c(0, 1)))
lam <- list(NULL, 5, 10)
p_zip <- list(NULL, 0.05, 0.1)
size <- list(NULL, 15, 20)
prob <- list(NULL, 0.4, 0.5)
mu <- mapply(function(x, y) x * (1 - y)/y, size, prob, SIMPLIFY = FALSE)
p_zinb <- list(NULL, 0.05, 0.1)
pois_eps <- nb_eps <- list()

K <- 4
corr.x <- list(NULL, list(NULL, matrix(0.4, K, K), matrix(0.4, K, K)),
               list(NULL, matrix(0.4, K, K), matrix(0.4, K, K)))
diag(corr.x[[2]][[2]]) <- 1
diag(corr.x[[3]][[3]]) <- 1

same.var <- NULL
subj.var <- matrix(c(2, 1), 1, 2, byrow = TRUE)
int.var <- NULL
tint.var <- matrix(c(2, 1), 1, 2, byrow = TRUE)
betas <- list(NULL, rep(0.5, K), rep(0.5, K))
betas.subj <- list(NULL, rep(0.5, K - 1), NULL)
betas.int <- list()
betas.tint <- list(NULL, 0.25, rep(0.25, K))

method <- "Polynomial"
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: Y1 no X, 1 ord, 1 cont, 1 Poisson, 1 NB,
          tint.var, subj.var, error mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    -2, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    -2, tolerance = tol, check.attributes = FALSE), TRUE)
})

method <- "Fleishman"
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, N$E_mix, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, N2$E_mix, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: Y1 no X, 1 ord, 1 cont, 1 Poisson, 1 NB,
          tint.var, subj.var, error mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    0, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    -2, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    -2, tolerance = tol, check.attributes = FALSE), TRUE)
})

skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

mix_pis <- list(list(c(0.4, 0.6)), NULL, list(c(0.4, 0.6)))
mix_mus <- list(list(c(-2, 2)), NULL, list(c(-2, 2)))
mix_sigmas <- list(list(c(1, 1)), NULL, list(c(1, 1)))
mix_skews <- mix_skurts <- mix_fifths <- mix_sixths <-
  list(list(c(0, 0)), NULL, list(c(0, 0)))
mix_Six <- list()
means <- list(c(Nstcum[1], B[1]), B[1], c(Nstcum[1], B[1]))
vars <- list(c(Nstcum[2]^2, B[2]^2), B[2]^2, c(Nstcum[2]^2, B[2]^2))

marginal <- list(list(0.4), NULL, list(0.4))
support <- list(list(c(0, 1)), NULL, NULL)
lam <- list(1, NULL, 10)
p_zip <- list(NULL, NULL, 0.1)
size <- list(10, NULL, 20)
prob <- list(0.3, NULL, 0.5)
mu <- mapply(function(x, y) x * (1 - y)/y, size, prob, SIMPLIFY = FALSE)
p_zinb <- list(NULL, NULL, 0.1)
pois_eps <- nb_eps <- list()

same.var <- matrix(c(1, 2, 3, 2, 1, 3, 3, 3), 2, 4, TRUE)

K <- 5
corr.x <- list(list(matrix(0.4, K, K), NULL, matrix(0.4, K, K)), NULL,
               list(matrix(0.4, K, K), NULL, matrix(0.4, K, K)))
diag(corr.x[[1]][[1]]) <- 1
corr.x[[1]][[1]][2:3, 2:3] <- diag(2)
corr.x[[1]][[3]][, 2:3] <- corr.x[[1]][[1]][, 2:3]
diag(corr.x[[3]][[3]]) <- 1
corr.x[[3]][[1]] <- t(corr.x[[1]][[3]])
corr.x[[3]][[3]][2:3, 2:3] <- diag(2)

error_type <- "non_mix"
corr.e <- matrix(c(1, 0.4, 0.4^2,
                   0.4, 1, 0.4,
                   0.4^2, 0.4, 1), M, M, byrow = TRUE)

K <- 4
subj.var <- matrix(c(1, 2, 1, 4), 2, 2, byrow = TRUE)
int.var <- tint.var <- NULL
betas <- list(rep(0.5, K), NULL, rep(0.5, K))
betas.subj <- list(c(0.25, 0.5, 0.75, 1), NULL, NULL)
betas.int <- list()
betas.tint <- list(c(0.25, 0.5), NULL, rep(0.5, K))

method <- "Polynomial"
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: Y2 no X, 1 ord, 1 mix, 1 Poisson, 1 NB,
          subj.var, error non_mix", {
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
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)
N2 <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE, errorloop = TRUE, epsilon = 0.01)
S2 <- summary_sys(N2$Y, N2$E, NULL, N2$X, N2$X_all, M, method, means,
  vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: Y2 no X, 1 ord, 1 mix, 1 Poisson, 1 NB,
          subj.var, error non_mix", {
  expect_equal(checkpar(M, method, error_type, means, vars, skews,
    skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas,
    mix_skews, mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal,
    support, lam, p_zip, pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x,
    corr.yx = list(), corr.e, same.var, subj.var, int.var, tint.var,
    betas.0, betas, betas.subj, betas.int, betas.t, betas.tint), TRUE)
  expect_equal(all.equal(N$constants[[1]][nrow(N$constants[[1]]), "c3"],
    -0.04565308, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(N2$constants[[1]][nrow(N2$constants[[1]]), "c3"],
    -0.04565308, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
  expect_equal(all.equal(S2$target_sum_e[1, "Mean"],
    0.7272727, tolerance = tol, check.attributes = FALSE), TRUE)
})

M <- 2
Time <- 1:M
skews <- lapply(seq_len(M), function(x) B[3])
skurts <- lapply(seq_len(M), function(x) B[4])
fifths <- lapply(seq_len(M), function(x) B[5])
sixths <- lapply(seq_len(M), function(x) B[6])
Six <- lapply(seq_len(M), function(x) list(0.03))

mix_pis <- list(NULL, list(c(0.4, 0.6)))
mix_mus <- list(NULL, list(c(-2, 2)))
mix_sigmas <- list(NULL, list(c(1, 1)))
mix_skews <- mix_skurts <- mix_fifths <- mix_sixths <-
  list(NULL, list(c(0, 0)))
mix_Six <- list()
means <- list(B[1], c(Nstcum[1], B[1]))
vars <- list(B[2]^2, c(Nstcum[2]^2, B[2]^2))

marginal <- list(NULL, list(0.4))
support <- list()
lam <- list(NULL, 10)
p_zip <- list()
size <- list(NULL, 20)
prob <- list(NULL, 0.5)
mu <- mapply(function(x, y) x * (1 - y)/y, size, prob, SIMPLIFY = FALSE)
p_zinb <- list()
pois_eps <- nb_eps <- list()

corr.x <- list(NULL, list(NULL, matrix(0.4, 5, 5)))
diag(corr.x[[2]][[2]]) <- 1
corr.x[[2]][[2]][2:3, 2:3] <- diag(2)
corr.e <- corr.e[1:2, 1:2]

K <- 4
same.var <- NULL
subj.var <- matrix(c(2, 2, 2, 3), 2, 2, byrow = TRUE)
int.var <- matrix(c(2, 2, 3), 1, 3)
tint.var <- NULL
betas.0 <- 0
betas.t <- 1
betas <- list(NULL, rep(0.5, K))
betas.subj <- list(NULL, c(0.25, 0.5, 0.75, 1, 1, 1))
betas.int <- list(NULL, 1)
betas.tint <- list(NULL, c(0.5, 0.5))

method <- "Polynomial"
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Polynomial method: X for Y2 only,
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
N <- corrsys2(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  pois_eps, size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means, vars,
  skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, marginal, support, lam, p_zip, size,
  prob, mu, p_zinb, corr.x, corr.e)

test_that("works for Fleishman method: X for Y2 only,
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


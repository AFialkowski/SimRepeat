<!-- README.md is generated from README.Rmd. Please edit that file -->
SimRepeat
=========

The goal of **SimRepeat** is to generate correlated systems of statistical equations which represent **repeated measurements** or clustered data. These systems contain either: **a)** continuous normal, non-normal, and mixture variables based on the techniques of Headrick and Beasley (2004, <doi:10.1081/SAC-120028431>) or *b)* continuous (normal, non-normal and mixture), ordinal, and count (regular or zero-inflated, Poisson and Negative Binomial) variables based on the hierarchical linear models (HLM) approach. Headrick and Beasley's method for continuous variables calculates the beta (slope) coefficients based on the target correlations between independent variables and between outcomes and independent variables. The package provides functions to calculate the expected correlations between outcomes, between outcomes and error terms, and between outcomes and independent variables, extending Headrick and Beasley's equations to include mixture variables. These theoretical values can be compared to the simulated correlations. The HLM approach requires specification of the beta coefficients, but permits group and subject-level independent variables, interactions among independent variables, and fixed and random effects, providing more flexibility in the system of equations. Both methods permit simulation of data sets that mimic real-world clinical or genetic data sets (i.e. plasmodes, as in Vaughan et al., 2009, <doi:10.1016/j.csda.2008.02.032>). The techniques extend those found in the **SimMultiCorrData** and **SimCorrMix** packages. Standard normal variables with an imposed intermediate correlation matrix are transformed to generate the desired distributions. Continuous variables are simulated using either Fleishman's third-order (<doi:10.1007/BF02293811>) or Headrick's fifth-order (<doi:10.1016/S0167-9473(02)00072-5>) power method transformation (PMT). Simulation occurs at the component-level for continuous mixture distributions. These components are transformed into the desired mixture variables using random multinomial variables based on the mixing probabilities. The target correlation matrices are specified in terms of correlations with components of continuous mixture variables. Binary and ordinal variables are simulated using a modification of **GenOrd**'s `ordsample` function. Count variables are simulated using the inverse CDF method. There are two simulation pathways for the multi-variable type systems which differ by intermediate correlations involving count variables. Correlation Method 1 adapts Yahav and Shmueli's 2012 method (<doi:10.1002/asmb.901>). Correlation Method 2 adapts Barbiero and Ferrari's 2015 modification of **GenOrd** (<doi:10.1002/asmb.2072>). The optional error loop may be used to improve the accuracy of the final correlation matrices. The package also provides function to check parameter inputs and summarize the generated systems of equations.

There are vignettes which accompany this package that may help the user understand the simulation and analysis methods.

1.  **Theory and Equations for Correlated Systems of Continuous Variables** describes the system of continuous variables generated with `nonnormsys` and derives the equations used in `calc_betas`, `calc_corr_y`, `calc_corr_ye`, and `calc_corr_yx`.

2.  **Correlated Systems of Statistical Equations with Non-Mixture and Mixture Continuous Variables** provides examples of using `nonnormsys`.

3.  **The Hierarchical Linear Models Approach for a System of Correlated Equations with Multiple Variable Types** describes the system of ordinal, continuous, and count variables generated with `corrsys` and `corrsys2`.

4.  **Correlated Systems of Statistical Equations with Multiple Variable Types** provides examples of using `corrsys` and `corrsys2`.

Installation instructions
-------------------------

**SimRepeat** can be installed using the following code:

``` r
## from GitHub
install.packages("devtools")
devtools::install_github("AFialkowski/SimRepeat", build_vignettes = TRUE)
```

Example: System of three equations for 5 independent variables with no random effects
-------------------------------------------------------------------------------------

### Description of Variables

1.  **Ordinal variable:** *X*<sub>*o**r**d*(1)</sub> has 3 categories (i.e., drug treatment) and is the same in each equation
2.  **Continuous variables:**

<!-- -->

1.  *X*<sub>*c**o**n**t*</sub> is a time-varying covariate (subject-level term) with an AR(1, p = 0.5) correlation structure

<!-- -->

1.  *X*<sub>*c**o**n**t*(11)</sub> has a Chisq(df = 2) distribution
2.  *X*<sub>*c**o**n**t*(21)</sub> has a Chisq(df = 4) distribution
3.  *X*<sub>*c**o**n**t*(31)</sub> has a Chisq(df = 8) distribution

<!-- -->

1.  *X*<sub>*m**i**x*</sub> is a normal mixture time-varying covariate (subject-level term), components have an AR(1, p = 0.4) correlation structure across Y

<!-- -->

1.  **Poisson variable:** *X*<sub>*p**o**i**s*(1)</sub> is a zero-inflated Poisson variable with *λ* = 15, the probability of a structural zero set at 0.10, and is the same in each equation
2.  **Negative Binomial variable:** *X*<sub>*n**b*(1)</sub> is a regular NB time-varying covariate (subject-level term) with an AR(1, p = 0.3) correlation structure and increasing mean and variance

<!-- -->

1.  *X*<sub>*n**b*(11)</sub> has a size of 10 and mean of 3
2.  *X*<sub>*n**b*(21)</sub> has a size of 10 and mean of 4
3.  *X*<sub>*n**b*(31)</sub> has a size of 10 and mean of 5

<!-- -->

1.  **Error terms** have a Beta(4, 1.5) distribution with an AR(1, p = 0.4) correlation structure. These require a sixth cumulant correction of 0.03.

There is an interaction between X\_ord(1) and X\_pois(1) for each Y. Since they are both group-level covariates, the interaction is also a group-level covariate that will interact with the subject-level covariates. However, only X\_ord(1) and X\_pois(1) interact with time. The ordering in the equations below reflects the ordering in the simulation process.

``` r
library("SimRepeat")
library("printr")
library("lme4")
library("reshape2")
options(scipen = 999)
```

### Step 1: Set up parameter inputs

This is the most time-consuming part of the simulation process. It is important to read the function documentation carefully to understand the formats for each parameter input. Incorrect formatting will lead to errors. Most of these can be prevented by using the `checkpar` function in **Step 2**.

``` r
seed <- 11
n <- 10000
M <- 3

# Ordinal variable
marginal <- lapply(seq_len(M), function(x) list(c(1/3, 2/3)))
support <- lapply(seq_len(M), function(x) list(c(0, 1, 2)))

# Non-mixture continuous variables
method <- "Polynomial"
Stcum1 <- calc_theory("Chisq", 2)
Stcum2 <- calc_theory("Chisq", 4)
Stcum3 <- calc_theory("Chisq", 8)

# Error terms
error_type <- "non_mix"
Stcum4 <- calc_theory("Beta", c(4, 1.5))
corr.e <- matrix(c(1, 0.4, 0.4^2, 0.4, 1, 0.4, 0.4^2, 0.4, 1), M, M, 
  byrow = TRUE)

skews <- list(c(Stcum1[3], Stcum4[3]), c(Stcum2[3], Stcum4[3]), 
  c(Stcum3[3], Stcum4[3]))
skurts <- list(c(Stcum1[4], Stcum4[4]), c(Stcum2[4], Stcum4[4]), 
  c(Stcum3[4], Stcum4[4]))
fifths <- list(c(Stcum1[5], Stcum4[5]), c(Stcum2[5], Stcum4[5]), 
  c(Stcum3[5], Stcum4[5]))
sixths <- list(c(Stcum1[6], Stcum4[6]), c(Stcum2[6], Stcum4[6]), 
  c(Stcum3[6], Stcum4[6]))
Six <- lapply(seq_len(M), function(x) list(NULL, 0.03))

# Mixture continuous variable
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

means <- list(c(Stcum1[1], Nstcum[1], Stcum4[1]),
              c(Stcum2[1], Nstcum[1], Stcum4[1]),
              c(Stcum3[1], Nstcum[1], Stcum4[1]))
vars <- list(c(Stcum1[2]^2, Nstcum[2]^2, Stcum4[2]^2),
             c(Stcum2[2]^2, Nstcum[2]^2, Stcum4[2]^2),
             c(Stcum3[2]^2, Nstcum[2]^2, Stcum4[2]^2))

# Poisson variable
lam <- list(15, 15, 15)
p_zip <- 0.10

# Negative Binomial variables
size <- list(10, 10, 10)
mu <- list(3, 4, 5)
prob <- list()
p_zinb <- 0

# X_ord(11) and X_pois(11) are the same across Y
same.var <- c(1, 5)

# set up X correlation matrix
corr.x <- list()
corr.x[[1]] <- list(matrix(0.4, 6, 6), matrix(0.35, 6, 6), matrix(0.25, 6, 6))
diag(corr.x[[1]][[1]]) <- 1
# set correlations between components of X_mix(11) to 0
corr.x[[1]][[1]][3:4, 3:4] <- diag(2)
# set correlations between time-varying covariates of Y1 and Y2
corr.x[[1]][[2]][2, 2] <- 0.5
corr.x[[1]][[2]][3:4, 3:4] <- matrix(0.4, 2, 2)
corr.x[[1]][[2]][6, 6] <- 0.3
# set correlations between time-varying covariates of Y1 and Y3
corr.x[[1]][[3]][2, 2] <- 0.5^2
corr.x[[1]][[3]][3:4, 3:4] <- matrix(0.4^2, 2, 2)
corr.x[[1]][[3]][6, 6] <- 0.3^2
# set correlations for the same variables equal across outcomes
corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
  corr.x[[1]][[1]][, same.var]

corr.x[[2]] <- list(t(corr.x[[1]][[2]]), matrix(0.35, 6, 6), 
  matrix(0.25, 6, 6))
diag(corr.x[[2]][[2]]) <- 1
# set correlations between components of X_mix(21) to 0
corr.x[[2]][[2]][3:4, 3:4] <- diag(2)
# set correlations between time-varying covariates of Y2 and Y3
corr.x[[2]][[3]][2, 2] <- 0.5
corr.x[[2]][[3]][3:4, 3:4] <- matrix(0.4, 2, 2)
corr.x[[2]][[3]][6, 6] <- 0.3
# set correlations for the same variables equal across outcomes
corr.x[[2]][[2]][same.var, ] <- corr.x[[1]][[2]][same.var, ]
corr.x[[2]][[2]][, same.var] <- corr.x[[2]][[3]][, same.var] <- 
  t(corr.x[[1]][[2]][same.var, ])
corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]

corr.x[[3]] <- list(t(corr.x[[1]][[3]]), t(corr.x[[2]][[3]]), 
  matrix(0.3, 6, 6))
diag(corr.x[[3]][[3]]) <- 1
# set correlations between components of X_mix(31) to 0
corr.x[[3]][[3]][3:4, 3:4] <- diag(2)
# set correlations for the same variables equal across outcomes
corr.x[[3]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
corr.x[[3]][[3]][, same.var] <- t(corr.x[[3]][[3]][same.var, ])

Time <- 1:M
betas.0 <- 0
betas.t <- 1
# use a list of length 1 so that betas are the same across Y
betas <- list(seq(0.5, 1.5, 0.25))
# interaction between ordinal and Poisson variable, becomes 
# another group-level variable
int.var <- matrix(c(1, 1, 4, 2, 1, 4, 3, 1, 4), 3, 3, byrow = TRUE)
betas.int <- list(0.5)
# continuous non-mixture, continuous mixture, and NB variables are 
# subject-level variables
subj.var <- matrix(c(1, 2, 1, 3, 1, 5, 2, 2, 2, 3, 2, 5, 3, 2, 3, 3, 3, 5), 
  nrow = 9, ncol = 2, byrow = TRUE)
# there are 3 subject-level variables and 3 group-level variables forming 
# 9 group-subject interactions
betas.subj <- list(seq(0.5, 0.5 + (9 - 1) * 0.1, 0.1))
# only ordinal and Poisson variable interact with time (excluding the 
# ordinal-Poisson interaction variable)
tint.var <- matrix(c(1, 1, 1, 4, 2, 1, 2, 4, 3, 1, 3, 4), 6, 2, byrow = TRUE)
betas.tint <- list(c(0.25, 0.5))
```

### Step 2: Check parameter inputs

``` r
checkpar(M, method, error_type, means, vars, skews, skurts, fifths, sixths, 
  Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, mix_fifths, 
  mix_sixths, mix_Six, marginal, support, lam, p_zip, pois_eps = list(), 
  size, prob, mu, p_zinb, nb_eps = list(), corr.x, corr.yx = list(), corr.e, 
  same.var, subj.var, int.var, tint.var, betas.0, betas, betas.subj, betas.int, 
  betas.t, betas.tint)
#> [1] TRUE
```

### Step 3: Generate system

Note that `use.nearPD = FALSE` so that negative eigen-values will be replaced with 0 instead of using the nearest positive-definite matrix (found with @Matrix's `Matrix::nearPD` function with @Higham's algorithm).

``` r
Sys1 <- corrsys(n, M, Time, method, error_type, means, vars,
  skews, skurts, fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews,
  mix_skurts, mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip,
  size, prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var,
  tint.var, betas.0, betas, betas.subj, betas.int, betas.t, betas.tint,
  seed = seed, use.nearPD = FALSE)
#> Total Simulation time: 0.212 minutes
```

``` r
knitr::kable(Sys1$constants[[1]], booktabs = TRUE, 
  caption = "PMT constants for Y_1")
```

|          c0|         c1|          c2|          c3|          c4|         c5|
|-----------:|----------:|-----------:|-----------:|-----------:|----------:|
|  -0.3077396|  0.8005605|   0.3187640|   0.0335001|  -0.0036748|  0.0001587|
|   0.0000000|  1.0000000|   0.0000000|   0.0000000|   0.0000000|  0.0000000|
|   0.0000000|  1.0000000|   0.0000000|   0.0000000|   0.0000000|  0.0000000|
|   0.1629657|  1.0899841|  -0.1873287|  -0.0449503|   0.0081210|  0.0014454|

``` r
Sys1$valid.pdf
#> [[1]]
#> [1] "TRUE" "TRUE" "TRUE" "TRUE"
#> 
#> [[2]]
#> [1] "TRUE" "TRUE" "TRUE" "TRUE"
#> 
#> [[3]]
#> [1] "TRUE" "TRUE" "TRUE" "TRUE"
```

### Step 4: Describe results

``` r
Sum1 <- summary_sys(Sys1$Y, Sys1$E, E_mix = NULL, Sys1$X, Sys1$X_all, M, 
  method, means, vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, 
  mix_sigmas, mix_skews, mix_skurts, mix_fifths, mix_sixths, marginal, 
  support, lam, p_zip, size, prob, mu, p_zinb, corr.x, corr.e)
names(Sum1)
#>  [1] "cont_sum_y"   "rho.y"        "cont_sum_e"   "target_sum_e"
#>  [5] "rho.e"        "rho.ye"       "ord_sum_x"    "cont_sum_x"  
#>  [9] "target_sum_x" "sum_xall"     "mix_sum_x"    "target_mix_x"
#> [13] "pois_sum_x"   "nb_sum_x"     "rho.x"        "rho.xall"    
#> [17] "rho.yx"       "rho.yxall"    "maxerr"
```

``` r
knitr::kable(Sum1$cont_sum_y, digits = 3, booktabs = TRUE, 
  caption = "Simulated Distributions of Outcomes")
```

|     |  Outcome|      N|     Mean|       SD|   Median|       Min|       Max|   Skew|  Skurtosis|  Fifth|   Sixth|
|-----|--------:|------:|--------:|--------:|--------:|---------:|---------:|------:|----------:|------:|-------:|
| Y1  |        1|  10000|  247.580|  260.176|  162.216|  -127.270|  1958.432|  1.702|      3.628|  8.056|  12.565|
| Y2  |        2|  10000|  338.192|  324.591|  238.280|   -57.635|  2437.983|  1.561|      3.044|  6.294|   9.708|
| Y3  |        3|  10000|  458.381|  397.247|  347.081|   -11.825|  3154.605|  1.376|      2.352|  4.561|   9.584|

``` r
knitr::kable(Sum1$target_sum_e, digits = 3, booktabs = TRUE, 
  caption = "Target Distributions of Error Terms")
```

|     |  Outcome|   Mean|     SD|    Skew|  Skurtosis|  Fifth|   Sixth|
|-----|--------:|------:|------:|-------:|----------:|------:|-------:|
| E1  |        1|  0.727|  0.175|  -0.694|     -0.069|  1.828|  -3.379|
| E2  |        2|  0.727|  0.175|  -0.694|     -0.069|  1.828|  -3.379|
| E3  |        3|  0.727|  0.175|  -0.694|     -0.069|  1.828|  -3.379|

``` r
knitr::kable(Sum1$cont_sum_e, digits = 3, booktabs = TRUE, 
  caption = "Simulated Distributions of Error Terms")
```

|     |  Outcome|      N|   Mean|     SD|  Median|    Min|    Max|    Skew|  Skurtosis|  Fifth|   Sixth|
|-----|--------:|------:|------:|------:|-------:|------:|------:|-------:|----------:|------:|-------:|
| E1  |        1|  10000|  0.727|  0.175|   0.756|  0.102|  1.016|  -0.702|     -0.080|  1.973|  -3.829|
| E2  |        2|  10000|  0.727|  0.175|   0.757|  0.083|  1.042|  -0.704|     -0.034|  1.754|  -3.450|
| E3  |        3|  10000|  0.727|  0.174|   0.756|  0.090|  1.054|  -0.675|     -0.136|  1.934|  -3.192|

``` r
knitr::kable(Sum1$target_sum_x, digits = 3, booktabs = TRUE, 
  caption = "Target Distributions of Continuous Non-Mixture and Components of 
  Mixture Variables")
```

|          |  Outcome|    X|  Mean|     SD|   Skew|  Skurtosis|   Fifth|  Sixth|
|----------|--------:|----:|-----:|------:|------:|----------:|-------:|------:|
| cont1\_1 |        1|    1|     2|  2.000|  2.000|        6.0|  24.000|  120.0|
| cont1\_2 |        1|    2|    -2|  1.000|  0.000|        0.0|   0.000|    0.0|
| cont1\_3 |        1|    3|     2|  1.000|  0.000|        0.0|   0.000|    0.0|
| cont2\_1 |        2|    1|     4|  2.828|  1.414|        3.0|   8.485|   30.0|
| cont2\_2 |        2|    2|    -2|  1.000|  0.000|        0.0|   0.000|    0.0|
| cont2\_3 |        2|    3|     2|  1.000|  0.000|        0.0|   0.000|    0.0|
| cont3\_1 |        3|    1|     8|  4.000|  1.000|        1.5|   3.000|    7.5|
| cont3\_2 |        3|    2|    -2|  1.000|  0.000|        0.0|   0.000|    0.0|
| cont3\_3 |        3|    3|     2|  1.000|  0.000|        0.0|   0.000|    0.0|

``` r
knitr::kable(Sum1$cont_sum_x, digits = 3, booktabs = TRUE, 
  caption = "Simulated Distributions of Continuous Non-Mixture and Components 
  of Mixture Variables")
```

|          |  Outcome|    X|      N|    Mean|     SD|  Median|     Min|     Max|    Skew|  Skurtosis|   Fifth|   Sixth|
|----------|--------:|----:|------:|-------:|------:|-------:|-------:|-------:|-------:|----------:|-------:|-------:|
| cont1\_1 |        1|    1|  10000|   2.002|  2.022|   1.385|  -0.697|  18.383|   1.996|      5.539|  17.959|  59.310|
| cont1\_2 |        1|    2|  10000|  -2.000|  1.001|  -2.003|  -6.346|   2.532|   0.032|     -0.017|   0.131|   0.633|
| cont1\_3 |        1|    3|  10000|   2.000|  1.001|   1.998|  -1.976|   6.014|   0.000|     -0.016|   0.107|   0.049|
| cont2\_1 |        2|    1|  10000|   4.000|  2.825|   3.396|  -0.150|  27.612|   1.424|      3.077|   8.786|  31.791|
| cont2\_2 |        2|    2|  10000|  -2.000|  1.002|  -2.001|  -5.383|   1.733|   0.040|     -0.003|   0.100|  -0.433|
| cont2\_3 |        2|    3|  10000|   2.000|  1.002|   1.993|  -1.753|   5.600|  -0.004|     -0.022|  -0.109|   0.102|
| cont3\_1 |        3|    1|  10000|   8.000|  3.999|   7.325|   0.254|  35.800|   0.983|      1.438|   3.030|   9.026|
| cont3\_2 |        3|    2|  10000|  -2.000|  1.001|  -1.998|  -5.857|   1.620|   0.005|     -0.040|  -0.061|   0.024|
| cont3\_3 |        3|    3|  10000|   2.000|  1.001|   2.003|  -2.311|   5.578|  -0.047|      0.120|  -0.204|  -0.157|

``` r
knitr::kable(Sum1$target_mix_x, digits = 3, booktabs = TRUE, 
  caption = "Target Distributions of Continuous Mixture Variables")
```

|         |  Outcome|    X|  Mean|   SD|    Skew|  Skurtosis|  Fifth|  Sixth|
|---------|--------:|----:|-----:|----:|-------:|----------:|------:|------:|
| mix1\_1 |        1|    1|   0.4|  2.2|  -0.289|     -1.154|  1.793|  6.173|
| mix2\_1 |        2|    1|   0.4|  2.2|  -0.289|     -1.154|  1.793|  6.173|
| mix3\_1 |        3|    1|   0.4|  2.2|  -0.289|     -1.154|  1.793|  6.173|

``` r
knitr::kable(Sum1$mix_sum_x, digits = 3, booktabs = TRUE, 
  caption = "Simulated Distributions of Continuous Mixture Variables")
```

|         |  Outcome|    X|      N|  Mean|   SD|  Median|     Min|    Max|    Skew|  Skurtosis|  Fifth|  Sixth|
|---------|--------:|----:|------:|-----:|----:|-------:|-------:|------:|-------:|----------:|------:|------:|
| mix1\_1 |        1|    1|  10000|   0.4|  2.2|   1.022|  -5.384|  6.046|  -0.273|     -1.172|  1.721|  6.394|
| mix2\_1 |        2|    1|  10000|   0.4|  2.2|   1.020|  -5.367|  5.602|  -0.281|     -1.179|  1.772|  6.380|
| mix3\_1 |        3|    1|  10000|   0.4|  2.2|   1.044|  -5.917|  5.579|  -0.293|     -1.147|  1.839|  6.062|

``` r
Nplot <- plot_simpdf_theory(sim_y = Sys1$X_all[[1]][, 3], ylower = -10, 
  yupper = 10, 
  title = "PDF of X_mix(11): Mixture of N(-5, 2) and N(3, 1) Distributions",
  fx = function(x) mix_pis[[1]][[1]][1] * dnorm(x, mix_mus[[1]][[1]][1], 
    mix_sigmas[[1]][[1]][1]) + mix_pis[[1]][[1]][2] * 
    dnorm(x, mix_mus[[1]][[1]][2], mix_sigmas[[1]][[1]][2]), 
  lower = -Inf, upper = Inf)
Nplot
```

![](man/figures/README-unnamed-chunk-17-1.png)

Summary of Ordinal Variable: (for *Y*<sub>1</sub>)

``` r
knitr::kable(Sum1$ord_sum_x[[1]][1:2, ], digits = 3, row.names = FALSE,
             booktabs = TRUE, caption = "Simulated Distribution of X_ord(1)")
```

|  Outcome|  Support|  Target|  Simulated|
|--------:|--------:|-------:|----------:|
|        1|        0|   0.333|      0.334|
|        1|        1|   0.667|      0.665|

Summary of Poisson Variable:

``` r
knitr::kable(Sum1$pois_sum_x, digits = 3, row.names = FALSE,
             booktabs = TRUE, caption = "Simulated Distribution of X_pois(1)")
```

|  Outcome|    X|      N|     P0|  Exp\_P0|    Mean|  Exp\_Mean|    Var|  Exp\_Var|  Median|  Min|  Max|    Skew|  Skurtosis|
|--------:|----:|------:|------:|--------:|-------:|----------:|------:|---------:|-------:|----:|----:|-------:|----------:|
|        1|    1|  10000|  0.102|      0.1|  13.494|       13.5|  33.94|        40|      14|    0|   32|  -0.876|      0.671|
|        2|    1|  10000|  0.102|      0.1|  13.494|       13.5|  33.94|        40|      14|    0|   32|  -0.876|      0.671|
|        3|    1|  10000|  0.102|      0.1|  13.494|       13.5|  33.94|        40|      14|    0|   32|  -0.876|      0.671|

``` r
Pplot <- plot_simpdf_theory(sim_y = Sys1$X_all[[1]][, 4], 
  title = "PMF of X_pois(1): Zero-Inflated Poisson Distribution", 
  Dist = "Poisson", params = c(lam[[1]][1], p_zip), cont_var = FALSE)
Pplot
```

![](man/figures/README-unnamed-chunk-19-1.png)

Summary of Negative Binomial Variables *X*<sub>*n**b*(11)</sub>, *X*<sub>*n**b*(21)</sub>, and *X*<sub>*n**b*(31)</sub>:

``` r
knitr::kable(Sum1$nb_sum_x, digits = 3, row.names = FALSE,
             booktabs = TRUE, caption = "Simulated Distributions")
```

|  Outcome|    X|      N|     P0|  Exp\_P0|   Prob|   Mean|  Exp\_Mean|    Var|  Exp\_Var|  Median|  Min|  Max|   Skew|  Skurtosis|
|--------:|----:|------:|------:|--------:|------:|------:|----------:|------:|---------:|-------:|----:|----:|------:|----------:|
|        1|    1|  10000|  0.074|    0.073|  0.769|  2.996|          3|  3.864|       3.9|       3|    0|   13|  0.803|      0.890|
|        2|    1|  10000|  0.034|    0.035|  0.714|  3.998|          4|  5.631|       5.6|       4|    0|   18|  0.781|      0.907|
|        3|    1|  10000|  0.018|    0.017|  0.667|  5.001|          5|  7.514|       7.5|       5|    0|   20|  0.720|      0.680|

``` r
NBplot <- plot_simtheory(sim_y = Sys1$X_all[[1]][, 5], 
  title = "Simulated Values for X_nb(11)", Dist = "Negative_Binomial", 
  params = c(size[[1]][1], mu[[1]][1], p_zinb), cont_var = FALSE)
NBplot
```

![](man/figures/README-unnamed-chunk-20-1.png)

Maximum Correlation Errors for X Variables by Outcome:

``` r
maxerr <- do.call(rbind, Sum1$maxerr)
rownames(maxerr) <- colnames(maxerr) <- paste("Y", 1:M, sep = "")
knitr::kable(as.data.frame(maxerr), digits = 5, booktabs = TRUE, 
  caption = "Maximum Correlation Errors for X Variables")
```

|     |       Y1|       Y2|       Y3|
|-----|--------:|--------:|--------:|
| Y1  |  0.02183|  0.02183|  0.02183|
| Y2  |  0.02183|  0.00920|  0.01444|
| Y3  |  0.02183|  0.01444|  0.00709|

### Linear model

A linear model will be fit to the data using `glm` in order to see if the slope coefficients can be recovered \[@Stats\]. First, the data is reshaped into long format using `reshape2::melt` \[@Reshape2\]. Note that since *X*<sub>*o**r**d*(1)</sub> and *X*<sub>*p**o**i**s*(1)</sub> are the same for each outcome, they will be used as factors (`id.vars`) and are only needed once.

``` r
data1 <- as.data.frame(cbind(factor(1:n), Sys1$Y, Sys1$X_all[[1]][, 1:5],
  Sys1$X_all[[2]][, c(2, 3, 5)], Sys1$X_all[[3]][, c(2, 3, 5)]))
colnames(data1)[1] <- "Subject"
data1.a <- melt(data1[, c("Subject", "ord1_1", "pois1_1", "Y1", "Y2", "Y3")], 
  id.vars = c("Subject", "ord1_1", "pois1_1"),
  measure.vars = c("Y1", "Y2", "Y3"), variable.name = "Time", value.name = "Y")
data1.b <- melt(data1[, c("Subject", "cont1_1", "cont2_1", "cont3_1")],
  id.vars = c("Subject"), variable.name = "Time", value.name = "cont1")
data1.c <- melt(data1[, c("Subject", "mix1_1", "mix2_1", "mix3_1")],
  id.vars = c("Subject"), variable.name = "Time", value.name = "mix1")
data1.d <- melt(data1[, c("Subject", "nb1_1", "nb2_1", "nb3_1")],
  id.vars = c("Subject"), variable.name = "Time", value.name = "nb1")
data1.a$Time <- data1.b$Time <- data1.c$Time <- data1.d$Time <- 
  c(rep(1, n), rep(2, n), rep(3, n))
data1 <- merge(merge(merge(data1.a, data1.b, by = c("Subject", "Time")), 
  data1.c, by = c("Subject", "Time")), data1.d, by = c("Subject", "Time"))
```

Errors *E*<sub>1</sub>, *E*<sub>2</sub>, and *E*<sub>3</sub> modeled as having Normal distributions:

``` r
fm1 <- glm(Y ~ ord1_1 + cont1 + mix1 + pois1_1 + nb1 + ord1_1:pois1_1 + 
  ord1_1:cont1 + pois1_1:cont1 + ord1_1:pois1_1:cont1 + 
  ord1_1:mix1 + pois1_1:mix1 + ord1_1:pois1_1:mix1 + 
  ord1_1:nb1 + pois1_1:nb1 + ord1_1:pois1_1:nb1 + 
  Time + ord1_1:Time + pois1_1:Time, data = data1)
summary(fm1)
#> 
#> Call:
#> glm(formula = Y ~ ord1_1 + cont1 + mix1 + pois1_1 + nb1 + ord1_1:pois1_1 + 
#>     ord1_1:cont1 + pois1_1:cont1 + ord1_1:pois1_1:cont1 + ord1_1:mix1 + 
#>     pois1_1:mix1 + ord1_1:pois1_1:mix1 + ord1_1:nb1 + pois1_1:nb1 + 
#>     ord1_1:pois1_1:nb1 + Time + ord1_1:Time + pois1_1:Time, data = data1)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -0.64658  -0.11230   0.02915   0.13954   0.32411  
#> 
#> Coefficients:
#>                        Estimate Std. Error  t value            Pr(>|t|)
#> (Intercept)          0.72192167 0.00812134    88.89 <0.0000000000000002
#> ord1_1               0.50952075 0.00712231    71.54 <0.0000000000000002
#> cont1                0.74968551 0.00141077   531.40 <0.0000000000000002
#> mix1                 1.00182443 0.00153563   652.39 <0.0000000000000002
#> pois1_1              1.25022093 0.00062176  2010.79 <0.0000000000000002
#> nb1                  1.49744169 0.00180595   829.17 <0.0000000000000002
#> Time                 1.00637627 0.00435544   231.06 <0.0000000000000002
#> ord1_1:pois1_1       0.49948318 0.00043819  1139.88 <0.0000000000000002
#> ord1_1:cont1         0.49827604 0.00110853   449.49 <0.0000000000000002
#> cont1:pois1_1        0.60002268 0.00010046  5973.03 <0.0000000000000002
#> ord1_1:mix1          0.79944395 0.00157757   506.76 <0.0000000000000002
#> mix1:pois1_1         0.89990435 0.00012033  7478.34 <0.0000000000000002
#> ord1_1:nb1           1.10077369 0.00160481   685.92 <0.0000000000000002
#> pois1_1:nb1          1.20018779 0.00013388  8964.59 <0.0000000000000002
#> ord1_1:Time          0.24969072 0.00219144   113.94 <0.0000000000000002
#> pois1_1:Time         0.49959373 0.00030845  1619.68 <0.0000000000000002
#> ord1_1:cont1:pois1_1 0.70010275 0.00006804 10289.74 <0.0000000000000002
#> ord1_1:mix1:pois1_1  1.00005780 0.00010343  9668.68 <0.0000000000000002
#> ord1_1:pois1_1:nb1   1.29993583 0.00010260 12670.37 <0.0000000000000002
#>                         
#> (Intercept)          ***
#> ord1_1               ***
#> cont1                ***
#> mix1                 ***
#> pois1_1              ***
#> nb1                  ***
#> Time                 ***
#> ord1_1:pois1_1       ***
#> ord1_1:cont1         ***
#> cont1:pois1_1        ***
#> ord1_1:mix1          ***
#> mix1:pois1_1         ***
#> ord1_1:nb1           ***
#> pois1_1:nb1          ***
#> ord1_1:Time          ***
#> pois1_1:Time         ***
#> ord1_1:cont1:pois1_1 ***
#> ord1_1:mix1:pois1_1  ***
#> ord1_1:pois1_1:nb1   ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 0.0305322)
#> 
#>     Null deviance: 3531879594.35  on 29999  degrees of freedom
#> Residual deviance:        915.39  on 29981  degrees of freedom
#> AIC: -19512
#> 
#> Number of Fisher Scoring iterations: 2
```

Each effect in the model was found to be statistically significant at the *α* = 0.001 level. Now, compare betas used in simulation to those returned by `glm`:

``` r
fm1.coef <- fm1$coefficients[c("(Intercept)", "ord1_1", "cont1", "mix1", 
  "pois1_1", "nb1", "ord1_1:pois1_1", "Time", "ord1_1:cont1", "cont1:pois1_1", 
  "ord1_1:cont1:pois1_1", "ord1_1:mix1", "mix1:pois1_1", 
  "ord1_1:mix1:pois1_1", "ord1_1:nb1", "pois1_1:nb1", 
  "ord1_1:pois1_1:nb1", "ord1_1:Time", "pois1_1:Time")]
coef <- rbind(c(betas.0, betas[[1]], betas.int[[1]], betas.t, 
  betas.subj[[1]], betas.tint[[1]]), fm1.coef)
colnames(coef) <- names(fm1.coef)
rownames(coef) <- c("Simulated", "Estimated")
knitr::kable(as.data.frame(coef[, 1:6]), digits = 3, booktabs = TRUE, 
  caption = "Beta Coefficients for Repeated Measures Model 1")
```

|           |  (Intercept)|  ord1\_1|  cont1|   mix1|  pois1\_1|    nb1|
|-----------|------------:|--------:|------:|------:|---------:|------:|
| Simulated |        0.000|     0.50|   0.75|  1.000|      1.25|  1.500|
| Estimated |        0.722|     0.51|   0.75|  1.002|      1.25|  1.497|

``` r
knitr::kable(as.data.frame(coef[, 7:12]), digits = 3, booktabs = TRUE)
```

|           |  ord1\_1:pois1\_1|   Time|  ord1\_1:cont1|  cont1:pois1\_1|  ord1\_1:cont1:pois1\_1|  ord1\_1:mix1|
|-----------|-----------------:|------:|--------------:|---------------:|-----------------------:|-------------:|
| Simulated |             0.500|  1.000|          0.500|             0.6|                     0.7|         0.800|
| Estimated |             0.499|  1.006|          0.498|             0.6|                     0.7|         0.799|

``` r
knitr::kable(as.data.frame(coef[, 13:19]), digits = 3, booktabs = TRUE)
```

|           |  mix1:pois1\_1|  ord1\_1:mix1:pois1\_1|  ord1\_1:nb1|  pois1\_1:nb1|  ord1\_1:pois1\_1:nb1|  ord1\_1:Time|  pois1\_1:Time|
|-----------|--------------:|----------------------:|------------:|-------------:|---------------------:|-------------:|--------------:|
| Simulated |            0.9|                      1|        1.100|           1.2|                   1.3|          0.25|            0.5|
| Estimated |            0.9|                      1|        1.101|           1.2|                   1.3|          0.25|            0.5|

All of the slope coefficients are estimated well except for the intercept. This could result from the non-normal error terms.

#' @title Generate Correlated Systems of Equations Containing Normal, Non-Normal, and Mixture Continuous Variables
#'
#' @description This function extends the techniques of Headrick and Beasley (2004, \doi{10.1081/SAC-120028431}) to
#'     create correlated systems of statistical equations containing continuous variables with normal,
#'     non-normal, or mixture distributions.  The method allows the user to control the distributions
#'     for the stochastic disturbance (error) terms \eqn{E} and independent variables \eqn{X}.  The user specifies the correlation structure
#'     between \eqn{X} terms within an outcome and across outcomes.  For a given equation, the user also specifies the correlation between
#'     the outcome \eqn{Y} and \eqn{X} terms.  These correlations are used to calculate the beta (slope) coefficients for the equations with
#'     \code{\link[SimRepeat]{calc_betas}}.  If the system contains mixture variables and \code{corr.yx} is specified in terms of non-mixture
#'     and mixture variables, the \code{betas} will be calculated in terms of non-mixture and mixture independent variables.  If \code{corr.yx}
#'     Finally, the user specifies the correlations across error terms.  The assumptions are that
#'     1) there are at least 2 equations and a total of at least 1 independent variable, 2) the independent variables are uncorrelated
#'     with the error terms, 3) each equation has an error term,
#'     and 4) all error terms have either a non-mixture or mixture distribution.  The outcomes \eqn{Y} are calculated as the \eqn{E} terms
#'     added to the products of the beta coefficients and the \eqn{X} terms.  There are no
#'     interactions between independent variables or distinction between subject and group-level terms (as in the hierarchical linear models approach
#'     of \code{\link[SimRepeat]{corrsys}} or \code{\link[SimRepeat]{corrsys2}}).  However, the user does not have to provide the beta coefficients
#'     (except for the intercepts).  Since the equations do not include random slopes (i.e. for the \eqn{X} terms), the effects of the independent
#'     variables are all considered "fixed."  However, a random intercept or a "time" effect with a random slope could be added by specifying
#'     them as independent variables.  There are no parameter input checks in order to decrease simulation time.  All inputs should be checked
#'     prior to simulation with \code{\link[SimRepeat]{checkpar}}.  Summaries of the simulation results can be found with \code{\link[SimRepeat]{summary_sys}}.
#'     The functions \code{\link[SimRepeat]{calc_corr_y}}, \code{\link[SimRepeat]{calc_corr_yx}}, and \code{\link[SimRepeat]{calc_corr_ye}} use equations
#'     from Headrick and Beasley (2004) to calculate the expected correlations for the outcomes, among a given outcome and covariates from
#'     the other outcomes, and for the error terms.  The vignette \strong{Theory and Equations for Correlated Systems of Continuous Variables}
#'     gives the equations, and the vignette \strong{Correlated Systems of Statistical Equations with Non-Mixture and Mixture Continuous
#'     Variables} gives examples.  There are also vignettes in \code{\link[SimCorrMix]{SimCorrMix}} which provide more details on continuous
#'     non-mixture and mixture variables.
#'
#' @section Generation of Continuous Non-Mixture and Mixture Variables:
#'     Mixture distributions describe random variables that are drawn from more than one component distribution.  For a random variable
#'     \eqn{X_{mix}} from a finite continuous mixture distribution with \eqn{k} components, the probability density function (PDF) can be
#'     described by:
#'
#'     \deqn{h_X(x) = \sum_{i=1}^{k} \pi_i f_{X_{comp_i}}(x), \sum_{i=1}^{k} \pi_i = 1.}
#'
#'     The \eqn{\pi_i} are mixing parameters which determine the weight of each component distribution \eqn{f_{X_{comp_i}}(x)} in the overall
#'     probability distribution.  As long as each component has a valid PDF, the overall distribution \eqn{h_X()} has a valid PDF.
#'     The main assumption is statistical independence between the process of randomly selecting the component distribution and the
#'     distributions themselves.  Simulation is done at the component-level for mixture variables.
#'
#'     All continuous variables are simulated using either Fleishman's third-order (\code{method} = "Fleishman",
#'     \doi{10.1007/BF02293811}) or Headrick's fifth-order (\code{method} = "Polynomial", \doi{10.1016/S0167-9473(02)00072-5})
#'     power method transformation (PMT).  It works by matching standardized cumulants -- the first four (mean,
#'     variance, skew, and standardized kurtosis) for Fleishman's method, or the first six (mean, variance,
#'     skew, standardized kurtosis, and standardized fifth and sixth cumulants) for Headrick's method.  The transformation is expressed
#'     as follows:
#'
#'     \deqn{Y = c_{0} + c_{1} * Z + c_{2} * Z^2 + c_{3} * Z^3 + c_{4} * Z^4 + c_{5} * Z^5, Z \sim N(0,1),}
#'
#'     where \eqn{c_{4}} and \eqn{c_{5}} both equal \eqn{0} for Fleishman's method.  The real constants are calculated by
#'     \code{\link[SimMultiCorrData]{find_constants}} for non-mixture and components of mixture variables.  Continuous mixture
#'     variables are generated componentwise and then transformed to the desired mixture variables using random multinomial variables
#'     generated based on mixing probabilities.  The correlation matrices are specified in terms of correlations with components of the
#'     mixture variables.
#'
#' @section Choice of Fleishman's third-order or Headrick's fifth-order method:
#'     Using the fifth-order approximation allows additional control over the fifth and sixth moments of the generated distribution, improving
#'     accuracy.  In addition, the range of feasible standardized kurtosis values, given skew and standardized fifth (\eqn{\gamma_3}) and sixth
#'     (\eqn{\gamma_4}) cumulants, is larger than with Fleishman's method (see \code{\link[SimMultiCorrData]{calc_lower_skurt}}).
#'     For example, the Fleishman method can not be used to generate a non-normal distribution with a ratio of
#'     \eqn{\gamma_3^2/\gamma_4 > 9/14} (see Headrick & Kowalchuk, 2007).  This eliminates the Chi-squared family of distributions, which has
#'     a constant ratio of \eqn{\gamma_3^2/\gamma_4 = 2/3}.  The fifth-order method also generates more distributions with valid PDF's.
#'     However, if the fifth and sixth cumulants are unknown or do not exist, the Fleishman approximation should be used.
#'
#' @section Reasons for Function Errors:
#'     1) The most likely cause for function errors is that the parameter inputs are mispecified.  Using \code{\link[SimRepeat]{checkpar}}
#'     prior to simulation can help decrease these errors.
#'
#'     2) No solutions to \code{\link[SimMultiCorrData]{fleish}} or
#'     \code{\link[SimMultiCorrData]{poly}} converged when using \code{\link[SimMultiCorrData]{find_constants}}.  If this happens,
#'     the simulation will stop.  It may help to first use \code{\link[SimMultiCorrData]{find_constants}} for each continuous variable to
#'     determine if a sixth cumulant correction value is needed.  If the standardized cumulants are obtained from \code{calc_theory}, the user
#'     may need to use rounded values as inputs (i.e. \code{skews = round(skews, 8)}).  For example, in order to ensure that skew is exactly 0 for symmetric distributions.
#'
#'     3) The kurtosis for a continuous variable may be outside the region of possible values.  There is an associated lower kurtosis boundary for
#'     associated with a given skew (for Fleishman's method) or skew and fifth and sixth cumulants (for Headrick's method).  Use
#'     \code{\link[SimMultiCorrData]{calc_lower_skurt}} to determine the boundary for a given set of cumulants.
#'
#'     4) No solutions to \code{\link[SimRepeat]{calc_betas}} converged when trying to find the beta coefficients.  Try different correlation
#'     matrices.
#'
#' @param n the sample size (i.e. the length of each simulated variable; default = 10000)
#' @param M the number of dependent variables \eqn{Y} (outcomes); equivalently, the number of equations in the system
#' @param method the PMT method used to generate all continuous variables, including independent variables (covariates) and error terms;
#'     "Fleishman" uses Fleishman's third-order polynomial transformation and "Polynomial" uses Headrick's fifth-order transformation
#' @param error_type "non_mix" if all error terms have continuous non-mixture distributions, "mix" if all error terms have continuous mixture distributions
#' @param means a list of length \code{M} of vectors of means for the non-mixture (\eqn{X_{cont}}) and mixture (\eqn{X_{mix}}) independent variables
#'     and for the error terms (\eqn{E}); the order in each vector should be: \eqn{X_{cont}, X_{mix}, E} so that the order for
#'     \eqn{X_{cont}, X_{mix}} is
#'     the same as in \code{corr.x} (considering the components of mixture variables)
#' @param vars a list of length \code{M} of vectors of variances for \eqn{X_{cont}, X_{mix}, E}; same order and dimension as \code{means}
#' @param skews a list of length \code{M} of vectors of skew values for \eqn{X_{cont}} and \eqn{E} (if \code{error_type = "non_mix"});
#'     same order as in \code{corr.x} and \code{means}
#' @param skurts a list of length \code{M} of vectors of standardized kurtoses (kurtosis - 3) for \eqn{X_{cont}} and \eqn{E} (if \code{error_type = "non_mix"});
#'     same order and dimension as \code{skews}
#' @param fifths a list of length \code{M} of vectors of standardized fifth cumulants for \eqn{X_{cont}} and \eqn{E} (if \code{error_type = "non_mix"});
#'     same order and dimension as \code{skews}; not necessary for \code{method = "Fleishman"}
#' @param sixths a list of length \code{M} of vectors of standardized sixth cumulants for \eqn{X_{cont}} and \eqn{E} (if \code{error_type = "non_mix"});
#'     same order and dimension as \code{skews}; not necessary for \code{method = "Fleishman"}
#' @param Six a list of length \code{M}, where \code{Six[[p]][[j]]} is a vector of sixth cumulant correction values to aid in finding a valid PDF for
#'     \eqn{X_{cont(pj)}}, the j-th continuous non-mixture covariate for outcome \eqn{Y_p}; the last element of \code{Six[[p]]} is for \eqn{E_p}
#'     (if \code{error_type = "non_mix"}); use \code{Six[[p]][[j]] = NULL} if no correction desired for \eqn{X_{cont(pj)}};
#'     use \code{Six[[p]] = NULL} if no correction desired for any non-mixture covariate or error term in equation p;
#'     keep \code{Six = list()} if no corrections desired for all covariates or if \code{method = "Fleishman"}
#' @param mix_pis a list of length \code{M}, where \code{mix_pis[[p]][[j]]} is a vector of mixing probabilities that sum to 1 for \eqn{X_{mix(pj)}},
#'     the j-th continuous mixture covariate for outcome \eqn{Y_p}; the last element of \code{mix_pis[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_pis[[p]] = NULL}; components should be ordered as in \code{corr.x}
#' @param mix_mus a list of length \code{M}, where \code{mix_mus[[p]][[j]]} is a vector of means of the component distributions for
#'     \eqn{X_{mix(pj)}}; the last element of \code{mix_mus[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_mus[[p]] = NULL}
#' @param mix_sigmas a list of length \code{M}, where \code{mix_sigmas[[p]][[j]]} is a vector of standard deviations of the component distributions for
#'     \eqn{X_{mix(pj)}}; the last element of \code{mix_sigmas[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_sigmas[[p]] = NULL}
#' @param mix_skews a list of length \code{M}, where \code{mix_skews[[p]][[j]]} is a vector of skew values of the component distributions for
#'     \eqn{X_{mix(pj)}}; the last element of \code{mix_skews[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_skews[[p]] = NULL}
#' @param mix_skurts a list of length \code{M}, where \code{mix_skurts[[p]][[j]]} is a vector of standardized kurtoses of the component distributions for
#'     \eqn{X_{mix(pj)}}; the last element of \code{mix_skurts[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_skurts[[p]] = NULL}
#' @param mix_fifths a list of length \code{M}, where \code{mix_fifths[[p]][[j]]} is a vector of standardized fifth cumulants of the component distributions for
#'     \eqn{X_{mix(pj)}}; the last element of \code{mix_fifths[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_fifths[[p]] = NULL}; not necessary for \code{method = "Fleishman"}
#' @param mix_sixths a list of length \code{M}, where \code{mix_sixths[[p]][[j]]} is a vector of standardized sixth cumulants of the component distributions for
#'     \eqn{X_{mix(pj)}}; the last element of \code{mix_sixths[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"});
#'     if \eqn{Y_p} has no mixture variables, use \code{mix_sixths[[p]] = NULL}; not necessary for \code{method = "Fleishman"}
#' @param mix_Six a list of length \code{M}, where \code{mix_Six[[p]]} is a list of length equal to the total number of component distributions for
#'     the \eqn{X_{mix(p)}} and \eqn{E_p} (if \code{error_type = "mix"});
#'     \code{mix_Six[[p]][[j]]} is a vector of sixth cumulant corrections for the j-th component distribution (i.e., if there are 2
#'     continuous mixture independent variables for \eqn{Y_p}, where \eqn{X_{mix(p1)}} has 2 components and \eqn{X_{mix(p2)}} has 3
#'     components, then \code{length(mix_Six[[p]]) = 5} and \code{mix_Six[[p]][[3]]} would correspond to the 1st component of
#'     \eqn{X_{mix(p2)}}); use \code{mix_Six[[p]][[j]] = NULL} if no correction desired for that component;
#'     use \code{mix_Six[[p]] = NULL} if no correction desired for any component of \eqn{X_{mix(p)}} and \eqn{E_p};
#'     keep \code{mix_Six = list()} if no corrections desired for all covariates or if \code{method = "Fleishman"}
#' @param corr.x list of length \code{M}, each component a list of length \code{M}; \code{corr.x[[p]][[q]]} is matrix of correlations
#'     for independent variables in equations p (\eqn{X_{(pj)}} for outcome \eqn{Y_p}) and q (\eqn{X_{(qj)}} for outcome \eqn{Y_q});
#'     order: 1st continuous non-mixture (same order as in \code{skews}) and 2nd components of continuous mixture (same order as in \code{mix_pis});
#'     if p = q, \code{corr.x[[p]][[q]]} is a correlation matrix with \code{nrow(corr.x[[p]][[q]])} = # of non-mixture + # of mixture components for outcome \eqn{Y_p};
#'     if p != q, \code{corr.x[[p]][[q]]} is a non-symmetric matrix of correlations where rows correspond to covariates for \eqn{Y_p}
#'     so that \code{nrow(corr.x[[p]][[q]])} = # of non-mixture + # of mixture components for outcome \eqn{Y_p} and
#'     columns correspond to covariates for \eqn{Y_q} so that \code{ncol(corr.x[[p]][[q]])} = # of non-mixture + # of mixture components for outcome \eqn{Y_q};
#'     use \code{corr.x[[p]][[q]] = NULL} if equation q has no \eqn{X_{(qj)}}; use \code{corr.x[[p]] = NULL} if equation p has no \eqn{X_{(pj)}}
#' @param corr.yx a list of length \code{M}, where the p-th component is a 1 row matrix of correlations between \eqn{Y_p} and \eqn{X_{(pj)}};
#'     if there are mixture variables and the \code{betas} are desired in terms of these (and not the components), then \code{corr.yx}
#'     should be specified in terms of correlations between outcomes and non-mixture or mixture variables, and the number of columns of the matrices
#'     of \code{corr.yx} should not match the dimensions of the matrices in \code{corr.x}; if the \code{betas} are desired in terms of
#'     the components, then \code{corr.yx} should be specified in terms of correlations between outcomes and non-mixture or components of
#'     mixture variables, and the number of columns of the matrices of \code{corr.yx} should match the dimensions of the matrices in \code{corr.x};
#'     use \code{corr.yx[[p]] = NULL} if equation p has no \eqn{X_{(pj)}}
#' @param corr.e correlation matrix for continuous non-mixture or components of mixture error terms
#' @param same.var either a vector or a matrix; if a vector, \code{same.var} includes column numbers of \code{corr.x[[1]][[1]]}
#'     corresponding to independent variables that should be identical across equations; these terms must have the same indices for all \code{p = 1, ..., M};
#'     i.e., if the 1st variable represents height which should be the same for each equation, then
#'     \code{same.var[1] = 1} and the 1st term for all other outcomes must also be height;
#'     if a matrix, columns 1 and 2 are outcome p and column index in \code{corr.x[[p]][[p]]} for 1st instance of variable,
#'     columns 3 and 4 are outcome q and column index in \code{corr.x[[q]][[q]]} for subsequent instances of variable; i.e., if
#'     1st term for all outcomes is height and \code{M = 3}, then \code{same.var = matrix(c(1, 1, 2, 1, 1, 1, 3, 1), 2, 4, byrow = TRUE)}; the
#'     independent variable index corresponds to continuous non-mixture and component of continuous mixture covariate
#' @param betas.0 vector of length \code{M} containing intercepts, if \code{NULL} all set equal to 0; if length 1, all intercepts set to \code{betas.0}
#' @param seed the seed value for random number generation (default = 1234)
#' @param use.nearPD TRUE to convert the overall intermediate correlation matrix formed by the \eqn{X} (for all outcomes and independent
#'     variables) or \eqn{E} to the nearest positive definite matrix with \code{Matrix::nearPD} if necessary; if FALSE the negative
#'     eigenvalues are replaced with 0 if necessary
#' @param errorloop if TRUE, uses \code{\link[SimCorrMix]{corr_error}} to attempt to correct the correlation of the independent
#'     variables within and across outcomes to be within \code{epsilon} of the target correlations \code{corr.x} until the number of iterations
#'     reaches \code{maxit} (default = FALSE)
#' @param epsilon the maximum acceptable error between the final and target correlation matrices (default = 0.001) in the error loop
#' @param maxit the maximum number of iterations to use (default = 1000) in the error loop
#' @importFrom psych describe
#' @import SimMultiCorrData
#' @import SimCorrMix
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @importFrom Matrix nearPD
#' @importFrom VGAM qzipois qzinegbin
#' @export
#' @return A list with the following components:
#'
#' @return \code{Y} matrix with \code{n} rows and \code{M} columns of outcomes
#'
#' @return \code{X} list of length \code{M} containing \eqn{X_{cont(pj)}, X_{comp(pj)}}
#'
#' @return \code{X_all} list of length \code{M} containing \eqn{X_{cont(pj)}, X_{mix(pj)}}
#'
#' @return \code{E} matrix with \code{n} rows containing continuous non-mixture or components of continuous mixture error terms
#'
#' @return \code{E_mix} matrix with \code{n} rows containing continuous mixture error terms
#'
#' @return \code{betas} a matrix of the slope coefficients calculated with \code{\link[SimRepeat]{calc_betas}}, rows represent the
#'     outcomes
#'
#' @return \code{constants} a list of length \code{M} with data.frames of the constants for the \eqn{X_{cont(pj)}},
#'     \eqn{X_comp(pj)} and \eqn{E_p}
#'
#' @return \code{SixCorr} a list of length \code{M} of lists of sixth cumulant correction values used to obtain valid
#'     \emph{pdf}'s for the \eqn{X_{cont(pj)}}, \eqn{X_comp(pj)}, and \eqn{E_p}
#'
#' @return \code{valid.pdf} a list of length \code{M} of vectors where the i-th element is "TRUE" if the constants for the i-th
#'         continuous variable generate a valid pdf, else "FALSE"
#'
#' @return \code{Sigma.X} matrix of intermediate correlations applied to generate \eqn{Z_{cont(pj)}, Z_{comp(pj)}};
#'     these are the normal variables transformed to get the desired distributions
#'
#' @return \code{Error_Time} the time in minutes required to use the error loop
#'
#' @return \code{Time} the total simulation time in minutes
#'
#' @return \code{niter} a matrix of the number of iterations used in the error loop
#'
#' @keywords continuous mixture Headrick Beasley
#' @seealso \code{\link[SimRepeat]{calc_betas}}, \code{\link[SimRepeat]{calc_corr_y}}, \code{\link[SimRepeat]{calc_corr_yx}},
#'     \code{\link[SimRepeat]{calc_corr_ye}}, \code{\link[SimRepeat]{checkpar}}, \code{\link[SimRepeat]{summary_sys}}
#' @references
#' Davenport JW, Bezder JC, & Hathaway RJ (1988). Parameter Estimation for Finite Mixture Distributions.
#'     Computers & Mathematics with Applications, 15(10):819-28.
#'
#' Everitt BS (1996). An Introduction to Finite Mixture Distributions. Statistical Methods in Medical Research, 5(2):107-127. \doi{10.1177/096228029600500202}.
#'
#' Fialkowski AC (2017). SimMultiCorrData: Simulation of Correlated Data with Multiple Variable Types. R package version 0.2.1.
#'     \url{https://CRAN.R-project.org/package=SimMultiCorrData}.
#'
#' Fialkowski AC (2018). SimCorrMix: Simulation of Correlated Data of Multiple Variable Types including Continuous and Count
#'     Mixture Distributions. R package version 0.1.0. \url{https://github.com/AFialkowski/SimCorrMix}
#'
#' Fleishman AI (1978). A Method for Simulating Non-normal Distributions. Psychometrika, 43:521-532. \doi{10.1007/BF02293811}.
#'
#' Headrick TC (2002). Fast Fifth-order Polynomial Transforms for Generating Univariate and Multivariate
#'     Non-normal Distributions. Computational Statistics & Data Analysis, 40(4):685-711. \doi{10.1016/S0167-9473(02)00072-5}.
#'     (\href{http://www.sciencedirect.com/science/article/pii/S0167947302000725}{ScienceDirect})
#'
#' Headrick TC (2004). On Polynomial Transformations for Simulating Multivariate Nonnormal Distributions.
#'     Journal of Modern Applied Statistical Methods, 3(1):65-71. \doi{10.22237/jmasm/1083370080}.
#'
#' Headrick TC, Beasley TM (2004).  A Method for Simulating Correlated Non-Normal Systems of Linear Statistical Equations.
#'     Communications in Statistics - Simulation and Computation, 33(1).  \doi{10.1081/SAC-120028431}
#'
#' Headrick TC, Kowalchuk RK (2007). The Power Method Transformation: Its Probability Density Function, Distribution
#'     Function, and Its Further Use for Fitting Data. Journal of Statistical Computation and Simulation, 77:229-249. \doi{10.1080/10629360600605065}.
#'
#' Headrick TC, Sawilowsky SS (1999). Simulating Correlated Non-normal Distributions: Extending the Fleishman Power
#'     Method. Psychometrika, 64:25-35. \doi{10.1007/BF02294317}.
#'
#' Headrick TC, Sheng Y, & Hodis FA (2007). Numerical Computing and Graphics for the Power Method Transformation Using
#'     Mathematica. Journal of Statistical Software, 19(3):1 - 17. \cr \doi{10.18637/jss.v019.i03}.
#'
#' Higham N (2002). Computing the nearest correlation matrix - a problem from finance; IMA Journal of Numerical Analysis 22:329-343.
#'
#' McCulloch CE, Searle SR, Neuhaus JM (2008). \emph{Generalized, Linear, and Mixed Models} (2nd ed.). Wiley Series in Probability and
#'     Statistics. Hoboken, New Jersey: John Wiley & Sons, Inc.
#'
#' Pearson RK (2011). Exploring Data in Engineering, the Sciences, and Medicine. In. New York: Oxford University Press.
#'
#' Schork NJ, Allison DB, & Thiel B (1996). Mixture Distributions in Human Genetics Research. Statistical Methods in Medical Research,
#'     5:155-178. \doi{10.1177/096228029600500204}.
#'
#' Vale CD & Maurelli VA (1983). Simulating Multivariate Nonnormal Distributions. Psychometrika, 48:465-471. \doi{10.1007/BF02293687}.
#'
#' @examples \dontrun{
#' # Example: system of three equations for 2 independent variables, where each
#' # error term has unit variance, from Headrick & Beasley (2002)
#' # Y_1 = beta_10 + beta_11 * X_11 + beta_12 * X_12 + sigma_1 * e_1
#' # Y_2 = beta_20 + beta_21 * X_21 + beta_22 * X_22 + sigma_2 * e_2
#' # Y_3 = beta_30 + beta_31 * X_31 + beta_32 * X_32 + sigma_3 * e_3
#'
#' # X_11 = X_21 = X_31 = Exponential(2)
#' # X_12 = X_22 = X_32 = Laplace(0, 1)
#' # e_1 = e_2 = e_3 = Cauchy(0, 1)
#'
#' seed <- 1234
#' M <- 3
#' Stcum1 <- calc_theory("Exponential", 2)
#' Stcum2 <- calc_theory("Laplace", c(0, 1))
#' Stcum3 <- c(0, 1, 0, 25, 0, 1500) # taken from paper
#' means <- lapply(seq_len(M), function(x) c(0, 0, 0))
#' vars <- lapply(seq_len(M), function(x) c(1, 1, 1))
#' skews <- lapply(seq_len(M), function(x) c(Stcum1[3], Stcum2[3], Stcum3[3]))
#' skurts <- lapply(seq_len(M), function(x) c(Stcum1[4], Stcum2[4], Stcum3[4]))
#' fifths <- lapply(seq_len(M), function(x) c(Stcum1[5], Stcum2[5], Stcum3[5]))
#' sixths <- lapply(seq_len(M), function(x) c(Stcum1[6], Stcum2[6], Stcum3[6]))
#'
#' # No sixth cumulant corrections will be used in order to match the results
#' # from the paper.  Otherwise, the following should be used in order to
#' # produce variables with valid PDF's:
#' # Six <- lapply(seq_len(M), function(x) list(NULL, 25.14, NULL))
#'
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
#'
#' # Check the parameter inputs
#' checkpar(M, "Polynomial", "non_mix", means, vars, skews,
#'   skurts, fifths, sixths, corr.x = corr.x, corr.yx = corr.yx,
#'   corr.e = corr.e)
#' # Generate the system
#' Sys1 <- nonnormsys(10000, M, "Polynomial", "non_mix", means, vars, skews,
#'   skurts, fifths, sixths, corr.x = corr.x, corr.yx = corr.yx,
#'   corr.e = corr.e, seed = seed)
#' # Summarize the results
#' Sum1 <- summary_sys(Sys1$Y, Sys1$E, E_mix = NULL, Sys1$X, X_all = list(), M,
#'   "Polynomial", means, vars, skews, skurts, fifths, sixths, corr.x = corr.x,
#'   corr.e = corr.e)
#'
#' # Calculate theoretical correlations for comparison to simulated values
#' calc_corr_y(Sys1$betas, corr.x, corr.e, vars)
#' Sum1$rho.y
#' calc_corr_ye(Sys1$betas, corr.x, corr.e, vars)
#' Sum1$rho.ye
#' calc_corr_yx(Sys1$betas, corr.x, vars)
#' Sum1$rho.yx
#' }
nonnormsys <- function(n = 10000, M = NULL,
                       method = c("Fleishman", "Polynomial"),
                       error_type = c("non_mix", "mix"),
                       means = list(), vars = list(), skews = list(),
                       skurts = list(), fifths = list(), sixths = list(),
                       Six = list(), mix_pis = list(), mix_mus = list(),
                       mix_sigmas = list(), mix_skews =  list(),
                       mix_skurts =  list(), mix_fifths =  list(),
                       mix_sixths =  list(), mix_Six = list(),
                       same.var = NULL, betas.0 = NULL, corr.x = list(),
                       corr.yx = list(), corr.e = NULL, seed = 1234,
                       use.nearPD = TRUE, errorloop = FALSE, epsilon = 0.001,
                       maxit = 1000) {
  start.time <- Sys.time()
  if (length(error_type) != 1)
    stop("Please choose one type of distribution for all of the error terms:
         mix if all errors have continuous mixture distributions,
         non_mix if all errors have continuous non-mixture distributions.")
  K.mix <- rep(0, M)
  K.comp <- rep(0, M)
  K.error <- rep(0, M)
  if (length(mix_pis) > 0) {
    K.mix <- lengths(mix_pis)
    K.comp <- sapply(lapply(mix_pis, unlist), length)
  }
  k.comp <- c(0, cumsum(K.comp))
  K.cont <- lengths(vars) - K.mix
  k.cont <- c(0, cumsum(K.cont))
  k.mix <- c(0, cumsum(K.mix))
  if (error_type == "mix") {
    K.error <- sapply(mix_pis, function(x) lengths(x[length(x)]))
    k.error <- c(0, cumsum(K.error))
    K.x <- K.cont + K.comp - K.error
    K.mix2 <- K.mix - 1
  } else {
    K.x <- K.cont + K.comp - 1
    K.mix2 <- K.mix
  }
  if (is.null(betas.0)) betas.0 <- rep(0, M)
  if (length(betas.0) == 1) betas.0 <- rep(betas.0, M)
  if (length(betas.0) != M)
    stop("Length of intercept vector betas.0 should be M!")
  same.cont <- NULL
  if (!is.null(same.var)) {
    if (class(same.var) == "numeric") {
      for (j in 1:length(same.var)) {
        same.cont <- rbind(same.cont, cbind(rep(1, M - 1),
          rep(same.var[j], M - 1), 2:M, rep(same.var[j], M - 1)))
      }
      same.var <- same.cont
    } else same.cont <- same.var
  }
  csame.dist <- NULL
  csame.dist2 <- NULL
  skews2 <- NULL
  skurts2 <- NULL
  fifths2 <- NULL
  sixths2 <- NULL
  Six2 <- list()
  if (length(skews) > 0) {
    skews2 <- unlist(skews)
    skurts2 <- unlist(skurts)
    if (method == "Polynomial") {
      fifths2 <- unlist(fifths)
      sixths2 <- unlist(sixths)
      Six2 <- unlist(Six, recursive = FALSE)
    }
    if (length(skews2) > 1) {
      for (i in 2:length(skews2)) {
        if (skews2[i] %in% skews2[1:(i - 1)]) {
          csame <- which(skews2[1:(i - 1)] == skews2[i])
          for (j in 1:length(csame)) {
            if (method == "Polynomial") {
              if ((skurts2[i] == skurts2[csame[j]]) &
                  (fifths2[i] == fifths2[csame[j]]) &
                  (sixths2[i] == sixths2[csame[j]])) {
                csame.dist <- rbind(csame.dist, c(csame[j], i))
                break
              }
            }
            if (method == "Fleishman") {
              if (skurts2[i] == skurts2[csame[j]]) {
                csame.dist <- rbind(csame.dist, c(csame[j], i))
                break
              }
            }
          }
        }
      }
      if (!is.null(csame.dist)) {
        ind.p <- sapply(csame.dist[, 1],
                        function(x) min(which(k.cont >= x)) - 1)
        ind.q <- sapply(csame.dist[, 2],
                        function(x) min(which(k.cont >= x)) - 1)
        ind.p2 <- numeric(length(ind.p))
        ind.q2 <- numeric(length(ind.q))
        for (i in 1:length(ind.p)) {
          ind.p2[i] <- csame.dist[i, 1] - k.cont[ind.p[i]]
          ind.q2[i] <- csame.dist[i, 2] - k.cont[ind.q[i]]
        }
        csame.dist2 <- cbind(ind.p, ind.p2, ind.q, ind.q2)
      }
    }
  }
  msame.dist2 <- NULL
  if (length(mix_pis) > 0) {
    mix_skews2 <- unlist(mix_skews)
    mix_skurts2 <- unlist(mix_skurts)
    if (method == "Polynomial") {
      mix_fifths2 <- unlist(mix_fifths)
      mix_sixths2 <- unlist(mix_sixths)
    }
    for (i in 1:length(mix_skews2)) {
      msame.d <- NULL
      if (length(skews) > 0) {
        if (mix_skews2[i] %in% skews2) {
          msame <- which(skews2 == mix_skews2[i])
          for (j in 1:length(msame)) {
            if (method == "Polynomial") {
              if ((mix_skurts2[i] == skurts2[msame[j]]) &
                  (mix_fifths2[i] == fifths2[msame[j]]) &
                  (mix_sixths2[i] == sixths2[msame[j]])) {
                msame.d <- c(msame[j], i)
                msame.d <- c(min(which(k.cont >= msame.d[1])) - 1,
                  msame.d[1] - k.cont[min(which(k.cont >= msame.d[1])) - 1],
                  min(which(k.comp >= msame.d[2])) - 1,
                  msame.d[2] - k.comp[min(which(k.comp >= msame.d[2])) - 1])
                msame.d[4] <- msame.d[4] + K.cont[msame.d[3]]
                break
              }
            }
            if (method == "Fleishman") {
              if (mix_skurts2[i] == skurts2[msame[j]]) {
                msame.d <- c(msame[j], i)
                msame.d <- c(min(which(k.cont >= msame.d[1])) - 1,
                  msame.d[1] - k.cont[min(which(k.cont >= msame.d[1])) - 1],
                  min(which(k.comp >= msame.d[2])) - 1,
                  msame.d[2] - k.comp[min(which(k.comp >= msame.d[2])) - 1])
                msame.d[4] <- msame.d[4] + K.cont[msame.d[3]]
                break
              }
            }
          }
        }
      }
      if (is.null(msame.d) & i >= 2) {
        if (mix_skews2[i] %in% mix_skews2[1:(i-1)]) {
          msame <- which(mix_skews2[1:(i - 1)] == mix_skews2[i])
          for (j in 1:length(msame)) {
            if (method == "Polynomial") {
              if ((mix_skurts2[i] == mix_skurts2[msame[j]]) &
                  (mix_fifths2[i] == mix_fifths2[msame[j]]) &
                  (mix_sixths2[i] == mix_sixths2[msame[j]])) {
                msame.d <- c(msame[j], i)
                msame.d <- c(min(which(k.comp >= msame.d[1])) - 1,
                  msame.d[1] - k.comp[min(which(k.comp >= msame.d[1])) - 1],
                  min(which(k.comp >= msame.d[2])) - 1,
                  msame.d[2] - k.comp[min(which(k.comp >= msame.d[2])) - 1])
                msame.d[2] <- msame.d[2] + K.cont[msame.d[1]]
                msame.d[4] <- msame.d[4] + K.cont[msame.d[3]]
                break
              }
            }
            if (method == "Fleishman") {
              if (mix_skurts2[i] == mix_skurts2[msame[j]]) {
                msame.d <- c(msame[j], i)
                msame.d <- c(min(which(k.comp >= msame.d[1])) - 1,
                  msame.d[1] - k.comp[min(which(k.comp >= msame.d[1])) - 1],
                  min(which(k.comp >= msame.d[2])) - 1,
                  msame.d[2] - k.comp[min(which(k.comp >= msame.d[2])) - 1])
                msame.d[2] <- msame.d[2] + K.cont[msame.d[1]]
                msame.d[4] <- msame.d[4] + K.cont[msame.d[3]]
                break
              }
            }
          }
        }
      }
      msame.dist2 <- rbind(msame.dist2, msame.d)
    }
    mix_skews2 <- lapply(mix_skews, unlist)
    mix_skurts2 <- lapply(mix_skurts, unlist)
    if (method == "Polynomial") {
      mix_fifths2 <- lapply(mix_fifths, unlist)
      mix_sixths2 <- lapply(mix_sixths, unlist)
    }
  }
  constants <- list()
  SixCorr <- list()
  Valid.PDF <- list()
  for (i in 1:M) {
    SixCorr[[i]] <- numeric(K.cont[i] + K.comp[i])
    Valid.PDF[[i]] <- numeric(K.cont[i] + K.comp[i])
    if (method == "Fleishman") {
      constants[[i]] <- matrix(NA, nrow = K.cont[i] + K.comp[i], ncol = 4)
      colnames(constants[[i]]) <- c("c0", "c1", "c2", "c3")
    }
    if (method == "Polynomial") {
      constants[[i]] <- matrix(NA, nrow = K.cont[i] + K.comp[i], ncol = 6)
      colnames(constants[[i]]) <- c("c0", "c1", "c2", "c3", "c4", "c5")
    }
    if (K.cont[i] > 0) {
      for (j in 1:K.cont[i]) {
        if (!is.null(csame.dist2)) {
          rind <- which(csame.dist2[, 3] == i & csame.dist2[, 4] == j)
          if (length(rind) > 0) {
            constants[[i]][j, ] <-
              constants[[csame.dist2[rind, 1]]][csame.dist2[rind, 2], ]
            SixCorr[[i]][j] <-
              SixCorr[[csame.dist2[rind, 1]]][csame.dist2[rind, 2]]
            Valid.PDF[[i]][j] <-
              Valid.PDF[[csame.dist2[rind, 1]]][csame.dist2[rind, 2]]
          }
        }
        if (sum(is.na(constants[[i]][j, ])) > 0) {
          if (length(Six) == 0) Six2 <- NULL else
            if (length(Six[[i]]) == 0) Six2 <- NULL else
              Six2 <- Six[[i]][[j]]
          cons <-
            suppressWarnings(find_constants(method, skews[[i]][j],
              skurts[[i]][j], fifths[[i]][j], sixths[[i]][j],
              Six = Six2, n = 25, seed = seed))
          if (length(cons) == 1 | is.null(cons)) {
            stop(paste("Constants can not be found for outcome ", i,
                       " non-mixture variable ", j, ".", sep = ""))
          }
          SixCorr[[i]][j] <- ifelse(is.null(cons$SixCorr1), NA, cons$SixCorr1)
          Valid.PDF[[i]][j] <- cons$valid
          constants[[i]][j, ] <- cons$constants
        }
      }
    }
  }
  for (i in 1:M) {
    if (K.mix[i] > 0) {
      for (j in (K.cont[i] + 1):(K.cont[i] + K.comp[i])) {
        if (!is.null(msame.dist2)) {
          rind <- which(msame.dist2[, 3] == i & msame.dist2[, 4] == j)
          if (length(rind) > 0) {
            constants[[i]][j, ] <-
              constants[[msame.dist2[rind, 1]]][msame.dist2[rind, 2], ]
            SixCorr[[i]][j] <-
              SixCorr[[msame.dist2[rind, 1]]][msame.dist2[rind, 2]]
            Valid.PDF[[i]][j] <-
              Valid.PDF[[msame.dist2[rind, 1]]][msame.dist2[rind, 2]]
          }
        }
        if (sum(is.na(constants[[i]][j, ])) > 0) {
          if (length(mix_Six) == 0) mix_Six2 <- NULL else
            if (length(mix_Six[[i]]) == 0) mix_Six2 <- NULL else
              mix_Six2 <- mix_Six[[i]][[j - K.cont[i]]]
          cons <-
            suppressWarnings(find_constants(method,
              mix_skews2[[i]][j - K.cont[i]], mix_skurts2[[i]][j - K.cont[i]],
              mix_fifths2[[i]][j - K.cont[i]], mix_sixths2[[i]][j - K.cont[i]],
              Six = mix_Six2, n = 25, seed = seed))
          if (length(cons) == 1 | is.null(cons)) {
            stop(paste("Constants can not be found for outcome ", i,
                 " mixture component ", j - K.cont[i], ".", sep = ""))
          }
          SixCorr[[i]][j] <- ifelse(is.null(cons$SixCorr1), NA, cons$SixCorr1)
          Valid.PDF[[i]][j] <- cons$valid
          constants[[i]][j, ] <- cons$constants
        }
      }
    }
  }
  if (error_type == "non_mix" & max(K.mix2) > 0) {
    for (i in 1:M) {
      constants[[i]] <- rbind(constants[[i]][-K.cont[i], ],
        constants[[i]][K.cont[i], ])
    }
  }
  e.constants <- NULL
  if (error_type == "non_mix") {
    e.means <- mapply('[[', means, lengths(means))
    e.vars <- mapply('[[', vars, lengths(vars))
    for (i in 1:M) {
      e.constants <- rbind(e.constants, constants[[i]][nrow(constants[[i]]), ])
    }
  } else {
    e.means <- unlist(mapply('[', mix_mus, lengths(mix_mus)))
    e.vars <- (unlist(mapply('[', mix_sigmas, lengths(mix_sigmas))))^2
    for (i in 1:M) {
      e.constants <- rbind(e.constants,
        constants[[i]][(nrow(constants[[i]]) - K.error[i] +
                          1):nrow(constants[[i]]), ])
    }
  }
  Sigma_E <- intercorr_cont(method, e.constants, corr.e)
  set.seed(seed)
  if (!is.null(same.var)) {
    Z <- matrix(rnorm((ncol(Sigma_E) + sum(K.x) - nrow(same.var)) * n), n)
  } else Z <- matrix(rnorm((ncol(Sigma_E) + sum(K.x)) * n), n)
  Z <- scale(Z, TRUE, FALSE)
  Z <- Z %*% svd(Z, nu = 0)$v
  Z <- scale(Z, FALSE, TRUE)
  if (min(eigen(Sigma_E, symmetric = TRUE)$values) < 0) {
    if (use.nearPD == TRUE) {
      message("Intermediate E correlation matrix is not positive definite.
Nearest positive definite matrix is used.")
      Sigma_E <- as.matrix(nearPD(Sigma_E, corr = T, keepDiag = T)$mat)
    } else {
      message("Intermediate E correlation matrix is not positive definite.
Negative eigenvalues are replaced with 0.  Set use.nearPD = TRUE to use nearest
positive-definite matrix instead.")
    }
  }
  eig <- eigen(Sigma_E, symmetric = TRUE)
  sqrteigval <- diag(sqrt(pmax(eig$values, 0)), nrow(Sigma_E))
  eigvec <- eig$vectors
  fry <- eigvec %*% sqrteigval
  E <- fry %*% t(Z[, (ncol(Z) - ncol(Sigma_E) + 1):ncol(Z)])
  E <- t(E)
  E2 <- matrix(1, nrow = n, ncol = ncol(Sigma_E))
  for (i in 1:ncol(Sigma_E)) {
    if (method == "Fleishman") {
      E2[, i] <- e.constants[i, 1] + e.constants[i, 2] * E[, i] +
        e.constants[i, 3] * E[, i]^2 + e.constants[i, 4] * E[, i]^3
    }
    if (method == "Polynomial") {
      E2[, i] <- e.constants[i, 1] + e.constants[i, 2] * E[, i] +
        e.constants[i, 3] * E[, i]^2 + e.constants[i, 4] * E[, i]^3 +
        e.constants[i, 5] * E[, i]^4 + e.constants[i, 6] * E[, i]^5
    }
    E2[, i] <- e.means[i] + sqrt(e.vars[i]) * E2[, i]
  }
  E_comp <- E2
  colnames(E_comp) <- paste("E", 1:ncol(E_comp), sep = "")
  if (error_type == "mix") {
    E2 <- matrix(1, n, M)
    for (i in 1:M) {
      seed <- seed + 1
      set.seed(seed)
      R <- rmultinom(n, size = 1, prob = mix_pis[[i]][[length(mix_pis[[i]])]])
      E2[, i] <- apply(t(R) * E_comp[, (k.error[i] + 1):k.error[i + 1]], 1,
                       sum)
      E2[, i] <- scale(E2[, i])
      E2[, i] <- means[[i]][length(means[[i]])] +
        sqrt(vars[[i]][length(vars[[i]])]) * E2[, i]
    }
  }
  colnames(E2) <- paste("E", 1:M, sep = "")
  constants0 <- constants
  constants2 <- NULL
  for (p in 1:M) {
    if (K.x[p] == 0) next
    rownames(corr.x[[p]][[p]]) <- paste("cont", p, "_", 1:K.x[p], sep = "")
    if (!is.null(same.cont)) {
      if (p %in% same.cont[, 3])
        constants[[p]] <-
          constants[[p]][-(same.cont[same.cont[, 3] == p, 4]), , drop = FALSE]
    }
    if (error_type == "non_mix") {
      constants2 <- rbind(constants2,
        constants[[p]][-nrow(constants[[p]]), , drop = FALSE])
    } else {
      constants2 <- rbind(constants2,
        constants[[p]][-c((nrow(constants[[p]]) - K.error[p] +
                             1):nrow(constants[[p]])), , drop = FALSE])
    }
    if (!is.null(same.var)) {
      same.var0 <- same.var[same.var[, 3] == p, , drop = FALSE]
      if (p != 1 & nrow(same.var0) != 0) {
        q <- unique(same.var0[, 1])
        for (i in 1:length(q)) {
          rownames(corr.x[[p]][[p]])[same.var0[same.var0[, 1] == q[i], 4]] <-
            rownames(corr.x[[q[i]]][[q[i]]])[same.var0[same.var0[, 1] == q[i],
                                                       2]]
        }
      }
    }
    colnames(corr.x[[p]][[p]]) <- rownames(corr.x[[p]][[p]])
  }
  names0 <- list()
  for (p in 1:M) {
    if (length(corr.x[[p]]) == 0) {
      names0 <- append(names0, list(NULL))
      next
    } else {
      names0[[p]] <- colnames(corr.x[[p]][[p]])
      for (q in 1:M) {
        if (length(corr.x[[p]][[q]]) == 0) next
        colnames(corr.x[[p]][[q]]) <- colnames(corr.x[[q]][[q]])
        rownames(corr.x[[p]][[q]]) <- rownames(corr.x[[p]][[p]])
      }
    }
  }
  corr.x0 <- corr.x
  Corr_X <- list()
  if (!is.null(same.var)) {
    for (p in 1:(M - 1)) {
      if (K.x[p] == 0) next
      for (q in (p + 1):M) {
        if (K.x[q] == 0) next
        same.var0 <- same.var[(same.var[, 1] == p & same.var[, 3] == q), ,
                              drop = FALSE]
        if (nrow(same.var0) != 0) {
          for (i in 1:M) {
            if (K.x[i] == 0) next
            for (j in 1:M) {
              if (K.x[j] == 0) next
              if (i == q & j != q)
                corr.x0[[i]][[j]] <- corr.x0[[i]][[j]][-same.var0[, 4], ,
                                                       drop = FALSE]
              if (i != q & j == q)
                corr.x0[[i]][[j]] <- corr.x0[[i]][[j]][, -same.var0[, 4],
                                                       drop = FALSE]
              if (i == q & j == q)
                corr.x0[[i]][[j]] <- corr.x0[[i]][[j]][-same.var0[, 4],
                  -same.var0[, 4], drop = FALSE]
            }
          }
        }
      }
    }
  }
  for (p in 1:M) {
    Corr_X <- append(Corr_X, list(NULL))
    if (length(corr.x0[[p]]) > 0) {
      for (q in 1:M) {
        if (length(corr.x0[[p]][[q]]) == 0) next else
          if (nrow(corr.x0[[p]][[q]]) == 0 |
              ncol(corr.x0[[p]][[q]]) == 0) next else
                Corr_X[[p]] <- cbind(Corr_X[[p]], corr.x0[[p]][[q]])
      }
    }
  }
  Corr_X <- do.call(rbind, Corr_X)
  names1 <- colnames(Corr_X)
  Sigma_X <- intercorr_cont(method = method, constants = constants2,
                            rho_cont = Corr_X)
  if (min(eigen(Sigma_X, symmetric = TRUE)$values) < 0) {
    if (use.nearPD == TRUE) {
      message("Intermediate correlation matrix is not positive definite.
Nearest positive definite matrix is used.")
      Sigma_X <- as.matrix(nearPD(Sigma_X, corr = T, keepDiag = T)$mat)
    } else {
      message("Intermediate correlation matrix is not positive definite.
Negative eigenvalues are replaced with 0.  Set use.nearPD = TRUE to use nearest
positive-definite matrix instead.")
    }
  }
  eig <- eigen(Sigma_X, symmetric = TRUE)
  sqrteigval <- diag(sqrt(pmax(eig$values, 0)), nrow(Sigma_X))
  eigvec <- eig$vectors
  fry <- eigvec %*% sqrteigval
  X <- fry %*% t(Z[, 1:ncol(Sigma_X), drop = FALSE])
  X <- t(X)
  colnames(X) <- colnames(Corr_X)
  colnames(Sigma_X) <- colnames(Corr_X)
  rownames(Sigma_X) <- rownames(Corr_X)
  X_cont <- list()
  for (i in 1:M) {
    if (K.x[i] != 0) {
      if (i == 1) {
        X_cont[[i]] <- X[, 1:K.x[i], drop = FALSE]
      } else {
        if (!is.null(same.var)) {
          if (K.x[i] == nrow(same.var[same.var[, 3] == i, , drop = FALSE])) {
            X_cont <- append(X_cont, list(NULL))
          } else {
            X_cont[[i]] <- X[, (cumsum(K.x)[i - 1] -
              nrow(same.var[same.var[, 3] <= (i - 1), , drop = FALSE]) +
              1):(cumsum(K.x)[i] - nrow(same.var[same.var[, 3] <= i, ,
                                               drop = FALSE])), drop = FALSE]
          }
        } else {
          X_cont[[i]] <- X[, (cumsum(K.x)[i - 1] + 1):cumsum(K.x)[i],
                       drop = FALSE]
        }
      }
    } else X_cont <- append(X_cont, list(NULL))
  }
  Y_cont <- list()
  means2 <- list()
  vars2 <- list()
  temp_means <- NULL
  temp_vars <- NULL
  for (m in 1:M) {
    Y_cont <- append(Y_cont, list(NULL))
    means2 <- append(means2, list(NULL))
    vars2 <- append(vars2, list(NULL))
    if (K.x[m] == 0 | length(grep("cont", colnames(X_cont[[m]]))) == 0)
      next else {
      if (error_type == "non_mix") {
        if ((K.cont[m] - 1) > 0) {
          means2[[m]] <- means[[m]][1:(K.cont[m] - 1)]
          vars2[[m]] <- vars[[m]][1:(K.cont[m] - 1)]
        }
        if (K.mix2[m] > 0) {
          means2[[m]] <- c(means2[[m]], unlist(mix_mus[[m]]))
          vars2[[m]] <- c(vars2[[m]], (unlist(mix_sigmas[[m]]))^2)
        }
      }
      if (error_type == "mix") {
        if (K.cont[m] > 0) {
          means2[[m]] <- means[[m]][1:K.cont[m]]
          vars2[[m]] <- vars[[m]][1:K.cont[m]]
        }
        if (K.mix2[m] > 0) {
          means2[[m]] <- c(means2[[m]],
            unlist(mix_mus[[m]][-length(mix_mus[[m]])]))
          vars2[[m]] <- c(vars2[[m]],
            (unlist(mix_sigmas[[m]][-length(mix_sigmas[[m]])]))^2)
        }
      }
      if (!is.null(same.cont)) {
        if (m %in% same.cont[, 3]) {
          means2[[m]] <- means2[[m]][-(same.cont[same.cont[, 3] == m, 4])]
          vars2[[m]] <- vars2[[m]][-(same.cont[same.cont[, 3] == m, 4])]
        }
      }
      Y_cont[[m]] <- matrix(1, nrow = n, ncol = ncol(X_cont[[m]]))
      for (i in 1:ncol(X_cont[[m]])) {
        if (method == "Fleishman") {
          Y_cont[[m]][, i] <- constants[[m]][i, 1] +
            constants[[m]][i, 2] * X_cont[[m]][, i] +
            constants[[m]][i, 3] * X_cont[[m]][, i]^2 +
            constants[[m]][i, 4] * X_cont[[m]][, i]^3
        }
        if (method == "Polynomial") {
          Y_cont[[m]][, i] <- constants[[m]][i, 1] +
            constants[[m]][i, 2] * X_cont[[m]][, i] +
            constants[[m]][i, 3] * X_cont[[m]][, i]^2 +
            constants[[m]][i, 4] * X_cont[[m]][, i]^3 +
            constants[[m]][i, 5] * X_cont[[m]][, i]^4 +
            constants[[m]][i, 6] * X_cont[[m]][, i]^5
        }
        Y_cont[[m]][, i] <- means2[[m]][i] +
          sqrt(vars2[[m]][i]) * Y_cont[[m]][, i]
      }
      colnames(Y_cont[[m]]) <- colnames(X_cont[[m]])
      temp_means <- c(temp_means, means2[[m]])
      temp_vars <- c(temp_vars, vars2[[m]])
    }
  }
  Y_all <- NULL
  for (p in 1:M) {
    Y_all <- cbind(Y_all, Y_cont[[p]])
  }
  rho.x0 <- cor(Y_all)
  emax0 <- max(abs(rho.x0 - Corr_X))
  niter <- diag(0, nrow(rho.x0), ncol(rho.x0))
  start.time.error <- Sys.time()
  if (emax0 > epsilon & errorloop == TRUE) {
    EL <- corr_error(n = n, k_cont = length(grep("cont", colnames(Y_all))),
      method = method, means = temp_means, vars = temp_vars,
      constants = constants2, seed = seed, epsilon = epsilon, maxit = maxit,
      rho0 = Corr_X, Sigma = Sigma_X, rho_calc = rho.x0)
    Y_all <- EL$Y_cont
    colnames(Y_all) <- colnames(Corr_X)
    niter <- EL$niter
    colnames(niter) <- colnames(Corr_X)
    rownames(niter) <- colnames(Corr_X)
    Sigma_X <- EL$Sigma
  }
  stop.time.error <- Sys.time()
  Y_cont <- list()
  for (i in 1:M) {
    if (K.x[i] != 0) {
      if (i == 1) {
        Y_cont[[i]] <- Y_all[, 1:K.x[i], drop = FALSE]
      } else {
        if (!is.null(same.var)) {
          if (K.x[i] == nrow(same.var[same.var[, 3] == i, , drop = FALSE])) {
            Y_cont <- append(Y_cont, list(NULL))
          } else {
            Y_cont[[i]] <- Y_all[, (cumsum(K.x)[i - 1] -
              nrow(same.var[same.var[, 3] <= (i - 1), , drop = FALSE]) +
              1):(cumsum(K.x)[i] - nrow(same.var[same.var[, 3] <= i, ,
              drop = FALSE])), drop = FALSE]
          }
        } else {
          Y_cont[[i]] <- Y_all[, (cumsum(K.x)[i - 1] + 1):cumsum(K.x)[i],
                           drop = FALSE]
        }
      }
    } else Y_cont <- append(Y_cont, list(NULL))
  }
  for (m in 1:M) {
    if (K.x[m] > 0) {
      if (!is.null(same.cont)) {
        if (nrow(same.cont[same.cont[, 3] == m, , drop = FALSE]) > 0) {
          temp.cont <- same.cont[same.cont[, 3] == m, , drop = FALSE]
          for (i in 1:nrow(temp.cont)) {
            Y_cont[[m]] <- cbind(Y_cont[[m]],
                                 Y_cont[[temp.cont[i, 1]]][, temp.cont[i, 2]])
            colnames(Y_cont[[m]])[ncol(Y_cont[[m]])] <-
              colnames(Y_cont[[temp.cont[i, 1]]])[temp.cont[i, 2]]
          }
        }
      }
      if (length(Y_cont[[m]]) > 0) {
        Y_cont[[m]] <- Y_cont[[m]][, names0[[m]], drop = FALSE]
      }
    }
  }
  Y_comp <- Y_cont
  if (max(K.mix2) > 0) {
    for (i in 1:M) {
      if (K.mix2[i] != 0) {
        if ((ncol(Y_cont[[i]]) - K.comp[i] + K.error[i]) > 0) {
          Y_cont[[i]] <- cbind(Y_cont[[i]][, 1:(ncol(Y_cont[[i]]) -
            K.comp[i] + K.error[i]), drop = FALSE], matrix(1, n, K.mix2[i]))
        } else {
          Y_cont[[i]] <- matrix(1, n, K.mix2[i])
        }
        for (j in 1:K.mix2[i]) {
          ind <- K.x[i] - K.comp[i] + K.error[i] + j
          ind2 <- c(0, cumsum(sapply(mix_pis[[i]], length)))
          if (!is.null(same.cont)) {
            if (nrow(same.cont[same.cont[, 3] == i, , drop = FALSE]) > 0) {
              temp.cont <- same.cont[same.cont[, 3] == i, , drop = FALSE]
              if (all((ind - j + 1 + ind2[j]):(ind - j + ind2[j + 1]) %in%
                      temp.cont[, 4])) {
                temp.cont <- temp.cont[temp.cont[, 4] %in%
                  (ind - j + 1 + ind2[j]):(ind - j + ind2[j + 1]), ,
                                       drop = FALSE]
                if (length(unique(temp.cont[, 1])) == 1 &
                    length(unique(temp.cont[, 2])) ==
                    length(mix_pis[[i]][[j]])) {
                  ind <- ncol(Y_comp[[i]]) - K.comp[i] + K.error[i] + j
                  ind3 <- cumsum(sapply(mix_pis[[temp.cont[1, 1]]], length))
                  ind4 <- ncol(Y_cont[[temp.cont[1, 1]]]) -
                    K.mix2[temp.cont[1, 1]] +
                    min(which(ind3 >= max(temp.cont[1, 2])))
                  Y_cont[[i]][, ind] <- Y_cont[[temp.cont[1, 1]]][, ind4,
                                                                  drop = FALSE]
                  colnames(Y_cont[[i]])[ind] <-
                    colnames(Y_cont[[temp.cont[1, 1]]])[ind4]
                  next
                }
              }
            }
          }
          seed <- seed + 1
          set.seed(seed)
          R <- rmultinom(n, size = 1, prob = mix_pis[[i]][[j]])
          ind <- ncol(Y_comp[[i]]) - K.comp[i] + K.error[i] + j
          Y_cont[[i]][, ind] <- apply(t(R) *
            Y_comp[[i]][, (ind - j + 1 + ind2[j]):(ind - j + ind2[j + 1]),
                                                    drop = FALSE], 1, sum)
          Y_cont[[i]][, ind] <- means[[i]][ind] +
            sqrt(vars[[i]][ind]) * scale(Y_cont[[i]][, ind])
          colnames(Y_cont[[i]])[ind] <- paste("mix", i, "_", j, sep = "")
        }
      }
    }
  }
  betas <- calc_betas(corr.yx, corr.x, vars, mix_pis, mix_mus, mix_sigmas,
    error_type, n = 25, seed = seed)
  if (length(betas) == 1 | is.null(betas)) {
    stop("The beta coefficients could not be found.")
  }
  for (p in 1:M) {
    if (K.x[p] == 0) next
    colnames(Y_comp[[p]]) <- paste("X", p, "_", 1:K.x[p], sep = "")
    colnames(Y_cont[[p]]) <- paste("X", p, "_", 1:ncol(Y_cont[[p]]), sep = "")
  }
  Y <- matrix(1, nrow = n, ncol = M)
  colnames(Y) <- paste("Y", 1:M, sep = "")
  for (i in 1:M) {
    if (K.x[i] == 0) {
      Y[, i] <- matrix(betas.0[i], n, 1) + E2[, i]
    } else if ((length(betas[i, ]) - sum(betas[i, ] == 0)) == K.x[i]) {
      Y[, i] <- matrix(betas.0[i], n, 1) +
        Y_comp[[i]] %*% matrix(betas[i, 1:K.x[i]], ncol = 1) + E2[, i]
    } else {
      Y[, i] <- matrix(betas.0[i], n, 1) + Y_cont[[i]] %*%
        matrix(betas[i, 1:ncol(Y_cont[[i]])], ncol = 1) + E2[, i]
    }
  }
  stop.time <- Sys.time()
  Time.error <- round(difftime(stop.time.error, start.time.error,
                               units = "min"), 3)
  Time <- round(difftime(stop.time, start.time, units = "min"), 3)
  cat("Total Simulation time:", Time, "minutes \n")
  result <- list(Y = Y, X = Y_comp, X_all = Y_cont, E = E_comp, betas = betas,
    Sigma.X = Sigma_X, constants = constants0, SixCorr = SixCorr,
    valid.pdf = Valid.PDF, niter = niter, Error_Time = Time.error,
    Time = Time)
  if (error_type == "mix") result <- append(result, list(E_mix = E2))
  result
}

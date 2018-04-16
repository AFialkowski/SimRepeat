#' @title Summary of Correlated Systems of Variables
#'
#' @description This function summarizes the results of \code{\link[SimRepeat]{nonnormsys}}, \code{\link[SimRepeat]{corrsys}}, or
#'     \code{\link[SimRepeat]{corrsys2}}.  The inputs are either the simulated variables or inputs for those functions.  See their
#'     documentation for more information.  If only selected descriptions are desired, keep the non-relevant parameter inputs at their
#'     defaults.  For example, if only a description of the error terms are desired, \code{error_type = "non_mix"}, and
#'     \code{method = "Polynomial"}, specify \code{E, M, method, means, vars, skews, skurts, fifths, sixths, corr.e}.
#'
#' @param Y the matrix of outcomes simulated with \code{corrsys} or \code{corrsys2}
#' @param E the matrix of continuous non-mixture or components of mixture error terms
#' @param E_mix the matrix of continuous mixture error terms
#' @param X a list of length \code{M} where \code{X[[p]] = cbind(X_cat(pj), X_cont(pj), X_comp(pj), X_pois(pj), X_nb(pj))}; keep
#'     \code{X[[p]] = NULL} if \eqn{Y_p} has no independent variables
#' @param X_all a list of length \code{M} where \code{X_all[[p]]} contains all independent variables, interactions, and time for \eqn{Y_p}; keep
#'     \code{X_all[[p]] = NULL} if \eqn{Y_p} has no independent variables
#' @param M the number of dependent variables \eqn{Y} (outcomes); equivalently, the number of equations in the system
#' @param method the PMT method used to generate all continuous variables, including independent variables (covariates), error terms, and random effects;
#'     "Fleishman" uses Fleishman's third-order polynomial transformation and "Polynomial" uses Headrick's fifth-order transformation
#' @param means if no random effects, a list of length \code{M} where \code{means[[p]]} contains a vector of means for the continuous independent variables
#'     in equation p with non-mixture (\eqn{X_{cont}}) or mixture (\eqn{X_{mix}}) distributions and for the error terms (\eqn{E});
#'     order in vector is \eqn{X_{cont}, X_{mix}, E}
#'
#'     if there are random effects, a list of length \code{M + 1} if the effects are the same across equations or \code{2 * M} if they differ;
#'     where \code{means[M + 1]} or \code{means[(M + 1):(2 * M)]} are vectors of means for all random effects with continuous non-mixture or mixture distributions;
#'     order in vector is 1st random intercept \eqn{U_0} (if \code{rand.int != "none"}), 2nd random time slope \eqn{U_1} (if \code{rand.tsl != "none"}),
#'     3rd other random slopes with non-mixture distributions \eqn{U_{cont}}, 4th other random slopes with mixture distributions \eqn{U_{mix}}
#' @param vars a list of same length and order as \code{means} containing vectors of variances for the continuous variables, error terms, and any random effects
#' @param skews if no random effects, a list of length \code{M} where \code{skews[[p]]} contains a vector of skew values for the continuous independent variables
#'     in equation p with non-mixture (\eqn{X_{cont}}) distributions and for \eqn{E} if \code{error_type = "non_mix"}; order in vector is \eqn{X_{cont}, E}
#'
#'     if there are random effects, a list of length \code{M + 1} if the effects are the same across equations or \code{2 * M} if they differ;
#'     where \code{skews[M + 1]} or \code{skews[(M + 1):(2 * M)]} are vectors of skew values for all random effects with continuous non-mixture distributions;
#'     order in vector is 1st random intercept \eqn{U_0} (if \code{rand.int = "non_mix"}), 2nd random time slope \eqn{U_1} (if \code{rand.tsl = "non_mix"}),
#'     3rd other random slopes with non-mixture distributions \eqn{U_{cont}}
#' @param skurts a list of same length and order as \code{skews} containing vectors of standardized kurtoses (kurtosis - 3) for the continuous variables,
#'     error terms, and any random effects with non-mixture distributions
#' @param fifths a list of same length and order as \code{skews} containing vectors of standardized fifth cumulants for the continuous variables,
#'     error terms, and any random effects with non-mixture distributions; not necessary for \code{method = "Fleishman"}
#' @param sixths a list of same length and order as \code{skews} containing vectors of standardized sixth cumulants for the continuous variables,
#'     error terms, and any random effects with non-mixture distributions; not necessary for \code{method = "Fleishman"}
#' @param mix_pis list of length \code{M}, \code{M + 1} or \code{2 * M}, where \code{mix_pis[1:M]} are for \eqn{X_{cont}, E} (if \code{error_type = "mix"}) and
#'     \code{mix_pis[M + 1]} or \code{mix_pis[(M + 1):(2 * M)]} are for mixture \eqn{U}; use \code{mix_pis[[p]] = NULL} if equation p has no continuous mixture terms
#'     if \code{error_type = "non_mix"} and there are only random effects (i.e., \code{length(corr.x) = 0}), use \code{mix_pis[1:M] = NULL} so that
#'     \code{mix_pis[M + 1]} or \code{mix_pis[(M + 1):(2 * M)]} describes the mixture \eqn{U};
#'
#'     \code{mix_pis[[p]][[j]]} is a vector of mixing probabilities of the component distributions for \eqn{X_{mix(pj)}}, the j-th mixture covariate for outcome \eqn{Y_p};
#'     the last vector in \code{mix_pis[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"}); components should be ordered as in \code{corr.x}
#'
#'     \code{mix_pis[[M + p]][[j]]} is a vector of mixing probabilities of the component distributions for \eqn{U_{(pj)}}, the j-th random effect with a mixture distribution
#'     for outcome \eqn{Y_p}; order is 1st random intercept (if \code{rand.int = "mix"}), 2nd random time slope (if \code{rand.tsl = "mix"}),
#'     3rd other random slopes with mixture distributions; components should be ordered as in \code{corr.u}
#' @param mix_mus list of same length and order as \code{mix_pis};
#'
#'     \code{mix_mus[[p]][[j]]} is a vector of means of the component distributions for \eqn{X_{mix(pj)}},
#'     the last vector in \code{mix_mus[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"})
#'
#'     \code{mix_mus[[p]][[j]]} is a vector of means of the component distributions for \eqn{U_{mix(pj)}}
#' @param mix_sigmas list of same length and order as \code{mix_pis};
#'
#'     \code{mix_sigmas[[p]][[j]]} is a vector of standard deviations of the component distributions for \eqn{X_{mix(pj)}},
#'     the last vector in \code{mix_sigmas[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"})
#'
#'     \code{mix_sigmas[[p]][[j]]} is a vector of standard deviations of the component distributions for \eqn{U_{mix(pj)}}
#' @param mix_skews list of same length and order as \code{mix_pis};
#'
#'     \code{mix_skews[[p]][[j]]} is a vector of skew values of the component distributions for \eqn{X_{mix(pj)}},
#'     the last vector in \code{mix_skews[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"})
#'
#'     \code{mix_skews[[p]][[j]]} is a vector of skew values of the component distributions for \eqn{U_{mix(pj)}}
#' @param mix_skurts list of same length and order as \code{mix_pis};
#'
#'     \code{mix_skurts[[p]][[j]]} is a vector of standardized kurtoses of the component distributions for \eqn{X_{mix(pj)}},
#'     the last vector in \code{mix_skurts[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"})
#'
#'     \code{mix_skurts[[p]][[j]]} is a vector of standardized kurtoses of the component distributions for \eqn{U_{mix(pj)}}
#' @param mix_fifths list of same length and order as \code{mix_pis}; not necessary for \code{method = "Fleishman"};
#'
#'     \code{mix_fifths[[p]][[j]]} is a vector of standardized fifth cumulants of the component distributions for \eqn{X_{mix(pj)}},
#'     the last vector in \code{mix_fifths[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"})
#'
#'     \code{mix_fifths[[p]][[j]]} is a vector of standardized fifth cumulants of the component distributions for \eqn{U_{mix(pj)}}
#' @param mix_sixths list of same length and order as \code{mix_pis}; not necessary for \code{method = "Fleishman"};
#'
#'     \code{mix_sixths[[p]][[j]]} is a vector of standardized sixth cumulants of the component distributions for \eqn{X_{mix(pj)}},
#'     the last vector in \code{mix_sixths[[p]]} is for \eqn{E_p} (if \code{error_type = "mix"})
#'
#'     \code{mix_sixths[[p]][[j]]} is a vector of standardized sixth cumulants of the component distributions for \eqn{U_{mix(pj)}}
#' @param marginal a list of length \code{M}, with the p-th component a list of cumulative probabilities for the ordinal variables
#'     associated with outcome \eqn{Y_p} (use \code{marginal[[p]] = NULL} if outcome \eqn{Y_p} has no ordinal variables);
#'     \code{marginal[[p]][[j]]} is a vector of the cumulative probabilities defining the marginal distribution of \eqn{X_{ord(pj)}},
#'     the j-th ordinal variable for outcome \eqn{Y_p}; if the variable can take r values, the vector will contain r - 1 probabilities
#'     (the r-th is assumed to be 1); for binary variables, the probability is the probability of the 1st category, which has the smaller support value;
#'     \code{length(marginal[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}
#' @param support a list of length \code{M}, with the p-th component a list of support values for the ordinal variables associated
#'     with outcome \eqn{Y_p}; use \code{support[[p]] = NULL} if outcome \eqn{Y_p} has no ordinal variables;
#'     \code{support[[p]][[j]]} is a vector of the support values defining the marginal distribution of \eqn{X_{ord(pj)}},
#'     the j-th ordinal variable for outcome \eqn{Y_p}; if not provided, the default for r categories is 1, ..., r
#' @param lam list of length \code{M}, p-th component a vector of lambda (means > 0) values for Poisson variables for outcome \eqn{Y_p}
#'     (see \code{stats::dpois}); order is 1st regular Poisson and 2nd zero-inflated Poisson; use \code{lam[[p]] = NULL} if outcome \eqn{Y_p} has no Poisson variables;
#'     \code{length(lam[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}
#' @param p_zip a list of vectors of probabilities of structural zeros (not including zeros from the Poisson distribution) for the
#'     zero-inflated Poisson variables (see \code{VGAM::dzipois}); if \code{p_zip} = 0, \eqn{Y_{pois}} has a regular Poisson
#'     distribution; if \code{p_zip} is in (0, 1), \eqn{Y_{pois}} has a zero-inflated Poisson distribution;
#'     if \code{p_zip} is in \code{(-(exp(lam) - 1)^(-1), 0)}, \eqn{Y_{pois}} has a zero-deflated Poisson distribution and \code{p_zip}
#'     is not a probability; if \code{p_zip = -(exp(lam) - 1)^(-1)}, \eqn{Y_{pois}} has a positive-Poisson distribution
#'     (see \code{VGAM::dpospois}); order is 1st regular Poisson and 2nd zero-inflated Poisson;
#'     if a single number, all Poisson variables given this value; if a vector of length \code{M}, all Poisson variables in equation p
#'     given \code{p_zip[p]}; otherwise, missing values are set to 0 and ordered 1st
#' @param size list of length \code{M}, p-th component a vector of size parameters for the Negative Binomial variables for outcome \eqn{Y_p}
#'     (see \code{stats::nbinom}); order is 1st regular NB and 2nd zero-inflated NB; use \code{size[[p]] = NULL} if outcome \eqn{Y_p} has no Negative Binomial variables;
#'     \code{length(size[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}
#' @param prob list of length \code{M}, p-th component a vector of success probabilities for the Negative Binomial variables for outcome \eqn{Y_p}
#'     (see \code{stats::nbinom}); order is 1st regular NB and 2nd zero-inflated NB; use \code{prob[[p]] = NULL} if outcome \eqn{Y_p} has no Negative Binomial variables;
#'     \code{length(prob[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}
#' @param mu list of length \code{M}, p-th component a vector of mean values for the Negative Binomial variables for outcome \eqn{Y_p}
#'     (see \code{stats::nbinom}); order is 1st regular NB and 2nd zero-inflated NB; use \code{mu[[p]] = NULL} if outcome \eqn{Y_p} has no Negative Binomial variables;
#'     \code{length(mu[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}; for zero-inflated NB variables,
#'     this refers to the mean of the NB distribution (see \code{VGAM::dzinegbin})
#'     (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables, not a mixture)
#' @param p_zinb a vector of probabilities of structural zeros (not including zeros from the NB distribution) for the zero-inflated NB variables
#'     (see \code{VGAM::dzinegbin}); if \code{p_zinb} = 0, \eqn{Y_{nb}} has a regular NB distribution;
#'     if \code{p_zinb} is in \code{(-prob^size/(1 - prob^size),} \code{0)}, \eqn{Y_{nb}} has a zero-deflated NB distribution and \code{p_zinb}
#'     is not a probability; if \code{p_zinb = -prob^size/(1 - prob^size)}, \eqn{Y_{nb}} has a positive-NB distribution (see
#'     \code{VGAM::dposnegbin}); order is 1st regular NB and 2nd zero-inflated NB;
#'     if a single number, all NB variables given this value; if a vector of length \code{M}, all NB variables in equation p
#'     given \code{p_zinb[p]}; otherwise, missing values are set to 0 and ordered 1st
#' @param corr.x list of length \code{M}, each component a list of length \code{M}; \code{corr.x[[p]][[q]]} is matrix of correlations
#'     for independent variables in equations p (\eqn{X_{(pj)}} for outcome \eqn{Y_p}) and q (\eqn{X_{(qj)}} for outcome \eqn{Y_q});
#'     order: 1st ordinal (same order as in \code{marginal}), 2nd continuous non-mixture (same order as in \code{skews}),
#'     3rd components of continuous mixture (same order as in \code{mix_pis}), 4th regular Poisson, 5th zero-inflated Poisson (same order as in \code{lam}),
#'     6th regular NB, and 7th zero-inflated NB (same order as in \code{size});
#'     if p = q, \code{corr.x[[p]][[q]]} is a correlation matrix with \code{nrow(corr.x[[p]][[q]])} = # \eqn{X_{(pj)}} for outcome \eqn{Y_p};
#'     if p != q, \code{corr.x[[p]][[q]]} is a non-symmetric matrix of correlations where rows correspond to covariates for \eqn{Y_p}
#'     so that \code{nrow(corr.x[[p]][[q]])} = # \eqn{X_{(pj)}} for outcome \eqn{Y_p} and
#'     columns correspond to covariates for \eqn{Y_q} so that \code{ncol(corr.x[[p]][[q]])} = # \eqn{X_{(qj)}} for outcome \eqn{Y_q};
#'     use \code{corr.x[[p]][[q]] = NULL} if equation q has no \eqn{X_{(qj)}}; use \code{corr.x[[p]] = NULL} if equation p has no \eqn{X_{(pj)}}
#' @param corr.e correlation matrix for continuous non-mixture or components of mixture error terms
#' @param U a list of length \code{M} of continuous non-mixture and components of mixture random effects
#' @param U_all a list of length \code{M} of continuous non-mixture and mixture random effects
#' @param rand.int "none" (default) if no random intercept term for all outcomes, "non_mix" if all random intercepts have a continuous
#'     non-mixture distribution, "mix" if all random intercepts have a continuous mixture distribution;
#'     also can be a vector of length \code{M} containing a combination (i.e., \code{c("non_mix", "mix", "none")} if the 1st has a non-mixture
#'     distribution, the 2nd has a mixture distribution, and 3rd outcome has no random intercept)
#' @param rand.tsl "none" (default) if no random slope for time for all outcomes, "non_mix" if all random time slopes have a
#'     continuous non-mixture distribution, "mix" if all random time slopes have a continuous mixture distribution; also can
#'     be a vector of length \code{M} as in \code{rand.int}
#' @param corr.u if the random effects are the same variables across equations, a matrix of correlations for \eqn{U};
#'     if the random effects are different variables across equations, a list of length \code{M}, each component a list of length \code{M};
#'     \code{corr.u[[p]][[q]]} is matrix of correlations for random effects in equations p (\eqn{U_{(pj)}} for outcome \eqn{Y_p}) and
#'     q (\eqn{U_{(qj)}} for outcome \eqn{Y_q});
#'     if p = q, \code{corr.u[[p]][[q]]} is a correlation matrix with \code{nrow(corr.u[[p]][[q]])} = # \eqn{U_{(pj)}} for outcome \eqn{Y_p};
#'     if p != q, \code{corr.u[[p]][[q]]} is a non-symmetric matrix of correlations where rows correspond to \eqn{U_{(pj)}} for \eqn{Y_p}
#'     so that \code{nrow(corr.u[[p]][[q]])} = # \eqn{U_{(pj)}} for outcome \eqn{Y_p} and
#'     columns correspond to \eqn{U_{(qj)}} for \eqn{Y_q} so that \code{ncol(corr.u[[p]][[q]])} = # \eqn{U_{(qj)}} for outcome \eqn{Y_q};
#'     the number of random effects for \eqn{Y_p} is taken from \code{nrow(corr.u[[p]][[1]])} so that if there should be random effects,
#'     there must be entries for \code{corr.u};
#'     use \code{corr.u[[p]][[q]] = NULL} if equation q has no \eqn{U_{(qj)}}; use \code{corr.u[[p]] = NULL} if equation p has no \eqn{U_{(pj)}};
#'
#'     correlations are specified in terms of components of mixture variables (if present);
#'     order is 1st random intercept (if \code{rand.int != "none"}), 2nd random time slope (if \code{rand.tsl != "none"}),
#'     3rd other random slopes with non-mixture distributions, 4th other random slopes with mixture distributions
#' @param rmeans2 a list returned from \code{corrsys} or \code{corrsys2} which has the non-mixture and component means ordered according to
#'     types of random intercept and time slope
#' @param rvars2 a list returned like \code{rmeans}
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
#' @keywords summary
#' @seealso \code{\link[SimRepeat]{nonnormsys}}, \code{\link[SimRepeat]{corrsys}}, \code{\link[SimRepeat]{corrsys2}}
#'
#' @return A list with the following components:
#'
#' @return \code{cont_sum_y} a data.frame summarizing the simulated distributions of the \eqn{Y_p},
#'
#' @return \code{cont_sum_e} a data.frame summarizing the simulated distributions of the non-mixture or components of mixture \eqn{E_p},
#'
#' @return \code{target_sum_e} a data.frame summarizing the target distributions of the non-mixture or components of mixture \eqn{E_p},
#'
#' @return \code{mix_sum_e} a data.frame summarizing the simulated distributions of the mixture \eqn{E_p},
#'
#' @return \code{target_mix_e} a data.frame summarizing the target distributions of the mixture \eqn{E_p},
#'
#' @return \code{rho.y} correlation matrix of dimension \code{M x M} for \eqn{Y_p}
#'
#' @return \code{rho.e} correlation matrix for the non-mixture or components of mixture \eqn{E_p}
#'
#' @return \code{rho.emix} correlation matrix for the mixture \eqn{E_p}
#'
#' @return \code{rho.ye} matrix with correlations between \eqn{Y_p} (rows) and the non-mixture or components of mixture \eqn{E_p} (columns)
#'
#' @return \code{rho.yemix} matrix with correlations between \eqn{Y_p} (rows) and the mixture \eqn{E_p} (columns)
#'
#' @return \code{sum_xall} a data.frame summarizing \code{X_all} without the Time variable,
#'
#' @return \code{rho.yx} a list of length \code{M}, where \code{rho.yx[[p]]} is matrix of correlations
#'     between \eqn{Y} (rows) and \code{X[[p]]} = \eqn{{X_ord(pj), X_cont(pj), X_comp(pj), X_pois(pj), X_nb(pj)}} (columns)
#'
#' @return \code{rho.yxall} a list of length \code{M}, where \code{rho.yx[[p]]} is matrix of correlations
#'     between \eqn{Y} (rows) and \code{X_all[[p]]} (columns) not including \code{Time}
#'
#' @return \code{rho.x} a list of length \code{M} of lists of length \code{M} where
#'     \code{rho.x[[p]][[q]] = cor(cbind(X[[p]], X[[q]]))} if p!= q or
#'     \code{rho.x[[p]][[q]] = cor(X[[p]]))} if p = q, where \code{X[[p]]} = \eqn{{X_ord(pj), X_cont(pj), X_comp(pj), X_pois(pj), X_nb(pj)}}
#'
#' @return \code{rho.xall} a list of length \code{M} of lists of length \code{M} where
#'     \code{rho.xall[[p]][[q]] = cor(cbind(X_all[[p]], X_all[[q]]))} if p!= q or
#'     \code{rho.xall[[p]][[q]] = cor(X_all[[p]]))} if p = q, not including \code{Time}
#'
#' @return \code{maxerr} a list of length \code{M} containing a vector of length \code{M} with the maximum correlation errors between outcomes,
#'     \code{maxerr[[p]]][q] = abs(max(corr.x[[p]][[q]] - rho.x[[p]][[q]]))}
#'
#' @return Additional components vary based on the type of simulated variables:
#'
#' @return If \bold{ordinal variables} are produced:
#'     \code{ord_sum_x} a list where \code{ord_sum_x[[j]]} is a data.frame summarizing \eqn{X_{ord(pj)}} for all p = 1, ..., \code{M}
#'
#' @return If \bold{continuous variables} are produced:
#'     \code{cont_sum_x} a data.frame summarizing the simulated distributions of the \eqn{X_{cont(pj)}} and \eqn{X_comp(pj)},
#'
#'     \code{target_sum_x} a data.frame summarizing the target distributions of the \eqn{X_{cont(pj)}} and \eqn{X_comp(pj)},
#'
#'     \code{mix_sum_x} a data.frame summarizing the simulated distributions of the \eqn{X_{mix(pj)}},
#'
#'     \code{target_mix_x} a data.frame summarizing the target distributions of the \eqn{X_{mix(pj)}}
#'
#' @return If \bold{Poisson variables} are produced:
#'     \code{pois_sum_x} a data.frame summarizing the simulated distributions of the \eqn{X_{pois(pj)}}
#'
#' @return If \bold{Negative Binomial variables} are produced:
#'     \code{nb_sum_x} a data.frame summarizing the simulated distributions of the \eqn{X_{nb(pj)}}
#'
#' @return If \bold{random effects} are produced:
#'     \code{cont_sum_u} a data.frame summarizing the simulated distributions of the \eqn{U_{cont(pj)}} and \eqn{U_{comp(pj)}},
#'
#'     \code{target_sum_u} a data.frame summarizing the target distributions of the \eqn{U_{cont(pj)}} and \eqn{U_{comp(pj)}},
#'
#'     \code{sum_uall} a data.frame summarizing the simulated distributions of \code{U_all},
#'
#'     \code{mix_sum_u} a data.frame summarizing the simulated distributions of the \eqn{U_{mix(pj)}},
#'
#'     \code{target_mix_u} a data.frame summarizing the target distributions of the \eqn{U_{mix(pj)}},
#'
#'     \code{rho.u} list of length \code{M}, each component a list of length \code{M};
#'     \code{rho.u[[p]][[q]] = cor(cbind(U[[p]], U[[q]]))} if p != q or \code{rho.u[[p]][[q]] = cor(U[[p]]))} if p = q
#'
#'     \code{rho.uall} list of length \code{M}, each component a list of length \code{M};
#'     \code{rho.uall[[p]][[q]] = cor(cbind(U_all[[p]], U_all[[q]]))} if p != q or \code{rho.uall[[p]][[q]] = cor(U_all[[p]]))} if p = q
#'
#'     \code{maxerr_u} list of length \code{M} containing a vector of length \code{M} with the maximum correlation errors for \eqn{U} between outcomes
#'     \code{maxerr_u[[p]]][q] = abs(max(corr.u[[p]][[q]] - rho.u[[p]][[q]]))}
#'
#' @examples
#' M <- 3
#' B <- calc_theory("Beta", c(4, 1.5))
#' skews <- lapply(seq_len(M), function(x) B[3])
#' skurts <- lapply(seq_len(M), function(x) B[4])
#' fifths <- lapply(seq_len(M), function(x) B[5])
#' sixths <- lapply(seq_len(M), function(x) B[6])
#' Six <- lapply(seq_len(M), function(x) list(0.03))
#' corr.e <- matrix(c(1, 0.4, 0.4^2, 0.4, 1, 0.4, 0.4^2, 0.4, 1), M, M,
#'   byrow = TRUE)
#' means <- lapply(seq_len(M), function(x) B[1])
#' vars <- lapply(seq_len(M), function(x) B[2]^2)
#' marginal <- list(0.3, 0.4, 0.5)
#' support <- lapply(seq_len(M), function(x) list(0:1))
#' corr.x <- list(list(matrix(1, 1, 1), matrix(0.4, 1, 1), matrix(0.4, 1, 1)),
#'   list(matrix(0.4, 1, 1), matrix(1, 1, 1), matrix(0.4, 1, 1)),
#'   list(matrix(0.4, 1, 1), matrix(0.4, 1, 1), matrix(1, 1, 1)))
#' betas <- list(0.5)
#' betas.t <- 1
#' betas.tint <- list(0.25)
#' Sys1 <- corrsys(10000, M, Time = 1:M, "Polynomial", "non_mix", means, vars,
#'   skews, skurts, fifths, sixths, Six, marginal = marginal, support = support,
#'   corr.x = corr.x, corr.e = corr.e, betas = betas, betas.t = betas.t,
#'   betas.tint = betas.tint, quiet = TRUE)
#' Sum1 <- summary_sys(Sys1$Y, Sys1$E, E_mix = NULL, Sys1$X, Sys1$X_all, M,
#'   "Polynomial", means, vars, skews, skurts, fifths, sixths,
#'   marginal = marginal, support = support, corr.x = corr.x, corr.e = corr.e)
#'
#' \dontrun{
#' seed <- 276
#' n <- 10000
#' M <- 3
#' Time <- 1:M
#'
#' # Error terms have a beta(4, 1.5) distribution with an AR(1, p = 0.4)
#' # correlation structure
#' B <- calc_theory("Beta", c(4, 1.5))
#' skews <- lapply(seq_len(M), function(x) B[3])
#' skurts <- lapply(seq_len(M), function(x) B[4])
#' fifths <- lapply(seq_len(M), function(x) B[5])
#' sixths <- lapply(seq_len(M), function(x) B[6])
#' Six <- lapply(seq_len(M), function(x) list(0.03))
#' error_type <- "non_mix"
#' corr.e <- matrix(c(1, 0.4, 0.4^2, 0.4, 1, 0.4, 0.4^2, 0.4, 1), M, M,
#'   byrow = TRUE)
#'
#' # 1 continuous mixture of Normal(-2, 1) and Normal(2, 1) for each Y
#' mix_pis <- lapply(seq_len(M), function(x) list(c(0.4, 0.6)))
#' mix_mus <- lapply(seq_len(M), function(x) list(c(-2, 2)))
#' mix_sigmas <- lapply(seq_len(M), function(x) list(c(1, 1)))
#' mix_skews <- lapply(seq_len(M), function(x) list(c(0, 0)))
#' mix_skurts <- lapply(seq_len(M), function(x) list(c(0, 0)))
#' mix_fifths <- lapply(seq_len(M), function(x) list(c(0, 0)))
#' mix_sixths <- lapply(seq_len(M), function(x) list(c(0, 0)))
#' mix_Six <- list()
#' Nstcum <- calc_mixmoments(mix_pis[[1]][[1]], mix_mus[[1]][[1]],
#'   mix_sigmas[[1]][[1]], mix_skews[[1]][[1]], mix_skurts[[1]][[1]],
#'   mix_fifths[[1]][[1]], mix_sixths[[1]][[1]])
#'
#' means <- lapply(seq_len(M), function(x) c(Nstcum[1], B[1]))
#' vars <- lapply(seq_len(M), function(x) c(Nstcum[2]^2, B[2]^2))
#'
#' # 1 binary variable for each Y
#' marginal <- lapply(seq_len(M), function(x) list(0.4))
#' support <- list(NULL, list(c(0, 1)), NULL)
#'
#' # 1 Poisson variable for each Y
#' lam <- list(1, 5, 10)
#' # Y2 and Y3 are zero-inflated Poisson variables
#' p_zip <- list(NULL, 0.05, 0.1)
#'
#' # 1 NB variable for each Y
#' size <- list(10, 15, 20)
#' prob <- list(0.3, 0.4, 0.5)
#' # either prob or mu is required (not both)
#' mu <- mapply(function(x, y) x * (1 - y)/y, size, prob, SIMPLIFY = FALSE)
#' # Y2 and Y3 are zero-inflated NB variables
#' p_zinb <- list(NULL, 0.05, 0.1)
#'
#' # The 2nd (the normal mixture) variable is the same across Y
#' same.var <- 2
#'
#' # Create the correlation matrix in terms of the components of the normal
#' # mixture
#' K <- 5
#' corr.x <- list()
#' corr.x[[1]] <- list(matrix(0.1, K, K), matrix(0.2, K, K), matrix(0.3, K, K))
#' diag(corr.x[[1]][[1]]) <- 1
#' # set correlation between components to 0
#' corr.x[[1]][[1]][2:3, 2:3] <- diag(2)
#' # set correlations with the same variable equal across outcomes
#' corr.x[[1]][[2]][, same.var] <- corr.x[[1]][[3]][, same.var] <-
#'   corr.x[[1]][[1]][, same.var]
#' corr.x[[2]] <- list(t(corr.x[[1]][[2]]), matrix(0.35, K, K),
#'   matrix(0.4, K, K))
#'   diag(corr.x[[2]][[2]]) <- 1
#'   corr.x[[2]][[2]][2:3, 2:3] <- diag(2)
#' corr.x[[2]][[2]][, same.var] <- corr.x[[2]][[3]][, same.var] <-
#'   t(corr.x[[1]][[2]][same.var, ])
#' corr.x[[2]][[3]][same.var, ] <- corr.x[[1]][[3]][same.var, ]
#' corr.x[[2]][[2]][same.var, ] <- t(corr.x[[2]][[2]][, same.var])
#' corr.x[[3]] <- list(t(corr.x[[1]][[3]]), t(corr.x[[2]][[3]]),
#'   matrix(0.5, K, K))
#' diag(corr.x[[3]][[3]]) <- 1
#' corr.x[[3]][[3]][2:3, 2:3] <- diag(2)
#' corr.x[[3]][[3]][, same.var] <- t(corr.x[[1]][[3]][same.var, ])
#' corr.x[[3]][[3]][same.var, ] <- t(corr.x[[3]][[3]][, same.var])
#'
#' # The 2nd and 3rd variables of each Y are subject-level variables
#' subj.var <- matrix(c(1, 2, 1, 3, 2, 2, 2, 3, 3, 2, 3, 3), 6, 2, byrow = TRUE)
#' int.var <- tint.var <- NULL
#' betas.0 <- 0
#' betas <- list(seq(0.5, 0.5 + (K - 2) * 0.25, 0.25))
#' betas.subj <- list(seq(0.5, 0.5 + (K - 2) * 0.1, 0.1))
#' betas.int <- list()
#' betas.t <- 1
#' betas.tint <- list(c(0.25, 0.5))
#'
#' method <- "Polynomial"
#'
#' # Check parameter inputs
#' checkpar(M, method, error_type, means, vars, skews, skurts, fifths, sixths,
#'   Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts, mix_fifths,
#'   mix_sixths, mix_Six, marginal, support, lam, p_zip, pois_eps = list(),
#'   size, prob, mu, p_zinb, nb_eps = list(), corr.x, corr.yx = list(),
#'   corr.e, same.var, subj.var, int.var, tint.var, betas.0, betas,
#'   betas.subj, betas.int, betas.t, betas.tint)
#'
#' # Simulated system using correlation method 1
#' N <- corrsys(n, M, Time, method, error_type, means, vars, skews, skurts,
#'   fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
#'   mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip, size,
#'   prob, mu, p_zinb, corr.x, corr.e, same.var, subj.var, int.var, tint.var,
#'   betas.0, betas, betas.subj, betas.int, betas.t, betas.tint, seed = seed,
#'   use.nearPD = FALSE)
#'
#' # Summarize the results
#' S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
#'   vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
#'   mix_skews, mix_skurts, mix_fifths, mix_sixths, marginal, support, lam,
#'   p_zip, size, prob, mu, p_zinb, corr.x, corr.e)
#' S$sum_xall
#' S$maxerr
#' }
#'
#' @references
#' See references for \code{\link[SimRepeat]{SimRepeat}}.
#'
summary_sys <- function(Y = NULL, E = NULL, E_mix = NULL, X = list(),
                        X_all = list(), M = NULL,
                        method = c("Fleishman", "Polynomial"),
                        means =  list(), vars =  list(), skews =  list(),
                        skurts =  list(), fifths =  list(), sixths =  list(),
                        mix_pis = list(), mix_mus = list(),
                        mix_sigmas = list(), mix_skews =  list(),
                        mix_skurts =  list(), mix_fifths =  list(),
                        mix_sixths =  list(), marginal = list(),
                        support = list(), lam  =  list(), p_zip = list(),
                        size = list(), prob = list(), mu = list(),
                        p_zinb = list(), corr.x = list(), corr.e = NULL,
                        U = list(), U_all = list(),
                        rand.int = c("none", "non_mix", "mix"),
                        rand.tsl = c("none", "non_mix", "mix"),
                        corr.u = list(), rmeans2 = list(), rvars2 = list()) {
  error_type <- ifelse(is.null(E_mix) & is.null(E), "none",
                       ifelse(is.null(E_mix), "non_mix", "mix"))
  if (length(X_all) == 0 & length(X) > 0) X_all <- X
  K.cat <- rep(0, M)
  K.pois <- rep(0, M)
  K.nb <- rep(0, M)
  if (length(marginal) > 0) K.cat <- lengths(marginal)
  if (length(lam) > 0) K.pois <- lengths(lam)
  if (length(size) > 0) K.nb <- lengths(size)
  if (max(K.pois) > 0) {
    p_zip0 <- p_zip
    p_zip <- list()
    if (class(p_zip0) == "numeric" & length(p_zip0) == 1) {
      for (i in 1:M) {
        if (K.pois[i] == 0) p_zip <- append(p_zip, list(NULL)) else
          p_zip[[i]] <- rep(p_zip0, K.pois[i])
      }
    } else if (class(p_zip0) == "numeric" & length(p_zip0) == M) {
      for (i in 1:M) {
        if (K.pois[i] == 0) p_zip <- append(p_zip, list(NULL)) else
          p_zip[[i]] <- rep(p_zip0[i], K.pois[i])
      }
    } else if (class(p_zip0) == "list" & length(p_zip0) == M) {
      for (i in 1:M) {
        if (K.pois[i] == 0) {
          p_zip <- append(p_zip, list(NULL))
          next
        }
        if (length(p_zip0[[i]]) == K.pois[i])
          p_zip[[i]] <- p_zip0[[i]] else
            p_zip[[i]] <- c(rep(0, K.pois[i] - length(p_zip0[[i]])),
                            p_zip0[[i]])
      }
    } else {
      p_zip <- lapply(K.pois, function(x) if (x == 0) NULL else rep(0, x))
    }
  }
  if (max(K.nb) > 0) {
    if (length(prob) > 0)
      mu <- mapply(function(x, y) x * (1 - y)/y, size, prob, SIMPLIFY = FALSE)
    p_zinb0 <- p_zinb
    p_zinb <- list()
    if (class(p_zinb0) == "numeric" & length(p_zinb0) == 1) {
      for (i in 1:M) {
        if (K.nb[i] == 0) p_zinb <- append(p_zinb, list(NULL)) else
          p_zinb[[i]] <- rep(p_zinb0, K.nb[i])
      }
    } else if (class(p_zinb0) == "numeric" & length(p_zinb0) == M) {
      for (i in 1:M) {
        if (K.nb[i] == 0) p_zinb <- append(p_zinb, list(NULL)) else
          p_zinb[[i]] <- rep(p_zinb0[i], K.nb[i])
      }
    } else if (class(p_zinb0) == "list" & length(p_zinb0) == M) {
      for (i in 1:M) {
        if (K.nb[i] == 0) {
          p_zinb <- append(p_zinb, list(NULL))
          next
        }
        if (length(p_zinb0[[i]]) == K.nb[i])
          p_zinb[[i]] <- p_zinb0[[i]] else
            p_zinb[[i]] <- c(rep(0, K.nb[i] - length(p_zinb0[[i]])),
                             p_zinb0[[i]])
      }
    } else {
      p_zinb <- lapply(K.nb, function(x) if (x == 0) NULL else rep(0, x))
    }
  }
  K.r <- rep(0, M)
  if (class(corr.u) == "list" & length(corr.u) > 0)
    K.r <- sapply(mapply('[', corr.u, lengths(corr.u), SIMPLIFY = FALSE),
                  function(x) if (is.null(x)) 0 else nrow(x[[1]]))
  if (class(corr.u) == "matrix") K.r <- rep(ncol(corr.u), M)
  rmix_pis <- list()
  if (length(means) %in% c(2 * M, M + 1)) {
    means <- means[1:M]
    vars <- vars[1:M]
  }
  if (length(skews) == (M + 1)) {
    rskews <- skews[M + 1]
    rskurts <- skurts[M + 1]
    if (method == "Polynomial") {
      rfifths <- fifths[M + 1]
      rsixths <- sixths[M + 1]
    }
  }
  if (length(skews) == (2 * M)) {
    rskews <- skews[(M + 1):(2 * M)]
    rskurts <- skurts[(M + 1):(2 * M)]
    if (method == "Polynomial") {
      rfifths <- fifths[(M + 1):(2 * M)]
      rsixths <- sixths[(M + 1):(2 * M)]
    }
  }
  if (length(skews) %in% c(2 * M, M + 1)) {
    skews <- skews[1:M]
    skurts <- skurts[1:M]
    if (method == "Polynomial") {
      fifths <- fifths[1:M]
      sixths <- sixths[1:M]
    }
  }
  if (length(mix_pis) == (M + 1)) {
    rmix_pis <- mix_pis[M + 1]
    rmix_mus <- mix_mus[M + 1]
    rmix_sigmas <- mix_sigmas[M + 1]
    rmix_skews <- mix_skews[M + 1]
    rmix_skurts <- mix_skurts[M + 1]
    if (method == "Polynomial") {
      rmix_fifths <- mix_fifths[M + 1]
      rmix_sixths <- mix_sixths[M + 1]
    }
  }
  if (length(mix_pis) == (2 * M)) {
    rmix_pis <- mix_pis[(M + 1):(2 * M)]
    rmix_mus <- mix_mus[(M + 1):(2 * M)]
    rmix_sigmas <- mix_sigmas[(M + 1):(2 * M)]
    rmix_skews <- mix_skews[(M + 1):(2 * M)]
    rmix_skurts <- mix_skurts[(M + 1):(2 * M)]
    if (method == "Polynomial") {
      rmix_fifths <- mix_fifths[(M + 1):(2 * M)]
      rmix_sixths <- mix_sixths[(M + 1):(2 * M)]
    }
  }
  if (length(mix_pis) %in% c(2 * M, M + 1)) {
    mix_pis <- mix_pis[1:M]
    mix_mus <- mix_mus[1:M]
    mix_sigmas <- mix_sigmas[1:M]
    mix_skews <- mix_skews[1:M]
    mix_skurts <- mix_skurts[1:M]
    if (method == "Polynomial") {
      mix_fifths <- mix_fifths[1:M]
      mix_sixths <- mix_sixths[1:M]
    }
  }
  K.mix <- rep(0, M)
  K.comp <- rep(0, M)
  K.error <- rep(0, M)
  if (length(mix_pis) > 0) {
    K.mix <- lengths(mix_pis)
    K.comp <- sapply(lapply(mix_pis, unlist), length)
    k.comp <- c(0, cumsum(K.comp))
  }
  K.cont <- lengths(vars) - K.mix
  k.cont <- c(0, cumsum(K.cont))
  k.mix <- c(0, cumsum(K.mix))
  if (error_type == "mix") {
    K.error <- sapply(mix_pis, function(x) lengths(x[length(x)]))
    k.error <- c(0, cumsum(K.error))
    K.x <- K.cont + K.comp - K.error + K.cat + K.pois + K.nb
    K.mix2 <- K.mix - 1
  } else {
    K.x <- K.cont + K.comp - 1 + K.cat + K.pois + K.nb
    K.mix2 <- K.mix
  }
  if (length(marginal) > 0) {
    if (length(support) == 0) {
      for (m in 1:M) {
        support <- append(support, list(NULL))
        if (K.cat[m] > 0) support[[m]] <- lapply(marginal[[m]],
          function(x) 1:(length(x) + 1))
      }
    } else {
      for (m in 1:M) {
        if (K.cat[m] > 0 & is.null(support[[m]])) {
          support[[m]] <- lapply(marginal[[m]], function(x) 1:(length(x) + 1))
        } else if (K.cat[m] > 0) {
          for (i in 1:K.cat[m]) {
            if (class(support[[m]]) != "list")
              support[[m]] <- append(support[[m]], list(NULL))
            if (length(support[[m]][[i]]) != (length(marginal[[m]][[i]]) + 1))
              support[[m]][[i]] <- 1:(length(marginal[[m]][[i]]) + 1)
          }
       } else support <- append(support, list(NULL))
      }
    }
  }
  if (!is.null(Y)) {
    rho.y <- cor(Y)
    mom <- apply(Y, 2, calc_moments)
    medians <- apply(Y, 2, median)
    mins <- apply(Y, 2, min)
    maxs <- apply(Y, 2, max)
    cont_sum_y <- as.data.frame(cbind(1:M, rep(nrow(Y), M), mom[1, ], mom[2, ],
      medians, mins, maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ]))
    colnames(cont_sum_y) <- c("Outcome", "N", "Mean", "SD", "Median",
      "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
    result <- list(cont_sum_y = cont_sum_y, rho.y = rho.y)
  }
  if (error_type != "none") {
    mom <- apply(E, 2, calc_moments)
    medians <- apply(E, 2, median)
    mins <- apply(E, 2, min)
    maxs <- apply(E, 2, max)
    if (error_type == "non_mix") {
      e.means <- mapply('[[', means, lengths(means))
      e.vars <- mapply('[[', vars, lengths(vars))
      cont_sum_e <- as.data.frame(cbind(1:M, rep(nrow(E), M), mom[1, ],
        mom[2, ], medians, mins, maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ]))
      colnames(cont_sum_e) <- c("Outcome", "N", "Mean", "SD", "Median",
        "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
      if (method == "Fleishman") {
        target_sum_e <- as.data.frame(cbind(1:M, e.means, sqrt(e.vars),
          mapply('[[', skews, lengths(skews)),
          mapply('[[', skurts, lengths(skurts))))
        colnames(target_sum_e) <- c("Outcome", "Mean", "SD", "Skew",
                                    "Skurtosis")
      } else {
        target_sum_e <- as.data.frame(cbind(1:M, e.means, sqrt(e.vars),
          mapply('[[', skews, lengths(skews)),
          mapply('[[', skurts, lengths(skurts)),
          mapply('[[', fifths, lengths(fifths)),
          mapply('[[', sixths, lengths(sixths))))
        colnames(target_sum_e) <- c("Outcome", "Mean", "SD", "Skew",
                                    "Skurtosis", "Fifth", "Sixth")
      }
    }
    if (error_type == "mix") {
      e.means <- unlist(mapply('[', mix_mus, lengths(mix_mus)))
      e.vars <- (unlist(mapply('[', mix_sigmas, lengths(mix_sigmas))))^2
      cont_sum_e <- as.data.frame(cbind(rep(1:M, K.error),
        unlist(lapply(K.error, seq)), rep(nrow(E), ncol(mom)), mom[1, ],
        mom[2, ], medians, mins, maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ]))
      colnames(cont_sum_e) <- c("Outcome", "Component", "N", "Mean", "SD",
        "Median", "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
      rho.emix <- cor(E_mix)
      mom <- apply(E_mix, 2, calc_moments)
      medians <- apply(E_mix, 2, median)
      mins <- apply(E_mix, 2, min)
      maxs <- apply(E_mix, 2, max)
      mix_sum_e <- as.data.frame(cbind(1:M, rep(nrow(E_mix), M), mom[1, ],
        mom[2, ], medians, mins, maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ]))
      colnames(mix_sum_e) <- c("Outcome", "N", "Mean", "SD", "Median",
        "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
      target_mix_e <- NULL
      if (method == "Fleishman") {
        target_sum_e <- as.data.frame(cbind(rep(1:M, K.error),
          unlist(lapply(K.error, seq)), e.means, sqrt(e.vars),
          unlist(mapply('[', mix_skews, lengths(mix_skews))),
          unlist(mapply('[', mix_skurts, lengths(mix_skurts)))))
        colnames(target_sum_e) <- c("Outcome", "Component", "Mean", "SD",
                                    "Skew", "Skurtosis")
        for (i in 1:M) {
          target_mix_e <- rbind(target_mix_e,
            calc_mixmoments(mix_pis[[i]][[length(mix_pis[[i]])]],
              mix_mus[[i]][[length(mix_pis[[i]])]],
              mix_sigmas[[i]][[length(mix_pis[[i]])]],
              mix_skews[[i]][[length(mix_pis[[i]])]],
              mix_skurts[[i]][[length(mix_pis[[i]])]]))
        }
        target_mix_e <- as.data.frame(cbind(1:M, target_mix_e))
        colnames(target_mix_e) <- c("Outcome", "Mean", "SD", "Skew",
                                    "Skurtosis")
      } else {
        target_sum_e <- as.data.frame(cbind(rep(1:M, K.error),
          unlist(lapply(K.error, seq)), e.means, sqrt(e.vars),
          unlist(mapply('[', mix_skews, lengths(mix_skews))),
          unlist(mapply('[', mix_skurts, lengths(mix_skurts))),
          unlist(mapply('[', mix_fifths, lengths(mix_fifths))),
          unlist(mapply('[', mix_sixths, lengths(mix_sixths)))))
        colnames(target_sum_e) <- c("Outcome", "Component", "Mean", "SD",
          "Skew", "Skurtosis", "Fifth", "Sixth")
        for (i in 1:M) {
          target_mix_e <- rbind(target_mix_e,
            calc_mixmoments(mix_pis[[i]][[length(mix_pis[[i]])]],
              mix_mus[[i]][[length(mix_pis[[i]])]],
              mix_sigmas[[i]][[length(mix_pis[[i]])]],
              mix_skews[[i]][[length(mix_pis[[i]])]],
              mix_skurts[[i]][[length(mix_pis[[i]])]],
              mix_fifths[[i]][[length(mix_pis[[i]])]],
              mix_sixths[[i]][[length(mix_pis[[i]])]]))
        }
        target_mix_e <- as.data.frame(cbind(1:M, target_mix_e))
        colnames(target_mix_e) <- c("Outcome", "Mean", "SD", "Skew",
          "Skurtosis", "Fifth", "Sixth")
      }
      rownames(target_mix_e) <- rownames(mix_sum_e)
    }
    rownames(target_sum_e) <- rownames(cont_sum_e)
    rho.e <- cor(E)
    result <- append(result, list(cont_sum_e = cont_sum_e,
      target_sum_e = target_sum_e, rho.e = rho.e))
    if (error_type == "mix") {
      result <- append(result, list(mix_sum_e = mix_sum_e,
        target_mix_e = target_mix_e, rho.emix = rho.emix))
    }
    if (!is.null(Y)) {
      rho.ye <- cor(cbind(Y, E))[1:M, (M + 1):(M + ncol(E)), drop = FALSE]
      result <- append(result, list(rho.ye = rho.ye))
      if (error_type == "mix") {
        rho.yemix <- cor(cbind(Y, E_mix))[1:M, (M + 1):(M + ncol(E_mix)),
                                          drop = FALSE]
        result <- append(result, list(rho.yemix = rho.yemix))
      }
    }
  }
  if (max(K.x) > 0) {
    n <- max(unlist(lapply(X, nrow)))
    for (i in 1:M) {
      if (!is.null(X_all[[i]])) {
        if (var(X_all[[i]][, ncol(X_all[[i]])]) == 0)
          X_all[[i]] <- X_all[[i]][, -ncol(X_all[[i]]), drop = FALSE]
      }
    }
    if (length(marginal) > 0) {
      ord_sum_x <- list()
      for (j in 1:max(K.cat)) {
        ord_sum_x <- append(ord_sum_x, list(NULL))
        for (i in 1:M) {
          if (K.cat[i] == 0 | K.cat[i] < j) next
          csum <- cumsum(table(X[[i]][, j]))/n
          ord_sum_x[[j]] <- rbind(ord_sum_x[[j]],
            cbind(rep(i, length(support[[i]][[j]])), support[[i]][[j]],
            append(marginal[[i]][[j]], 1),
            c(csum, rep(0, length(support[[i]][[j]]) - length(csum)))))
        }
        rownames(ord_sum_x[[j]]) <- 1:nrow(ord_sum_x[[j]])
        ord_sum_x[[j]] <- as.data.frame(ord_sum_x[[j]])
        colnames(ord_sum_x[[j]]) <- c("Outcome", "Support", "Target",
                                      "Simulated")
        names(ord_sum_x)[j] <- paste("O", j, sep = "")
      }
      result <- append(result, list(ord_sum_x = ord_sum_x))
    }
    sum_xall <- list()
    for (i in 1:M) {
      sum_xall <- append(sum_xall, list(NULL))
      if (K.x[i] > 0) {
        mom <- apply(X_all[[i]], 2, calc_moments)
        medians <- apply(X_all[[i]], 2, median)
        mins <- apply(X_all[[i]], 2, min)
        maxs <- apply(X_all[[i]], 2, max)
        sum_xall[[i]] <- cbind(rep(i, ncol(X_all[[i]])),
          1:ncol(X_all[[i]]), rep(n, ncol(X_all[[i]])), mom[1, ], mom[2, ],
          medians, mins, maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ])
        rownames(sum_xall[[i]]) <- paste(colnames(X_all[[i]]), ":", i, sep = "")
      }
    }
    sum_xall <- as.data.frame(do.call(rbind, sum_xall))
    colnames(sum_xall) <- c("Outcome", "X", "N", "Mean", "SD", "Median",
      "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
    if (max(K.x - K.cat - K.pois - K.nb) > 0) {
      cont_sum_x <- list()
      target_sum_x <- list()
      mix_sum_x <- list()
      target_mix_x <- list()
      Y_comp <- list()
      for (i in 1:M) {
        cont_sum_x <- append(cont_sum_x, list(NULL))
        target_sum_x <- append(target_sum_x, list(NULL))
        mix_sum_x <- append(mix_sum_x, list(NULL))
        target_mix_x <- append(target_mix_x, list(NULL))
        Y_comp <- append(Y_comp, list(NULL))
        if ((K.x[i] - K.cat[i] - K.pois[i] - K.nb[i]) > 0) {
          Y_comp[[i]] <- X[[i]][, (K.cat[i] + 1):(ncol(X[[i]]) - K.pois[i] -
                                                    K.nb[i]), drop = FALSE]
          mom <- apply(Y_comp[[i]], 2, calc_moments)
          medians <- apply(Y_comp[[i]], 2, median)
          mins <- apply(Y_comp[[i]], 2, min)
          maxs <- apply(Y_comp[[i]], 2, max)
          cont_sum_x[[i]] <- cbind(rep(i, ncol(Y_comp[[i]])),
            1:ncol(Y_comp[[i]]), rep(n, ncol(Y_comp[[i]])), mom[1, ], mom[2, ],
            medians, mins, maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ])
          rownames(cont_sum_x[[i]]) <- paste("cont", rep(i, ncol(Y_comp[[i]])),
            "_", 1:ncol(Y_comp[[i]]), sep = "")
          ind <- ncol(Y_comp[[i]]) - K.comp[i] + K.error[i]
          if (method == "Fleishman") {
            if (ind > 0) {
              target_sum_x[[i]] <- cbind(means[[i]][1:ind],
                sqrt(vars[[i]][1:ind]), skews[[i]][1:ind],
                skurts[[i]][1:ind])
            }
            if (K.mix2[i] > 0) {
              target_sum_x[[i]] <- rbind(target_sum_x[[i]],
                cbind(unlist(mix_mus[[i]][1:K.mix2[i]]),
                unlist(mix_sigmas[[i]][1:K.mix2[i]]),
                unlist(mix_skews[[i]][1:K.mix2[i]]),
                unlist(mix_skurts[[i]][1:K.mix2[i]])))
            }
          } else {
            if (ind > 0) {
              target_sum_x[[i]] <- cbind(means[[i]][1:ind],
                sqrt(vars[[i]][1:ind]), skews[[i]][1:ind], skurts[[i]][1:ind],
                fifths[[i]][1:ind], sixths[[i]][1:ind])
            }
            if (K.mix2[i] > 0) {
              target_sum_x[[i]] <- rbind(target_sum_x[[i]],
                cbind(unlist(mix_mus[[i]][1:K.mix2[i]]),
                unlist(mix_sigmas[[i]][1:K.mix2[i]]),
                unlist(mix_skews[[i]][1:K.mix2[i]]),
                unlist(mix_skurts[[i]][1:K.mix2[i]]),
                unlist(mix_fifths[[i]][1:K.mix2[i]]),
                unlist(mix_sixths[[i]][1:K.mix2[i]])))
            }
          }
          target_sum_x[[i]] <- cbind(rep(i, ncol(Y_comp[[i]])),
            1:ncol(Y_comp[[i]]), target_sum_x[[i]])
          rownames(target_sum_x[[i]]) <- rownames(cont_sum_x[[i]])
          if (K.mix2[i] > 0) {
            mom <- apply(X_all[[i]][, (K.cat[i] + ind + 1):(K.cat[i] + ind +
                                  K.mix2[i]), drop = FALSE], 2, calc_moments)
            medians <- apply(X_all[[i]][, (K.cat[i] + ind + 1):(K.cat[i] + ind +
                                  K.mix2[i]), drop = FALSE], 2, median)
            mins <- apply(X_all[[i]][, (K.cat[i] + ind + 1):(K.cat[i] + ind +
                                  K.mix2[i]), drop = FALSE], 2, min)
            maxs <- apply(X_all[[i]][, (K.cat[i] + ind + 1):(K.cat[i] + ind +
                                  K.mix2[i]), drop = FALSE], 2, max)
            mix_sum_x[[i]] <- cbind(rep(i, K.mix2[i]), 1:K.mix2[i],
              rep(n, K.mix2[i]), mom[1, ], mom[2, ],
              medians, mins, maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ])
            rownames(mix_sum_x[[i]]) <- paste("mix", rep(i, K.mix2[i]), "_",
                                              1:K.mix2[i], sep = "")
            if (method == "Fleishman") {
              for (j in 1:K.mix2[i]) {
                target_mix_x[[i]] <- rbind(target_mix_x[[i]],
                  calc_mixmoments(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
                  mix_sigmas[[i]][[j]], mix_skews[[i]][[j]],
                  mix_skurts[[i]][[j]]))
              }
            } else {
              for (j in 1:K.mix2[i]) {
                target_mix_x[[i]] <- rbind(target_mix_x[[i]],
                  calc_mixmoments(mix_pis[[i]][[j]], mix_mus[[i]][[j]],
                  mix_sigmas[[i]][[j]], mix_skews[[i]][[j]],
                  mix_skurts[[i]][[j]], mix_fifths[[i]][[j]],
                  mix_sixths[[i]][[j]]))
              }
            }
            target_mix_x[[i]] <- cbind(rep(i, K.mix2[i]), 1:K.mix2[i],
                                       target_mix_x[[i]])
            rownames(target_mix_x[[i]]) <- rownames(mix_sum_x[[i]])
          }
        }
      }
      cont_sum_x <- as.data.frame(do.call(rbind, cont_sum_x))
      colnames(cont_sum_x) <- c("Outcome", "X", "N", "Mean", "SD", "Median",
        "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
      target_sum_x <- as.data.frame(do.call(rbind, target_sum_x))
      if (method == "Fleishman") {
        colnames(target_sum_x) <- c("Outcome", "X", "Mean", "SD", "Skew",
                                    "Skurtosis")
      } else {
        colnames(target_sum_x) <- c("Outcome", "X", "Mean", "SD", "Skew",
                                    "Skurtosis", "Fifth", "Sixth")
      }
      result <- append(result, list(cont_sum_x = cont_sum_x,
        target_sum_x = target_sum_x, sum_xall = sum_xall))
      if (max(K.mix2) > 0) {
        target_mix_x <- as.data.frame(do.call(rbind, target_mix_x))
        colnames(target_mix_x) <- colnames(target_sum_x)
        mix_sum_x <- as.data.frame(do.call(rbind, mix_sum_x))
        colnames(mix_sum_x) <- colnames(cont_sum_x)
        result <- append(result, list(mix_sum_x = mix_sum_x,
                                      target_mix_x = target_mix_x))
      }
    }
    if (length(lam) > 0) {
      pois_sum_x <- list()
      for (i in 1:M) {
        if (K.pois[i] == 0) {
          pois_sum_x <- append(pois_sum_x, list(NULL))
        } else {
          mom <- apply(X[[i]][, (ncol(X[[i]]) - K.pois[i] - K.nb[i] +
            1):(ncol(X[[i]]) - K.nb[i]), drop = FALSE], 2, calc_moments)
          medians <- apply(X[[i]][, (ncol(X[[i]]) - K.pois[i] - K.nb[i] +
            1):(ncol(X[[i]]) - K.nb[i]), drop = FALSE], 2, median)
          mins <- apply(X[[i]][, (ncol(X[[i]]) - K.pois[i] - K.nb[i] +
            1):(ncol(X[[i]]) - K.nb[i]), drop = FALSE], 2, min)
          maxs <- apply(X[[i]][, (ncol(X[[i]]) - K.pois[i] - K.nb[i] +
            1):(ncol(X[[i]]) - K.nb[i]), drop = FALSE], 2, max)
          p_0 <- apply(X[[i]][, (ncol(X[[i]]) - K.pois[i] - K.nb[i] +
            1):(ncol(X[[i]]) - K.nb[i]), drop = FALSE], 2,
            function(x) sum(x == 0)/n)
          pois_sum_x[[i]] <- cbind(rep(i, K.pois[i]), 1:K.pois[i],
            rep(n, K.pois[i]), p_0,
            mapply(function(x, y) dzipois(0, x, y), lam[[i]], p_zip[[i]]),
            mom[1, ], (1 - p_zip[[i]]) * lam[[i]], mom[2, ]^2,
            lam[[i]] + (lam[[i]]^2) * p_zip[[i]]/(1 - p_zip[[i]]),
            medians, mins, maxs, mom[3, ], mom[4, ])
          rownames(pois_sum_x[[i]]) <- paste("pois", rep(i, K.pois[i]), "_",
            1:K.pois[i], sep = "")
        }
      }
      pois_sum_x <- as.data.frame(do.call(rbind, pois_sum_x))
      colnames(pois_sum_x) <- c("Outcome", "X", "N", "P0", "Exp_P0", "Mean",
        "Exp_Mean", "Var", "Exp_Var", "Median", "Min", "Max", "Skew",
        "Skurtosis")
      result <- append(result, list(pois_sum_x = pois_sum_x))
    }
    if (length(size) > 0) {
      nb_sum_x <- list()
      for (i in 1:M) {
        if (K.nb[i] == 0) {
          nb_sum_x <- append(nb_sum_x, list(NULL))
        } else {
          prob <- size[[i]]/(mu[[i]] + size[[i]])
          p_0 <- apply(X[[i]][, (ncol(X[[i]]) - K.nb[i] + 1):ncol(X[[i]]),
            drop = FALSE], 2, function(x) sum(x == 0)/n)
          mom <- apply(X[[i]][, (ncol(X[[i]]) - K.nb[i] + 1):ncol(X[[i]]),
            drop = FALSE], 2, calc_moments)
          medians <- apply(X[[i]][, (ncol(X[[i]]) - K.nb[i] + 1):ncol(X[[i]]),
            drop = FALSE], 2, median)
          mins <- apply(X[[i]][, (ncol(X[[i]]) - K.nb[i] + 1):ncol(X[[i]]),
            drop = FALSE], 2, min)
          maxs <- apply(X[[i]][, (ncol(X[[i]]) - K.nb[i] + 1):ncol(X[[i]]),
            drop = FALSE], 2, max)
          nb_sum_x[[i]] <- cbind(rep(i, K.nb[i]), 1:K.nb[i], rep(n, K.nb[i]),
            p_0, mapply(function(x, y, z) dzinegbin(0, x, munb = y, pstr0 = z),
              size[[i]], mu[[i]], p_zinb[[i]]),
            prob, mom[1, ], (1 - p_zinb[[i]]) * mu[[i]], mom[2, ]^2,
            (1 - p_zinb[[i]]) * mu[[i]] * (1 + mu[[i]] *
              (p_zinb[[i]] + 1/size[[i]])), medians, mins, maxs, mom[3, ],
            mom[4, ])
          rownames(nb_sum_x[[i]]) <- paste("nb", rep(i, K.nb[i]), "_",
            1:K.nb[i], sep = "")
        }
      }
      nb_sum_x <- as.data.frame(do.call(rbind, nb_sum_x))
      colnames(nb_sum_x) <- c("Outcome", "X", "N", "P0", "Exp_P0", "Prob",
        "Mean", "Exp_Mean", "Var", "Exp_Var", "Median", "Min", "Max", "Skew",
        "Skurtosis")
      result <- append(result, list(nb_sum_x = nb_sum_x))
    }
    rho.x <- list()
    rho.xall <- list()
    for (p in 1:M) {
      if (K.x[p] == 0) {
        rho.x <- append(rho.x, list(NULL))
        rho.xall <- append(rho.xall, list(NULL))
        next
      }
      rho.x[[p]] <- list()
      rho.xall[[p]] <- list()
      for (q in 1:M) {
        if (K.x[q] == 0) {
          rho.x[[p]] <- append(rho.x[[p]], list(NULL))
          rho.xall[[p]] <- append(rho.xall[[p]], list(NULL))
          next
        }
        if (q > p) {
          rho.x[[p]][[q]] <- cor(cbind(X[[p]], X[[q]]))
          rho.xall[[p]][[q]] <- cor(cbind(X_all[[p]], X_all[[q]]))
        } else if (q == p) {
          rho.x[[p]][[q]] <- cor(X[[p]])
          rho.xall[[p]][[q]] <- cor(X_all[[p]])
        } else {
          rho.x[[p]][[q]] <- t(rho.x[[q]][[p]])
          rho.xall[[p]][[q]] <- t(rho.xall[[q]][[p]])
        }
      }
    }
    result <- append(result, list(rho.x = rho.x, rho.xall = rho.xall))
    if (!is.null(Y)) {
      rho.yx <- list()
      rho.yxall <- list()
      for (p in 1:M) {
        if (K.x[p] == 0) {
          rho.yx <- append(rho.yx, list(NULL))
          rho.yxall <- append(rho.yxall, list(NULL))
          next
        }
        rho.yx[[p]] <- cor(cbind(Y, X[[p]]))[1:M, (M + 1):(M + ncol(X[[p]])),
                                             drop = FALSE]
        rho.yxall[[p]] <-
          cor(cbind(Y, X_all[[p]]))[1:M, (M + 1):(M + ncol(X_all[[p]])),
                                    drop = FALSE]
      }
      result <- append(result, list(rho.yx = rho.yx, rho.yxall = rho.yxall))
    }
  }
  if (length(corr.x) > 0) {
    emax2 <- list()
    for (p in 1:M) {
      if (K.x[p] == 0) {
        emax2 <- append(emax2, list(NULL))
        next
      }
      emax2[[p]] <- numeric(M)
      for (q in 1:M) {
        if (K.x[q] == 0) {
          emax2[[p]][q] <- NA
          next
        }
        if (q > p) {
          emax2[[p]][q] <- max(abs(rho.x[[p]][[q]][1:K.x[p],
            (K.x[p] + 1):(K.x[p] + K.x[q])] - corr.x[[p]][[q]]))
        } else if (q == p) {
          emax2[[p]][q] <- max(abs(rho.x[[p]][[q]] - corr.x[[p]][[q]]))
        } else {
          emax2[[p]][q] <- emax2[[q]][p]
        }
      }
    }
    result <- append(result, list(maxerr = emax2))
  }
  if (length(rmeans2) > 0) {
    if (length(rand.int) != M) rand.int <- rep(rand.int[1], M)
    if (length(rand.tsl) != M) rand.tsl <- rep(rand.tsl[1], M)
    M0 <- ifelse(class(corr.u) == "list", M, 1)
    K.rmix <- rep(0, M0)
    K.rcomp <- rep(0, M0)
    if (length(rmix_pis) > 0) {
      K.rmix <- lengths(rmix_pis)
      K.rcmix <- unlist(lapply(rmix_pis, lengths))
      K.rcomp <- sapply(lapply(rmix_pis, unlist), length)
      k.rcomp <- c(0, cumsum(K.rcomp))
    }
    K.rcont <- sapply(rvars2,
      function(x) if (is.null(x[[1]])) 0 else length(x)) - K.rcomp
    k.rcont <- c(0, cumsum(K.rcont))
    k.rmix <- c(0, cumsum(K.rmix))
    rskews2 <- list()
    rskurts2 <- list()
    if (method == "Polynomial") {
      rfifths2 <- list()
      rsixths2 <- list()
    }
    for (i in 1:M0) {
      rskews2[[i]] <- append(rskews2, list(NULL))
      rskurts2[[i]] <- append(rskurts2, list(NULL))
      if (method == "Polynomial") {
        rfifths2[[i]] <- append(rfifths2, list(NULL))
        rsixths2[[i]] <- append(rsixths2, list(NULL))
      }
      if (K.r[i] == 0) next
      if (K.rmix[i] > 0) ind <- c(0, cumsum(lengths(rmix_pis[[i]])))
      if ((rand.int[i] == "none" & rand.tsl[i] == "mix") |
          (rand.int[i] == "mix" & rand.tsl[i] == "none")) {
        if (K.rcont[i] > 0) {
          rskews2[[i]] <- c(rmix_skews[[i]][[1]], rskews[[i]],
            unlist(rmix_skews[[i]][-1]))
          rskurts2[[i]] <- c(rmix_skurts[[i]][[1]], rskurts[[i]],
            unlist(rmix_skurts[[i]][-1]))
          if (method == "Polynomial") {
            rfifths2[[i]] <- c(rmix_fifths[[i]][[1]], rfifths[[i]],
              unlist(rmix_fifths[[i]][-1]))
            rsixths2[[i]] <- c(rmix_sixths[[i]][[1]], rsixths[[i]],
              unlist(rmix_sixths[[i]][-1]))
          }
        } else {
          rskews2[[i]] <- unlist(rmix_skews[[i]])
          rskurts2[[i]] <- unlist(rmix_skurts[[i]])
          if (method == "Polynomial") {
            rfifths2[[i]] <- unlist(rmix_fifths[[i]])
            rsixths2[[i]] <- unlist(rmix_sixths[[i]])
          }
        }
      } else if (rand.int[i] == "non_mix" & rand.tsl[i] == "mix") {
        if (K.rcont[i] > 1) {
          rskews2[[i]] <- c(rskews[[i]][1], rmix_skews[[i]][[1]],
            rskews[[i]][-1], unlist(rmix_skews[[i]][-1]))
          rskurts2[[i]] <- c(rskurts[[i]][1], rmix_skurts[[i]][[1]],
            rskurts[[i]][-1], unlist(rmix_skurts[[i]][-1]))
          if (method == "Polynomial") {
            rfifths2[[i]] <- c(rfifths[[i]][1], rmix_fifths[[i]][[1]],
              rfifths[[i]][-1], unlist(rmix_fifths[[i]][-1]))
            rsixths2[[i]] <- c(rsixths[[i]][1], rmix_sixths[[i]][[1]],
              rsixths[[i]][-1], unlist(rmix_sixths[[i]][-1]))
          }
        } else {
          rskews2[[i]] <- c(rskews[[i]][1], unlist(rmix_skews[[i]]))
          rskurts2[[i]] <- c(rskurts[[i]][1], unlist(rmix_skurts[[i]]))
          if (method == "Polynomial") {
            rfifths2[[i]] <- c(rfifths[[i]][1], unlist(rmix_fifths[[i]]))
            rsixths2[[i]] <- c(rsixths[[i]][1], unlist(rmix_sixths[[i]]))
          }
        }
      } else if (rand.int[i] == "mix" & rand.tsl[i] == "non_mix") {
        rskews2[[i]] <- c(rmix_skews[[i]][[1]], rskews[[i]],
          unlist(rmix_skews[[i]][-1]))
        rskurts2[[i]] <- c(rmix_skurts[[i]][[1]], rskurts[[i]],
          unlist(rmix_skurts[[i]][-1]))
        if (method == "Polynomial") {
          rfifths2[[i]] <- c(rmix_fifths[[i]][[1]], rfifths[[i]],
            unlist(rmix_fifths[[i]][-1]))
          rsixths2[[i]] <- c(rmix_sixths[[i]][[1]], rsixths[[i]],
            unlist(rmix_sixths[[i]][-1]))
        }
      } else if (rand.int[i] == "mix" & rand.tsl[i] == "mix") {
        if (K.rcont[i] > 0) {
          rskews2[[i]] <- c(unlist(rmix_skews[[i]][1:2]),
            rskews[[i]], unlist(rmix_skews[[i]][-c(1:2)]))
          rskurts2[[i]] <- c(unlist(rmix_skurts[[i]][1:2]),
            rskurts[[i]], unlist(rmix_skurts[[i]][-c(1:2)]))
          if (method == "Polynomial") {
            rfifths2[[i]] <- c(unlist(rmix_fifths[[i]][1:2]),
              rfifths[[i]], unlist(rmix_fifths[[i]][-c(1:2)]))
            rsixths2[[i]] <- c(unlist(rmix_sixths[[i]][1:2]),
              rsixths[[i]], unlist(rmix_sixths[[i]][-c(1:2)]))
          }
        } else {
          rskews2[[i]] <- unlist(rmix_skews[[i]])
          rskurts2[[i]] <- unlist(rmix_skurts[[i]])
          if (method == "Polynomial") {
            rfifths2[[i]] <- unlist(rmix_fifths[[i]])
            rsixths2[[i]] <- unlist(rmix_sixths[[i]])
          }
        }
      } else {
        if (K.rcont[i] > 0) {
          rskews2[[i]] <- rskews[[i]]
          rskurts2[[i]] <- rskurts[[i]]
          if (method == "Polynomial") {
            rfifths2[[i]] <- rfifths[[i]]
            rsixths2[[i]] <- rsixths[[i]]
          }
        }
        if (K.rmix[i] > 0) {
          rskews2[[i]] <- c(rskews2[[i]], unlist(rmix_skews[[i]]))
          rskurts2[[i]] <- c(rskurts2[[i]], unlist(rmix_skurts[[i]]))
          if (method == "Polynomial") {
            rfifths2[[i]] <- c(rfifths2[[i]], unlist(rmix_fifths[[i]]))
            rsixths2[[i]] <- c(rsixths2[[i]], unlist(rmix_sixths[[i]]))
          }
        }
      }
    }
    if (max(K.rmix) > 0) {
      U_mix <- list()
      for (i in 1:M0) {
        if (K.rmix[i] == 0) {
          U_mix <- append(U_mix, list(NULL))
          next
        }
        if ((rand.int[i] == "none" & rand.tsl[i] == "mix") |
           (rand.int[i] == "mix" & rand.tsl[i] == "none") |
           (rand.int[i] == "mix" & rand.tsl[i] == "non_mix")) {
          U_mix[[i]] <- U_all[[i]][, 1, drop = FALSE]
          if (K.rmix[i] > 1)
            U_mix[[i]] <- cbind(U_mix[[i]],
              U_all[[i]][, (ncol(U_all[[i]]) - K.rmix[i] +
                               2):ncol(U_all[[i]]), drop = FALSE])
        } else if (rand.int[i] == "non_mix" & rand.tsl[i] == "mix") {
          U_mix[[i]] <- U_all[[i]][, 2, drop = FALSE]
          if (K.rmix[i] > 1)
            U_mix[[i]] <- cbind(U_mix[[i]],
              U_all[[i]][, (ncol(U_all[[i]]) - K.rmix[i] +
                               2):ncol(U_all[[i]])])
        } else if (rand.int[i] == "mix" & rand.tsl[i] == "mix") {
          U_mix[[i]] <- U_all[[i]][, 1:2]
          if (K.rmix[i] > 2)
            U_mix[[i]] <- cbind(U_mix[[i]],
              U_all[[i]][, (ncol(U_all[[i]]) - K.rmix[i] +
                               3):ncol(U_all[[i]]), drop = FALSE])
        } else {
          U_mix[[i]] <- U_all[[i]][, (ncol(U_all[[i]]) - K.rmix[i] +
                               1):ncol(U_all[[i]]), drop = FALSE]
        }
      }
    }
    rho.u <- list()
    for (i in 1:M0) {
      if (K.r[i] == 0) {
        rho.u <- append(rho.u, list(NULL))
        next
      }
      rho.u[[i]] <- list()
      for (j in 1:M0) {
        if (K.r[j] == 0) {
          rho.u[[i]] <- append(rho.u[[i]], list(NULL))
          next
        }
        if (j > i) {
          rhou <- cor(cbind(U[[i]], U[[j]]))
          rho.u[[i]][[j]] <- rhou[1:K.r[i], (K.r[i] + 1):ncol(rhou),
                                  drop = FALSE]
        } else if (j == i) {
          rho.u[[i]][[j]] <- cor(U[[i]])
        } else {
          rho.u[[i]][[j]] <- t(rho.u[[j]][[i]])
        }
      }
    }
    cont_sum_u <- list()
    target_sum_u <- list()
    mix_sum_u <- list()
    target_mix_u <- list()
    sum_uall <- list()
    for (i in 1:M0) {
      sum_uall <- append(sum_uall, list(NULL))
      cont_sum_u <- append(cont_sum_u, list(NULL))
      target_sum_u <- append(target_sum_u, list(NULL))
      mix_sum_u <- append(mix_sum_u, list(NULL))
      target_mix_u <- append(target_mix_u, list(NULL))
      if (K.r[i] > 0) {
        mom <- apply(U_all[[i]], 2, calc_moments)
        medians <- apply(U_all[[i]], 2, median)
        mins <- apply(U_all[[i]], 2, min)
        maxs <- apply(U_all[[i]], 2, max)
        sum_uall[[i]] <- cbind(rep(i, ncol(U_all[[i]])), 1:ncol(U_all[[i]]),
          rep(nrow(U_all[[i]]), ncol(U_all[[i]])), mom[1, ], mom[2, ], medians,
          mins, maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ])
        rownames(sum_uall[[i]]) <- colnames(U_all[[i]])
        mom <- apply(U[[i]], 2, calc_moments)
        medians <- apply(U[[i]], 2, median)
        mins <- apply(U[[i]], 2, min)
        maxs <- apply(U[[i]], 2, max)
        cont_sum_u[[i]] <- cbind(rep(i, K.r[i]), 1:K.r[i],
          rep(nrow(U[[i]]), ncol(U[[i]])), mom[1, ], mom[2, ], medians, mins,
          maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ])
        rownames(cont_sum_u[[i]]) <- paste("cont", rep(i, ncol(U[[i]])),
          "_", 1:ncol(U[[i]]), sep = "")
        if (method == "Fleishman") {
          target_sum_u[[i]] <- cbind(rmeans2[[i]], sqrt(rvars2[[i]]),
                                     rskews2[[i]], rskurts2[[i]])
        } else {
          target_sum_u[[i]] <- cbind(rmeans2[[i]], sqrt(rvars2[[i]]),
            rskews2[[i]], rskurts2[[i]], rfifths2[[i]], rsixths2[[i]])
        }
        target_sum_u[[i]] <- cbind(rep(i, K.r[i]), 1:K.r[i], target_sum_u[[i]])
        rownames(target_sum_u[[i]]) <- rownames(cont_sum_u[[i]])
        if (K.rmix[i] > 0) {
          mom <- apply(U_mix[[i]], 2, calc_moments)
          medians <- apply(U_mix[[i]], 2, median)
          mins <- apply(U_mix[[i]], 2, min)
          maxs <- apply(U_mix[[i]], 2, max)
          mix_sum_u[[i]] <- cbind(rep(i, K.rmix[i]), 1:K.rmix[i],
            rep(nrow(U_mix[[i]]), ncol(U_mix[[i]])), mom[1, ], mom[2, ],
            medians, mins, maxs, mom[3, ], mom[4, ], mom[5, ], mom[6, ])
          rownames(mix_sum_u[[i]]) <- paste("mix", rep(i, K.rmix[i]),
            "_", 1:K.rmix[i], sep = "")
          if (method == "Fleishman") {
            for (j in 1:K.rmix[i]) {
              target_mix_u[[i]] <- rbind(target_mix_u[[i]],
                calc_mixmoments(rmix_pis[[i]][[j]], rmix_mus[[i]][[j]],
                rmix_sigmas[[i]][[j]], rmix_skews[[i]][[j]],
                rmix_skurts[[i]][[j]]))
            }
          } else {
            for (j in 1:K.rmix[i]) {
              target_mix_u[[i]] <- rbind(target_mix_u[[i]],
                calc_mixmoments(rmix_pis[[i]][[j]], rmix_mus[[i]][[j]],
                rmix_sigmas[[i]][[j]], rmix_skews[[i]][[j]],
                rmix_skurts[[i]][[j]], rmix_fifths[[i]][[j]],
                rmix_sixths[[i]][[j]]))
            }
          }
          target_mix_u[[i]] <- cbind(rep(i, K.rmix[i]), 1:K.rmix[i],
                                     target_mix_u[[i]])
          rownames(target_mix_u[[i]]) <- rownames(mix_sum_u[[i]])
        }
      }
    }
    cont_sum_u <- as.data.frame(do.call(rbind, cont_sum_u))
    colnames(cont_sum_u) <- c("Outcome", "U", "N", "Mean", "SD", "Median",
      "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
    sum_uall <- as.data.frame(do.call(rbind, sum_uall))
    colnames(sum_uall) <- c("Outcome", "U", "N", "Mean", "SD", "Median",
      "Min", "Max", "Skew", "Skurtosis", "Fifth", "Sixth")
    target_sum_u <- as.data.frame(do.call(rbind, target_sum_u))
    if (method == "Fleishman") {
      colnames(target_sum_u) <- c("Outcome", "U", "Mean", "SD", "Skew",
                                  "Skurtosis")
    } else {
      colnames(target_sum_u) <- c("Outcome", "U", "Mean", "SD", "Skew",
                                  "Skurtosis", "Fifth", "Sixth")
    }
    result <- append(result, list(target_sum_u = target_sum_u,
      cont_sum_u = cont_sum_u, sum_uall = sum_uall, rho.u = rho.u))
    if (max(K.rmix) > 0) {
      rho.uall <- list()
      for (i in 1:M0) {
        if (K.rmix[i] == 0) {
          rho.uall <- append(rho.uall, list(NULL))
          next
        }
        rho.uall[[i]] <- list()
        for (j in 1:M0) {
          if (K.rmix[j] == 0) {
            rho.uall[[i]] <- append(rho.uall[[i]], list(NULL))
            next
          }
          if (j >= i) {
            rhou <- cor(cbind(U_all[[i]], U_all[[j]]))
            if (j == i) rho.uall[[i]][[j]] <-
                rhou[1:ncol(U_all[[i]]), 1:ncol(U_all[[i]]), drop = FALSE]
            if (j > i) rho.uall[[i]][[j]] <-
                rhou[1:ncol(U_all[[i]]), (ncol(U_all[[i]]) + 1):ncol(rhou),
                     drop = FALSE]
          }
          if (j < i) rho.uall[[i]][[j]] <- t(rho.uall[[j]][[i]])
        }
      }
      target_mix_u <- as.data.frame(do.call(rbind, target_mix_u))
      colnames(target_mix_u) <- colnames(target_sum_u)
      mix_sum_u <- as.data.frame(do.call(rbind, mix_sum_u))
      colnames(mix_sum_u) <- colnames(cont_sum_u)
      result <- append(result, list(rho.uall = rho.uall,
        target_mix_u = target_mix_u, mix_sum_u = mix_sum_u))
    }
  }
  if (class(corr.u) == "matrix") {
    emax3 <- max(abs(rho.u[[1]][[1]] - corr.u))
    result <- append(result, list(maxerr_u = emax3))
  }
  if (class(corr.u) == "list" & length(corr.u) > 0) {
    emax3 <- list()
    for (i in 1:M) {
      if (K.r[i] == 0) {
        emax3 <- append(emax3, list(NULL))
        next
      }
      emax3[[i]] <- numeric(M)
      for (j in 1:M) {
        if (K.r[j] == 0) {
          emax3[[i]][j] <- NA
          next
        }
        if (j > i) {
          emax3[[i]][j] <- max(abs(rho.u[[i]][[j]] - corr.u[[i]][[j]]))
        } else if (j == i) {
          emax3[[i]][j] <- max(abs(rho.u[[i]][[j]] - corr.u[[i]][[j]]))
        } else {
          emax3[[i]][j] <- emax3[[j]][i]
        }
      }
    }
    result <- append(result, list(maxerr_u = emax3))
  }
  result
}

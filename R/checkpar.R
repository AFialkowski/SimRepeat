#' @title Parameter Check for Simulation Functions
#'
#' @description This function checks the parameter inputs to the simulation functions \code{\link[SimRepeat]{nonnormsys}},
#'     \code{\link[SimRepeat]{corrsys}}, and \code{\link[SimRepeat]{corrsys2}}.  It should be used prior to execution of these
#'     functions to ensure all inputs are of the correct format.  Those functions do not contain parameter checks in order to decrease
#'     simulation time.  This would be important if the user is running several simulation repetitions so that the inputs only have to
#'     be checked once.  Note that the inputs do not include all of the inputs to the simulation functions.  See the appropriate function
#'     documentation for more details about parameter inputs.  Since the parameter input list is extensive and this function does not check
#'     for all possible errors, if simulation gives an error, the user should still check the parameter inputs.
#'
#' @param M the number of dependent variables \eqn{Y} (outcomes); equivalently, the number of equations in the system
#' @param method the PMT method used to generate all continuous variables, including independent variables (covariates), error terms, and random effects;
#'     "Fleishman" uses Fleishman's third-order polynomial transformation and "Polynomial" uses Headrick's fifth-order transformation
#' @param error_type "non_mix" if all error terms have continuous non-mixture distributions, "mix" if all error terms have continuous mixture distributions
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
#' @param Six a list of length \code{M}, \code{M + 1}, or \code{2 * M}, where \code{Six[1:M]} are for \eqn{X_{cont}, E} (if \code{error_type = "non_mix"}) and
#'     \code{Six[M + 1]} or \code{Six[(M + 1):(2 * M)]} are for non-mixture \eqn{U};
#'     if \code{error_type = "mix"} and there are only random effects (i.e., \code{length(corr.x) = 0}), use \code{Six[1:M] = rep(list(NULL), M)} so that
#'     \code{Six[M + 1]} or \code{Six[(M + 1):(2 * M)]} describes the non-mixture \eqn{U};
#'
#'     \code{Six[[p]][[j]]} is a vector of sixth cumulant correction values to aid in finding a valid PDF for \eqn{X_{cont(pj)}}, the
#'     j-th continuous non-mixture covariate for outcome \eqn{Y_p}; the last vector in \code{Six[[p]]} is for \eqn{E_p} (if \code{error_type = "non_mix"});
#'     use \code{Six[[p]][[j]] = NULL} if no correction desired for \eqn{X_{cont(pj)}};
#'     use \code{Six[[p]] = NULL} if no correction desired for any continuous non-mixture covariate or error term in equation p
#'
#'     \code{Six[[M + p]][[j]]} is a vector of sixth cumulant correction values to aid in finding a valid PDF for \eqn{U_{(pj)}}, the
#'     j-th non-mixture random effect for outcome \eqn{Y_p}; use \code{Six[[M + p]][[j]] = NULL} if no correction desired for \eqn{U_{(pj)}};
#'     use \code{Six[[M + p]] = NULL} if no correction desired for any continuous non-mixture random effect in equation p
#'
#'     keep \code{Six = list()} if no corrections desired for all equations or if \code{method = "Fleishman"}
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
#' @param mix_Six a list of same length and order as \code{mix_pis};
#'     keep \code{mix_Six = list()} if no corrections desired for all equations or if \code{method = "Fleishman"}
#'
#'     p-th component of \code{mix_Six[1:M]} is a list of length equal to the total number of component distributions for the \eqn{X_{mix(p)}}
#'     and \eqn{E_p} (if \code{error_type = "mix"});
#'     \code{mix_Six[[p]][[j]]} is a vector of sixth cumulant corrections for the j-th component distribution (i.e., if there are 2
#'     continuous mixture independent variables for \eqn{Y_p}, where \eqn{X_{mix(p1)}} has 2 components and \eqn{X_{mix(p2)}} has 3
#'     components, then \code{length(mix_Six[[p]]) = 5} and \code{mix_Six[[p]][[3]]} would correspond to the 1st component of
#'     \eqn{X_{mix(p2)}}); use \code{mix_Six[[p]][[j]] = NULL} if no correction desired for that component;
#'     use \code{mix_Six[[p]] = NULL} if no correction desired for any component of \eqn{X_{mix(p)}} and \eqn{E_p}
#'
#'     q-th component of \code{mix_Six[M + 1]} or \code{mix_Six[(M + 1):(2 * M)]} is a list of length equal to the total number of component distributions for
#'     the \eqn{U_{mix(q)}}; \code{mix_Six[[q]][[j]]} is a vector of sixth cumulant corrections for the j-th component distribution; use
#'     \code{mix_Six[[q]][[j]] = NULL} if no correction desired for that component;
#'     use \code{mix_Six[[q]] = NULL} if no correction desired for any component of \eqn{U_{mix(q)}}
#'
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
#'     (see \code{\link[stats;Poisson]{dpois}}); order is 1st regular Poisson and 2nd zero-inflated Poisson; use \code{lam[[p]] = NULL} if outcome \eqn{Y_p} has no Poisson variables;
#'     \code{length(lam[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}
#' @param p_zip a list of vectors of probabilities of structural zeros (not including zeros from the Poisson distribution) for the
#'     zero-inflated Poisson variables (see \code{\link[VGAM;Zipois]{dzipois}}); if \code{p_zip} = 0, \eqn{Y_{pois}} has a regular Poisson
#'     distribution; if \code{p_zip} is in (0, 1), \eqn{Y_{pois}} has a zero-inflated Poisson distribution;
#'     if \code{p_zip} is in \code{(-(exp(lam) - 1)^(-1), 0)}, \eqn{Y_{pois}} has a zero-deflated Poisson distribution and \code{p_zip}
#'     is not a probability; if \code{p_zip = -(exp(lam) - 1)^(-1)}, \eqn{Y_{pois}} has a positive-Poisson distribution
#'     (see \code{\link[VGAM;Pospois]{dpospois}}); order is 1st regular Poisson and 2nd zero-inflated Poisson;
#'     if a single number, all Poisson variables given this value; if a vector of length \code{M}, all Poisson variables in equation p
#'     given \code{p_zip[p]}; otherwise, missing values are set to 0 and ordered 1st
#' @param pois_eps list of length \code{M}, p-th component a vector of length \code{lam[[p]]} containing cumulative probability truncation values
#'     used to calculate intermediate correlations involving Poisson variables; order is 1st regular Poisson and 2nd zero-inflated Poisson;
#'     if a single number, all Poisson variables given this value; if a vector of length \code{M}, all Poisson variables in equation p
#'     given \code{pois_eps[p]}; otherwise, missing values are set to 0.0001 and ordered 1st
#' @param size list of length \code{M}, p-th component a vector of size parameters for the Negative Binomial variables for outcome \eqn{Y_p}
#'     (see \code{\link[stats;NegBinomial]{dnbinom}}); order is 1st regular NB and 2nd zero-inflated NB; use \code{size[[p]] = NULL} if outcome \eqn{Y_p} has no Negative Binomial variables;
#'     \code{length(size[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}
#' @param prob list of length \code{M}, p-th component a vector of success probabilities for the Negative Binomial variables for outcome \eqn{Y_p}
#'     (see \code{\link[stats;NegBinomial]{dnbinom}}); order is 1st regular NB and 2nd zero-inflated NB; use \code{prob[[p]] = NULL} if outcome \eqn{Y_p} has no Negative Binomial variables;
#'     \code{length(prob[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}
#' @param mu list of length \code{M}, p-th component a vector of mean values for the Negative Binomial variables for outcome \eqn{Y_p}
#'     (see \code{\link[stats;NegBinomial]{dnbinom}}); order is 1st regular NB and 2nd zero-inflated NB; use \code{mu[[p]] = NULL} if outcome \eqn{Y_p} has no Negative Binomial variables;
#'     \code{length(mu[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}; for zero-inflated NB variables,
#'     this refers to the mean of the NB distribution (see \code{\link[VGAM;Zinegbin]{dzinegbin}})
#'     (*Note: either \code{prob} or \code{mu} should be supplied for all Negative Binomial variables, not a mixture)
#' @param p_zinb a vector of probabilities of structural zeros (not including zeros from the NB distribution) for the zero-inflated NB variables
#'     (see \code{\link[VGAM;Zinegbin]{dzinegbin}}); if \code{p_zinb} = 0, \eqn{Y_{nb}} has a regular NB distribution;
#'     if \code{p_zinb} is in \code{(-prob^size/(1 - prob^size),} \code{0)}, \eqn{Y_{nb}} has a zero-deflated NB distribution and \code{p_zinb}
#'     is not a probability; if \code{p_zinb = -prob^size/(1 - prob^size)}, \eqn{Y_{nb}} has a positive-NB distribution (see
#'     \code{\link[VGAM;Posnegbin]{dposnegbin}}); order is 1st regular NB and 2nd zero-inflated NB;
#'     if a single number, all NB variables given this value; if a vector of length \code{M}, all NB variables in equation p
#'     given \code{p_zinb[p]}; otherwise, missing values are set to 0 and ordered 1st
#' @param nb_eps list of length \code{M}, p-th component a vector of length \code{size[[p]]} containing cumulative probability truncation values
#'     used to calculate intermediate correlations involving Negative Binomial variables; order is 1st regular NB and 2nd zero-inflated NB;
#'     if a single number, all NB variables given this value; if a vector of length \code{M}, all NB variables in equation p
#'     given \code{nb_eps[p]}; otherwise, missing values are set to 0.0001 and ordered 1st
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
#' @param corr.yx input for \code{nonnormsys} only; a list of length \code{M}, where the p-th component is a 1 row matrix of correlations between \eqn{Y_p} and \eqn{X_{(pj)}};
#'     if there are mixture variables and the \code{betas} are desired in terms of these (and not the components), then \code{corr.yx}
#'     should be specified in terms of correlations between outcomes and non-mixture or mixture variables, and the number of columns of the matrices
#'     of \code{corr.yx} should not match the dimensions of the matrices in \code{corr.x}; if the \code{betas} are desired in terms of
#'     the components, then \code{corr.yx} should be specified in terms of correlations between outcomes and non-mixture or components of
#'     mixture variables, and the number of columns of the matrices of \code{corr.yx} should match the dimensions of the matrices in \code{corr.x};
#'     use \code{corr.yx[[p]] = NULL} if equation p has no \eqn{X_{(pj)}}
#' @param corr.e correlation matrix for continuous non-mixture or components of mixture error terms
#' @param same.var either a vector or a matrix; if a vector, \code{same.var} includes column numbers of \code{corr.x[[1]][[1]]}
#'     corresponding to independent variables that should be identical across equations; these terms must have the same indices for all
#'     \code{p = 1, ..., M}; i.e., if the 1st ordinal variable represents sex which should be the same for each equation, then
#'     \code{same.var[1] = 1} since ordinal variables are 1st in \code{corr.x[[1]][[1]]} and sex is the 1st ordinal variable, and
#'     the 1st term for all other outcomes must also be sex;
#'     if a matrix, columns 1 and 2 are outcome p and column index in \code{corr.x[[p]][[p]]} for 1st instance of variable,
#'     columns 3 and 4 are outcome q and column index in \code{corr.x[[q]][[q]]} for subsequent instances of variable; i.e., if
#'     1st term for all outcomes is sex and \code{M = 3}, then \code{same.var = matrix(c(1, 1, 2, 1, 1, 1, 3, 1), 2, 4, byrow = TRUE)}; the
#'     independent variable index corresponds to ordinal, continuous non-mixture, component of continuous mixture, Poisson, or
#'     NB variable
#' @param subj.var matrix where 1st column is outcome index (\code{p = 1, ..., M}), 2nd column is independent variable index
#'     corresponding to covariate which is a a subject-level term (not including time), including time-varying covariates;
#'     the independent variable index corresponds to ordinal, continuous non-mixture, continuous mixture (not mixture component), Poisson, or
#'     NB variable; assumes all other variables are group-level terms; these subject-level terms are used to form interactions with the group
#'     level terms
#' @param int.var matrix where 1st column is outcome index (\code{p = 1, ..., M}), 2nd and 3rd columns are indices
#'     corresponding to independent variables to form interactions between; this includes all interactions that are not accounted for by
#'     a subject-group level interaction (as indicated by \code{subj.var}) or by a time-covariate interaction (as indicated by
#'     \code{tint.var}); ex: 1, 2, 3 indicates that for outcome 1, the 2nd and 3rd independent variables form an interaction term;
#'     the independent variable index corresponds to ordinal, continuous non-mixture, continuous mixture (not mixture component), Poisson, or
#'     NB variable
#' @param tint.var matrix where 1st column is outcome index (\code{p = 1, ..., M}), 2nd column is index of
#'     independent variable to form interaction with time; if \code{tint.var = NULL} or no \eqn{X_{(pj)}} are indicated for outcome \eqn{Y_p},
#'     this includes all group-level variables (variables not indicated as subject-level variables in \code{subj.var}), else includes only
#'     terms indicated by 2nd column (i.e., in order to include subject-level variables);
#'     ex: 1, 1 indicates that for outcome 1, the 1st independent variable has an interaction with time;
#'     the independent variable index corresponds to ordinal, continuous non-mixture, continuous mixture (not mixture component), Poisson, or
#'     NB variable
#' @param betas.0 vector of length \code{M} containing intercepts, if \code{NULL} all set equal to 0; if length 1, all intercepts set to
#'     \code{betas.0}
#' @param betas list of length \code{M}, p-th component a vector of coefficients for outcome \eqn{Y_p}, including group and subject-level terms;
#'     order is order of variables in \code{corr.x[[p]][[p]]}; if \code{betas = list()}, all set to 0 so that all \eqn{Y} only have intercept
#'     and/or interaction terms plus error terms; if all outcomes have the same betas, use list of length 1; if \eqn{Y_p} only has intercept
#'     and/or interaction terms, set \code{betas[[p]] = NULL}; if there are continuous mixture variables, beta is for mixture variable
#'     (not for components)
#' @param betas.subj list of length \code{M}, p-th component a vector of coefficients for interaction terms between group-level terms and
#'     subject-level terms given in \code{subj.var}; order is the same order as given in \code{subj.var}; if all outcomes have the same betas,
#'     use list of length 1; if \eqn{Y_p} only has group-level terms, set \code{betas.subj[[p]] = NULL};
#'     if there are continuous mixture variables, beta is for mixture variable (not for components)
#' @param betas.int list of length \code{M}, p-th component a vector of coefficients for interaction terms indicated in \code{int.var};
#'     order is the same order as given in \code{int.var}; if all outcomes have the same betas, use list of length 1;
#'     if \eqn{Y_p} has none, set \code{betas.int[[p]] = NULL};
#'     if there are continuous mixture variables, beta is for mixture variable (not for components)
#' @param betas.t vector of length \code{M} of coefficients for time terms, if \code{NULL} all set equal to 1;
#'     if length 1, all intercepts set to \code{betas.t}
#' @param betas.tint list of length \code{M}, p-th component a vector of coefficients for interaction terms indicated in \code{tint.var};
#'     order is the same order as given in \code{tint.var}; if all outcomes have the same betas, use list of length 1;
#'     if \eqn{Y_p} has none, set \code{betas.tint[[p]] = NULL};
#'     if there are continuous mixture variables, beta is for mixture variable (not for components)
#' @param rand.int "none" (default) if no random intercept term for all outcomes, "non_mix" if all random intercepts have a continuous
#'     non-mixture distribution, "mix" if all random intercepts have a continuous mixture distribution;
#'     also can be a vector of length \code{M} containing a combination (i.e., \code{c("non_mix", "mix", "none")} if the 1st has a non-mixture
#'     distribution, the 2nd has a mixture distribution, and 3rd outcome has no random intercept)
#' @param rand.tsl "none" (default) if no random slope for time for all outcomes, "non_mix" if all random time slopes have a
#'     continuous non-mixture distribution, "mix" if all random time slopes have a continuous mixture distribution; also can
#'     be a vector of length \code{M} as in \code{rand.int}
#' @param rand.var matrix where 1st column is outcome index (\code{p = 1, ..., M}), 2nd column is independent variable index
#'     corresponding to covariate to assign random effect to (not including the random intercept or time slope if present);
#'     the independent variable index corresponds to ordinal, continuous non-mixture, continuous mixture (not mixture component), Poisson, or
#'     NB variable; order is 1st continuous non-mixture and 2nd continuous mixture random effects; note that the order of the rows corresponds
#'     to the order of the random effects in \code{corr.u} not the order of the
#'     independent variable so that a continuous mixture covariate with a non-mixture random effect would be ordered before a
#'     continuous non-mixture covariate with a mixture random effect (the 2nd column of \code{rand.var} indicates the specific covariate)
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
#' @param quiet if FALSE prints messages, if TRUE suppresses messages
#'
#' @import utils
#' @export
#' @return TRUE if all inputs are correct, else it will stop with a correction message
#' @keywords ParameterCheck
#' @seealso \code{\link[SimRepeat]{nonnormsys}}, \code{\link[SimRepeat]{corrsys}}, \code{\link[SimRepeat]{corrsys2}}
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
#'
#' # Examples for the corrsys and corrsys2 functions can be found in the
#' # function documentation.
#' }
#'
#' @references
#' Headrick TC, Beasley TM (2004).  A Method for Simulating Correlated Non-Normal Systems of Linear Statistical Equations.
#'     Communications in Statistics - Simulation and Computation, 33(1).  \doi{10.1081/SAC-120028431}
#'
checkpar <- function(M = NULL, method = c("Fleishman", "Polynomial"),
                     error_type = c("non_mix", "mix"), means =  list(),
                     vars =  list(), skews =  list(), skurts =  list(),
                     fifths =  list(), sixths =  list(), Six = list(),
                     mix_pis = list(), mix_mus = list(), mix_sigmas = list(),
                     mix_skews =  list(), mix_skurts =  list(),
                     mix_fifths =  list(), mix_sixths =  list(),
                     mix_Six = list(), marginal = list(), support = list(),
                     lam  =  list(), p_zip = list(), pois_eps = list(),
                     size = list(), prob = list(), mu = list(),
                     p_zinb = list(), nb_eps = list(), corr.x = list(),
                     corr.yx = list(), corr.e = NULL, same.var = NULL,
                     subj.var = NULL, int.var = NULL, tint.var = NULL,
                     betas.0 = NULL, betas = list(), betas.subj = list(),
                     betas.int = list(), betas.t = NULL, betas.tint = list(),
                     rand.int = c("none", "non_mix", "mix"),
                     rand.tsl = c("none", "non_mix", "mix"), rand.var = NULL,
                     corr.u = list(), quiet = FALSE) {
  if (length(error_type) != 1)
    stop("Please choose one type of distribution for all of the error terms:
         mix if all errors have continuous mixture distributions,
         non_mix if all errors have continuous non-mixture distributions.")
  if (length(method) != 1)
    stop("Choose a PMT method for the continuous variables.")
  K.cat <- rep(0, M)
  K.pois <- rep(0, M)
  K.nb <- rep(0, M)
  if (length(marginal) > 0) K.cat <- lengths(marginal)
  if (max(K.cat) > 0) {
    if (!all(unlist(lapply(marginal, function(x) unlist(lapply(x,
      function(y) (sort(y) == y & min(y) > 0 & max(y) < 1)))))))
      stop("Error in given marginal distributions.")
    if (length(support) > 0) {
      if (!all(K.cat %in% lengths(support)))
        stop("Lengths of support do not match lengths of marginal.")
      k.cat <- unlist(lapply(marginal, function(x) lengths(x) + 1))
      k.cat2 <- unlist(lapply(support, function(x) lengths(x)))
      if (!all(mapply(function (x, y) (if (y == 0 | y == x) TRUE else FALSE),
                      k.cat, k.cat2)))
        stop("Error in given support values.")
    }
  }
  if (length(lam) > 0) K.pois <- lengths(lam)
  if (max(K.pois) > 0) {
    if (sum(sapply(lam, function(x) sum(x < 0)) > 0) > 0)
      stop("Lambda values cannot be negative.")
    if (class(p_zip) == "numeric" & length(p_zip) == 1 & quiet == FALSE) {
      message("All p_zip values will be set at p_zip.")
    } else if (class(p_zip) == "numeric" & length(p_zip) == M & quiet == FALSE) {
      message("All p_zip[[p]] values will be set at p_zip[p].")
    } else if (class(p_zip) == "list" & length(p_zip) == M & quiet == FALSE) {
      for (p in 1:M) {
        if (length(p_zip[[p]]) != K.pois[p])
          message(paste("Missing value of p_zip for equation ", p,
            " will be set at 0.", sep = ""))
      }
    } else if (quiet == FALSE) {
      message("All p_zip values will be set at 0.")
    }
    if (class(pois_eps) == "numeric" & length(pois_eps) == 1 &
        quiet == FALSE) {
      message("All pois_eps values will be set at pois_eps.")
    } else if (class(pois_eps) == "numeric" & length(pois_eps) == M &
               quiet == FALSE) {
      message("All pois_eps[[p]] values will be set at pois_eps[p].")
    } else if (class(pois_eps) == "list" & length(pois_eps) == M &
               quiet == FALSE) {
      for (p in 1:M) {
        if (length(pois_eps[[p]]) != K.pois[p])
          message(paste("Missing value of pois_eps for equation ", p,
            " will be set at 0.0001.", sep = ""))
      }
    } else if (quiet == FALSE) {
      message("All pois_eps values will be set at 0.0001.")
    }
  }
  if (length(size) > 0) K.nb <- lengths(size)
  if (max(K.nb) > 0) {
    if (sum(sapply(size, function(x) sum(x < 0)) > 0) > 0)
      stop("Size values cannot be negative.")
    if (length(prob) > 0) {
      if (sum(sapply(prob, function(x) sum(x < 0)) > 0) > 0)
        stop("Prob values cannot be negative.")
      if (!all(lengths(prob) %in% K.nb))
        stop("Check lengths of prob vectors.")
    }
    if (length(mu) > 0) {
      if (sum(sapply(mu, function(x) sum(x < 0)) > 0) > 0)
        stop("mu values cannot be negative.")
      if (!all(lengths(mu) %in% K.nb))
        stop("Check lengths of mu vectors.")
    }
    if (class(p_zinb) == "numeric" & length(p_zinb) == 1 & quiet == FALSE) {
      message("All p_zinb values will be set at p_zinb.")
    } else if (class(p_zinb) == "numeric" & length(p_zinb) == M &
               quiet == FALSE) {
      message("All p_zinb[[p]] values will be set at p_zinb[p].")
    } else if (class(p_zinb) == "list" & length(p_zinb) == M &
               quiet == FALSE) {
      for (p in 1:M) {
        if (length(p_zinb[[p]]) != K.nb[p])
          message(paste("Missing value of p_zinb for equation ", p,
            " will be set at 0.", sep = ""))
      }
    } else if (quiet == FALSE) {
      message("All p_zinb values will be set at 0.")
    }
    if (class(nb_eps) == "numeric" & length(nb_eps) == 1 & quiet == FALSE) {
      message("All nb_eps values will be set at nb_eps.")
    } else if (class(nb_eps) == "numeric" & length(nb_eps) == M &
               quiet == FALSE) {
      message("All nb_eps[[p]] values will be set at nb_eps[p].")
    } else if (class(nb_eps) == "list" & length(nb_eps) == M &
               quiet == FALSE) {
      for (p in 1:M) {
        if (length(nb_eps[[p]]) != K.nb[p])
          message(paste("Missing value of nb_eps for equation ", p,
            " will be set at 0.0001.", sep = ""))
      }
    } else if (quiet == FALSE) {
      message("All nb_eps values will be set at 0.0001.")
    }
  }
  if (!(length(means) %in% c(M, 1 + M, 2 * M)) |
      !(length(vars) %in% c(M, 1 + M, 2 * M)))
    stop("Each equation must have a continuous error term with means given in
         means and variances in vars.  Length of means and vars
         should be M, 1 + M, or 2 * M.")
  M0 <- length(means)
  K.mix <- rep(0, M0)
  K.comp <- rep(0, M0)
  if (length(mix_pis) > 0) {
    if (length(mix_pis) == M0) {
      K.mix <- lengths(mix_pis)
      K.comp <- lengths(lapply(mix_pis, unlist))
    } else {
      K.mix <- c(lengths(mix_pis), rep(0, M0 - M))
      K.comp <- c(lengths(lapply(mix_pis, unlist)), rep(0, M0 - M))
    }
  }
  K.cont <- lengths(vars) - K.mix
  if (!all(lengths(means) %in% (K.cont + K.mix)))
    stop("Lengths of mean vectors should be same as lengths of vars vectors.")
  if (max(K.cont) > 0) {
    if (!all(lengths(skews) %in% K.cont) | !all(lengths(skurts) %in% K.cont))
      stop("Lengths of continuous non-mixture lists should equal lengths of
           vars minus the number of continuous mixture variables.")
    if (method == "Polynomial") {
      if (!all(lengths(fifths) %in% K.cont) |
          !all(lengths(sixths) %in% K.cont))
        stop("Lengths of continuous non-mixture lists should equal lengths of
             vars minus the number of continuous mixture variables.")
      if (!(length(Six) %in% c(0, M, 1 + M, 2 * M)))
        stop("Six should be either list() or a list of length M, M + 1,
             or 2 * M.")
      if (length(Six) != 0) {
        for (i in 1:length(Six)) {
          if (!(length(Six[[i]]) %in% c(0, K.cont[i])))
            stop("The i-th element of Six should be either NULL or a list of
                 the same length as the i-th element of sixths.")
        }
      }
    }
  }
  if (length(mix_pis) > 0) {
    if (!all(lengths(mix_mus) %in% K.mix) |
        !all(lengths(mix_sigmas) %in% K.mix) |
        !all(lengths(mix_skews) %in% K.mix) |
        !all(lengths(mix_skurts) %in% K.mix))
      stop("Lengths of mixing parameter lists should be equal to the number
           of mixture variables.")
    if (method == "Polynomial") {
      if (!all(lengths(mix_fifths) %in% K.mix) |
          !all(lengths(mix_sixths) %in% K.mix))
        stop("Lengths of mixing parameter lists should be equal to the number
             of mixture variables.")
      if (!(length(mix_Six) %in% c(0, M, M + 1, 2 * M)))
        stop("mix_Six should be either list() or a list of length M, M + 1,
             or 2 * M.")
      if (length(mix_Six) != 0) {
        for (i in 1:length(mix_Six)) {
          if (!(length(mix_Six[[i]]) %in% c(0, K.comp[i])))
            stop("The i-th element of mix_Six should be either NULL or a list
                 of the same length as the i-th element of mix_sixths.")
        }
      }
    }
  }
  if (length(means) %in% c(2 * M, M + 1)) {
    means <- means[1:M]
    vars <- vars[1:M]
  }
  if (length(mix_pis) %in% c(2 * M, M + 1)) mix_pis <- mix_pis[1:M]
  K.mix <- rep(0, M)
  K.comp <- rep(0, M)
  K.error <- rep(0, M)
  if (length(mix_pis) > 0) {
    K.mix <- lengths(mix_pis)
    K.comp <- sapply(lapply(mix_pis, unlist), length)
  }
  K.cont <- lengths(vars) - K.mix
  k.mix <- c(0, cumsum(K.mix))
  if (error_type == "mix") {
    K.error <- sapply(mix_pis, function(x) lengths(x[length(x)]))
    K.x <- K.cont + K.comp - K.error + K.cat + K.pois + K.nb
    K.cont2 <- K.cont
    K.mix2 <- K.mix - 1
  } else {
    K.x <- K.cont + K.comp - 1 + K.cat + K.pois + K.nb
    K.cont2 <- K.cont - 1
    K.mix2 <- K.mix
  }
  if (sum(lengths(means) == 0) > 0 | sum(lengths(vars) == 0) > 0)
    stop("Each equation must have a continuous error term with means given in
         means and variances in vars.")
  if (length(corr.x) > 0) {
    if (length(corr.x) != M) stop("corr.x should be list of length M.")
    if (!all(lengths(corr.x) %in% c(0, M)))
      stop("corr.x should be list of lists of length M.")
    for (i in 1:M) {
      if (K.x[i] == 0) next
      if (!all(sapply(corr.x[[i]], function(x) if (is.null(x)) 0
                      else nrow(x)) %in% c(0, K.x[i])))
        stop("corr.x[[i]] should be list of matrices with number of rows equal
             to number of independent variables for equation i.")
      for (j in 1:M) {
        if (K.x[j] == 0) next
        if (ncol(corr.x[[i]][[j]]) != K.x[j])
          stop("corr.x[[i]][[j]] should be matrix with number of cols equal
               to number of independent variables for equation j.")
      }
    }
  }
  if (length(corr.yx) > 0) {
    if (length(corr.yx) != M)
      stop("corr.yx should be a list of length M.")
    if (!all(sapply(corr.yx, function(x) if (is.null(x)) 0
                    else nrow(x)) %in% c(0, 1)))
      stop("corr.yx should be a list of matrices with 1 row each.")
    if (!all(sapply(corr.yx, function(x) if (is.null(x)) 0
                    else ncol(x)) %in% K.x) &
        !all(sapply(corr.yx, function(x) if (is.null(x)) 0
                    else ncol(x)) %in% (K.cont2 + K.mix2)))
      stop("corr.yx should be a list of matrices with 1 row each.
           The number of columns should either equal the number of non-mixture
           plus components of mixture variables or the number of non-mixture
           plus mixture variables.")
  }
  K.r <- rep(0, M)
  if (class(corr.u) == "list" & length(corr.u) > 0) {
    K.r <- sapply(mapply('[', corr.u, lengths(corr.u), SIMPLIFY = FALSE),
                  function(x) if (is.null(x)) 0 else nrow(x[[1]]))
    if (length(corr.u) != M) stop("corr.u should be list of length M.")
    if (!all(lengths(corr.u) %in% c(0, M)))
      stop("corr.u should be list of lists of length M.")
    for (i in 1:M) {
      if (K.r[i] == 0) next
      if (!all(sapply(corr.u[[i]],
        function(x) if (is.null(x)) 0 else nrow(x)) %in% c(0, K.r[i])))
        stop("corr.u[[i]] should be list of matrices with number of rows equal
             to number of random effects for equation i.")
      for (j in 1:M) {
        if (K.r[j] == 0) next
        if (ncol(corr.u[[i]][[j]]) != K.r[j])
          stop("corr.u[[i]][[j]] should be matrix with number of cols equal
               to number of random effects for equation j.")
      }
    }
  }
  if (max(K.r) > 0) {
    if (length(rand.int) != M & quiet == FALSE)
      message("The random intercept for all equations will be set at
              rand.int[1].")
    if (length(rand.tsl) != M & quiet == FALSE)
      message("The random time slope for all equations will be set at
              rand.tsl[1].")
    if (!is.null(rand.var)) {
      if (ncol(rand.var) != 2)
        stop("rand.var should be matrix with 2 columns.")
    }
  }
  if (quiet == FALSE) {
    if (is.null(betas.0)) message("The intercepts will all be set at 0.")
    if (length(betas.0) < M & length(betas.0) > 0)
      message("The intercepts will all be set at betas.0[1].")
    if (length(corr.yx) == 0 & is.null(betas.t))
      message("The time slopes will all be set at 1.")
    if (length(betas.t) < M & length(betas.t) > 0)
      message("The time slopes will all be set at betas.t[1].")
  }
  if (!is.null(same.var)) {
    if (class(same.var) == "matrix") {
      if (ncol(same.var) != 4)
        stop("same.var should be matrix with 4 columns.")
    } else if (class(same.var) == "numeric") {
      if (length(same.var) > ncol(corr.x[[1]][[1]]))
        stop("same.var should be vector with at most ncol(corr.x[[1]][[1]])
             elements.")
    } else stop("same.var should be vector or matrix.")
  }
  if (!is.null(subj.var)) {
    if (ncol(subj.var) != 2) stop("subj.var should be matrix with 2 columns.")
  }
  if (!is.null(int.var)) {
    if (ncol(int.var) != 3) stop("int.var should be matrix with 3 columns.")
  }
  if (!is.null(tint.var)) {
    if (ncol(tint.var) != 2) stop("tint.var should be matrix with 2 columns.")
  }
  if (length(corr.yx) == 0 & length(betas) == 0 & quiet == FALSE)
    message("All betas will be set to 0 if using corrsys or corrsys2.")
  if (class(betas) == "list" & length(betas) == 1 & quiet == FALSE)
    message("All betas will be set to betas[[1]].")
  if (class(betas) != "list" |
      (class(betas) == "list" & !(length(betas) %in% c(0, 1, M))))
    stop("Betas should be list of length 0, 1 or M.")
  if (length(betas.subj) == 0 & !is.null(subj.var) & quiet == FALSE)
    message("All betas.subj will be set to 0 if using corrsys or corrsys2.")
  if (class(betas.subj) == "list" & length(betas.subj) == 1 & quiet == FALSE)
    message("All betas.subj will be set to betas.subj[[1]] if using corrsys or
            corrsys2.")
  if (class(betas.subj) != "list" |
      (class(betas.subj) == "list" & !(length(betas.subj) %in% c(0, 1, M))))
    stop("betas.subj should be list of length 0, 1 or M.")
  if (length(betas.int) == 0 & !is.null(int.var) & quiet == FALSE)
    message("All betas.int will be set to 0 if using corrsys or corrsys2.")
  if (class(betas.int) == "list" & length(betas.int) == 1 & quiet == FALSE)
    message("All betas.int will be set to betas.int[[1]] if using corrsys or
            corrsys2.")
  if (class(betas.int) != "list" |
      (class(betas.int) == "list" & !(length(betas.int) %in% c(0, 1, M))))
    stop("betas.int should be list of length 0, 1 or M.")
  if (length(betas.tint) == 0  & !is.null(tint.var) & quiet == FALSE)
    message("All betas.tint will be set to 0 if using corrsys or corrsys2.")
  if (class(betas.tint) == "list" & length(betas.tint) == 1 & quiet == FALSE)
    message("All betas.tint will be set to betas.tint[[1]] if using corrsys or
            corrsys2.")
  if (class(betas.tint) != "list" |
      (class(betas.tint) == "list" & !(length(betas.tint) %in% c(0, 1, M))))
    stop("betas.tint should be list of length 0, 1 or M.")
  return(TRUE)
}

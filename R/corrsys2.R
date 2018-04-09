#' @title Generate Correlated Systems of Equations with Ordinal, Continuous, and/or Count Variables: Correlation Method 2
#'
#' @description This function generates a correlated system of \code{M} equations representing a \strong{system of repeated measures} at
#'     \code{M} time points.  The equations may contain 1) ordinal (\eqn{r \geq 2} categories), continuous (normal, non-normal, and mixture distributions),
#'     count (regular and zero-inflated, Poisson and Negative Binomial) independent variables \eqn{X};
#'     2) continuous error terms \eqn{E}; 3) a discrete time variable \eqn{Time}; and 4) random effects \eqn{U}.  The assumptions are that
#'     1) there are at least 2 equations, 2) the independent variables, random effect terms, and error terms are uncorrelated,
#'     3) each equation has an error term, 4) all error terms have a continuous non-mixture distribution or all have a continuous mixture
#'     distribution, 5) all random effects are continuous, and 6) growth is linear (with respect to time).  The random effects may be a
#'     random intercept, a random slope for time, or a random slope for any of the \eqn{X} variables.  Continuous variables are simulated
#'     using either Fleishman's third-order (\code{method = "Fleishman"}, \doi{10.1007/BF02293811})
#'     or Headrick's fifth-order (\code{method = "Polynomial"}, \doi{10.1016/S0167-9473(02)00072-5}) power method transformation (PMT).
#'     Simulation occurs at the component-level for continuous mixture distributions.  The target correlation matrix is specified in terms of
#'     correlations with components of continuous mixture variables.  These components are transformed into
#'     the desired mixture variables using random multinomial variables based on the mixing probabilities.  The \eqn{X} terms can be the same across equations (i.e., modeling sex or height) or may be time-varying covariates.  The
#'     equations may contain different numbers of \eqn{X} terms (i.e., a covariate could be missing for a given equation).
#'
#'     The outcomes \eqn{Y} are generated using a hierarchical linear models (HLM) approach, which allows the data to be structured in at least two levels.
#'     \strong{Level-1} is the repeated measure (time or condition) and other subject-level variables.  \strong{Level-1} is nested
#'     within \strong{Level-2}, which describes the average of the outcome (the intercept) and growth (slope for time) as a function
#'     of group-level variables.  The first level captures the within-subject variation, while the second level describes the between-subjects variability.
#'     Using a HLM provides a way to determine if: a) subjects differ at a specific time point with respect to the dependent
#'     variable, b) growth rates differ across conditions, or c) growth rates differ across subjects.  Random effects
#'     describe deviation at the subject-level from the average (fixed) effect described by the slope coefficients (betas).  See the
#'     \strong{The Hierarchical Linear Models Approach for a System of Correlated Equations with Multiple Variable Types} vignette for a
#'     description of the HLM model.  The user can specify subject-level \eqn{X} terms, and each subject-level \eqn{X} term is crossed with
#'     all group-level \eqn{X} terms.  The equations may also contain interactions between \eqn{X} variables.  Interactions specified in
#'     \code{int.var} between two group-level covariates are themselves considered group-level covariates and will be crossed with subject-level
#'     covariates.  Interactions between two subject-level covariates are considered subject-level covariates and will be crossed with
#'     group-level covariates.  Since \code{Time} is a subject-level variable, each group-level term is crossed with \code{Time}
#'     unless otherwise specified.
#'
#'     \strong{Random effects} may be added for the intercept, time slope, or effects of any of the covariates.  The type of random intercept and
#'     time slope (i.e., non-mixture or mixture) is specified in \code{rand.int} and \code{rand.tsl}.  This type may vary by equation.
#'     The random effects for independent variables are specified in \code{rand.var} and may also contain a combination of non-mixture and
#'     mixture continuous distributions.  If the parameter lists are of length \code{M + 1}, the random effects are the same variables across equations and the
#'     correlation for the effects \code{corr.u} is a matrix.  If the parameter lists are of length \code{2 * M}, the random effects are different variables
#'     across equations and the correlation for the effects \code{corr.u} is a list.
#'
#'     The independent variables, interactions, \code{Time} effect, random effects, and error terms are summed
#'     together to produce the outcomes \eqn{Y}. The beta coefficients may be the same or differ across equations.  The user specifies the betas for
#'     the independent variables in \code{betas}, for the interactions between two group-level or two subject-level covariates in \code{betas.int},
#'     for the group-subject level interactions in \code{betas.subj}, and for the \code{Time} interactions in \code{betas.tint}.
#'     Setting a coefficient to 0 will eliminate that term.  The user also provides the correlations 1) between \eqn{E} terms; 2) between \eqn{X}
#'     variables within each outcome \eqn{Y_p}, \code{p = 1, ..., M}, and between outcome pairs; and 3) between \eqn{U} variables within each
#'     outcome \eqn{Y_p}, \code{p = 1, ..., M}, and between outcome pairs.  The order of the independent variables in \code{corr.x} must be
#'     1st ordinal (same order as in \code{marginal}), 2nd continuous non-mixture (same order as in \code{skews}), 3rd components of continuous
#'     mixture (same order as in \code{mix_pis}), 4th regular Poisson, 5th zero-inflated Poisson (same order as in \code{lam}), 6th regular NB, and
#'     7th zero-inflated NB (same order as in \code{size}).  The order of the random effects in \code{corr.u} must be 1st random intercept, 2nd random
#'     time slope, 3rd continuous non-mixture random effects, and 4th components of continuous mixture random effects.
#'
#'     The variables are generated from multivariate normal variables with intermediate correlations calculated using
#'     \code{\link[SimCorrMix]{intercorr2}}, which employs \strong{correlation method 2}.  See \code{\link[SimCorrMix]{SimCorrMix}} for a
#'     description of the correlation method and the techniques used to generate each variable type.  The order of the variables returned is
#'     1st covariates \eqn{X} (as specified in \code{corr.x}), 2nd group-group or subject-subject interactions (ordered as in \code{int.var}),
#'     3rd subject-group interactions (1st by subject-level variable as specified in \code{subj.var}, 2nd by covariate as specified in \code{corr.x}),
#'     and 4th time interactions (either as specified in \code{corr.x} with group-level covariates or in \code{tint.var}).
#'
#'     This function contains no parameter checks in order to decrease simulation time.  That should be done first using
#'     \code{\link[SimRepeat]{checkpar}}.  Summaries of the system can be obtained using \code{\link[SimRepeat]{summary_sys}}.  The
#'     \strong{Correlated Systems of Statistical Equations with Multiple Variable Types} demonstrates examples.
#'
#' @section Reasons for Function Errors:
#'     1) The most likely cause for function errors is that the parameter inputs are mispecified.  Using \code{\link[SimRepeat]{checkpar}}
#'     prior to simulation can help decrease these errors.
#'
#'     2) Another reason for error is that no solutions to \code{\link[SimMultiCorrData]{fleish}} or
#'     \code{\link[SimMultiCorrData]{poly}} converged when using \code{\link[SimMultiCorrData]{find_constants}}.  If this happens,
#'     the simulation will stop.  It may help to first use \code{\link[SimMultiCorrData]{find_constants}} for each continuous variable to
#'     determine if a sixth cumulant correction value is needed.  If the standardized cumulants are obtained from \code{calc_theory}, the user
#'     may need to use rounded values as inputs (i.e. \code{skews = round(skews, 8)}).  For example, in order to ensure that skew is exactly 0 for symmetric distributions.
#'
#'     3) The kurtosis for a continuous variable may be outside the region of possible values.  There is an associated lower kurtosis boundary for
#'     associated with a given skew (for Fleishman's method) or skew and fifth and sixth cumulants (for Headrick's method).  Use
#'     \code{\link[SimMultiCorrData]{calc_lower_skurt}} to determine the boundary for a given set of cumulants.
#'
#' @param n the sample size (i.e. the length of each simulated variable; default = 10000)
#' @param M the number of dependent variables \eqn{Y} (outcomes); equivalently, the number of equations in the system
#' @param Time a vector of values to use for time; each subject receives the same time value; if \code{NULL}, \code{Time = 1:M}
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
#' @param pois_eps list of length \code{M}, p-th component a vector of length \code{lam[[p]]} containing cumulative probability truncation values
#'     used to calculate intermediate correlations involving Poisson variables; order is 1st regular Poisson and 2nd zero-inflated Poisson;
#'     if a single number, all Poisson variables given this value; if a vector of length \code{M}, all Poisson variables in equation p
#'     given \code{pois_eps[p]}; otherwise, missing values are set to 0.0001 and ordered 1st
#' @param size list of length \code{M}, p-th component a vector of size parameters for the Negative Binomial variables for outcome \eqn{Y_p}
#'     (see \code{stats::dnbinom}); order is 1st regular NB and 2nd zero-inflated NB; use \code{size[[p]] = NULL} if outcome \eqn{Y_p} has no Negative Binomial variables;
#'     \code{length(size[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}
#' @param prob list of length \code{M}, p-th component a vector of success probabilities for the Negative Binomial variables for outcome \eqn{Y_p}
#'     (see \code{stats::dnbinom}); order is 1st regular NB and 2nd zero-inflated NB; use \code{prob[[p]] = NULL} if outcome \eqn{Y_p} has no Negative Binomial variables;
#'     \code{length(prob[[p]])} can differ across outcomes; the order should be the same as in \code{corr.x}
#' @param mu list of length \code{M}, p-th component a vector of mean values for the Negative Binomial variables for outcome \eqn{Y_p}
#'     (see \code{stats::dnbinom}); order is 1st regular NB and 2nd zero-inflated NB; use \code{mu[[p]] = NULL} if outcome \eqn{Y_p} has no Negative Binomial variables;
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
#' @param corr.e correlation matrix for continuous non-mixture or components of mixture error terms
#' @param same.var either a vector or a matrix; if a vector, \code{same.var} includes column numbers of \code{corr.x[[1]][[1]]}
#'     corresponding to independent variables that should be identical across equations; these terms must have the same indices for all
#'     \code{p = 1, ..., M}; i.e., if the 1st ordinal variable represents sex which should be the same for each equation, then
#'     \code{same.var[1] = 1} since ordinal variables are 1st in \code{corr.x[[1]][[1]]} and sex is the 1st ordinal variable, and
#'     the 1st term for all other outcomes must also be sex;
#'     if a matrix, columns 1 and 2 are outcome p and column index in \code{corr.x[[p]][[p]]} for 1st instance of variable,
#'     columns 3 and 4 are outcome q and column index in \code{corr.x[[q]][[q]]} for subsequent instances of variable; i.e., if
#'     1st term for all outcomes is sex and \code{M = 3}, then \code{same.var = matrix(c(1,} \code{1, 2, 1, 1, 1, 3, 1), 2, 4, byrow = TRUE)}; the
#'     independent variable index corresponds to ordinal, continuous non-mixture, component of continuous mixture, Poisson, or
#'     NB variable
#' @param subj.var matrix where 1st column is outcome index (\code{p = 1, ..., M}), 2nd column is independent variable index
#'     corresponding to covariate which is a a subject-level term (not including time), including time-varying covariates;
#'     the independent variable index corresponds to ordinal, continuous non-mixture, continuous mixture (not mixture component), Poisson, or
#'     NB variable; assumes all other variables are group-level terms; these subject-level terms are used to form interactions with the group
#'     level terms
#' @param int.var matrix where 1st column is outcome index (\code{p = 1, ..., M}), 2nd and 3rd columns are indices
#'     corresponding to two group-level or two subject-level independent variables to form interactions between; this includes all
#'     interactions that are not accounted for by a subject-group level interaction (as indicated by \code{subj.var}) or by a
#'     time-covariate interaction (as indicated by \code{tint.var}); ex: 1, 2, 3 indicates that for outcome 1, the 2nd and 3rd
#'     independent variables form an interaction term; the independent variable index corresponds to ordinal, continuous non-mixture,
#'     continuous mixture (not mixture component), Poisson, or NB variable
#' @param tint.var matrix where 1st column is outcome index (\code{p = 1, ..., M}), 2nd column is index of
#'     independent variable to form interaction with time; if \code{tint.var = NULL} or no \eqn{X_{(pj)}} are indicated for outcome \eqn{Y_p},
#'     all group-level variables (variables not indicated as subject-level variables in \code{subj.var}) will be crossed with time, else includes only
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
#'     subject-level terms given in \code{subj.var}; order is 1st by subject-level covariate as given in \code{subj.var} and 2nd by group-level covariate
#'     as given in \code{corr.x} or an interaction between group-level terms; if all outcomes have the same betas,
#'     use list of length 1; if \eqn{Y_p} only has group-level terms, set \code{betas.subj[[p]] = NULL}; since subject-subject interactions are
#'     treated as subject-level variables, these will also be crossed with all group-level variables and require coefficients;
#'     if there are continuous mixture variables, beta is for mixture variable (not for components)
#' @param betas.int list of length \code{M}, p-th component a vector of coefficients for interaction terms indicated in \code{int.var};
#'     order is the same order as given in \code{int.var}; if all outcomes have the same betas, use list of length 1;
#'     if \eqn{Y_p} has none, set \code{betas.int[[p]] = NULL};
#'     if there are continuous mixture variables, beta is for mixture variable (not for components)
#' @param betas.t vector of length \code{M} of coefficients for time terms, if \code{NULL} all set equal to 1;
#'     if length 1, all intercepts set to \code{betas.t}
#' @param betas.tint list of length \code{M}, p-th component a vector of coefficients for all interactions with time; this includes interactions
#'     with group-level covariates or terms indicated in \code{tint.var};
#'     order is the same order as given in \code{corr.x} or \code{tint.var}; if all outcomes have the same betas, use list of length 1;
#'     if \eqn{Y_p} has none, set \code{betas.tint[[p]] = NULL}; since group-group interactions are treated as group-level variables, these will
#'     also be crossed with time (unless otherwise specified for that outcome in \code{tint.var}) and require coefficients;
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
#' @param seed the seed value for random number generation (default = 1234)
#' @param use.nearPD TRUE to convert the overall intermediate correlation matrix formed by the \eqn{X} (for all outcomes and independent
#'     variables), \eqn{E}, or the random effects to the nearest positive definite matrix with \code{Matrix::nearPD} if necessary;
#'     if FALSE and \code{adjgrad = FALSE} the negative eigenvalues are replaced with \code{eigmin} if necessary
#' @param eigmin minimum replacement eigenvalue if overall intermediate correlation matrix is not positive-definite (default = 0)
#' @param adjgrad TRUE to use \code{adj_grad} to convert overall intermediate correlation matrix to a positive-definite matrix and next
#'     5 inputs can be used
#' @param B1 the initial matrix for algorithm; if NULL, uses a scaled initial matrix with diagonal elements \code{sqrt(nrow(Sigma))/2}
#' @param tau parameter used to calculate theta (default = 0.5)
#' @param tol maximum error for Frobenius norm distance between new matrix and original matrix (default = 0.1)
#' @param steps maximum number of steps for k (default = 100)
#' @param msteps maximum number of steps for m (default = 10)
#' @param errorloop if TRUE, uses \code{\link[SimCorrMix]{corr_error}} to attempt to correct the correlation of the independent
#'     variables within and across outcomes to be within \code{epsilon} of the target correlations \code{corr.x} until the number of iterations
#'     reaches \code{maxit} (default = FALSE)
#' @param epsilon the maximum acceptable error between the final and target correlation matrices (default = 0.001)
#'     in the calculation of ordinal intermediate correlations with \code{\link[SimCorrMix]{ord_norm}} or in the error loop
#' @param maxit the maximum number of iterations to use (default = 1000) in the calculation of ordinal
#'     intermediate correlations with \code{\link[SimCorrMix]{ord_norm}} or in the error loop
#' @param quiet if FALSE prints messages, if TRUE suppresses messages
#' @import SimMultiCorrData
#' @import SimCorrMix
#' @importFrom stats cor dbeta dbinom dchisq density dexp df dgamma dlnorm dlogis dmultinom dnbinom dnorm dpois dt dunif dweibull ecdf
#'     median pbeta pbinom pchisq pexp pf pgamma plnorm plogis pnbinom pnorm ppois pt punif pweibull qbeta qbinom qchisq qexp qf qgamma
#'     qlnorm qlogis qnbinom qnorm qpois qt quantile qunif qweibull rbeta rbinom rchisq rexp rf rgamma rlnorm rlogis rmultinom rnbinom
#'     rnorm rpois rt runif rweibull sd uniroot var
#' @import utils
#' @importFrom Matrix nearPD
#' @importFrom VGAM qzipois qzinegbin dzipois dzinegbin rzipois rzinegbin
#' @export
#' @keywords simulation continuous mixture ordinal Poisson NegativeBinomial Fleishman Headrick method2
#' @seealso \code{\link[SimMultiCorrData]{find_constants}}, \code{\link[SimCorrMix]{intercorr2}},
#'     \code{\link[SimRepeat]{checkpar}}, \code{\link[SimRepeat]{summary_sys}}
#'
#' @return A list with the following components:
#'
#' @return \code{Y} matrix with \code{n} rows and \code{M} columns of outcomes
#'
#' @return \code{X} list of length \code{M} containing \eqn{X_{ord(pj)}, X_{cont(pj)}, X_{comp(pj)}, X_{pois(pj)}, X_{nb(pj)}}
#'
#' @return \code{X_all} list of length \code{M} containing \eqn{X_{ord(pj)}, X_{cont(pj)}, X_{mix(pj)}, X_{pois(pj)}, X_{nb(pj)},} \eqn{X} interactions as indicated by
#'     \code{int.var}, subject-group level term interactions as indicated by \code{subj.var}, \eqn{Time_p}, and \eqn{Time} interactions as indicated by
#'     \code{tint.var}; order is 1st covariates \eqn{X} (as specified in \code{corr.x}), 2nd group-group or subject-subject interactions (ordered as in \code{int.var}),
#'     3rd subject-group interactions (1st by subject-level variable as specified in \code{subj.var}, 2nd by covariate as specified in
#'     \code{corr.x}), and 4th time interactions (either as specified in \code{corr.x} with group-level covariates or in \code{tint.var})
#'
#' @return \code{E} matrix with \code{n} rows containing continuous non-mixture or components of continuous mixture error terms
#'
#' @return \code{E_mix} matrix with \code{n} rows containing continuous mixture error terms
#'
#' @return \code{Sigma_X0} matrix of intermediate correlations calculated by \code{intercorr2}
#'
#' @return \code{Sigma_X} matrix of intermediate correlations after \code{nearPD} or \code{adj_grad} function has been used;
#'     applied to generate the normal variables transformed to get the desired distributions
#'
#' @return \code{Error_Time} the time in minutes required to use the error loop
#'
#' @return \code{Time} the total simulation time in minutes
#'
#' @return \code{niter} a matrix of the number of iterations used in the error loop
#'
#' @return If \bold{continuous variables} are produced:
#'     \code{constants} a list of maximum length \code{2 * M}, the 1st \code{M} components are data.frames of the constants for the
#'         \eqn{X_{cont(pj)}}, \eqn{X_comp(pj)} and \eqn{E_p}, the 2nd \code{M} components are for random effects (if present),
#'
#'     \code{SixCorr} a list of maximum length \code{2 * M}, the 1st \code{M} components are lists of sixth cumulant correction
#'         values used to obtain valid \emph{pdf}'s for the \eqn{X_{cont(pj)}}, \eqn{X_comp(pj)}, and \eqn{E_p}, the 2nd \code{M} components are for random effects (if present),
#'
#'     \code{valid.pdf} a list of maximum length \code{2 * M} of vectors where the i-th element is "TRUE" if the constants for the i-th
#'         continuous variable generate a valid pdf, else "FALSE"; the 1st \code{M} components are for the \eqn{X_{cont(pj)}},
#'         \eqn{X_comp(pj)}, and \eqn{E_p}, the 2nd \code{M} components are for random effects (if present)
#'
#' @return If \bold{random effects} are produced:
#'     \code{U} a list of length \code{M} containing matrices of continuous non-mixture and components of mixture random effects,
#'
#'     \code{U_all} a list of length \code{M} containing matrices of continuous non-mixture and mixture random effects,
#'
#'     \code{V} a list of length \code{M} containing matrices of design matrices for random effects,
#'
#'     \code{rmeans2} and \code{rvars2} the means and variances of the non-mixture and components reordered in accordance with the
#'     random intercept and time slope types (input for \code{summary_sys})
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
#' Sys2 <- corrsys2(10000, M, Time = 1:M, "Polynomial", "non_mix", means, vars,
#'   skews, skurts, fifths, sixths, Six, marginal = marginal, support = support,
#'   corr.x = corr.x, corr.e = corr.e, betas = betas, betas.t = betas.t,
#'   betas.tint = betas.tint, quiet = TRUE)
#'
#' \dontrun{
#' seed <- 276
#' n <- 10000
#' M <- 3
#' Time <- 1:M
#'
#' # Error terms have a beta(4, 1.5) distribution with an AR(1, p = 0.4)
#' correlation structure
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
#' 1 continuous mixture of Normal(-2, 1) and Normal(2, 1) for each Y
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
#' # Y2 and Y3 have zero-inflated Poisson variables
#' p_zip <- list(NULL, 0.05, 0.1)
#'
#' # 1 NB variable for each Y
#' size <- list(10, 15, 20)
#' prob <- list(0.3, 0.4, 0.5)
#' # either prob or mu is required (not both)
#' mu <- mapply(function(x, y) x * (1 - y)/y, size, prob, SIMPLIFY = FALSE)
#' # Y2 and Y3 have zero-inflated NB variables
#' p_zinb <- list(NULL, 0.05, 0.1)
#' pois_eps <- nb_eps <- list()
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
#'   mix_sixths, mix_Six, marginal, support, lam, p_zip, pois_eps, size, prob,
#'   mu, p_zinb, nb_eps, corr.x, corr.yx = list(), corr.e, same.var, subj.var,
#'   int.var, tint.var, betas.0, betas, betas.subj, betas.int, betas.t,
#'   betas.tint)
#'
#' # Simulated system using correlation method 2
#' N <- corrsys2(n, M, Time, method, error_type, means, vars, skews, skurts,
#'   fifths, sixths, Six, mix_pis, mix_mus, mix_sigmas, mix_skews, mix_skurts,
#'   mix_fifths, mix_sixths, mix_Six, marginal, support, lam, p_zip, pois_eps,
#'   size, prob, mu, p_zinb, nb_eps, corr.x, corr.e, same.var, subj.var,
#'   int.var, tint.var, betas.0, betas, betas.subj, betas.int, betas.t,
#'   betas.tint, seed = seed, use.nearPD = FALSE)
#'
#' # Summarize the results
#' S <- summary_sys(N$Y, N$E, E_mix = NULL, N$X, N$X_all, M, method, means,
#'   vars, skews, skurts, fifths, sixths, mix_pis, mix_mus, mix_sigmas,
#'   mix_skews, mix_skurts, mix_fifths, mix_sixths, marginal, support, lam,
#'   p_zip, size, prob, mu, p_zinb, corr.x, corr.e)
#' }
#'
#' @references
#' Barbiero A & Ferrari PA (2015). Simulation of correlated Poisson variables. Applied Stochastic Models in
#'     Business and Industry, 31:669-80. \doi{10.1002/asmb.2072}.
#'
#' Barbiero A & Ferrari PA (2015). GenOrd: Simulation of Discrete Random Variables with Given
#'     Correlation Matrix and Marginal Distributions. R package version 1.4.0. \cr \url{https://CRAN.R-project.org/package=GenOrd}
#'
#' Davenport JW, Bezder JC, & Hathaway RJ (1988). Parameter Estimation for Finite Mixture Distributions.
#'     Computers & Mathematics with Applications, 15(10):819-28.
#'
#' Demirtas H (2006). A method for multivariate ordinal data generation given marginal distributions and correlations. Journal of Statistical
#'     Computation and Simulation, 76(11):1017-1025. \cr \doi{10.1080/10629360600569246}.
#'
#' Demirtas H (2014). Joint Generation of Binary and Nonnormal Continuous Data. Biometrics & Biostatistics, S12.
#'
#' Demirtas H, Hedeker D, & Mermelstein RJ (2012). Simulation of massive public health data by power polynomials.
#'     Statistics in Medicine, 31(27):3337-3346. \doi{10.1002/sim.5362}.
#'
#' Everitt BS (1996). An Introduction to Finite Mixture Distributions. Statistical Methods in Medical Research, 5(2):107-127. \doi{10.1177/096228029600500202}.
#'
#' Ferrari PA & Barbiero A (2012). Simulating ordinal data. Multivariate Behavioral Research, 47(4): 566-589.
#'     \doi{10.1080/00273171.2012.692630}.
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
#' Ismail N & Zamani H (2013). Estimation of Claim Count Data Using Negative Binomial, Generalized Poisson, Zero-Inflated Negative Binomial and
#'     Zero-Inflated Generalized Poisson Regression Models. Casualty Actuarial Society E-Forum 41(20):1-28.
#'
#' Kincaid C (2005). Guidelines for Selecting the Covariance Structure in Mixed Model Analysis. Computational Statistics and
#'     Data Analysis, 198(30):1-8.
#'
#' Lambert D (1992). Zero-Inflated Poisson Regression, with an Application to Defects in Manufacturing. Technometrics 34(1):1-14.
#'
#' Lininger M, Spybrook J, & Cheatham CC (2015). Hierarchical Linear Model: Thinking Outside the Traditional Repeated-Measures
#'     Analysis-of-Variance Box. Journal of Athletic Training, 50(4):438-441. \doi{10.4085/1062-6050-49.5.09}.
#'
#' McCulloch CE, Searle SR, Neuhaus JM (2008). \emph{Generalized, Linear, and Mixed Models} (2nd ed.). Wiley Series in Probability and
#'     Statistics. Hoboken, New Jersey: John Wiley & Sons, Inc.
#'
#' Olsson U, Drasgow F, & Dorans NJ (1982). The Polyserial Correlation Coefficient. Psychometrika, 47(3):337-47.
#'     \doi{10.1007/BF02294164}.
#'
#' Pearson RK (2011). Exploring Data in Engineering, the Sciences, and Medicine. In. New York: Oxford University Press.
#'
#' Schork NJ, Allison DB, & Thiel B (1996). Mixture Distributions in Human Genetics Research. Statistical Methods in Medical Research,
#'     5:155-178. \doi{10.1177/096228029600500204}.
#'
#' Vale CD & Maurelli VA (1983). Simulating Multivariate Nonnormal Distributions. Psychometrika, 48:465-471. \doi{10.1007/BF02293687}.
#'
#' Van Der Leeden R (1998). Multilevel Analysis of Repeated Measures Data. Quality & Quantity, 32(1):15-29.
#'
#' Yee TW (2017). VGAM: Vector Generalized Linear and Additive Models. \cr \url{https://CRAN.R-project.org/package=VGAM}.
#'
#' Zhang X, Mallick H, & Yi N (2016). Zero-Inflated Negative Binomial Regression for Differential Abundance Testing in Microbiome
#'     Studies. Journal of Bioinformatics and Genomics 2(2):1-9. \doi{10.18454/jbg.2016.2.2.1}.
#'
corrsys2 <- function(n = 10000, M = NULL, Time = NULL,
                     method = c("Fleishman", "Polynomial"),
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
                     corr.e = NULL, same.var = NULL, subj.var = NULL,
                     int.var = NULL, tint.var = NULL, betas.0 = NULL,
                     betas = list(), betas.subj = list(), betas.int = list(),
                     betas.t = NULL, betas.tint = list(),
                     rand.int = c("none", "non_mix", "mix"),
                     rand.tsl = c("none", "non_mix", "mix"), rand.var = NULL,
                     corr.u = list(), seed = 1234, use.nearPD = TRUE,
                     eigmin = 0, adjgrad = FALSE, B1 = NULL, tau = 0.5,
                     tol = 0.1, steps = 100, msteps = 10, errorloop = FALSE,
                     epsilon = 0.001, maxit = 1000, quiet = FALSE) {
  start.time <- Sys.time()
  if (length(error_type) != 1)
    stop("Please choose one type of distribution for all of the error terms:
         mix if all errors have continuous mixture distributions,
         non_mix if all errors have continuous non-mixture distributions.")
  K.cat <- rep(0, M)
  K.pois <- rep(0, M)
  K.nb <- rep(0, M)
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
  k.comp <- c(0, cumsum(K.comp))
  K.cont <- lengths(vars) - K.mix
  k.cont <- c(0, cumsum(K.cont))
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
    pois_eps0 <- pois_eps
    pois_eps <- list()
    if (class(pois_eps0) == "numeric" & length(pois_eps0) == 1) {
      for (i in 1:M) {
        if (K.pois[i] == 0) pois_eps <- append(pois_eps, list(NULL)) else
          pois_eps[[i]] <- rep(pois_eps0, K.pois[i])
      }
    } else if (class(pois_eps0) == "numeric" & length(pois_eps0) == M) {
      for (i in 1:M) {
        if (K.pois[i] == 0) pois_eps <- append(pois_eps, list(NULL)) else
          pois_eps[[i]] <- rep(pois_eps0[i], K.pois[i])
      }
    } else if (class(pois_eps0) == "list" & length(pois_eps0) == M) {
      for (i in 1:M) {
        if (K.pois[i] == 0) {
          pois_eps <- append(pois_eps, list(NULL))
          next
        }
        if (length(pois_eps0[[i]]) == K.pois[i])
          pois_eps[[i]] <- pois_eps0[[i]] else
            pois_eps[[i]] <- c(rep(0.0001, K.pois[i] - length(pois_eps0[[i]])),
                               pois_eps0[[i]])
      }
    } else {
      pois_eps <- lapply(K.pois,
                         function(x) if (x == 0) NULL else rep(0.0001, x))
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
    nb_eps0 <- nb_eps
    nb_eps <- list()
    if (class(nb_eps0) == "numeric" & length(nb_eps0) == 1) {
      for (i in 1:M) {
        if (K.nb[i] == 0) nb_eps <- append(nb_eps, list(NULL)) else
          nb_eps[[i]] <- rep(nb_eps0, K.nb[i])
      }
    } else if (class(nb_eps0) == "numeric" & length(nb_eps0) == M) {
      for (i in 1:M) {
        if (K.nb[i] == 0) nb_eps <- append(nb_eps, list(NULL)) else
          nb_eps[[i]] <- rep(nb_eps0[i], K.nb[i])
      }
    } else if (class(nb_eps0) == "list" & length(nb_eps0) == M) {
      for (i in 1:M) {
        if (K.nb[i] == 0) {
          nb_eps <- append(nb_eps, list(NULL))
          next
        }
        if (length(nb_eps0[[i]]) == K.nb[i])
          nb_eps[[i]] <- nb_eps0[[i]] else
            nb_eps[[i]] <- c(rep(0.0001, K.nb[i] - length(nb_eps0[[i]])),
                               nb_eps0[[i]])
      }
    } else {
      nb_eps <- lapply(K.nb,
                         function(x) if (x == 0) NULL else rep(0.0001, x))
    }
  }
  K.r <- rep(0, M)
  if (class(corr.u) == "list" & length(corr.u) > 0)
    K.r <- sapply(mapply('[', corr.u, lengths(corr.u), SIMPLIFY = FALSE),
                  function(x) if (is.null(x)) 0 else nrow(x[[1]]))
  if (class(corr.u) == "matrix") K.r <- rep(ncol(corr.u), M)
  if (is.null(betas.0)) betas.0 <- rep(0, M)
  if (length(betas.0) < M) betas.0 <- rep(betas.0[1], M)
  if (is.null(betas.t)) betas.t <- rep(1, M)
  if (length(betas.t) < M) betas.t <- rep(betas.t[1], M)
  if (is.null(Time)) Time <- 1:M
  csame.dist <- NULL
  csame.dist2 <- NULL
  skews2 <- NULL
  skurts2 <- NULL
  fifths2 <- NULL
  sixths2 <- NULL
  if (length(skews) > 0) {
    skews2 <- unlist(skews)
    skurts2 <- unlist(skurts)
    if (method == "Polynomial") {
      fifths2 <- unlist(fifths)
      sixths2 <- unlist(sixths)
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
  for (i in 1:M0) {
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
  for (i in 1:M0) {
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
  rmix_pis <- list()
  rmeans <- list()
  if (length(means) == (M + 1)) {
    rmeans <- means[M + 1]
    means <- means[1:M]
    rvars <- vars[M + 1]
    vars <- vars[1:M]
    rconstants <- constants[M + 1]
    constants <- constants[1:M]
    names(constants) <- paste("XE", 1:M, sep = "")
    names(rconstants) <- "U"
  }
  if (length(mix_pis) == (M + 1)) {
    rmix_pis <- mix_pis[M + 1]
    mix_pis <- mix_pis[1:M]
    rmix_mus <- mix_mus[M + 1]
    mix_mus <- mix_mus[1:M]
    rmix_sigmas <- mix_sigmas[M + 1]
    mix_sigmas <- mix_sigmas[1:M]
  }
  if (length(means) == (2 * M)) {
    rmeans <- means[(M + 1):(2 * M)]
    means <- means[1:M]
    rvars <- vars[(M + 1):(2 * M)]
    vars <- vars[1:M]
    rconstants <- constants[(M + 1):(2 * M)]
    constants <- constants[1:M]
    names(constants) <- paste("XE", 1:M, sep = "")
    names(rconstants) <- paste("U", 1:M, sep = "")
  }
  if (length(mix_pis) == (2 * M)) {
    rmix_pis <- mix_pis[(M + 1):(2 * M)]
    mix_pis <- mix_pis[1:M]
    rmix_mus <- mix_mus[(M + 1):(2 * M)]
    mix_mus <- mix_mus[1:M]
    rmix_sigmas <- mix_sigmas[(M + 1):(2 * M)]
    mix_sigmas <- mix_sigmas[1:M]
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
  if (error_type == "non_mix" & max(K.mix2) > 0) {
    for (i in 1:M) {
      constants[[i]] <- rbind(constants[[i]][-K.cont[i], ],
                              constants[[i]][K.cont[i], ])
    }
  }
  same.ord <- NULL
  same.cont <- NULL
  same.pois <- NULL
  same.nb <- NULL
  if (!is.null(same.var)) {
    if (class(same.var) == "numeric") {
      temp.var <- NULL
      for (j in 1:length(same.var)) {
        temp.var <- rbind(temp.var, cbind(rep(1, M - 1),
          rep(same.var[j], M - 1), 2:M, rep(same.var[j], M - 1)))
      }
      same.var <- temp.var
    }
    for (i in 1:nrow(same.var)) {
      if (same.var[i, 2] <= K.cat[same.var[i, 1]])
        same.ord <- rbind(same.ord, same.var[i, ])
      if (same.var[i, 2] <= (K.x[same.var[i, 1]] - K.pois[same.var[i, 1]] -
                             K.nb[same.var[i, 1]]) &
          same.var[i, 2] > K.cat[same.var[i, 1]])
        same.cont <- rbind(same.cont, same.var[i, ])
      if (same.var[i, 2] <= (K.x[same.var[i, 1]] - K.nb[same.var[i, 1]]) &
          same.var[i, 2] > (K.x[same.var[i, 1]] - K.pois[same.var[i, 1]] -
                            K.nb[same.var[i, 1]]))
        same.pois <- rbind(same.pois, same.var[i, ])
      if (same.var[i, 2] > (K.x[same.var[i, 1]] - K.nb[same.var[i, 1]]))
        same.nb <- rbind(same.nb, same.var[i, ])
    }
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
        if (K.cat[m] > 0 & length(support[[m]]) == 0) {
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
    Z <- matrix(rnorm((ncol(Sigma_E) + sum(K.x) - nrow(same.var) +
                         sum(K.r)) * n), n)
  } else Z <- matrix(rnorm((ncol(Sigma_E) + sum(K.x) + sum(K.r)) * n), n)
  Z <- scale(Z, TRUE, FALSE)
  Z <- Z %*% svd(Z, nu = 0)$v
  Z <- scale(Z, FALSE, TRUE)
  if (min(eigen(Sigma_E, symmetric = TRUE)$values) < 0) {
    if (quiet == FALSE)
      message("Intermediate E correlation matrix is not positive definite.")
    if (use.nearPD == TRUE) {
      Sigma_E <- as.matrix(nearPD(Sigma_E, corr = T, keepDiag = T)$mat)
    } else if (adjgrad == TRUE) {
      sadj <- adj_grad(Sigma_E, B1, tau, tol, steps, msteps)
      Sigma_E <- sadj$Sigma2
    }
  }
  eig <- eigen(Sigma_E, symmetric = TRUE)
  sqrteigval <- diag(sqrt(pmax(eig$values, 0)), nrow(Sigma_E))
  eigvec <- eig$vectors
  fry <- eigvec %*% sqrteigval
  E <- fry %*% t(Z[, 1:ncol(Sigma_E)])
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
  if (length(corr.x) > 0) {
    constants0 <- constants
    marginal0 <- marginal
    support0 <- support
    lam0 <- lam
    p_zip0 <- p_zip
    size0 <- size
    mu0 <- mu
    p_zinb0 <- p_zinb
    constants2 <- NULL
    marginal2 <- list()
    support2 <- list()
    lam2 <- NULL
    pois_eps2 <- NULL
    p_zip2 <- NULL
    size2 <- NULL
    mu2 <- NULL
    nb_eps2 <- NULL
    p_zinb2 <- NULL
    for (p in 1:M) {
      if (K.x[p] == 0) next
      rownames(corr.x[[p]][[p]]) <- 1:K.x[p]
      if (K.cat[p] > 0) {
        rownames(corr.x[[p]][[p]])[1:K.cat[p]] <-
          paste("cat", p, "_", 1:K.cat[p], sep = "")
        if (!is.null(same.ord)) {
          if (p %in% same.ord[, 3]) {
            marginal[[p]] <- marginal[[p]][-same.ord[same.ord[, 3] == p, 4]]
            support[[p]] <- support[[p]][-same.ord[same.ord[, 3] == p, 4]]
          }
        }
        marginal2 <- append(marginal2, marginal[[p]])
        support2 <- append(support2, support[[p]])
      }
      if ((K.x[p] - K.cat[p] - K.pois[p] - K.nb[p]) > 0) {
        rownames(corr.x[[p]][[p]])[(K.cat[p] + 1):(K.x[p] -
                                                     K.pois[p] - K.nb[p])] <-
          paste("cont", p, "_", 1:(K.x[p] - K.cat[p] - K.pois[p] - K.nb[p]),
                sep = "")
        if (!is.null(same.cont)) {
          if (p %in% same.cont[, 3])
            constants[[p]] <-
              constants[[p]][-(same.cont[same.cont[, 3] == p, 4] - K.cat[p]),
                             , drop = FALSE]
        }
        if (error_type == "non_mix") {
          constants2 <- rbind(constants2,
            constants[[p]][-nrow(constants[[p]]), , drop = FALSE])
        } else {
          constants2 <- rbind(constants2,
            constants[[p]][-c((nrow(constants[[p]]) - K.error[p] +
                                 1):nrow(constants[[p]])), , drop = FALSE])
        }
      }
      if (K.pois[p] > 0) {
        rownames(corr.x[[p]][[p]])[(K.x[p] - K.pois[p] - K.nb[p] +
                                      1):(K.x[p] - K.nb[p])] <-
          paste("pois", p, "_", 1:K.pois[p], sep = "")
        if (!is.null(same.pois)) {
          if (p %in% same.pois[, 3]) {
            lam[[p]] <- lam[[p]][-(same.pois[same.pois[, 3] == p, 4] -
                                     K.x[p] + K.pois[p] + K.nb[p])]
            p_zip[[p]] <- p_zip[[p]][-(same.pois[same.pois[, 3] == p, 4] -
                                     K.x[p] + K.pois[p] + K.nb[p])]
            pois_eps[[p]] <- pois_eps[[p]][-(same.pois[same.pois[, 3] == p, 4] -
                                         K.x[p] + K.pois[p] + K.nb[p])]
          }
        }
        lam2 <- c(lam2, lam[[p]])
        p_zip2 <- c(p_zip2, p_zip[[p]])
        pois_eps2 <- c(pois_eps2, pois_eps[[p]])
      }
      if (K.nb[p] > 0) {
        rownames(corr.x[[p]][[p]])[(K.x[p] - K.nb[p] + 1):K.x[p]] <-
          paste("xnb", p, "_", 1:K.nb[p], sep = "")
        if (!is.null(same.nb)) {
          if (p %in% same.nb[, 3]) {
            size[[p]] <- size[[p]][-(same.nb[same.nb[, 3] == p, 4] - K.x[p] +
                                       K.nb[p])]
            mu[[p]] <- mu[[p]][-(same.nb[same.nb[, 3] == p, 4] - K.x[p] +
                                   K.nb[p])]
            p_zinb[[p]] <- p_zinb[[p]][-(same.nb[same.nb[, 3] == p, 4] -
                                           K.x[p] + K.nb[p])]
            nb_eps[[p]] <- nb_eps[[p]][-(same.nb[same.nb[, 3] == p, 4] -
                                           K.x[p] + K.nb[p])]
          }
        }
        size2 <- c(size2, size[[p]])
        p_zinb2 <- c(p_zinb2, p_zinb[[p]])
        mu2 <- c(mu2, mu[[p]])
        nb_eps2 <- c(nb_eps2, nb_eps[[p]])
      }
      if (!is.null(same.var)) {
        same.var0 <- same.var[same.var[, 3] == p, , drop = FALSE]
        if (p != 1 & nrow(same.var0) != 0) {
          q <- unique(same.var0[, 1])
          for (i in 1:length(q)) {
            rownames(corr.x[[p]][[p]])[same.var0[same.var0[, 1] == q[i], 4]] <-
              rownames(corr.x[[q[i]]][[q[i]]])[same.var0[same.var0[, 1] ==
                                                           q[i], 2]]
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
    for (p in 1:(M - 1)) {
      if (K.x[p] == 0) next
      if (!is.null(same.var)) {
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
                                                         -same.var0[, 4],
                                                         drop = FALSE]
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
    Corr_X <- Corr_X[order(rownames(Corr_X)), order(colnames(Corr_X)),
                     drop = FALSE]
    Sigma_X <- intercorr2(k_cat = length(grep("cat", colnames(Corr_X))),
      k_cont = length(grep("cont", colnames(Corr_X))),
      k_pois = length(grep("pois", colnames(Corr_X))),
      k_nb = length(grep("xnb", colnames(Corr_X))), method = method,
      constants = constants2, marginal = marginal2, support = support2,
      lam = lam2, p_zip = p_zip2, size = size2, mu = mu2, p_zinb = p_zinb2,
      pois_eps = pois_eps2, nb_eps = nb_eps2, rho = Corr_X, epsilon = epsilon,
      maxit = maxit, quiet = quiet)
    colnames(Sigma_X) <- colnames(Corr_X)
    rownames(Sigma_X) <- rownames(Corr_X)
    Sigma_X0 <- Sigma_X[names1, names1, drop = FALSE]
    if (min(eigen(Sigma_X, symmetric = TRUE)$values) < 0) {
      if (quiet == FALSE)
        message("Intermediate correlation matrix is not positive definite.")
      if (use.nearPD == TRUE) {
        Sigma_X <- as.matrix(nearPD(Sigma_X, corr = T, keepDiag = T)$mat)
      } else if (adjgrad == TRUE) {
        sadj <- adj_grad(Sigma_X, B1, tau, tol, steps, msteps)
        Sigma_X <- sadj$Sigma2
      }
    }
    eig <- eigen(Sigma_X, symmetric = TRUE)
    sqrteigval <- diag(sqrt(pmax(eig$values, 0)), nrow(Sigma_X))
    eigvec <- eig$vectors
    fry <- eigvec %*% sqrteigval
    X <- fry %*% t(Z[, (ncol(Sigma_E) + 1):(ncol(Sigma_E) + ncol(Sigma_X)),
                     drop = FALSE])
    X <- t(X)
    colnames(X) <- colnames(Corr_X)
    colnames(Sigma_X) <- colnames(Corr_X)
    rownames(Sigma_X) <- rownames(Corr_X)
    X <- X[, names1, drop = FALSE]
    X0 <- list()
    for (i in 1:M) {
      if (K.x[i] != 0) {
        if (i == 1) {
          X0[[i]] <- X[, 1:K.x[i], drop = FALSE]
        } else {
          if (!is.null(same.var)) {
            if (K.x[i] == nrow(same.var[same.var[, 3] == i, , drop = FALSE])) {
              X0 <- append(X0, list(NULL))
            } else {
              X0[[i]] <- X[, (cumsum(K.x)[i - 1] -
                nrow(same.var[same.var[, 3] <= (i - 1), , drop = FALSE]) +
                1):(cumsum(K.x)[i] - nrow(same.var[same.var[, 3] <= i, ,
                drop = FALSE])), drop = FALSE]
            }
          } else {
            X0[[i]] <- X[, (cumsum(K.x)[i - 1] + 1):cumsum(K.x)[i],
                         drop = FALSE]
          }
        }
      } else X0 <- append(X0, list(NULL))
    }
    X_cat <- list()
    X_cont <- list()
    X_pois <- list()
    X_nb <- list()
    Y_cat <- list()
    Y_cont <- list()
    Y_pois <- list()
    Y_nb <- list()
    means2 <- list()
    vars2 <- list()
    temp_means <- NULL
    temp_vars <- NULL
    for (m in 1:M) {
      X_cat <- append(X_cat, list(NULL))
      Y_cat <- append(Y_cat, list(NULL))
      X_cont <- append(X_cont, list(NULL))
      Y_cont <- append(Y_cont, list(NULL))
      means2 <- append(means2, list(NULL))
      vars2 <- append(vars2, list(NULL))
      X_pois <- append(X_pois, list(NULL))
      Y_pois <- append(Y_pois, list(NULL))
      X_nb <- append(X_nb, list(NULL))
      Y_nb <- append(Y_nb, list(NULL))
      if (K.cat[m] != 0 & length(grep("cat", colnames(X0[[m]]))) != 0) {
        X_cat[[m]] <- X0[[m]][, grep("cat", colnames(X0[[m]])), drop = FALSE]
        Y_cat[[m]] <- matrix(1, nrow = n, ncol = ncol(X_cat[[m]]))
        for (i in 1:ncol(X_cat[[m]])) {
          Y_cat[[m]][, i] <- as.integer(cut(X_cat[[m]][, i],
            breaks = c(min(X_cat[[m]][, i]) - 1, qnorm(marginal[[m]][[i]]),
                       max(X_cat[[m]][, i])  +  1)))
          Y_cat[[m]][, i] <- support[[m]][[i]][Y_cat[[m]][, i]]
        }
        colnames(Y_cat[[m]]) <- colnames(X_cat[[m]])
      }
      if ((K.x[m] - K.cat[m] - K.pois[m] - K.nb[m]) != 0 &
          length(grep("cont", colnames(X0[[m]]))) != 0) {
        X_cont[[m]] <- X0[[m]][, grep("cont", colnames(X0[[m]])), drop = FALSE]
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
            means2[[m]] <- means2[[m]][-(same.cont[same.cont[, 3] == m, 4] -
                                           K.cat[m])]
            vars2[[m]] <- vars2[[m]][-(same.cont[same.cont[, 3] == m, 4] -
                                         K.cat[m])]
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
      if (K.pois[m] != 0 & length(grep("pois", colnames(X0[[m]]))) != 0) {
        X_pois[[m]] <- X0[[m]][, grep("pois", colnames(X0[[m]])), drop = FALSE]
        Y_pois[[m]] <- matrix(1, nrow = n, ncol = ncol(X_pois[[m]]))
        for (i in 1:ncol(X_pois[[m]])) {
          Y_pois[[m]][, i] <- qzipois(pnorm(X_pois[[m]][, i]), lam[[m]][i],
            pstr0 = p_zip[[m]][i])
        }
        colnames(Y_pois[[m]]) <- colnames(X_pois[[m]])
      }
      if (K.nb[m] != 0 & length(grep("xnb", colnames(X0[[m]]))) != 0) {
        X_nb[[m]] <- X0[[m]][, grep("xnb", colnames(X0[[m]])), drop = FALSE]
        Y_nb[[m]] <- matrix(1, nrow = n, ncol = ncol(X_nb[[m]]))
        for (i in 1:ncol(X_nb[[m]])) {
          Y_nb[[m]][, i] <- qzinegbin(pnorm(X_nb[[m]][, i]),
            size = size[[m]][i], munb = mu[[m]][i], pstr0 = p_zinb[[m]][i])
        }
        colnames(Y_nb[[m]]) <- colnames(X_nb[[m]])
      }
    }
    Y_all <- NULL
    for (p in 1:M) {
      Y_all <- cbind(Y_all, Y_cat[[p]], Y_cont[[p]], Y_pois[[p]], Y_nb[[p]])
    }
    Y_all <- Y_all[, order(colnames(Y_all)), drop = FALSE]
    rho.x0 <- cor(Y_all)
    emax0 <- max(abs(rho.x0 - Corr_X))
    niter <- diag(0, nrow(rho.x0), ncol(rho.x0))
    Sigma_X2 <- Sigma_X[names1, names1, drop = FALSE]
    start.time.error <- Sys.time()
    if (emax0 > epsilon & errorloop == TRUE) {
      EL <- corr_error(n = n, k_cat = length(grep("cat", colnames(Y_all))),
        k_cont = length(grep("cont", colnames(Y_all))),
        k_pois = length(grep("pois", colnames(Y_all))),
        k_nb = length(grep("xnb", colnames(Y_all))), method = method,
        means = temp_means, vars = temp_vars, constants = constants2,
        marginal = marginal2, support = support2, lam = lam2, p_zip = p_zip2,
        size = size2, mu = mu2, p_zinb = p_zinb2, seed = seed,
        epsilon = epsilon, maxit = maxit, rho0 = Corr_X, Sigma = Sigma_X,
        rho_calc = rho.x0)
      Y_all <- cbind(EL$Y_cat, EL$Y_cont, EL$Y_pois, EL$Y_nb)
      colnames(Y_all) <- colnames(Corr_X)
      niter <- EL$niter
      colnames(niter) <- colnames(Corr_X)
      rownames(niter) <- colnames(Corr_X)
      Sigma_X2 <- EL$Sigma
      Sigma_X2 <- Sigma_X2[names1, names1, drop = FALSE]
      niter <- niter[names1, names1, drop = FALSE]
    }
    stop.time.error <- Sys.time()
    Y_all <- Y_all[, names1, drop = FALSE]
    Y0 <- list()
    for (i in 1:M) {
      if (K.x[i] != 0) {
        if (i == 1) {
          Y0[[i]] <- Y_all[, 1:K.x[i], drop = FALSE]
        } else {
          if (!is.null(same.var)) {
            if (K.x[i] == nrow(same.var[same.var[, 3] == i, , drop = FALSE])) {
              Y0 <- append(Y0, list(NULL))
            } else {
              Y0[[i]] <- Y_all[, (cumsum(K.x)[i - 1] -
                nrow(same.var[same.var[, 3] <= (i - 1), , drop = FALSE]) +
                1):(cumsum(K.x)[i] - nrow(same.var[same.var[, 3] <= i, ,
                drop = FALSE])), drop = FALSE]
            }
          } else {
            Y0[[i]] <- Y_all[, (cumsum(K.x)[i - 1] + 1):cumsum(K.x)[i],
                         drop = FALSE]
          }
        }
      } else Y0 <- append(Y0, list(NULL))
    }
    Y_cat <- list()
    Y_cont <- list()
    Y_pois <- list()
    Y_nb <- list()
    for (m in 1:M) {
      Y_cat <- append(Y_cat, list(NULL))
      Y_cont <- append(Y_cont, list(NULL))
      Y_pois <- append(Y_pois, list(NULL))
      Y_nb <- append(Y_nb, list(NULL))
      if (K.cat[m] > 0) {
        if (length(grep("cat", colnames(Y0[[m]]))) > 0)
          Y_cat[[m]] <- Y0[[m]][, grep("cat", colnames(Y0[[m]])), drop = FALSE]
        if (!is.null(same.ord)) {
          if (nrow(same.ord[same.ord[, 3] == m, , drop = FALSE]) > 0) {
            temp.ord <- same.ord[same.ord[, 3] == m, , drop = FALSE]
            for (i in 1:nrow(temp.ord)) {
              Y_cat[[m]] <- cbind(Y_cat[[m]],
                Y_cat[[temp.ord[i, 1]]][, temp.ord[i, 2]])
              colnames(Y_cat[[m]])[ncol(Y_cat[[m]])] <-
                colnames(Y_cat[[temp.ord[i, 1]]])[temp.ord[i, 2]]
            }
          }
        }
        if (length(Y_cat[[m]]) != 0) {
          Y_cat[[m]] <- Y_cat[[m]][, names0[[m]][grep("cat", names0[[m]])],
                                   drop = FALSE]
        }
      }
      if ((K.x[m] - K.cat[m] - K.pois[m] - K.nb[m]) > 0) {
        if (length(grep("cont", colnames(Y0[[m]]))) > 0)
          Y_cont[[m]] <- Y0[[m]][, grep("cont", colnames(Y0[[m]])),
                                 drop = FALSE]
        if (!is.null(same.cont)) {
          if (nrow(same.cont[same.cont[, 3] == m, , drop = FALSE]) > 0) {
            temp.cont <- same.cont[same.cont[, 3] == m, , drop = FALSE]
            for (i in 1:nrow(temp.cont)) {
              temp.cont[i, 2] <- temp.cont[i, 2] - K.cat[temp.cont[i, 1]]
              Y_cont[[m]] <- cbind(Y_cont[[m]],
                Y_cont[[temp.cont[i, 1]]][, temp.cont[i, 2]])
              colnames(Y_cont[[m]])[ncol(Y_cont[[m]])] <-
                colnames(Y_cont[[temp.cont[i, 1]]])[temp.cont[i, 2]]
            }
          }
        }
        if (length(Y_cont[[m]]) != 0) {
          Y_cont[[m]] <- Y_cont[[m]][, names0[[m]][grep("cont", names0[[m]])],
                                     drop = FALSE]
        }
      }
      if (K.pois[m] > 0) {
        if (length(grep("pois", colnames(Y0[[m]]))) > 0)
          Y_pois[[m]] <- Y0[[m]][, grep("pois", colnames(Y0[[m]])),
                                 drop = FALSE]
        if (!is.null(same.pois)) {
          if (nrow(same.pois[same.pois[, 3] == m, , drop = FALSE]) > 0) {
            temp.pois <- same.pois[same.pois[, 3] == m, , drop = FALSE]
            for (i in 1:nrow(temp.pois)) {
              temp.pois[i, 2] <- temp.pois[i, 2] - K.x[temp.pois[i, 1]] +
                K.pois[temp.pois[i, 1]] + K.nb[temp.pois[i, 1]]
              Y_pois[[m]] <- cbind(Y_pois[[m]],
                Y_pois[[temp.pois[i, 1]]][, temp.pois[i, 2]])
              colnames(Y_pois[[m]])[ncol(Y_pois[[m]])] <-
                colnames(Y_pois[[temp.pois[i, 1]]])[temp.pois[i, 2]]
            }
          }
        }
        if (length(Y_pois[[m]]) != 0) {
          Y_pois[[m]] <- Y_pois[[m]][, names0[[m]][grep("pois", names0[[m]])],
                                     drop = FALSE]
        }
      }
      if (K.nb[m] > 0) {
        if (length(grep("xnb", colnames(Y0[[m]]))) > 0)
          Y_nb[[m]] <- Y0[[m]][, grep("xnb", colnames(Y0[[m]])), drop = FALSE]
        if (!is.null(same.nb)) {
          if (nrow(same.nb[same.nb[, 3] == m, , drop = FALSE]) > 0) {
            temp.nb <- same.nb[same.nb[, 3] == m, , drop = FALSE]
            for (i in 1:nrow(temp.nb)) {
              temp.nb[i, 2] <- temp.nb[i, 2] - K.x[temp.nb[i, 1]] +
                K.nb[temp.nb[i, 1]]
              Y_nb[[m]] <- cbind(Y_nb[[m]],
                Y_nb[[temp.nb[i, 1]]][, temp.nb[i, 2]])
              colnames(Y_nb[[m]])[ncol(Y_nb[[m]])] <-
                colnames(Y_nb[[temp.nb[i, 1]]])[temp.nb[i, 2]]
            }
          }
        }
        if (length(Y_nb[[m]]) != 0) {
          Y_nb[[m]] <- Y_nb[[m]][, names0[[m]][grep("xnb", names0[[m]])],
                                 drop = FALSE]
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
            ind <- K.x[i] - K.comp[i] + K.error[i] - K.pois[i] - K.nb[i] + j
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
    Y_all <- list()
    for (p in 1:M) {
      if (K.x[p] == 0) {
        Y_all <- append(Y_all, list(NULL))
        next
      }
      Y_all[[p]] <- cbind(Y_cat[[p]], Y_comp[[p]], Y_pois[[p]], Y_nb[[p]])
    }
  }
  Y <- matrix(1, nrow = n, ncol = M)
  X_all2 <- list()
  Time <- matrix(Time, n, M, byrow = TRUE)
  colnames(Time) <- colnames(Time, do.NULL = FALSE, prefix = "T")
  int.var2 <- NULL
  subj.var2 <- NULL
  tint.var2 <- NULL
  group.var2 <- NULL
  betas.tint2 <- NULL
  for (i in 1:M) {
    if (K.x[i] == 0) {
      Y[, i] <- matrix(betas.0[i], n, 1) +
        Time[, i] %*% matrix(betas.t[i], 1, 1) + E2[, i]
      X_all2 <- append(X_all2, list(NULL))
      next
    }
    X_all <- cbind(Y_cat[[i]], Y_cont[[i]], Y_pois[[i]], Y_nb[[i]])
    colnames(X_all) <- gsub("cat", "ord", colnames(X_all))
    colnames(X_all) <- gsub("xnb", "nb", colnames(X_all))
    if (length(betas) == 0) betas2 <- matrix(0, ncol(X_all), ncol = 1)
    if (length(betas) == 1) {
      if (length(betas[[1]]) == 1)
        betas2 <- matrix(rep(betas[[1]], ncol(X_all)), ncol = 1) else
          if (length(betas[[1]]) < ncol(X_all))
            betas2 <- matrix(c(betas[[1]],
              rep(0, ncol(X_all) - length(betas[[1]]))), ncol = 1) else
                betas2 <- matrix(betas[[1]], ncol = 1)
    }
    if (length(betas) == M) {
      if (length(betas[[i]]) == 1)
        betas2 <- matrix(rep(betas[[i]], ncol(X_all)), ncol = 1) else
          if (length(betas[[i]]) < ncol(X_all))
            betas2 <- matrix(c(betas[[i]],
              rep(0, ncol(X_all) - length(betas[[i]]))), ncol = 1) else
                betas2 <- matrix(betas[[i]], ncol = 1)
    }
    if (!(length(betas) %in% c(1, M)))
      stop("Length of betas should be 1 or M.")
    if (!is.null(int.var)) {
      int.var2 <- int.var[int.var[, 1] == i, , drop = FALSE]
      if (length(betas.int) == 0)
        betas.int2 <- matrix(0, nrow(int.var2), ncol = 1)
      if (length(betas.int) == 1) {
        if (length(betas.int[[1]]) == 1)
          betas.int2 <- matrix(rep(betas.int[[1]], nrow(int.var2)),
                               ncol = 1) else
            if (length(betas.int[[1]]) < nrow(int.var2))
              betas.int2 <- matrix(c(betas.int[[1]],
                rep(0, nrow(int.var2) - length(betas.int[[1]]))),
                ncol = 1) else
                  betas.int2 <- matrix(betas.int[[1]], ncol = 1)
      }
      if (length(betas.int) == M) {
        if (length(betas.int[[i]]) == 1)
          betas.int2 <- matrix(rep(betas.int[[i]], nrow(int.var2)),
                               ncol = 1) else
            if (length(betas.int[[i]]) < nrow(int.var2))
              betas.int2 <- matrix(c(betas.int[[i]],
                rep(0, nrow(int.var2) - length(betas.int[[i]]))),
                ncol = 1) else
                  if (length(betas.int[[i]]) == 0)
                    betas.int2 <- NULL else
                      betas.int2 <- matrix(betas.int[[i]], ncol = 1)
      }
      if (!(length(betas.int) %in% c(0, 1, M)))
        stop("Length of betas.int should be 0, 1, or M.")
    }
    if (!is.null(subj.var)) {
      subj.var2 <- subj.var[subj.var[, 1] == i, , drop = FALSE]
      if (length(betas.subj) == 0)
        betas.subj2 <- matrix(0, nrow(subj.var2), ncol = 1)
      if (length(betas.subj) == 1) {
        if (length(betas.subj[[1]]) == 1)
          betas.subj2 <- matrix(rep(betas.subj[[1]], nrow(subj.var2)),
                                ncol = 1) else
            if (length(betas.subj[[1]]) < nrow(subj.var2))
              betas.subj2 <- matrix(c(betas.subj[[1]],
                rep(0, nrow(subj.var2) - length(betas.subj[[1]]))),
                ncol = 1) else
                  betas.subj2 <- matrix(betas.subj[[1]], ncol = 1)
      }
      if (length(betas.subj) == M) {
        if (length(betas.subj[[i]]) == 1)
          betas.subj2 <- matrix(rep(betas.subj[[i]], nrow(subj.var2)),
                                ncol = 1) else
            if (length(betas.subj[[i]]) < nrow(subj.var2))
              betas.subj2 <- matrix(c(betas.subj[[i]],
                rep(0, nrow(subj.var2) - length(betas.subj[[i]]))),
                ncol = 1) else
                  if (length(betas.subj[[i]]) == 0)
                    betas.subj2 <- NULL else
                      betas.subj2 <- matrix(betas.subj[[i]], ncol = 1)
      }
      if (!(length(betas.subj) %in% c(0, 1, M)))
        stop("Length of betas.subj should be 0, 1, or M.")
    }
    X_all2[[i]] <- X_all
    Y[, i] <- matrix(betas.0[i], n, 1) + X_all %*% betas2 +
      Time[, i] %*% matrix(betas.t[i], 1, 1)
    if (!is.null(int.var2)) {
      if (nrow(int.var2) != 0) {
        X_int <- matrix(1, n, nrow(int.var2),
                        dimnames = list(1:n, 1:nrow(int.var2)))
        for (j in 1:nrow(int.var2)) {
          X_int[, j] <- X_all[, int.var2[j, 2]] * X_all[, int.var2[j, 3]]
          colnames(X_int)[j] <- paste(colnames(X_all)[int.var2[j, 2]],
            colnames(X_all)[int.var2[j, 3]], sep = "_")
          if (!is.null(subj.var2)) {
            if (int.var2[j, 2] %in% subj.var2[, 2] &
                int.var2[j, 3] %in% subj.var2[, 2])
              subj.var2 <- rbind(subj.var2, c(i, ncol(X_all) + j))
          }
        }
        Y[, i] <- Y[, i] + X_int %*% betas.int2
        X_all2[[i]] <- cbind(X_all2[[i]], X_int)
      }
    }
    if (is.null(subj.var2)) {
      group.var2 <- 1:ncol(X_all2[[i]])
    } else group.var2 <- which(!(col(X_all2[[i]])[1, ] %in% subj.var2[, 2]))
    if (!is.null(subj.var2)) {
      if (nrow(subj.var2) != 0 & length(group.var2) != 0) {
        X_subj <- matrix(0, n, 1, dimnames = list(1:n, 1))
        for (j in 1:nrow(subj.var2)) {
          for (k in 1:ncol(X_all2[[i]])) {
            if (!(k %in% subj.var2[, 2])) {
              X_subj <- cbind(X_subj, X_all2[[i]][, k] *
                                X_all2[[i]][, subj.var2[j, 2]])
              colnames(X_subj)[ncol(X_subj)] <- paste(colnames(X_all2[[i]])[k],
                colnames(X_all2[[i]])[subj.var2[j, 2]], sep = "_")
            }
          }
        }
        X_subj <- X_subj[, -1, drop = FALSE]
        Y[, i] <- Y[, i] + X_subj %*% betas.subj2
        X_all2[[i]] <- cbind(X_all2[[i]], X_subj)
      }
    }
    if (!is.null(tint.var)) {
      tint.var2 <- tint.var[tint.var[, 1] == i, , drop = FALSE]
      if (nrow(tint.var2) > 0) group.var2 <- tint.var2[, 2]
    }
    if (length(group.var2) > 0) {
      if (length(betas.tint) == 0)
        betas.tint2 <- matrix(0, length(group.var2), ncol = 1)
      if (length(betas.tint) == 1) {
        if (length(betas.tint[[1]]) == 1)
          betas.tint2 <- matrix(rep(betas.tint[[1]], length(group.var2)),
                                ncol = 1) else
            if (length(betas.tint[[1]]) < length(group.var2))
              betas.tint2 <- matrix(c(betas.tint[[1]],
                rep(0, length(group.var2) - length(betas.tint[[1]]))),
                ncol = 1) else
                  betas.tint2 <- matrix(betas.tint[[1]], ncol = 1)
      }
      if (length(betas.tint) == M) {
        if (length(betas.tint[[i]]) == 1)
          betas.tint2 <- matrix(rep(betas.tint[[i]], length(group.var2)),
                                ncol = 1) else
            if (length(betas.tint[[i]]) < length(group.var2))
              betas.tint2 <- matrix(c(betas.tint[[i]],
                rep(0, length(group.var2) - length(betas.tint[[i]]))),
                ncol = 1) else
                  if (length(betas.tint[[i]]) == 0)
                    betas.tint2 <- NULL else
                      betas.tint2 <- matrix(betas.tint[[i]], ncol = 1)
      }
      if (!(length(betas.tint) %in% c(0, 1, M)))
        stop("Length of betas.tint should be 0, 1, or M.")
    }
    if (length(betas.tint2) != 0 & length(group.var2) != 0) {
      X_tint <- matrix(1, n, length(group.var2),
                       dimnames = list(1:n, 1:length(group.var2)))
      for (j in 1:length(group.var2)) {
        X_tint[, j] <- X_all2[[i]][, group.var2[j]] * Time[, i]
        colnames(X_tint)[j] <- paste(colnames(X_all2[[i]])[group.var2[j]],
                                     colnames(Time)[i], sep = "_")
      }
      Y[, i] <- Y[, i] + X_tint %*% betas.tint2
      X_all2[[i]] <- cbind(X_all2[[i]], X_tint)
    }
    Y[, i] <- Y[, i] + E2[, i]
  }
  if ((class(corr.u) == "list" & length(corr.u) > 0) |
      class(corr.u) == "matrix") {
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
    K.rcont <- lengths(rvars) - K.rmix
    k.rcont <- c(0, cumsum(K.rcont))
    k.rmix <- c(0, cumsum(K.rmix))
    rmeans2 <- list()
    rvars2 <- list()
    for (i in 1:M0) {
      rmeans2 <- append(rmeans2, list(NULL))
      rvars2 <- append(rvars2, list(NULL))
      if (K.r[i] == 0) next
      if (K.rmix[i] > 0) ind <- c(0, cumsum(lengths(rmix_pis[[i]])))
      if ((rand.int[i] == "none" & rand.tsl[i] == "mix") |
          (rand.int[i] == "mix" & rand.tsl[i] == "none")) {
        rconstants[[i]] <- rbind(rconstants[[i]][(K.rcont[i] + 1):(K.rcont[i] +
                                 ind[2]), ],
          rconstants[[i]][-c((K.rcont[i] + 1):(K.rcont[i] + ind[2])), ])
        if (K.rcont[i] > 0) {
          rmeans2[[i]] <- c(rmix_mus[[i]][[1]],
            rmeans[[i]][2:(2 + K.rcont[i] - 1)], unlist(rmix_mus[[i]][-1]))
          rvars2[[i]] <- c((rmix_sigmas[[i]][[1]])^2,
            rvars[[i]][2:(2 + K.rcont[i] - 1)],
            (unlist(rmix_sigmas[[i]][-1]))^2)
        } else {
          rmeans2[[i]] <- unlist(rmix_mus[[i]])
          rvars2[[i]] <- (unlist(rmix_sigmas[[i]]))^2
        }
      } else if (rand.int[i] == "non_mix" & rand.tsl[i] == "mix") {
        rconstants[[i]] <- rbind(rconstants[[i]][1, ],
          rconstants[[i]][(K.rcont[i] + 1):(K.rcont[i] + ind[2]), ],
          rconstants[[i]][-c(1, (K.rcont[i] + 1):(K.rcont[i] + ind[2])), ])
        if (K.rcont[i] > 1) {
          rmeans2[[i]] <- c(rmeans[[i]][1], rmix_mus[[i]][[1]],
            rmeans[[i]][3:(3 + K.rcont[i] - 2)], unlist(rmix_mus[[i]][-1]))
          rvars2[[i]] <- c(rvars[[i]][1], (rmix_sigmas[[i]][[1]])^2,
            rvars[[i]][3:(3 + K.rcont[i] - 2)],
            (unlist(rmix_sigmas[[i]][-1]))^2)
        } else {
          rmeans2[[i]] <- c(rmeans[[i]][1], unlist(rmix_mus[[i]]))
          rvars2[[i]] <- c(rvars[[i]][1], (unlist(rmix_sigmas[[i]]))^2)
        }
      } else if (rand.int[i] == "mix" & rand.tsl[i] == "non_mix") {
        rconstants[[i]] <- rbind(rconstants[[i]][(K.rcont[i] + 1):(K.rcont[i] +
                                ind[2]), ], rconstants[[i]][1, ],
          rconstants[[i]][-c(1, (K.rcont[i] + 1):(K.rcont[i] + ind[2])), ])
        if (K.rcont[i] > 1) {
          rmeans2[[i]] <- c(rmix_mus[[i]][[1]], rmeans[[i]][2],
            rmeans[[i]][3:(3 + K.rcont[i] - 2)], unlist(rmix_mus[[i]][-1]))
          rvars2[[i]] <- c((rmix_sigmas[[i]][[1]])^2, rvars[[i]][2],
            rvars[[i]][3:(3 + K.rcont[i] - 2)],
            (unlist(rmix_sigmas[[i]][-1]))^2)
        } else {
          rmeans2[[i]] <- c(rmix_mus[[i]][[1]], rmeans[[i]][2],
                            unlist(rmix_mus[[i]][-1]))
          rvars2[[i]] <- c((rmix_sigmas[[i]][[1]])^2, rvars[[i]][2],
                           (unlist(rmix_sigmas[[i]][-1]))^2)
        }
      } else if (rand.int[i] == "mix" & rand.tsl[i] == "mix") {
        rconstants[[i]] <- rbind(rconstants[[i]][(K.rcont[i] + 1):(K.rcont[i] +
                                  ind[3]), ],
          rconstants[[i]][-c((K.rcont[i] + 1):(K.rcont[i] + ind[3])), ])
        if (K.rcont[i] > 0) {
          rmeans2[[i]] <- c(unlist(rmix_mus[[i]][1:2]),
            rmeans[[i]][3:(3 + K.rcont[i] - 1)],
            unlist(rmix_mus[[i]][-c(1:2)]))
          rvars2[[i]] <- c((unlist(rmix_sigmas[[i]][1:2]))^2,
            rvars[[i]][3:(3 + K.rcont[i] - 1)],
            (unlist(rmix_sigmas[[i]][-c(1:2)]))^2)
        } else {
          rmeans2[[i]] <- unlist(rmix_mus[[i]])
          rvars2[[i]] <- (unlist(rmix_sigmas[[i]]))^2
        }
      } else {
        if (K.rcont[i] > 0) {
          rmeans2[[i]] <- rmeans[[i]][1:K.rcont[i]]
          rvars2[[i]] <- rvars[[i]][1:K.rcont[i]]
        }
        if (K.rmix[i] > 0) {
          rmeans2[[i]] <- c(rmeans2[[i]], unlist(rmix_mus[[i]]))
          rvars2[[i]] <- c(rvars2[[i]], (unlist(rmix_sigmas[[i]]))^2)
        }
      }
    }
    rconstants2 <- NULL
    for (p in 1:M0) {
      if (K.r[p] > 0) rconstants2 <- rbind(rconstants2, rconstants[[p]])
    }
    if (class(corr.u) == "list") {
      Corr_U <- list()
      for (p in 1:M0) {
        if (K.r[p] > 0) {
          Corr_U[[p]] <- do.call(cbind, corr.u[[p]])
        } else Corr_U <- append(Corr_U, list(NULL))
      }
      Corr_U <- do.call(rbind, Corr_U)
    } else Corr_U <- corr.u
    Sigma_U <- intercorr2(k_cont = ncol(Corr_U), method = method,
      constants = rconstants2, rho = Corr_U, quiet = quiet)
    if (min(eigen(Sigma_U, symmetric = TRUE)$values) < 0) {
      if (quiet == FALSE)
        message("Intermediate random correlation matrix is not positive
                definite.")
      if (use.nearPD == TRUE) {
        Sigma_U <- as.matrix(nearPD(Sigma_U, corr = T, keepDiag = T)$mat)
      } else if (adjgrad == TRUE) {
        sadj <- adj_grad(Sigma_U, B1, tau, tol, steps, msteps)
        Sigma_U <- sadj$Sigma2
      }
    }
    eig <- eigen(Sigma_U, symmetric = TRUE)
    sqrteigval <- diag(sqrt(pmax(eig$values, 0)), nrow(Sigma_U))
    eigvec <- eig$vectors
    fry <- eigvec %*% sqrteigval
    U <- Z[, (ncol(Z) - ncol(Sigma_U) + 1):ncol(Z), drop = FALSE]
    U <- fry %*% t(U)
    U <- t(U)
    U0 <- list()
    U_comp <- list()
    for (m in 1:M0) {
      U0 <- append(U0, list(NULL))
      U_comp <- append(U_comp, list(NULL))
      if (K.r[m] > 0) {
        if (m == 1) {
          U0[[m]] <- U[, 1:K.r[m], drop = FALSE]
        } else {
          U0[[m]] <- U[, (cumsum(K.r)[m - 1] + 1):cumsum(K.r)[m], drop = FALSE]
        }
        U_comp[[m]] <- matrix(1, nrow = n, ncol = K.r[m])
        for (i in 1:K.r[m]) {
          if (method == "Fleishman") {
            U_comp[[m]][, i] <- rconstants[[m]][i, 1] +
              rconstants[[m]][i, 2] * U0[[m]][, i] +
              rconstants[[m]][i, 3] * U0[[m]][, i]^2 +
              rconstants[[m]][i, 4] * U0[[m]][, i]^3
          }
          if (method == "Polynomial") {
            U_comp[[m]][, i] <- rconstants[[m]][i, 1] +
              rconstants[[m]][i, 2] * U0[[m]][, i] +
              rconstants[[m]][i, 3] * U0[[m]][, i]^2 +
              rconstants[[m]][i, 4] * U0[[m]][, i]^3 +
              rconstants[[m]][i, 5] * U0[[m]][, i]^4 +
              rconstants[[m]][i, 6] * U0[[m]][, i]^5
          }
          U_comp[[m]][, i] <- rmeans2[[m]][i] + sqrt(rvars2[[m]][i]) *
            U_comp[[m]][, i]
        }
      }
    }
    U_cont <- U_comp
    if (max(K.rmix) > 0) {
      for (i in 1:M0) {
        if (K.rmix[i] == 0) next
        ind <- c(0, cumsum(lengths(rmix_pis[[i]])))
        if ((rand.int[i] == "none" & rand.tsl[i] == "mix") |
           (rand.int[i] == "mix" & rand.tsl[i] == "none")) {
          seed <- seed + 1
          set.seed(seed)
          R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[1]])
          U_cont[[i]] <- apply(t(R) *
            U_comp[[i]][, 1:ind[2], drop = FALSE], 1, sum)
          U_cont[[i]] <- matrix(rmeans[[i]][1] + sqrt(rvars[[i]][1]) *
            scale(U_cont[[i]]), n)
          if (K.rcont[i] > 0)
            U_cont[[i]] <- cbind(U_cont[[i]],
              U_comp[[i]][, (ind[2] + 1):(ind[2] + K.rcont[i])])
          if (K.rmix[i] > 1) {
            for (j in 2:K.rmix[i]) {
              seed <- seed + 1
              set.seed(seed)
              R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[j]])
              U_cont[[i]] <- cbind(U_cont[[i]], apply(t(R) *
                U_comp[[i]][, (ind[j] + K.rcont[i] + 1):(ind[j + 1] +
                  K.rcont[i]), drop = FALSE], 1, sum))
              U_cont[[i]][, K.rcont[i] + j] <- rmeans[[i]][K.rcont[i] + j] +
                sqrt(rvars[[i]][K.rcont[i] + j]) *
                scale(U_cont[[i]][, K.rcont[i] + j])
            }
          }
        } else if (rand.int[i] == "non_mix" & rand.tsl[i] == "mix") {
          U_cont[[i]] <- U_cont[[i]][, 1:2]
          seed <- seed + 1
          set.seed(seed)
          R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[1]])
          U_cont[[i]][, 2] <- apply(t(R) *
            U_comp[[i]][, 2:(ind[2] + 1), drop = FALSE], 1, sum)
          U_cont[[i]][, 2] <- rmeans[[i]][2] + sqrt(rvars[[i]][2]) *
            scale(U_cont[[i]][, 2])
          if (K.rcont[i] > 1)
            U_cont[[i]] <- cbind(U_cont[[i]][, 1:2],
              U_comp[[i]][, (ind[2] + 2):(ind[2] + K.rcont[i]), drop = FALSE])
          if (K.rmix[i] > 1) {
            for (j in 2:K.rmix[i]) {
              seed <- seed + 1
              set.seed(seed)
              R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[j]])
              U_cont[[i]] <- cbind(U_cont[[i]], apply(t(R) *
                U_comp[[i]][, (ind[j] + K.rcont[i] + 1):(ind[j + 1] +
                  K.rcont[i]), drop = FALSE], 1, sum))
              U_cont[[i]][, K.rcont[i] + j] <- rmeans[[i]][K.rcont[i] + j] +
                sqrt(rvars[[i]][K.rcont[i] + j]) *
                scale(U_cont[[i]][, K.rcont[i] + j])
            }
          }
        } else if (rand.int[i] == "mix" & rand.tsl[i] == "non_mix") {
          seed <- seed + 1
          set.seed(seed)
          R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[1]])
          U_cont[[i]] <- apply(t(R) *
            U_comp[[i]][, 1:ind[2], drop = FALSE], 1, sum)
          U_cont[[i]] <- rmeans[[i]][1] + sqrt(rvars[[i]][1]) *
            scale(U_cont[[i]])
          U_cont[[i]] <- cbind(U_cont[[i]],
              U_comp[[i]][, (ind[2] + 1):(ind[2] + K.rcont[i])])
          if (K.rmix[i] > 1) {
            for (j in 2:K.rmix[i]) {
              seed <- seed + 1
              set.seed(seed)
              R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[j]])
              U_cont[[i]] <- cbind(U_cont[[i]], apply(t(R) *
                U_comp[[i]][, (ind[j] + K.rcont[i] + 1):(ind[j + 1] +
                  K.rcont[i]), drop = FALSE], 1, sum))
              U_cont[[i]][, K.rcont[i] + j] <- rmeans[[i]][K.rcont[i] + j] +
                sqrt(rvars[[i]][K.rcont[i] + j]) *
                scale(U_cont[[i]][, K.rcont[i] + j])
            }
          }
        } else if (rand.int[i] == "mix" & rand.tsl[i] == "mix") {
          seed <- seed + 1
          set.seed(seed)
          R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[1]])
          U_cont[[i]] <- apply(t(R) *
            U_comp[[i]][, 1:ind[2], drop = FALSE], 1, sum)
          U_cont[[i]] <- rmeans[[i]][1] + sqrt(rvars[[i]][1]) *
            scale(U_cont[[i]])
          seed <- seed + 1
          set.seed(seed)
          R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[2]])
          U_cont[[i]] <- cbind(U_cont[[i]], apply(t(R) *
            U_comp[[i]][, (ind[2] + 1):ind[3], drop = FALSE], 1, sum))
          U_cont[[i]][, 2] <- rmeans[[i]][2] + sqrt(rvars[[i]][2]) *
            scale(U_cont[[i]][, 2])
          if (K.rcont[i] > 0)
            U_cont[[i]] <- cbind(U_cont[[i]],
              U_comp[[i]][, (ind[3] + 1):(ind[3] + K.rcont[i]), drop = FALSE])
          if (K.rmix[i] > 2) {
            for (j in 3:K.rmix[i]) {
              seed <- seed + 1
              set.seed(seed)
              R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[j]])
              U_cont[[i]] <- cbind(U_cont[[i]], apply(t(R) *
                U_comp[[i]][, (ind[j] + K.rcont[i] + 1):(ind[j + 1] +
                  K.rcont[i]), drop = FALSE], 1, sum))
              U_cont[[i]][, K.rcont[i] + j] <- rmeans[[i]][K.rcont[i] + j] +
                sqrt(rvars[[i]][K.rcont[i] + j]) *
                scale(U_cont[[i]][, K.rcont[i] + j])
            }
          }
        } else {
          U_cont[[i]] <- matrix(1, n, K.rcont[i] + K.rmix[i])
          if (K.rcont[i] > 0)
            U_cont[[i]][, 1:K.rcont[i]] <- U_comp[[i]][, 1:K.rcont[i]]
          for (j in 1:K.rmix[i]) {
            seed <- seed + 1
            set.seed(seed)
            R <- rmultinom(n, size = 1, prob = rmix_pis[[i]][[j]])
            U_cont[[i]][, K.rcont[i] + j] <- apply(t(R) *
              U_comp[[i]][, (ind[j] + K.rcont[i] + 1):(ind[j + 1] +
                K.rcont[i]), drop = FALSE], 1, sum)
            U_cont[[i]][, K.rcont[i] + j] <- rmeans[[i]][K.rcont[i] + j] +
              sqrt(rvars[[i]][K.rcont[i] + j]) *
              scale(U_cont[[i]][, K.rcont[i] + j])
          }
        }
      }
    }
    if (class(corr.u) == "matrix") {
      U_cont <- lapply(seq_len(M), function(x) U_cont[[1]])
      U_comp <- lapply(seq_len(M), function(x) U_comp[[1]])
    }
    V <- list()
    for (i in 1:M) {
      if (K.r[i] == 0) {
        V <- append(V, list(NULL))
        next
      }
      V[[i]] <- matrix(0, n, 1)
      if (rand.int[i] != "none") {
        V[[i]] <- cbind(V[[i]], matrix(1, n, 1, dimnames = list(1:n, "int")))
      }
      if (rand.tsl[i] != "none") {
        V[[i]] <- cbind(V[[i]], Time[, i, drop = FALSE])
      }
      if (!is.null(rand.var)) {
        rand.var2 <- rand.var[rand.var[, 1] == i, , drop = FALSE]
        if (nrow(rand.var2) > 0) {
          for (j in 1:nrow(rand.var2)) {
            V[[i]] <- cbind(V[[i]], X_all2[[i]][, rand.var2[j, 2],
                                                drop = FALSE])
            colnames(V[[i]])[ncol(V[[i]])] <-
              colnames(X_all2[[i]])[rand.var2[j, 2]]
          }
        }
      }
      V[[i]] <- V[[i]][, -1, drop = FALSE]
      Y[, i] <- Y[, i] + apply(V[[i]] * U_cont[[i]], 1, sum)
      colnames(U_comp[[i]]) <- paste("U", 1:ncol(U_comp[[i]]), sep = "")
      colnames(U_cont[[i]]) <- paste("U", colnames(V[[i]]), sep = "_")
    }
  }
  colnames(Y) <- paste("Y", 1:M, sep = "")
  result <- list(Y = Y, E = E_comp)
  if (error_type == "mix") result <- append(result, list(E_mix = E2))
  Time.error <- 0
  if (length(corr.x) > 0) {
    for (i in 1:M) {
      X_all2[[i]] <- cbind(X_all2[[i]], Time[, i, drop = FALSE])
    }
    result <- append(result, list(X = Y_all, X_all = X_all2,
      Sigma_X0 = Sigma_X0, Sigma_X = Sigma_X2, niter = niter))
    Time.error <- round(difftime(stop.time.error, start.time.error,
                                 units = "min"), 3)
  }
  if (max(K.r) > 0) {
    constants0 <- append(constants0, rconstants)
    result <- append(result, list(U = U_comp, U_all = U_cont, V = V,
      rmeans2 = rmeans2, rvars2 = rvars2))
  }
  stop.time <- Sys.time()
  Time <- round(difftime(stop.time, start.time, units = "min"), 3)
  cat("Total Simulation time:", Time, "minutes \n")
  result <- append(result, list(constants = constants0, SixCorr = SixCorr,
    valid.pdf = Valid.PDF, Error_Time = Time.error, Time = Time))
  result
}

#' Causal Mediation Analysis for Heterogeneous Effect
#'
#' This function is for the estimation of various causal heterogeneous effects for causal mediation analysis, including natural direct and indirect effect heterogeneity, pure direct and indirect effect heterogeneity, mediated interactive effect heterogeneity, controlled direct effect heterogeneity, and total effect heterogeneity.
#' @param model.m A fitted model object for mediator.
#' @param model.y A fitted model object for outcome.
#' @param modifier A character string indicating the name of the modifier used in the models.
#' @param treat A character string indicating the name of the treatment indicator used in the models.
#' @param mediator A character string indicating the name of the mediator used in the models.
#' @param n Sample size.
#' @param conf.level Level of the returned two-sided confidence intervals. Default is set to 0.95, which returns the 2.5 and 97.5 percentiles of the quantities.
#' @param mreg Type of mediator model. Can be of class 'linear', 'logistic', 'loglinear'. Default is 'linear'.
#' @param yreg Type of outcome model. Can be of class 'linear', 'logistic', 'loglinear'. Default is 'linear'.
#' @return  The estimate, standard error, p value and confidence interval of natural direct and indirect effect heterogeneity, pure direct and indirect effect heterogeneity, mediated interactive effect heterogeneity, controlled direct effect heterogeneity, and total effect heterogeneity.
#' @export


hetero.mediate <- function(model.m, model.y, modifier, treat, mediator, n, conf.level = 0.95, mreg='linear', yreg='linear'){
  if (mreg=='linear' & yreg=='linear'){
    return(hetero.mediate.mcontiyconti(model.m, model.y, modifier, treat, mediator, n, conf.level))
  } else if (mreg=='linear' & yreg=='logistic'){
    return(hetero.mediate.mcontiybinary(model.m, model.y, modifier, treat, mediator, n, conf.level))
  } else if (mreg=='linear' & yreg=='loglinear'){
    return(hetero.mediate.mcontiybinary(model.m, model.y, modifier, treat, mediator, n, conf.level))
  } else if (mreg=='logistic' & yreg=='linear'){
    return(hetero.mediate.mbinaryyconti(model.m, model.y, modifier, treat, mediator, n, conf.level))
  } else if (mreg=='logistic' & yreg=='logistic'){
    return(hetero.mediate.mbinaryybinary(model.m, model.y, modifier, treat, mediator, n, conf.level))
  } else if (mreg=='logistic' & yreg=='loglinear'){
    return(hetero.mediate.mbinaryybinary(model.m, model.y, modifier, treat, mediator, n, conf.level))
  }else {
    return(print('Error, Please check the input.'))
  }
}

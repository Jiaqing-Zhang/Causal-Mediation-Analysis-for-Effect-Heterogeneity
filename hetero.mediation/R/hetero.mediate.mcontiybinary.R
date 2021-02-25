#' Causal Mediation Analysis for Heterogeneous Effect (Continuous Mediator and Binary Outcome, on the odds ratio scale)
#'
#' This function is for the estimation of various causal heterogeneous effects for causal mediation analysis, including natural direct and indirect effect heterogeneity, pure direct and indirect effect heterogeneity, mediated interactive effect heterogeneity, controlled direct effect heterogeneity, and total effect heterogeneity. The mediator is continuous and the outcome is binary, and estimates are on odds ratio scale.
#' @param model.m A fitted model object for mediator.
#' @param model.y A fitted model object for outcome.
#' @param modifier A character string indicating the name of the modifier used in the models.
#' @param treat A character string indicating the name of the treatment indicator used in the models.
#' @param mediator A character string indicating the name of the mediator used in the models.
#' @param n Sample size.
#' @param conf.level Level of the returned two-sided confidence intervals. Default is set to 0.95, which returns the 2.5 and 97.5 percentiles of the quantities.
#' @return  The estimate, standard error, p value and confidence interval of natural direct and indirect effect heterogeneity, pure direct and indirect effect heterogeneity, mediated interactive effect heterogeneity, controlled direct effect heterogeneity, and total effect heterogeneity on the odds ratio scale.
#' @export



###continuous outcome and continuous mediator
hetero.mediate.mcontiybinary <- function (model.m, model.y, modifier, treat, mediator, n, conf.level = 0.95){
  ##with the modifier and the treatment with two levels and coded as 0, 1. Then the constant of (a-a*)(q-q*) =1
  treatment.modifier <- paste(treat, ':', modifier, sep="")
  treatment.mediator <- paste(treat, ':', mediator, sep="")
  modifier.mediator <- paste(mediator, ':', modifier, sep="")

  model.m.coef<- summary(model.m)$coefficient
  model.y.coef <- summary(model.y)$coefficient

  model.m.cov <- vcov(model.m)
  model.y.cov <- vcov(model.y)
  model.m.sigma.sq <- sigma(model.m)^2 ##this is the squared residual standard error; sigma^2
  model.m.sigma.var <-(2 * (model.m.sigma.sq^2))/summary(model.m)$df[2] ##the variance of sigma^2

  ##The residual standard error is the square root of the residual sum of squares divided by the residual degrees of freedom

  ###variance-covariance matrix for all coefficients and the error term of the mediator model, dimension: 12*12
  cov <- matrix(, nrow = dim(model.m.cov)[1] + dim(model.y.cov)[1] + 1, ncol = dim(model.m.cov)[1] + dim(model.y.cov)[1] + 1)
  cov[1:dim(model.m.cov)[1], 1:dim(model.m.cov)[1]] <- model.m.cov
  cov[(dim(model.m.cov)[1] + 1): (dim(model.m.cov)[1] + dim(model.y.cov)[1]), (dim(model.m.cov)[1] + 1):(dim(model.m.cov)[1] + dim(model.y.cov)[1])] <- model.y.cov
  cov[(dim(model.m.cov)[1] + dim(model.y.cov)[1] + 1), (dim(model.m.cov)[1] + dim(model.y.cov)[1] + 1)] <- model.m.sigma.var
  cov[is.na(cov)]=0


  c.l <- 0.5 + conf.level/2

  ######################NDEH: estimand and sd##############################
  NDEH <- exp(model.y.coef[treatment.mediator, 1] * model.m.coef[modifier, 1] + model.y.coef[treatment.mediator, 1] * model.y.coef[modifier.mediator, 1] * model.m.sigma.sq + model.y.coef[treatment.modifier, 1]) ##exp(theta4*beta2 + theta4*theta6*sigma^2 + theta5)

  NDEH.gamma <- matrix(c(0, ##beta0
                         0, ##beta1
                         model.y.coef[treatment.mediator, 1], ##beta2
                         0, ##beta3
                         0, ##theta0
                         0, ##theta1
                         0, ##theta2
                         0, ##theta3
                         (model.m.coef[modifier, 1] + model.y.coef[modifier.mediator, 1] * model.m.sigma.sq), ##theta4
                         1, ##theta5
                         model.y.coef[treatment.mediator, 1] * model.m.sigma.sq, ##theta6
                         model.y.coef[treatment.mediator, 1] * model.y.coef[modifier.mediator, 1] ##sigma
  ), nrow=1)

  NDEH.sd <- sqrt(NDEH.gamma %*% cov %*% t(NDEH.gamma)) ##on the logarithm scale

  NDEH.z <- log(NDEH)/NDEH.sd
  NDEH.p <- 2*pnorm(-abs(NDEH.z))

  NDEH.upper <- exp(log(NDEH) + qnorm(c.l) * (NDEH.sd/sqrt(n)))
  NDEH.lower <- exp(log(NDEH) - qnorm(c.l) * (NDEH.sd/sqrt(n)))

  NDEH.set <- c(sprintf("%10.5f", c(NDEH, NDEH.sd, NDEH.lower, NDEH.upper)), star.pval(round(NDEH.p, 3)))


  ######################NIEH: estimand and sd##############################
  NIEH <- exp(model.m.coef[treatment.modifier, 1] * (model.y.coef[mediator, 1] + model.y.coef[treatment.mediator, 1]) + model.y.coef[modifier.mediator, 1] * (model.m.coef[treat, 1] + model.m.coef[treatment.modifier, 1])) ##exp(beta3*(theta2+theta4) + theta6*(beta1+beta3))

  NIEH.gamma <- matrix(c(0, ##beat0
                         model.y.coef[modifier.mediator, 1], ##beta1
                         0, ##beta2
                         model.y.coef[mediator, 1] + model.y.coef[treatment.mediator, 1] + model.y.coef[modifier.mediator, 1], ##beta3
                         0, ##theta0
                         0, ##theta1
                         model.m.coef[treatment.modifier, 1], ##theta2
                         0, ##theta3
                         model.m.coef[treatment.modifier, 1], ##theta4
                         0, ##theta5
                         model.m.coef[treat, 1] + model.m.coef[treatment.modifier, 1],  ##theta6
                         0 ##sigma
  ), nrow=1)

  NIEH.sd <- sqrt(NIEH.gamma %*% cov %*% t(NIEH.gamma))

  NIEH.z <- log(NIEH)/NIEH.sd
  NIEH.p <- 2*pnorm(-abs(NIEH.z))

  NIEH.upper <- exp(log(NIEH) + qnorm(c.l) * (NIEH.sd/sqrt(n)))
  NIEH.lower <- exp(log(NIEH) - qnorm(c.l) * (NIEH.sd/sqrt(n)))

  NIEH.set <- c(sprintf("%10.5f", c(NIEH, NIEH.sd, NIEH.lower, NIEH.upper)), star.pval(round(NIEH.p, 3)))


  ######################PIEH: estimand and sd##############################
  PIEH <- exp(model.m.coef[treatment.modifier, 1] * model.y.coef[mediator, 1] + model.y.coef[modifier.mediator, 1] * (model.m.coef[treat, 1] + model.m.coef[treatment.modifier, 1])) ##exp(beta3*theta2 + theta6*(beta1+beta3))


  PIEH.gamma <- matrix(c(0, ##beat0
                         model.y.coef[modifier.mediator, 1], ##beta1
                         0, ##beta2
                         model.y.coef[mediator, 1] + model.y.coef[modifier.mediator, 1], ##beta3
                         0, ##theta0
                         0, ##theta1
                         model.m.coef[treatment.modifier, 1], ##theta2
                         0, ##theta3
                         0, ##theta4
                         0, ##theta5
                         model.m.coef[treat, 1] + model.m.coef[treatment.modifier, 1],  ##theta6
                         0 ##sigma
  ), nrow=1)

  PIEH.sd <- sqrt(PIEH.gamma %*% cov %*% t(PIEH.gamma))

  PIEH.z <- log(PIEH)/PIEH.sd
  PIEH.p <- 2*pnorm(-abs(PIEH.z))

  PIEH.upper <- exp(log(PIEH) + qnorm(c.l) * (PIEH.sd/sqrt(n)))
  PIEH.lower <- exp(log(PIEH) - qnorm(c.l) * (PIEH.sd/sqrt(n)))

  PIEH.set <- c(sprintf("%10.5f", c(PIEH, PIEH.sd, PIEH.lower, PIEH.upper)), star.pval(round(PIEH.p, 3)))


  ######################MIEH: estimand and sd##############################
  MIEH <- exp(model.y.coef[treatment.mediator, 1] * model.m.coef[treatment.modifier, 1])
  ####exp(theta4 * beta3)

  MIEH.gamma <- matrix(c(0, ##beta0
                         0, ##beta1
                         0, ##beta2
                         model.y.coef[treatment.mediator, 1], ##beta3
                         0, ##theta0
                         0, ##theta1
                         0, ##theta2
                         0, ##theta3
                         model.m.coef[treatment.modifier, 1], ##theta4
                         0, ##theta5
                         0, ##theta6
                         0 ##sigma
  ), nrow=1)

  MIEH.sd <- sqrt(MIEH.gamma %*% cov %*% t(MIEH.gamma))

  MIEH.z <- log(MIEH)/MIEH.sd
  MIEH.p <- 2*pnorm(-abs(MIEH.z))

  MIEH.upper <- exp(log(MIEH) + qnorm(c.l) * (MIEH.sd/sqrt(n)))
  MIEH.lower <- exp(log(MIEH) - qnorm(c.l) * (MIEH.sd/sqrt(n)))

  MIEH.set <- c(sprintf("%10.5f", c(MIEH, MIEH.sd, MIEH.lower, MIEH.upper)), star.pval(round(MIEH.p, 3)))

  ######################CDEH: estimand and sd##############################
  CDEH <- exp(model.y.coef[treatment.modifier, 1]) ##exp(theta5)


  CDEH.gamma <- matrix(c(0,
                         0,
                         0,
                         0,
                         0,
                         0,
                         0,
                         0,
                         0,
                         1,
                         0,
                         0), nrow=1)
  CDEH.sd <- sqrt(CDEH.gamma %*% cov %*% t(CDEH.gamma))

  CDEH.z <- log(CDEH)/CDEH.sd
  CDEH.p <- 2*pnorm(-abs(CDEH.z))

  CDEH.upper <- exp(log(CDEH) + qnorm(c.l) * (CDEH.sd/sqrt(n)))
  CDEH.lower <- exp(log(CDEH) - qnorm(c.l) * (CDEH.sd/sqrt(n)))

  CDEH.set <- c(sprintf("%10.5f", c(CDEH, CDEH.sd, CDEH.lower, CDEH.upper)), star.pval(round(CDEH.p, 3)))

  ######################TEH: estimand and sd##############################
  TEH <- NDEH * NIEH

  TEH.gamma <- matrix(c(0, ##beta0
                        model.y.coef[modifier.mediator, 1], ##beta1
                        model.y.coef[treatment.mediator, 1], ##beta2
                        model.y.coef[mediator, 1] + model.y.coef[treatment.mediator, 1] + model.y.coef[modifier.mediator, 1], ##beta3
                        0, ##theta0
                        0, ##theta1
                        model.m.coef[treatment.modifier, 1], ##theta2
                        0, ##theta3
                        model.m.coef[modifier, 1] + model.m.coef[treatment.modifier, 1] + model.y.coef[modifier.mediator, 1] * model.m.sigma.sq, ##theta4
                        1, ##theta5
                        model.y.coef[treatment.mediator, 1] * model.m.sigma.sq + model.m.coef[treat, 1] + model.m.coef[treatment.modifier, 1], ##theta6
                        model.y.coef[treatment.mediator, 1] * model.y.coef[modifier.mediator, 1] ##sigma
  ), nrow=1)

  TEH.sd <- sqrt(TEH.gamma %*% cov %*% t(TEH.gamma))

  TEH.z <- log(TEH)/TEH.sd
  TEH.p <- 2*pnorm(-abs(TEH.z))

  TEH.upper <- exp(log(TEH) + qnorm(c.l) * (TEH.sd/sqrt(n)))
  TEH.lower <- exp(log(TEH) - qnorm(c.l) * (TEH.sd/sqrt(n)))

  TEH.set <- c(sprintf("%10.5f", c(TEH, TEH.sd, TEH.lower, TEH.upper)), star.pval(round(TEH.p, 3)))

  out <- matrix(NA, nrow=6, ncol=5)
  row.names(out) <- c('NDEH', 'NIEH', 'PIEH', 'MIEH', 'CDEH', 'TEH')

  colnames(out) <- c('Estimate', 'Std. Error', paste(as.character(conf.level * 100), '% CI Lower', sep=''), paste(as.character(conf.level * 100), '% CI Upper', sep=''), 'p-value')


  out[1, ] <- NDEH.set
  out[2, ] <- NIEH.set
  out[3, ] <- PIEH.set
  out[4, ] <- MIEH.set
  out[5, ] <- CDEH.set
  out[6, ] <- TEH.set

  cat('Causal Heterogeneity Mediation Analysis')
  cat("\n")
  cat("\n")
  print(out, quote=FALSE)
  cat('-----------------------------------------------------------')
  cat('\n')
  cat('S.Es are for the logarithm of estimates')
  cat('\n')
  cat("\n")
  cat('Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1')
}

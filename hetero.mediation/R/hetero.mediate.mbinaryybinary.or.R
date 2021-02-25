#' Causal Mediation Analysis for Heterogeneous Effect (Binary Mediator and Binary Outcome, on the odds ratio scale)
#'
#' This function is for the estimation of various causal heterogeneous effects for causal mediation analysis, including natural direct and indirect effect heterogeneity, pure direct and indirect effect heterogeneity, mediated interactive effect heterogeneity, controlled direct effect heterogeneity, and total effect heterogeneity. Both the mediator and the outcome are binary, and estimates are on odds ratio scale.
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
hetero.mediate.mbinaryybinary <- function (model.m, model.y, modifier, treat, mediator, n, conf.level = 0.95){
  ##with the modifier and the treatment with two levels and coded as 0, 1. Then the constant of (a-a*)(q-q*) =1
  treatment.modifier <- paste(treat, ':', modifier, sep="")
  treatment.mediator <- paste(treat, ':', mediator, sep="")
  modifier.mediator <- paste(mediator, ':', modifier, sep="")

  model.m.coef<- summary(model.m)$coefficient
  model.y.coef <- summary(model.y)$coefficient

  model.m.cov <- vcov(model.m)
  model.y.cov <- vcov(model.y)

  ###variance-covariance matrix for all coefficients
  cov <- matrix(, nrow = dim(model.m.cov)[1] + dim(model.y.cov)[1], ncol = dim(model.m.cov)[1] + dim(model.y.cov)[1])
  cov[1:dim(model.m.cov)[1], 1:dim(model.m.cov)[1]] <- model.m.cov
  cov[(dim(model.m.cov)[1] + 1): nrow(cov), (dim(model.m.cov)[1] + 1):ncol(cov)] <- model.y.cov
  cov[is.na(cov)]=0

  c.l <- 0.5 + conf.level/2

  beta0 <- summary(model.m)$coefficients[1, 1]
  beta1 <- summary(model.m)$coefficients[treat, 1]
  beta2 <- summary(model.m)$coefficients[modifier, 1]
  beta3 <- summary(model.m)$coefficients[treatment.modifier, 1]
  theta0 <- summary(model.y)$coefficients[1, 1]
  theta1 <- summary(model.y)$coefficients[treat, 1]
  theta2 <- summary(model.y)$coefficients[mediator, 1]
  theta3 <- summary(model.y)$coefficients[modifier, 1]
  theta4 <- summary(model.y)$coefficients[treatment.mediator, 1]
  theta5 <- summary(model.y)$coefficients[treatment.modifier, 1]
  theta6 <- summary(model.y)$coefficients[modifier.mediator, 1]

  a1 <- exp(theta2 + theta4 + theta6 + beta0 + beta2) / (1 + exp(theta2 + theta4 + theta6 + beta0 + beta2))
  a2 <- exp(theta2 + beta0)/ (1 + exp(theta2 + beta0))
  a3 <- exp(theta2 + theta4 + beta0) / (1 + exp(theta2 + theta4 + beta0))
  a4 <- exp(theta2 + theta6 + beta0 + beta2) / (1 + exp(theta2 + theta6 + beta0 + beta2))
  a5 <- exp(theta2 + theta4 + theta6 + beta0 + beta1 + beta2 + beta3) / (1 + exp(theta2 + theta4 + theta6 + beta0 + beta1 + beta2 + beta3))
  a6 <- exp(beta0 + beta2) / (1 + exp(beta0 + beta2))
  a7 <- exp(beta0 + beta1) / (1 + exp(beta0 + beta1))
  a8 <- exp(beta0 + beta1 + beta2 + beta3) / (1 + exp(beta0 + beta1 + beta2 + beta3))
  a9 <- exp(theta2 + theta4  + beta0 + beta1) / (1 + exp(theta2 + theta4  + beta0 + beta1))
  a10 <- exp(beta0) / (1 + exp(beta0))
  a11 <- exp(theta2 + theta6 + beta0 + beta1 + beta2 + beta3) / (1 + exp(theta2 + theta6 + beta0 + beta1 + beta2 + beta3))
  a14 <- exp(theta2 + beta0 + beta1) / (1 + exp(theta2 + beta0 + beta1))


  ######################NDEH: estimand and sd##############################
  NDEH <- exp(theta5) * (1 + exp(theta2 + theta4 + theta6 + beta0 + beta2)) * (1 + exp(theta2 + beta0)) / ((1 + exp(theta2 + theta4 + beta0)) * (1 + exp(theta2 + theta6 + beta0 + beta2)))

  NDEH.gamma <- matrix(c(a1 + a2 - a3 - a4, ##beta0
                         0, ##beta1
                         a1 - a4, ##beta2
                         0, ##beta3
                         0, ##theta0
                         0, ##theta1
                         a1 + a2 - a3 - a4, ##theta2
                         0, ##theta3
                         a1 - a3, ##theta4
                         1, ##theta5
                         a1 - a4 ##theta6
  ), nrow=1)

  NDEH.sd <- sqrt(NDEH.gamma %*% cov %*% t(NDEH.gamma)) ##on the logarithm scale

  NDEH.z <- log(NDEH)/NDEH.sd
  NDEH.p <- 2*pnorm(-abs(NDEH.z))

  NDEH.upper <- exp(log(NDEH) + qnorm(c.l) * (NDEH.sd/sqrt(n)))
  NDEH.lower <- exp(log(NDEH) - qnorm(c.l) * (NDEH.sd/sqrt(n)))

  NDEH.set <- c(sprintf("%10.5f", c(NDEH, NDEH.sd, NDEH.lower, NDEH.upper)), star.pval(round(NDEH.p, 3)))


  ######################NIEH: estimand and sd##############################
  NIEH <- (1+exp(theta2 + theta4 + theta6 + beta0 + beta1 + beta2 + beta3)) * (1+exp(theta2 + theta4 + beta0)) * (1+exp(beta0 + beta1)) * (1+exp(beta0 + beta2)) / ((1+exp(theta2 + theta4 + theta6 + beta0 + beta2)) * (1+exp(theta2 + theta4 + beta0 + beta1)) * (1+exp(beta0)) * (1+exp(beta0 + beta1 + beta2 + beta3)))

  NIEH.gamma <- matrix(c(a5 + a6 + a3 + a7 - a1 - a8 - a9 - a10, ##beat0
                         a5 + a7 - a8 - a9, ##beta1
                         a5 + a6 - a1 - a8, ##beta2
                         a5 - a8, ##beta3
                         0, ##theta0
                         0, ##theta1
                         a5 + a3 - a1 - a9, ##theta2
                         0, ##theta3
                         a5 + a3 - a1 - a9, ##theta4
                         0, ##theta5
                         a5 - a1 ##theta6
  ), nrow=1)

  NIEH.sd <- sqrt(NIEH.gamma %*% cov %*% t(NIEH.gamma)) ##on the logarithm scale

  NIEH.z <- log(NIEH)/NIEH.sd
  NIEH.p <- 2*pnorm(-abs(NIEH.z))

  NIEH.upper <- exp(log(NIEH) + qnorm(c.l) * (NIEH.sd/sqrt(n)))
  NIEH.lower <- exp(log(NIEH) - qnorm(c.l) * (NIEH.sd/sqrt(n)))

  NIEH.set <- c(sprintf("%10.5f", c(NIEH, NIEH.sd, NIEH.lower, NIEH.upper)), star.pval(round(NIEH.p, 3)))

  ######################PIEH: estimand and sd##############################
  PIEH <- (1+exp(theta2 + theta6 + beta0 + beta1 + beta2 + beta3)) * (1+exp(theta2 + beta0)) * (1+exp(beta0 + beta1)) * (1+exp(beta0 + beta2)) / ((1+exp(theta2 + theta6 + beta0 + beta2)) * (1+exp(theta2 + beta0 + beta1)) * (1+exp(beta0)) * (1+exp(beta0 + beta1 + beta2 + beta3)))

  PIEH.gamma <- matrix(c(a11 + a6 + a2 + a7 - a4 - a8 - a14 - a10, ##beat0
                         a11 + a7 - a8 - a14, ##beta1
                         a11 + a6 - a4 - a8, ##beta2
                         a11 - a8, ##beta3
                         0, ##theta0
                         0, ##theta1
                         a11 + a2 - a4 - a14, ##theta2
                         0, ##theta3
                         0, ##theta4
                         0, ##theta5
                         a11- a4 ##theta6
  ), nrow=1)

  PIEH.sd <- sqrt(PIEH.gamma %*% cov %*% t(PIEH.gamma)) ##on the logarithm scale

  PIEH.z <- log(PIEH)/PIEH.sd
  PIEH.p <- 2*pnorm(-abs(PIEH.z))

  PIEH.upper <- exp(log(PIEH) + qnorm(c.l) * (PIEH.sd/sqrt(n)))
  PIEH.lower <- exp(log(PIEH) - qnorm(c.l) * (PIEH.sd/sqrt(n)))

  PIEH.set <- c(sprintf("%10.5f", c(PIEH, PIEH.sd, PIEH.lower, PIEH.upper)), star.pval(round(PIEH.p, 3)))

  ######################MIEH: estimand and sd##############################
  MIEH <- (1+exp(theta2 + theta4 + theta6 + beta0 + beta1 + beta2 + beta3)) * (1+exp(theta2 + theta6 + beta0 + beta2)) * (1+exp(theta2 + theta4 + beta0)) * (1+exp(theta2 + beta0 + beta1)) / ((1+exp(theta2 + theta4 + theta6 + beta0 + beta2)) * (1+exp(theta2 + theta6 + beta0 + beta1 + beta2 + beta3)) * (1+exp(theta2 + theta4 + beta0 + beta1)) * (1+exp(theta2 + beta0)))

  MIEH.gamma <- matrix(c(a5 + a4 + a3 + a14 - a1 - a11 - a9 - a2, ##beta0
                         a5 + a14 - a11 - a9, ##beta1
                         a5 + a4 - a1 - a11, ##beta2
                         a5 - a11, ##beta3
                         0, ##theta0
                         0, ##theta1
                         a5 + a4 + a3 + a14 - a1 - a11 - a9 - a2, ##theta2
                         0, ##theta3
                         a5 + a2 - a1 - a9, ##theta4
                         0, ##theta5
                         a5 + a4 - a1 - a11 ## theta6
  ), nrow=1)

  MIEH.sd <- sqrt(MIEH.gamma %*% cov %*% t(MIEH.gamma)) ##on the logarithm scale

  MIEH.z <- log(MIEH)/MIEH.sd
  MIEH.p <- 2*pnorm(-abs(MIEH.z))

  MIEH.upper <- exp(log(MIEH) + qnorm(c.l) * (MIEH.sd/sqrt(n)))
  MIEH.lower <- exp(log(MIEH) - qnorm(c.l) * (MIEH.sd/sqrt(n)))

  MIEH.set <- c(sprintf("%10.5f", c(MIEH, MIEH.sd, MIEH.lower, MIEH.upper)), star.pval(round(MIEH.p, 3)))

  ######################CDEH: estimand and sd##############################
  CDEH <- exp(theta5) ##theta5

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
                         0), nrow=1)

  CDEH.sd <- sqrt(CDEH.gamma %*% cov %*% t(CDEH.gamma)) ##on the logarithm scale

  CDEH.z <- log(CDEH)/CDEH.sd
  CDEH.p <- 2*pnorm(-abs(CDEH.z))

  CDEH.upper <- exp(log(CDEH) + qnorm(c.l) * (CDEH.sd/sqrt(n)))
  CDEH.lower <- exp(log(CDEH) - qnorm(c.l) * (CDEH.sd/sqrt(n)))

  CDEH.set <- c(sprintf("%10.5f", c(CDEH, CDEH.sd, CDEH.lower, CDEH.upper)), star.pval(round(CDEH.p, 3)))

  ######################TEH: estimand and sd##############################
  TEH <- NDEH * NIEH

  TEH.gamma <- NDEH.gamma + NIEH.gamma

  TEH.sd <- sqrt(TEH.gamma %*% cov %*% t(TEH.gamma)) ##on the logarithm scale

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

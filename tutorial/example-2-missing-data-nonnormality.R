

# -------------------------------------------------------------------------
#
# Example 2: Bifactor Model with Non-Normal Indicators and Missing Data
#
# -------------------------------------------------------------------------

library(simsem)

source('summaryECV.R')


## define population model

populationModel <- '

  # define loadings
  G =~ .6*y1 + .6*y2 + .6*y3 + .6*y4 + .6*y5 + .6*y6 + .6*y7 + .6*y8 + .6*y9
  S1 =~ .4*y1 + .4*y2 + .4*y3
  S2 =~ .4*y4 + .4*y5 + .4*y6
  S3 =~ .4*y7 + .4*y8 + .4*y9
  
  # define residual variances of indicators (= 1 - .6^2 - .4^2)
  y1 ~~ .48*y1
  y2 ~~ .48*y2
  y3 ~~ .48*y3
  y4 ~~ .48*y4
  y5 ~~ .48*y5
  y6 ~~ .48*y6
  y7 ~~ .48*y7
  y8 ~~ .48*y8
  y9 ~~ .48*y9
  
  # set factor variances to 1 (so that loadings are in a standardized metric)
  G ~~ 1*G
  S1 ~~ 1*S1
  S2 ~~ 1*S2
  S3 ~~ 1*S3
  
  # define factors to be orthogonal
  S1 + S2 + S3 ~~ 0*G
  S1 + S2 ~~ 0*S3
  S1 ~~ 0*S2
'


## define analysis model

analysisModel <- '
  G =~ NA*y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9
  S1 =~ NA*y1 + y2 + y3
  S2 =~ NA*y4 + y5 + y6
  S3 =~ NA*y7 + y8 + y9
  
  S1 + S2 + S3 ~~ 0*G
  S1 + S2 ~~ 0*S3
  S1 ~~ 0*S2
  
  G ~~ 1*G
  S1 ~~ 1*S1
  S2 ~~ 1*S2
  S3 ~~ 1*S3
'

## define amount of missing data
missing <- miss(pmMCAR = 0.1, m = 0)

## define marginal distribution
distribution <- bindDist(skewness = 0.5, kurtosis = 4, p = 9)




## do simulation

## N = 300
results <- sim(n = 300, generate=populationModel, model = analysisModel, lavaanfun = "sem", nRep = 500, completeRep = 5000, miss = missing, indDist = distribution, estimator = 'MLR', seed = 2021)


# convergence
summaryConverge(results)

## parameter and standard error estimates
param <- summaryParam(results, improper = FALSE)
param <- param[1:18, ] # only loadings
param$`Average Param` <- t(results@paramValue[1:18])
param$`Average Bias` <- param$`Estimate Average` - param$`Average Param`
param$`Rel Bias` <- (param$`Estimate Average` - param$`Average Param`)/param$`Average Param`
param$`Rel SE Bias` <- (param$`Average SE` - param$`Estimate SD`)/param$`Estimate SD`
names(param) <- c('Estimate Average', 'Estimate SD', 'Average SE', 'Power (Not Equal 0)', 'Std Est', 'Std Est SD', 'Std Ave SE', 'Average Param', 'Average Bias', 'Rel Bias', 'Rel SE Bias')
param

## estimated ECV
summaryECV(results, improper = FALSE) 




## N = 350
results.N350 <- sim(n = 350, generate=populationModel, model = analysisModel, lavaanfun = "sem", nRep = 500, completeRep = 5000, miss = missing, indDist = distribution, estimator = 'MLR', seed = 2021)


# convergence
summaryConverge(results.N350)

## parameter and standard error estimates
param <- summaryParam(results.N350, improper = FALSE)
param <- param[1:18, ] # only loadings and residual variances
param$`Average Param` <- t(results.N350@paramValue[1:18])
param$`Average Bias` <- param$`Estimate Average` - param$`Average Param`
param$`Rel Bias` <- (param$`Estimate Average` - param$`Average Param`)/param$`Average Param`
param$`Rel SE Bias` <- (param$`Average SE` - param$`Estimate SD`)/param$`Estimate SD`
names(param) <- c('Estimate Average', 'Estimate SD', 'Average SE', 'Power (Not Equal 0)', 'Std Est', 'Std Est SD', 'Std Ave SE', 'Average Param', 'Average Bias', 'Rel Bias', 'Rel SE Bias')
param

## estimated ECV
summaryECV(results.N350, improper = FALSE) 

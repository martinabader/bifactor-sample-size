

# -------------------------------------------------------------------------
#
# Example 3: Bifactor Model with Ordinal Indicators
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
  
  # set factor variances to 1 (so that loadinGS are in a standardized metric)
  G ~~ 1*G
  S1 ~~ 1*S1
  S2 ~~ 1*S2
  S3 ~~ 1*S3
  
  # define factors to be orthoGonal
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


## function to discretize data into 5 categories

discretize <- function(data){
  data[data > 1.5] <- 5
  data[data > 0.5 & data <= 1.5] <- 4
  data[data > -0.5 & data <= 0.5] <- 3
  data[data > -1.5 & data <= -0.5] <- 2
  data[data <= -1.5] <- 1
  data
}


## do simulation

# N = 300
results <- sim(n = 300, generate = populationModel, model = analysisModel, lavaanfun = "sem", nRep = 500, completeRep = 5000, datafun = discretize, estimator = 'WLSMV', seed = 2021) 

summaryConverge(results)
summaryParam(results, improper = FALSE, detail = TRUE)
summaryECV(results, improper = FALSE)


# N = 350
results.N350 <- sim(n = 350, generate = populationModel, model = analysisModel, lavaanfun = "sem", nRep = 500, completeRep = 5000, datafun = discretize, estimator = 'WLSMV', seed = 2021) 

summaryConverge(results.N350)
summaryParam(results.N350, improper = FALSE, detail = TRUE)
summaryECV(results.N350, improper = FALSE)


# N = 400
results.N400 <- sim(n = 400, generate = populationModel, model = analysisModel, lavaanfun = "sem", nRep = 500, completeRep = 5000, datafun = discretize, estimator = 'WLSMV', seed = 2021) 

summaryConverge(results.N400)
summaryParam(results.N400, improper = FALSE, detail = TRUE)
summaryECV(results.N400, improper = FALSE)




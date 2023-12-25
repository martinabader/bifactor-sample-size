
# -------------------------------------------------------------------------
#
# aggregate simulation results
#
# -------------------------------------------------------------------------


options(scipen = 999)
options(max.print = 999999)

library(psych)
library(data.table)


# get sim results ---------------------------------------------------------

simresults <- as.data.frame(fread("simresults.csv", sep=";"))

errors <- read.csv("convergence.csv", sep=";")



# conditions --------------------------------------------------------------

proportionality <- c(1, 0) # yes, no
meanLoadg <- c(.3, .5, .7) # mean loading on g
meanECV <- c(.25, .50, .75) # mean explained common variance
numIndicator <- c(3, 5, 7, 10) # number of indicators per specific factor
numSpecFac <- c(3, 5, 7, 10) # number of specific factors
Ncon <- c(150, 300, 500, 1000, 2000) # sample size conditions



# outcomes ----------------------------------------------------------------

# 1) convergence rates
# 2) bias in loadings on g factor
# 3) bias in loadings on specific factors
# 4) bias in ecv
# 5) bias in standard error of loadings on g factor
# 6) bias in standard errors of loadings on specific factors 




# results object ----------------------------------------------------------



outcomes <- simresults[, c('pop.lambda.g', 'pop.ecv', 'n.ind', 'n.spec', 'N', 'proportionality')]

outcomes[, c('convergence.rate',
             'bias.lambda.g', 'bias.lambda.s',
             'bias.ecv',
             'bias.se.lambda.g', 'bias.se.lambda.s')] <- NA



# loop through conditions -------------------------------------------------


row <- 1


for( a in 1:length(proportionality)){
  
  proportional <- proportionality[a]

for( b in 1:length(meanLoadg)){
  
  mean.gen <- meanLoadg[b]

  # depending on the mean loading on g, there are only some feasible ecvs
  if(mean.gen == .7){meanECV <- c(.75)
  } else if(mean.gen == .5){meanECV <- c(.50, .75) 
  } else {meanECV <- c(.25, .50, .75)}
  
  
  for( c in 1:length(meanECV)){
  ecv <- meanECV[c]

    for( d in 1:length(numIndicator)){
    n.ind <- numIndicator[d]

      for( e in 1:length(numSpecFac)){
      n.spec <- numSpecFac[e] 
        
        for(f in 1:length(Ncon)){
        N <- Ncon[f]
        
        p <- n.ind * n.spec
        
          
        # get results for current condition
        
        results <- simresults[simresults$pop.lambda.g == mean.gen & 
                              simresults$pop.ecv == ecv &
                              simresults$n.ind == n.ind &
                              simresults$n.spec == n.spec & 
                              simresults$N == N &
                              simresults$proportionality == proportional, ]
        
        
        c.errors <- errors[errors$pop.lambda.g == mean.gen &
                             errors$pop.ecv == ecv &
                             errors$n.ind == n.ind &
                             errors$n.spec == n.spec &
                             errors$N == N &
                             errors$proportionality == proportional, ]

        
        # 1) convergence rates ----------------------------------------------------
        
        convergence.rate <- c.errors$convergence.rate
        
        convergence.rate.rep <- rep(convergence.rate, nrow(results))
        
        if(c.errors$num.converged >= 1){
        
          
        # 2) bias in loadings on g factor -----------------------------------------
        
        # pop loadings on g
        pop.lambda.g <- results[, grep('pop.gload_x', names(results))]

        # std loadings on g
        std.lambda.g <- results[, grep('std.gload_x', names(results))]

        # relative percentage bias
        bias.lambda.g.spar <- ((std.lambda.g-pop.lambda.g)/pop.lambda.g)
        bias.lambda.g <- rowMeans(as.matrix(bias.lambda.g.spar), na.rm = T)
        
        
        # 3) bias in loadings on specific factors ---------------------------------
        
        # pop loadings on spec
        pop.lambda.s <- results[, grep('pop.sload_x', names(results))]

        # std loadings on spec
        std.lambda.s <- results[, grep('std.sload_x', names(results))]

        # bias
        bias.lambda.s.spar <- ((std.lambda.s-pop.lambda.s)/pop.lambda.s)
        bias.lambda.s <- rowMeans(as.matrix(bias.lambda.s.spar), na.rm = T)
        
        
        # 4) bias in ecv ----------------------------------------------------------
        
        # pop ecv
        pop.ecv <- results$pop.ecv_g
        
        # estimated ecv
        est.ecv <- results$est.ecv
        
        # bias
        bias.ecv <- ((est.ecv - pop.ecv)/pop.ecv)
        
        
        
        # 5) bias in standard error of loadings on g factor -----------------------
        
        # standard deviation of est. parameters (population) 
        est.lambda.g <- results[, grep('est.gload_x', names(results))]
        pop.se.lambda.g <- apply(est.lambda.g, 2, sd)
        pop.se.lambda.g <- matrix(rep(pop.se.lambda.g, nrow(results)), nrow = nrow(results), byrow = TRUE)
        
        # estimated standard errors
        est.se.lambda.g <- results[, grep('est.se.gload_x', names(results))]

        # bias
        bias.se.lambda.g.spar <- ((est.se.lambda.g - pop.se.lambda.g)/pop.se.lambda.g)
        bias.se.lambda.g <- rowMeans(as.matrix(bias.se.lambda.g.spar), na.rm = T)
        
        
        
        # 6) bias in standard errors of loadings on specific factors --------------
        
        # standard deviation of est. parameters (population)
        est.lambda.s <- results[, grep('est.sload_x', names(results))]
        pop.se.lambda.s <- apply(est.lambda.s, 2, sd)
        pop.se.lambda.s <- matrix(rep(pop.se.lambda.s, nrow(results)), nrow = nrow(results), byrow = TRUE)
        
        # estimated standard errors
        est.se.lambda.s <- results[, grep('est.se.sload_x', names(results))]

        # bias
        bias.se.lambda.s.spar <- ((est.se.lambda.s - pop.se.lambda.s)/pop.se.lambda.s)
        bias.se.lambda.s <- rowMeans(as.matrix(bias.se.lambda.s.spar), na.rm = T)
        
    

        # save results --------------------------------------------------------

        outcomes[outcomes$pop.lambda.g == mean.gen &
                  outcomes$pop.ecv == ecv &
                  outcomes$n.ind == n.ind &
                  outcomes$n.spec == n.spec &
                  outcomes$N == N &
                  outcomes$proportionality == proportional,
                 c('convergence.rate',
                   'bias.lambda.g', 'bias.lambda.s',
                   'bias.ecv',
                   'bias.se.lambda.g', 'bias.se.lambda.s')] <- cbind(convergence.rate.rep, 
                                              bias.lambda.g, bias.lambda.s,
                                              bias.ecv,
                                              bias.se.lambda.g, bias.se.lambda.s)
        

        
        } else {
          
          outcomes[outcomes$pop.lambda.g == mean.gen &
                     outcomes$pop.ecv == ecv &
                     outcomes$n.ind == n.ind &
                     outcomes$n.spec == n.spec &
                     outcomes$N == N &
                     outcomes$proportionality == proportional,
                   'convergence.rate'] <- convergence.rate
          
        }
        
        
        
        
        

        row <- row + 1


        } # end f loop
      } # end e loop
    } # end d loop
  } # end c loop
} # end b loop
} # end a loop


write.table(outcomes, 'outcomes.csv', sep = ';', row.names = FALSE)

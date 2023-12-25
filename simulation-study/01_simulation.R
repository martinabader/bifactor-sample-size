
# -------------------------------------------------------------------------
#
# simulation study using simsem
#
# -------------------------------------------------------------------------


library(lavaan)
library(simsem)


# conditions --------------------------------------------------------------

proportionality <- c(0, 1) # validity of proportionality condition

meanLoadg <- c(.3, .5, .7) # mean loading on g

meanECV <- c(.25, .50, .75) # mean explained common variance

numIndicator <- c(3, 5, 7, 10) # number of indicators per specific factor

numSpecFac <- c(3, 5, 7, 10) # number of specific factors

Ncon <- c(150, 300, 500, 1000, 2000) # sample size conditions

numRep <- 500

maxReps <- 5000 # maximum number of repetitions




# container ---------------------------------------------------------------

## results

results.colnames <- c('pop.lambda.g', 'pop.ecv', 'n.ind', 'n.spec', 'N', 'proportionality',
                           paste0('pop.gload_x', seq(1, (max(numIndicator)*max(numSpecFac)))),
                           paste0('pop.sload_x', seq(1, (max(numIndicator)*max(numSpecFac)))),
                           paste0('est.gload_x', seq(1, (max(numIndicator)*max(numSpecFac)))),
                           paste0('std.gload_x', seq(1, (max(numIndicator)*max(numSpecFac)))),
                           paste0('est.sload_x', seq(1, (max(numIndicator)*max(numSpecFac)))),
                           paste0('std.sload_x', seq(1, (max(numIndicator)*max(numSpecFac)))),
                           paste0('est.se.gload_x', seq(1, (max(numIndicator)*max(numSpecFac)))),
                           paste0('est.se.sload_x', seq(1, (max(numIndicator)*max(numSpecFac)))),
                           'pop.ecv_g',
                           paste0('pop.ecv_f', seq(1, max(numSpecFac))),
                           'est.ecv')


results <- as.data.frame(matrix(NA, 
                                   ncol = length(results.colnames),
                                   nrow = length(meanLoadg)*length(meanECV)*length(numIndicator)*length(numSpecFac)*length(Ncon)*length(proportionality)*numRep))

colnames(results) <- results.colnames

intermediate.storage <- rbind(results.colnames, results.colnames)

write.table(intermediate.storage, 'intermediate.storage.csv', sep = ';', col.names = FALSE, quote = FALSE, row.names = FALSE)



## convergence & errors

error.colnames <- c('condition','pop.lambda.g', 'pop.ecv', 'n.ind', 'n.spec', 'N', 'proportionality',
                    'num.converged', 
                    'num.nonconverged', 
                    'total.reps',
                    'convergence.rate', 
                    'nonconvergent', 
                    'improperSE', 
                    'nonpositivedefinite', 
                    'nooptimalestimates')

errors <- data.frame(matrix(NA, 
                               ncol = length(error.colnames),
                               nrow = 960))

colnames(errors) <- error.colnames





# simulation --------------------------------------------------------------

cnt.rep.all <- 1

cnt.cond <- 1


for( a in 1:length(proportionality) ){
  
  proportional <- proportionality[a]


  for( b in 1:length(meanLoadg) ){
    
    mean.gen <- meanLoadg[b]
    
    # depending on the mean loading on g, there are only some feasible ecvs
    if(mean.gen == .7){meanECV <- c(.75)
    } else if(mean.gen == .5){meanECV <- c(.50, .75)
    } else {meanECV <- c(.25, .50, .75)}
    
    
    for( c in 1:length(meanECV) ){
      
      ecv <- meanECV[c]
      
      for( d in 1:length(numIndicator) ){
        
        n.ind <- numIndicator[d]
      
        for( e in 1:length(numSpecFac) ){
          
          n.spec <- numSpecFac[e] 
      
  

          # set model parameters ----------------------------------------------------

          p <- n.ind * n.spec 
  
          fac.ind <- matrix(1:p, ncol = n.spec, nrow = n.ind)
  
          
          ## general factor loadings
          
          set.seed(472)
          load.g <- round(runif(p, min = mean.gen-0.1, max = mean.gen+0.1), 2)
          
          
          ## specific factor loadings (defined depending on validity of proportionality condition)
          
          if(proportional == 1){
            
            denum.g.f <- sapply(1:n.spec, function(x) sum(load.g[fac.ind[,x]]^2) )
            
            if(n.spec == 3){
              ecv.s.w <- c(.40, .33, .27)
            } else if (n.spec == 5){
              ecv.s.w <- c(.25, .16, .20, .21, .18)
            } else if (n.spec == 7){
              ecv.s.w <- c(.14, .175, .115, .15, .12, .16, .14)
            } else if (n.spec == 10){
              ecv.s.w <- c(.12, .101, .08, .099, .105, .095, .09, .1, .11, .1)
            }
            
            denum.s.f <- ecv.s.w * sum(load.g^2*(1-ecv)/ecv)
            
            ecv.s <- denum.g.f / (denum.g.f + denum.s.f) # ecvs for subscales
            
            load.s <- c( sapply(1:n.spec, function(x) sqrt(load.g[fac.ind[,x]]^2 * (1-ecv.s[x]) / ecv.s[x])) )
            
          } else{
            
            mean.spe <- sqrt(mean.gen^2*(1-ecv)/ecv)
            
            set.seed(2)
            load.s <- round(runif(p, min = mean.spe-0.1, max = mean.spe+0.1), 2)
            
          }
          

          ## loadings
  
          lambda <- matrix(0, ncol = (1+n.spec), nrow = p)
          lambda[ ,1] <- load.g
          for(ii in 1:n.spec){lambda[fac.ind[,ii], (1+ii)] <- load.s[fac.ind[,ii]]}
          if(any(rowSums(lambda^2) >= 1)) stop('cor too large')
          
          lambda.var <- lambda
          lambda.var[lambda != 0] <- NA
      
          LY <- bind(lambda.var, lambda)
          
          
          ## factor correlations
          phi <- diag(nrow = n.spec+1, ncol = n.spec+1)
          RPS <- binds(phi, 0)
          
          
          ## residual correlations
          theta <- diag(1, nrow=p, ncol=p)
          RTE <- binds(theta,0)
          
         
          # sanity check ecv total
          pop.ecv.g <- sum(lambda[,1]^2) / (sum(lambda[,1]^2)  + sum(lambda[,-1]^2) ) 
          
          # sanity check ecv subscales
          pop.ecv.s <- sapply(1:n.spec, function(x) # ecv subscale
            sum(lambda[fac.ind[,x],1]^2) / (sum(lambda[fac.ind[,x],1]^2)  + sum(lambda[fac.ind[,x],(x+1)]^2) )
          )
          
          # save lambda
          c.filename.l <- paste0('lambda_sf', n.spec, '_i', n.ind, '_g', mean.gen, '_ecv', ecv, '.csv')
          write.table(lambda, c.filename.l, sep = ';', row.names = FALSE, col.names = FALSE)
          
          

          
          # define model ------------------------------------------------------------

          mod <- model.cfa(LY = LY, RPS = RPS, RTE = RTE, indLab = paste0('x',1:p), facLab = c('g',paste0('s',1:n.spec)))
          

          
          
            # define sample size -------------------------------------------
            
            for(f in 1:length(Ncon)){
              
              N <- Ncon[f]
              
            
            
              # run simulation ---------------------------------------------
              
              print(paste('condition ', cnt.cond,
                          '| proportionality ', a, '/', length(proportionality),                          
                          '| gload ', b, '/', length(meanLoadg),
                          '| ecv ', c, '/', length(meanECV), 
                          '| ind ', d, '/', length(numIndicator), 
                          '| spec ', e, '/', length(numSpecFac),
                          '| N', f, '/', length(Ncon)))
              
            
              output <- sim(model = mod, n = N, nRep = numRep, completeRep = maxReps,
                            multicore = TRUE, numProc = 46, seed = 2021,
                            std.lv = TRUE, orthogonal = T)
             

              # convergence rate -------------------------------------------------------

              if(!is.null(output@converged) & !any(output@converged !=0)){
                
                counts <- sum(output@converged == 0)
                reps <- length(output@converged)
                rate <- counts/reps
                reasons <- rep(0,4)
                
                c.errors <- c(cnt.cond, mean.gen, ecv, n.ind, n.spec, N, proportional,
                              counts, 0, reps, rate, reasons)
                
                errors[cnt.cond,] <- c.errors
                
                
              } else if(any(output@converged !=0)){
                
                counts <- summaryConverge(output)$Converged
                reps <- sum(counts['num.nonconverged'], counts['num.converged'])
                rate <- counts['num.converged']/reps
                reasons <- summaryConverge(output)[['Nonconvergent Reasons']][,1]
                
                c.errors <- c(cnt.cond, mean.gen, ecv, n.ind, n.spec, N, proportional,
                              counts, reps, rate, reasons)
                
                errors[cnt.cond,] <- c.errors            
              }
              

              

              # parameters --------------------------------------------------------------

              param <- output@coef[output@converged == 0, ]
              if(nrow(param) > numRep){param <- param[1:numRep, ]}
              
              stdParam <- output@stdCoef[output@converged == 0, ]
              if(nrow(stdParam) > numRep){stdParam <- stdParam[1:numRep, ]}
              
              creps <- nrow(param)
              
              if(creps > 0){
                
                ## est loadings on g
                est.lambda.g <- param[, grep('g=~', names(param))]
                
                ## std loadings on g
                std.lambda.g <- stdParam[, grep('g=~', names(stdParam))]
                
                ## est loadings on specifics
                est.lambda.s <- param[, grep('s[[:digit:]]=~', names(param))]
                
                ## std loadings on specifics
                std.lambda.s <- stdParam[, grep('s[[:digit:]]=~', names(stdParam))]              
                
                ## implied ecv
                est.ecv <- rowSums(std.lambda.g^2)/( rowSums(std.lambda.g^2) +  rowSums(std.lambda.s^2))
                
                
                
                # standard errors ---------------------------------------------------------
                
                stdErrors <- output@se[output@converged == 0, ]
                if(nrow(stdErrors) > numRep){stdErrors <- stdErrors[1:numRep, ]}
                
                ## SE g loadings
                se.lambda.g <- stdErrors[, grep('g=~', names(stdErrors))]
                
                ## SE s loadings
                se.lambda.s <- stdErrors[, grep('s[[:digit:]]=~', names(stdErrors))]              
                
                
                
                # save results ------------------------------------------------------------
                
                
                ## condition
                results[cnt.rep.all:(cnt.rep.all+creps-1), c('pop.lambda.g', 'pop.ecv', 'n.ind', 'n.spec', 'N', 'proportionality')] <- matrix(rep(c(mean.gen, ecv, n.ind, n.spec, N, proportional), creps), nrow=creps, byrow = TRUE)
                
                ## pop loadings on g
                p1 <- which(colnames(results)=="pop.gload_x1")
                results[cnt.rep.all:(cnt.rep.all+creps-1), p1:(p1+length(load.g)-1)] <- matrix(rep(load.g, creps), nrow=creps, byrow=T)
                
                ## pop loadings on spec
                p1 <- which(colnames(results)=="pop.sload_x1")
                results[cnt.rep.all:(cnt.rep.all+creps-1), p1:(p1+length(load.s)-1)] <- matrix(rep(load.s, creps), nrow=creps, byrow=T)
                
                ## est loadings on g
                p1 <- which(colnames(results)=="est.gload_x1")
                results[cnt.rep.all:(cnt.rep.all+creps-1), p1:(p1+length(est.lambda.g)-1)] <- est.lambda.g
                
                ## std loadings on g
                p1 <- which(colnames(results)=="std.gload_x1")
                results[cnt.rep.all:(cnt.rep.all+creps-1), p1:(p1+length(std.lambda.g)-1)] <- std.lambda.g              
                
                ## est loadings on spec
                p1 <- which(colnames(results)=="est.sload_x1")
                results[cnt.rep.all:(cnt.rep.all+creps-1), p1:(p1+length(est.lambda.s)-1)] <- est.lambda.s
                
                ## std loadings on spec
                p1 <- which(colnames(results)=="std.sload_x1")
                results[cnt.rep.all:(cnt.rep.all+creps-1), p1:(p1+length(std.lambda.s)-1)] <- std.lambda.s          
                
                ## est se of loading on g
                p1 <- which(colnames(results)=="est.se.gload_x1")
                results[cnt.rep.all:(cnt.rep.all+creps-1), p1:(p1+length(se.lambda.g)-1)] <- se.lambda.g
                
                ## est se of loading on spec
                p1 <- which(colnames(results)=="est.se.sload_x1")
                results[cnt.rep.all:(cnt.rep.all+creps-1), p1:(p1+length(se.lambda.s)-1)] <- se.lambda.s
                
                # pop ecv of g
                results[cnt.rep.all:(cnt.rep.all+creps-1), 'pop.ecv_g'] <- pop.ecv.g
                
                # pop ecv of specific factors
                p1 <- which(colnames(results)=="pop.ecv_f1")
                results[cnt.rep.all:(cnt.rep.all+creps-1), p1:(p1+length(pop.ecv.s)-1)] <- matrix(rep(pop.ecv.s, creps), nrow=creps, byrow = T)
                
                ## est ecv
                results[cnt.rep.all:(cnt.rep.all+creps-1), 'est.ecv'] <- est.ecv
              
            
              }
            
                
              ## save to disk
              write.table(errors, 'convergence.csv', sep = ';', row.names = FALSE)
              
              if(creps > 0){
                c.res <- results[cnt.rep.all:(cnt.rep.all+creps-1), ]
                write.table(c.res, file = 'intermediate.storage.csv', sep = ";", col.names = FALSE, quote = FALSE, row.names = FALSE, append = T)                
              }

                
  
              # increase counter --------------------------------------------------------
  
              cnt.rep.all <- cnt.rep.all + creps
              cnt.cond <- cnt.cond + 1
                

          } # end f loop
        } # end e loop
      } # end d loop
    } # end c loop
  } # end b loop
} # end a loop


write.table(results, file = 'simresults.csv', sep = ";", row.names = FALSE)



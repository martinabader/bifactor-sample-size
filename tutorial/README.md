# How to conduct a Monte Carlo simulation to determine the sample size requirements for bifactor models?

This tutorial is based on the paper:

Bader, M., Jobst, L. J., & Moshagen, M. (2022). Sample size requirements for bifactor models. *Structural Equation Modeling: A Multidisciplinary Journal, 29*(5), 772–783. [https://doi.org/10.1080/10705511.2021.2019587]

Please cite this paper if you follow this tutorial.


## Preliminaries

This tutorial relies on the R packages `simsem` and `lavaan`. 

Before continuing, make sure that these packages have been installed and activated with the command `library(simsem)`.

In addition, a costum function labeled `summaryECV` is used which is not part of the `simsem` package. Please download this function from the tutorial folder and load the function into your R project using the command `source('summaryECV.R')`


## Model Specification

### Defining the population model

First, you need to define the population model for which you would like to determine the required sample size. The population model will be used by `simsem` for data generation. Our example model is a bifactor model with one general factor and three specific factors with three indicators each. You also need to specify the values of all parameters of the population model by relying on theoretical assumptions or prior research findings. In our example, all general factors loadings are .60 and all specific factor loadings are .40. The residual variances of each indicator can be found by subtracting the squared factor loadings from 1. The variances of the latent factors were set to 1 so that all parameters are in a standardized metric. In line with the standard specification of bifactor models, orthogonality between the specific factors and the general factor as well as among the specific factors is assumed. The model is specified via the syntax style of lavaan as follows: 

``` r
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
```


### Defining the analysis model

Next, we define the analysis model that will be used by simsem to analyse the generated data. The analysis model defines the general model structure but does not impose specific parameter values. Our analysis model correctly represents the structure of the population model

``` r
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
```

## Run Simulation

Now, we can run a Monte Carlo simulation by using the `sim` command:

``` r
results <- sim(n = 150, generate = populationModel, model = analysisModel, lavaanfun = "sem", nRep = 500, completeRep = 5000, seed = 2021)
```

In the `n` argument, we specify that a sample size of 150 observations should be used. The population model and analysis model are supplied to the `generate` and `model` arguments, respectively. The `lavaanfun` defines the name of the function that is used for data analysis (here, the `sem` function of lavaan). The number of simulation replications (`nRep`) is set to 500, which tells `simsem` that we aim to have a minimum of 500 properly converged solutions. This number can be increased, but I would recommend to not use fewer than 500 simulation replications. We also set that the maximum number of simulation replications (`completeRep`) to 5,000, meaning that the simulation should stop after 5000 replications even if our target number of 500 proper solutions has not been reached. The `seed` function can be used for reproducibility. 


## Check Results

### Convergence

First, we check the convergence of our model, that is, the proportion of simulation replications in which our model converged properly. The convergence rate is considered acceptable when being .90 or higher.

``` r
summaryConverge(results) 
```
This gives us the following output:

``` r
$Converged
   num.converged num.nonconverged 
             511              147 

$`Nonconvergent Reasons`
                                      count
Nonconvergent                            22
Improper SE                               0
Nonpositive definite matrix             125
Optimal estimates were not guaranteed     0

```
The upper part of the output summarizes that 511 out of 658 replications converged properly. This amounts to a convergence rate of .78 which is unacceptable. The reasons for nonconvergence were: nonconvergence of the estimation function (in 22 replications) and non-positive definiteness of latent or manifest (residual) covariance matrices (in 125 replications). Note that in the case of bifactor models with orthogonal factors, only the manifest residual covariances can be non-positive definite.

### Parameter and standard error estimates

To check how accurately parameters and standard errors were estimated, we use the following command (indicating that we wish to obtained detailed summary statistics while excluding improper solutions):

``` r
summaryParam(results, detail = TRUE, improper = FALSE) 
```

Parts of the output are shown below:

``` r

|            | Estimate Average | Estimate SD | Average SE | Power (Not equal 0) | Average Bias | Coverage | Rel Bias | Rel SE Bias |
|------------|------------------|-------------|------------|---------------------|--------------|----------|----------|-------------|
| G=~y1      | 0.59             | 0.09        | 0.09       | 1.00                | -0.01        | 0.93     | -0.01    | -0.05       |
| G=~y2      | 0.60             | 0.09        | 0.09       | 1.00                | 0.00         | 0.95     | -0.01    | 0.00        |
| G=~y3      | 0.60             | 0.09        | 0.09       | 1.00                | 0.00         | 0.95     | 0.00     | 0.01        |
| G=~y4      | 0.60             | 0.09        | 0.09       | 1.00                | 0.00         | 0.93     | -0.01    | -0.06       |
| G=~y5      | 0.60             | 0.09        | 0.09       | 1.00                | 0.00         | 0.93     | 0.00     | -0.03       |
| G=~y6      | 0.60             | 0.09        | 0.09       | 1.00                | 0.00         | 0.95     | 0.00     | -0.02       |
| G=~y7      | 0.60             | 0.09        | 0.09       | 1.00                | 0.00         | 0.96     | 0.01     | 0.00        |
| G=~y8      | 0.60             | 0.09        | 0.09       | 1.00                | 0.00         | 0.96     | 0.00     | 0.00        |
| G=~y9      | 0.60             | 0.09        | 0.09       | 1.00                | 0.00         | 0.94     | 0.00     | 0.01        |
| S1=~y1     | 0.40             | 0.13        | 0.16       | 0.79                | 0.00         | 0.98     | 0.01     | 0.21        |
| S1=~y2     | 0.40             | 0.13        | 0.16       | 0.78                | 0.00         | 0.97     | 0.00     | 0.22        |
| S1=~y3     | 0.40             | 0.14        | 0.16       | 0.79                | 0.00         | 0.97     | 0.01     | 0.16        |
| S2=~y4     | 0.39             | 0.13        | 0.16       | 0.75                | -0.01        | 0.98     | -0.02    | 0.21        |
| S2=~y5     | 0.41             | 0.15        | 0.16       | 0.76                | 0.01         | 0.96     | 0.01     | 0.12        |
| S2=~y6     | 0.40             | 0.14        | 0.17       | 0.76                | 0.00         | 0.97     | 0.00     | 0.19        |
| S3=~y7     | 0.40             | 0.13        | 0.17       | 0.77                | 0.00         | 0.98     | 0.00     | 0.29        |
| S3=~y8     | 0.40             | 0.14        | 0.16       | 0.78                | 0.00         | 0.98     | 0.01     | 0.21        |
| S3=~y9     | 0.40             | 0.13        | 0.17       | 0.78                | 0.00         | 0.98     | -0.01    | 0.30        |


```

The rows refer to the various parameters of the model and the columns provide different summary statistics. 
The initial two columns (`Estimate Average` and `Estimate SD`) offer a summary of the average and standard deviation of estimated parameters from various simulation solutions. When considering the standard deviation of parameter estimates across numerous replications, it serves as an approximation of the corresponding standard errors in the population, allowing a comparison against empirically calculated standard errors found in the `Average SE` column. 

Next, the `Power (Not Equal 0)` column indicates the proportion of instances where parameters significantly deviated from zero (at a default $\alpha$-level of 0.05). In our analysis, we observed notably high power for general factor loadings but comparatively lower power for specific factor loadings.

In the `Average Bias` column, the summary captures the disparity between parameter estimates and population parameter values. The `Coverage` column reveals the percentage of solutions where the $(1 - \alpha)$% confidence interval covered the true population value, ideally nearing .95 with an $\alpha$-level of 0.05. Our findings generally exhibit acceptable coverage, ranging between 0.93 and 0.98, in line with previous studies.

The last two columns (`Rel Bias` and `Rel SE Bias`) display the average relative biases of parameter estimates and standard errors. Typically, absolute RB values below .05 are considered negligible, values between .05 and .10 moderate, and values surpassing .10 strongly (and hence unacceptably) biased.

In our analysis, the relative bias for all factor loadings proved negligible. However, the estimated standard errors for general factor loadings displayed minor to moderate biases, while those for specific factor loadings exhibited strong biases.



## Next steps






## install official version of packages from CRAN

install.packages('lavaan')
install.packages('simsem')


## alternatively, install development version of packages (recommended)

install.packages('devtools') # install devtools, in case it is not yet installed
devtools::install_github("simsem/simsem/simsem")

install.packages("lavaan", repos = "https://www.da.ugent.be", type = "source")


## activate packages
library(lavaan)
library(simsem)

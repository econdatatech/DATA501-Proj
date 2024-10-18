## -----------------------------------------------------------------------------
library(baselineARPLss)
data("Abelsonite")
head(Abelsonite,10)

## -----------------------------------------------------------------------------
result<-baseline_estimation(Abelsonite$measurement)
#due to default values that is the equivalent for the following function call
#result <- -baseline_estimation(Abelsonite$measurement, lambda = 1e6, 
#ratio = 1e-6, max_iter = 50,verbose=FALSE, algo="banded")

## -----------------------------------------------------------------------------
plot(result)

## -----------------------------------------------------------------------------
summary(result)


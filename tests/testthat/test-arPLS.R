#Generate testdata
wavenumbers <- seq(0, 1000, length.out = 1000)
raman_spectra <- 100*exp(-((wavenumbers-300)/15)^2) + 200*exp(-((wavenumbers-750)/30)^2) #+ 100*exp(-((wavenumbers-800)/15)^2)



test_that("we test wrong data type of the first parameter of baseline_correction()", {
  expect_error(
    baseline_correction(y = "This is a sting"),
    "The input data must be of type numeric."
  )
})

test_that("we test if the first parameter of baseline_correction() is missing", {
  expect_error(
    baseline_correction(lambda = 1e6, ratio = 1e-6, max_iter = 50,verbose=FALSE),
    "Missing input data. You need to provide an input vector for y"
  )
})

test_that("we test if the first parameter of baseline_correction() has NAs", {
  expect_error(
    baseline_correction(c(1,2,3,4,NA)),
    "The input spectrum contains NA values. You need to correct your input data so it doesn't contain NAs"
  )
})

test_that("we test if the first parameter of baseline_correction() has Infs", {
  expect_error(
    baseline_correction(c(1,2,3,4,Inf)),
    "The input spectrum contains Inf values. You need to correct your input data so it doesn't contain Infs"
  )
})

test_that("we test if the first parameter of baseline_correction() has less than 3 values", {
  expect_error(
    baseline_correction(c(1,2)),
    "The input spectrum must contain at least 3 values"
  )
})

test_that("we test if the first parameter of baseline_correction() has only positive values", {
  expect_error(
    baseline_correction(c(-1,2,0)),
    "The input spectrum must contain only positive values"
  )
})

test_that("we test if the lambda parameter of baseline_correction() has correct format and only positive values", {
  expect_error(
    baseline_correction(c(1,2,0),lambda=c(-1,2)),
    "The parameter 'lambda' must be a single valid numeric value.
         NA, Inf or negative values are not valid."
  )
})

test_that("we test if the ratio parameter of baseline_correction() has correct format and only positive values", {
  expect_error(
    baseline_correction(c(1,2,0),ratio=c(-1,2)),
    "The parameter 'ratio' must be a single valid numeric value.
         NA, Inf or negative values are not valid"
  )
})


test_that("we test if the max_iter parameter of baseline_correction() has correct format and only positive values", {
  expect_error(
    baseline_correction(c(1,2,0),max_iter=c(-1,2)),
    "The parameter 'max_iter' must be a single valid integer value.
         NA, Inf or negative values or fractions and decimals are not valid."
  )
})

test_that("we test if the verbose parameter of baseline_correction() has correct format and only boolean values", {
  expect_error(
    baseline_correction(c(1,2,0),verbose=c(-1,2)),
    "The parameter 'verbose' must be a single boolean value. TRUE or FALSE are the only valid inputs."
  )
})


test_that("We test baseline_correction() with correct input and check if the output is the right class", {
  expect_no_error(if(!inherits(baseline_correction(raman_spectra,max_iter=5,verbose=FALSE),"arPLSresult")){stop("Wrong output class")})
})

test_that("We test baseline_correction() with correct input and check if the output baseline is as expected", {
  expect_no_error(if(sum(baseline_correction(wavenumbers+raman_spectra,verbose=FALSE)$baseline-wavenumbers)>=0.001)
    {stop("Algorithm doesn't seem to work as expected based on synthetic data")})
})

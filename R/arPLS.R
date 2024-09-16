#' Constructor for arPLSresult object
#'
#' @description
#' This is an constructor for the S3 object arPLSresult.
#'
#' @param rawinput The original spectrum fed into the algorithm.
#' @param lambda The lambda parameter fed into the algorithm.
#' @param ratio The ratio stopping parameter fed into the algorithm.
#' @param max_iter The maximum iteration stopping parameter fed into the algorithm.
#' @param baseline The fitted spectral baseline.
#' @param last_iter The number of iterations the algorithm did before stopping.
#' @param last_ratio The last value of the ratio stopping criterium before stopping.
#' @export
new_arPLSresult <- function(rawinput = numeric(), lambda = 1e6, ratio=1e-6, max_iter=50,
                            baseline=numeric(),last_iter=integer(),last_ratio=double()) {
  object<-list(rawinput=rawinput, lambda=lambda,ratio=ratio,max_iter=max_iter,baseline=baseline,
            last_iter=last_iter,last_ratio=last_ratio)
  attr(object,"class") <- "arPLSresult"
  object
}

#' Take an object of class arPLSresult and plot some results
#'
#' @description
#' This is an S3 generic. To plot an input spectrum and an estimated baseline spectrum.
#'
#' @param x A result object of class arPLSresult (mainly a list).
#' @param ... placeholder for arbitrary additional parameters (to stay in line with other generic plot functions)
#' @export
  plot.arPLSresult<- function(x,...){
  plot(x$rawinput,main="arPLS baseline estimation",ylab="Measurements")
  graphics::lines(x$baseline,col='red')
}

#' Take an object of class arPLSresult and summarize (print) some facts about it
#'
#' @description
#' This is an S3 generic.
#' To summarize (print) some facts about the arPLS baseline estimation that led to it.
#'
#'
#' @param object A result object of class arPLSresult (mainly a list).
#' @param ... placeholder for arbitrary additional parameters (to stay in line with other generic summary functions)
#' @export
summary.arPLSresult<- function(object, ...){
  print(paste("The lambda parmeter value used was: ",object$lambda))
  print(paste("The ratio parameter value used was: ",object$ratio))
  print(paste("The max_iter parameter value used was: ",object$max_iter))
  print(paste("The alogrithm stopped after the following number or iterations: ",object$last_iter))
  print(paste("The last weight ratio value was: ",object$last_ratio))
  if(object$last_iter==object$max_iter){
    print("It appears that the algorithm stopped because the maximum number of iterations was reached")
  }
  if(object$last_ratio<object$ratio){
    print("It appears that the algorithm stopped because the change in weights per iteration fell below the ratio threshold")
  }
}

#' @title asymmetrically reweighted penalized least squares
#'
#' @author Corvin Idler
#' @description
#' Baseline correction using asymmetrically reweighted penalized least squares smoothing (Baek et al. 2015).
#'
#' @details
#' The algorithm iteratively estimates a spectral baseline curve by updating a
#' weight vector by means of a generalized logistic function that focuses the
#' estimation efforts on regions where the baseline and the signal are close to
#' each other
#'
#' @references
#' Baek, S.-J., Park, A., Ahn, Y.-J., and Choo, J. (2015). Baseline correction
#' using asymmetrically reweighted penalized least squares smoothing.
#' Analyst, 140:250–257.
#'
#' @param y Numeric vector representing the spectrum.
#' @param lambda Smoothing parameter. The smaller the more curvature (wiggliness). (default: 1e6).
#' @param ratio Stopping criterion based on changes in weight vector per iteration (default: 1e-6).
#' @param max_iter Maximum number of iterations as fall back criterion if no conversion happens (default: 50).
#' @param verbose Boolean to print intermediary outputs (default: FALSE).
#' @param cpp Boolean to change between Armadillo CPP armaInv and native solver function  (default: TRUE).
#' @return object of class arPLSresult:
#' \itemize{
#'   \item \code{rawinput}: The original spectrum fed into the algorithm.
#'   \item \code{lambda}: The lambda parameter fed into the algorithm.
#'   \item \code{ratio}: The ratio stopping parameter fed into the algorithm.
#'   \item \code{max_iter}: The maximum iteration stopping parameter fed into the algorithm.
#'   \item \code{baseline}: The fitted spectral baseline.
#'   \item \code{last_iter}: The number of iterations the algorithm did before stopping.
#'   \item \code{last_ratio}: The last value of the ratio stopping criterium before stopping.
#' }
#'
#' @examples{
#' y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' baseline <- baseline_correction(y)
#' }
#' @export baseline_correction
baseline_correction <- function(y, lambda = 1e6, ratio = 1e-6, max_iter = 50,verbose=FALSE,cpp=TRUE) {
# default values from page 253 of original publication

#input validation
  if (missing(y)) {
    stop("Missing input data. You need to provide an input vector for y")
  }
  if (!is.numeric(y)) {
    stop("The input data must be of type numeric.")
  }
  if (anyNA(y)) {
    stop("The input spectrum contains NA values. You need to correct your input data so it doesn't contain NAs")
  }
  if (sum(is.infinite(y)) >= 1) {
    stop("The input spectrum contains Inf values. You need to correct your input data so it doesn't contain Infs")
  }
  if (length(y) <= 2) {
    stop("The input spectrum must contain at least 3 values")
  }
  if (length(y[y<0]) > 0) {
    stop("The input spectrum must contain only positive values")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || sum(is.infinite(lambda)) >= 1 || lambda <0 ) {
    stop("The parameter 'lambda' must be a single valid numeric value.
         NA, Inf or negative values are not valid.")
  }
  if (!is.numeric(ratio) || length(ratio) != 1 || sum(is.infinite(ratio)) >= 1   || lambda <0 ) {
    stop("The parameter 'ratio' must be a single valid numeric value.
         NA, Inf or negative values are not valid.")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 || sum(is.infinite(max_iter)) >= 1   || max_iter <0 || max_iter%%1!=0 ) {
    stop("The parameter 'max_iter' must be a single valid integer value.
         NA, Inf or negative values or fractions and decimals are not valid.")
  }
  if (!is.logical(verbose) || length(verbose) != 1  ) {
    stop("The parameter 'verbose' must be a single boolean value. TRUE or FALSE are the only valid inputs.")
  }

  #get ready with plotting setup
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(2, 1))

  n <- length(y)
  D <- diff(diag(n), differences = 2)
  H<-lambda*(t(D) %*% D)
  w <- rep(1, n)
  wold <- w
  arPLSresult<-new_arPLSresult(rawinput=y,lambda,ratio,max_iter,baseline=y,0,1.0)
  for (i in 1:max_iter) {
    if(verbose){print(paste("Iteration: ",i))}
    arPLSresult$last_iter<-i

    W <- diag(as.numeric(w))
    #z <- solve(W + H) %*% (W %*% y)
    #z<-limSolve::Solve.banded(W + H) %*% (W %*% y)
    z <- limSolve::Solve.banded(W+H,2,2,W %*%y)

    #z <- rcpparma_armaInv(W + H) %*% (W %*% y)
    arPLSresult$baseline<-z
    d=y-z
    dneg=d[d<0]

    #Not sure what to do when there are no or only 1 element with d<0
    #Doesn't allow for calculation of sd and new weight vector
    #So I abort at that point.
    if (length(dneg) <=1){
      return(arPLSresult)
      break;}

    m=mean(dneg)
    s=stats::sd(dneg)
    w <- 1/(1+exp(2*(d-(-m+2*s))/s))

    arPLSresult$last_ratio<-norm(wold-w, type = "2")/ norm(wold, type = "2")
    if(verbose){
      print(paste("Weight vector change ratio: ",arPLSresult$last_ratio))
      # store original settings and then divide frame in 2X1 grid
      plot(y,main=paste("Itreation: ",i))
      graphics::lines(z,col='red')
      plot(w,main="Weights")
    }
    if (arPLSresult$last_ratio < ratio){
      return(arPLSresult)
      break;}
    wold<-w
  }
  return(arPLSresult)
}

#' Raw Raman spectrum for Abelsonite
#' A data frame containing 3315 rows and 2 variables (wavenumber and measurement)
#' @author Bob Downs \email{rdowns@u.arizona.edu}
#' @references
#' Lafuente B, Downs R T, Yang H, Stone N (2015) The power of databases: the RRUFF project.
#' In: Highlights in Mineralogical Crystallography,
#' T Armbruster and R M Danisi, eds. Berlin, Germany, W. De Gruyter, pp 1-30
#' @source \url{https://rruff.info/repository/sample_child_record_raman_full/by_minerals/Abelsonite__R070007__Broad_Scan__532__0__unoriented__Raman_Data_RAW__13756.txt}
"Abelsonite"


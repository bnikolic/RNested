
##' Create a Gaussian Likelihood Function for multiple observations
##'
##' This function is appropriate when the uncertainty on the
##' measurements yobs are independent and normally distributed with
##' s.d. given by err
##' @title nested.GaussLkl
##' @param x Points at which the system is observed (these are not-the
##' fitted for parameters). For example, for fitting of spectra this
##' would be an array of frequencies at which observations have been
##' made
##' @param yobs The measurements at points x. 
##' @param model The function which models the observations yobs at points x given some parameters p
##' @param err The standard deviation of uncertainty at each observation
##' @return The generated likelihood function (of one variable only,
##' the point in parameter space) is retured 
##' @author bnikolic
nested.GaussLkl <- function(x, yobs, model, err)
  {
    lf <- function(p) {
      ym <- model(x, p)
      res <- 0.5 * sum(((ym -yobs)/err)**2) +0.5 * sum(log(2*pi* err**2))
      return (-res)
    }
    return (lf);
  }


##' Polynomial model that is sometimes used in literature for fitting
##' of spectra
##' @title nested.PolyModel
##' @param x 
##' @param p 
##' @return 
##' @author bnikolic
nested.PolyModel <- function(x, p)
  {
    r <- p[[1]]
    for (i in 2:length(p))
         {
           r <- r+p[[i]]*x**(i-1)
         }
    10**r
  }

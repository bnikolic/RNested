# (C) 2012 Bojan Nikolic <bojan@bnikolic.co.uk>
#
# GPL v2 License
#
# Output of results from nested sampling

library(plyr)
library(foreach)
library(plotrix)
library(ks)

##' A brief summary of output of the nested sampling
##'
##' @title nested.summary
##' @param r The data frame with output of  nested sampler
##' @return No return value
##' @author bnikolic
nested.summary <- function(r)
  {
    evidencec <- cumsum(exp(r$ll)*r$w)
    plot(stepfun(1:(length(evidencec)-1),evidencec ),
         main="Evidence growth curve",
         ylab="Evidence",
         xlab="Sample number",
         )
    cat(" *** Evidence: "  , tail(evidencec,1) , "\n")
    
  }


##' Histogram the marginal probability  distributions of each
##' parameter (diagonal) and the 2d-contour-plot the joint probability
##' distribution of each pair of model parameters
##'
##' @title nested.hist2
##' @param r Output of nested sampling
##' @return 
##' @author bnikolic
nested.hist2 <- function(r)
  {
    N <- dim(r$p)[2]
    par(mfrow=c(N,N))
    for (i in 1:N)
      for (j in 1:N)
        {
          if (i==j)
            {
                      weighted.hist(r$p[,i],
                                    exp(r$ll)*r$w,
                                    main=paste("Marginal probability distribution of parameter", i),
                                    breaks="FD")
                    }
          else
            {
              w <- exp(r$ll)*r$w
              x2 <- kde(r$p,
                        H=Hpi(r$p, binned=TRUE),
                        w=exp(r$ll)*r$w)
              plot(x2)
            }
        }
  }


##' Draw a fan-diagram for a one-dimensional function model, e.g., a
##' spectrum.  
##'
##' Each vertical slice through the fan-diagram represents the
##' probability distribution of the model at that abcissa value. 
##' 
##' @title 
##' @param r Results of nested sampling
##' @param m The model 
##' @param xmin 
##' @param xmax 
##' @param ymin 
##' @param ymax 
##' @param nbins Number of bins to use (same number is used in the obcisa and ordinate directions) 
##' @param logy  Transform the ordinate axis into a log-axis (useful
##' for power-law and exponential  models and their combinations)
##' @return 
##' @author bnikolic
nested.fan <- function(r, m, xmin, xmax, ymin, ymax,
                       nbins=100,
                       logy=FALSE)
  {
    xvals <- seq(xmin, xmax,  length.out=nbins)
    if (logy)
      {
        ylabs <- seq(log10(ymin), log10(ymax), length.out=nbins+1)
        ycuts <- 10**ylabs
      }
    else
      {
        ycuts <- seq(ymin, ymax, length.out=nbins+1)
        ylabs <- ycuts
      }

    midpoints <- function(x) (x[-1]+x[-length(x)])/2
    
    pw <-cbind(xvals, foreach(i=1:nrow(r), .combine=c) %do% m(xvals, r$p[i,]))
        
    m <- tapply(rep(exp(r$ll)*r$w,nbins),
                list(pw[,1],
                     cut( pw[,2], ycuts, include.lowest=TRUE)),
                sum)
    image(xvals,
          midpoints(ylabs),
          m,
          ylim=c(min(ylabs),max(ylabs)));
  }

##' Given a model and the result of nested sampling, compute the
##' probability distribution of the value of the model at some
##' evaluation point x
##'
##' @title 
##' @param r Output of the nested sampler
##' @param m Model
##' @return  Output of histogram routine
##' @author bnikolic
nested.mhist <- function(r, m, x)
  {
    w <- exp(r$ll)*r$w
    pw <-foreach(i=1:nrow(r), .combine=c) %do% m(x, r$p[i,])
    m <- w!=0
    w <- w[m]
    pw <- pw[m]
    weighted.hist(pw,
                  w,
                  main=paste("Probability distribution of model at ", x),
                  breaks="FD")
  }

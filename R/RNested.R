# Bojan Nikolic <bojan@bnikolic.co.uk>
# Nested sampling, implemented in R

library(MASS)


##' Create a starting set based on a box prior
##'
##' @title sset.box
##' @param box Matrix with boundaries of the box
##' @param nss Number of elements of starting set
##' @param llfn Function to initilalise the log-likelihod
##' @return The starting set as a frame
##' @author bnikolic
sset.box <- function(box, nss,
                     llfn)
  {
    p <- apply(box,
               2,
               function(r) {
                 runif(nss, min=r[[1]], max=r[[2]]) 
               })
    # Also calculate the likelihood and set prior probability for all
    # the points in the starting set
    ll <- apply(p,
                1,
                llfn)
    lpr <- rep(0.0,
               nss)
    data.frame(p=I(p), ll=ll, lpr=lpr)
  }

##' Create a function that will accept or reject proposals for the new sample
##'
##' @title mkPriorSamplePred
##' @param s Row of the live set which is being replaced
##' @return The predicate function
##' @author bnikolic
mkPriorSamplePred <- function(s)
{
  function(ll, lp)
    {
      if(ll < s$ll)
        {
          return (FALSE)
        }
      else
        {
          if(lp >= s$lpr)
            {
              return (TRUE);
            }
          else
            {
              if ( exp(lp-s$lpr) > runif(1))
                {
                  return (TRUE);
                }
              else
                {
                  return (FALSE);
                }
            }
        }
    }
}

##' Rectangular offseting with user supplied scales, pretty much the
##' simplest offsetter
##'
##' @title rectOffseter
##' @param scale  Array of scales (one for each parameter of the
##' problem) for generation of the offsets
##' @return Function that generates the offset
##' @author bnikolic
rectOffseter <- function(scale)
  {
    function()
      {
        rnorm(length(scale), sd=scale)
      }
  }


##' Return a random element of the live set
##'
##' @title randomEl
##' @param cs The live (or current) set
##' @return The selected element
##' @author bnikolic
randomEl <- function(cs)
  {
    N <- dim(cs)[1]
    return (cs[as.integer(runif(1, min=1, max=N)),]    )
  }

mkFixedRectProp <- function(scales)
  {
    off  <- rectOffseter(scales)
    function (x)
      {
        x+off()
      }
  }

##' Constrained Prior Sampler using modified MCMC
##'
##' Note that the performance of the overal nested sampler is close to
##' linearly dependent on the value of the parameter n. When  a good
##' proposer (in the sense that it proposes that are uncorrelated with
##' the starting or existing points but tend to satisfy the likelihood
##' constraint) for new points is used this value of n can be reduced
##' leading to better performance.
##' 
##' @title CPChain
##' @param s The starting point
##' @param proposer Function to propose new points
##' @param n Number of steps to make 
##' @param llf Log-likelihood function
##' @param lpf Log-prior function
##' @param cs Live (or current) set
##' @return Either the new point or if new point was not found FALSE
##' @author bnikolic
CPChain <- function(s,
                    proposer,
                    n,
                    llf, lpf,
                    cs)
  {
    pred <- mkPriorSamplePred(s)
    pcurr <- randomEl(cs)$p
    ll <- 0
    lp <- 0
    r <- sapply(1:n, function(x) {
      pnew <- proposer(pcurr)
      llnew <- llf(pnew)
      lpnew <- lpf(pnew)
      if (pred(llnew, lpnew))
        {
          pcurr <<- pnew
          ll <<- llnew
          lp <<- lpnew
          return (TRUE);
        }
      else
        {
          return (FALSE);
        }
    })
    if( any(r) )
      {
        return (list(p=pcurr, ll=ll, lpr=lp));
      }
    else
      {
      return (FALSE);
    }
  }

mkSimplestPSampler <- function(s)
  {
    proposer <- mkFixedRectProp(c(s, s));
    function(worst,llf,lpf,cs) { CPChain(worst,
                                         proposer,
                                         100,
                                         llf,
                                         lpf,
                                         cs)
                               }
  }

##' Use the covariance of the points in the live set to determine the
##' scales for sampling the prior space
##'
##' @title mkCovarianceSampler
##' @param s Scale the covariances of the live set by this factor
##' before using them for sampling
##' @return Constrained prior sampler 
##' @author bnikolic
mkCovarianceSampler <- function(s=1.0)
  {
    cvm <- 0
    proposer <- function(x)
      {
        N <- dim(cvm)[[1]]
        x+mvrnorm(n=1,
                  mu=rep(0,N),
                  Sigma=cvm*s)
      }
    function(worst,llf,lpf,cs) {
      cvm <<- cov(cs$p)
      CPChain(worst,
              proposer,
              100,
              llf,
              lpf,
              cs)
    }
  }


##' Create a box prior function
##'
##' The prior function returned takes position in parameter space as
##' an argument and returns log-prior probability. Value of -998 is
##' returned outside the box and 0 inside the box
##' @title boxp
##' @param box The box inside which the prior probability is finite
##' @return Box prior function
##' @author bnikolic
boxp <- function(box)
  {
    ff <- function (x)
      {
        for (i in 1:length(x))
            if (x[[i]] < box[[2*i-1]] || x[[i]] > box[[2*i]] )
              {
                return( -998)
              }
        return (0);
      }
    return(ff)
  }

##' Take one step of the nested sampler
##'
##' @title nested.step
##' @param cs The current (live) set
##' @param llf Log-likelihood function
##' @param lpf Log-prior function 
##' @param psampler The prior space sampler
##' @return list (new current set, eliminated row)
##' @author bnikolic
nested.step <- function(cs,
                        llf, lpf,
                        psampler)
  {
    worsti <- which.min(cs$ll)
    worst <- cs[worsti,]
    newp <- psampler(worst, llf, lpf, cs)
    if (identical(newp,FALSE))
      {
        return (newp);
      }
    else
      {
      cs[worsti,] <- newp
      return (list(cs, worst));
      }
  }

##' Do nested sampling of the suplied pair of log-likelihood and prior
##' probability functions
##'
##' This function is the top-level driver for the nested sampling
##' algorithm. It requires preparation of a starting set (using for
##' example sset.box) and a constrained prior space sampler (using for
##' example mkCovarianceSampler). It allows restart of the sampling if
##' inspection of outputs shows that the sampler may not have
##' converged
##' @title nested.sample
##' @param cs The current/live/starting set of poings
##' @param llf  The log-likelihood function (which should take one
##' argument which is the point in parameter space at which to evalue
##' the likelihood)
##' @param lpf the log-prior probability function (which should take
##' one argument, which is the point in parameter sample)
##' @param psampler The prior sampler to use to advance find sample a
##' new point  from the prior space 
##' @param cout Optional parameter which is the output of a previous
##' run of the nested sampler. Allows easy restart of the algorithm
##' @param N Number of samples to make
##' @return List of the live set  and the output of nested
##' sampling. See nested.summary and nested.hist2 on how to interpret
##' the output
##' @author bnikolic
nested.sample <- function(cs,
                          llf, lpf,
                          psampler,
                          cout=rbind(),
                          N=1)
  {
    for (i in 1:N)
      {
        r <- nested.step(cs, llf, lpf, psampler)
        if (identical(r,FALSE))
          break;
        cs <- r[[1]]
        nsamples <- dim(cout)[[1]]
        if (is.null(nsamples))
          nsamples <- 0;
        nlive <-    dim(cs)[[1]]
        weight <- exp(-1.0*nsamples/nlive) - exp(-1.0*(nsamples+1)/nlive)
        cout <- rbind(cout, data.frame(c(r[[2]][1,], w=weight)))
      }
    return (list(cs=cs, cout=data.frame(cout,
                          row.names=NULL)));    
  }



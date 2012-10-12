# Bojan Nikolic <bojan@bnikolic.co.uk>
# Nested sampling, implemented in R

##' Create a starting set based on a box prior
##'
##' .. content for \details{} ..
##' @title 
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
    ll <- apply(p,
                1,
                llfn)
    lpr <- rep(0.0,
               nss)
    data.frame(p=I(p), ll=ll, lpr=lpr)
  }

##' Create a function that will accept or reject proposals for the new sample
##'
##' .. content for \details{} ..
##' @title 
##' @param s Row of the live set which is being replaced
##' @return 
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
          if(lp > s$lpr)
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
##' .. content for \details{} ..
##' @title 
##' @param scale 
##' @return 
##' @author bnikolic
rectOffseter <- function(scale)
  {
    function()
      {
        rnorm(length(scale), sd=scale)
      }
  }

CPChain <- function(s, scale, n,
                    llf, lpf)
  {
    pred <- mkPriorSamplePred(s)
    off  <- rectOffseter(scale)
    pcurr <- s$p
    ll <- 0
    lp <- 0
    r <- sapply(1:n, function(x) {
      pnew <- pcurr+off()
      ll <<- llf(pnew)
      lp <<- lpf(pnew)
      if (pred(ll, lp))
        {
          pcurr <<- pnew
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

##' Create a box prior function (old)
##'
##' @title 
##' @param box 
##' @return Box prior function
##' @author bnikolic
boxp <- function(box)
  {
    ff <- function (x)
      {
        for (i in 1:length(x))
            if (x[i] < box[2*i-1] || x[i] > box[2*i] )
              {
                return( -998)
              }
        return (0);
      }
    return(ff)
  }


nested.sample <- function(ss, llf, lpf,
                        psampler)
  {
    #inclomplete
    psampler(ss[which.min(ss$ll)],
             lpf)
  }
              

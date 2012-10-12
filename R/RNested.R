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
    data.frame(p=p, ll=ll, lpr=lpr)
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


nested.step <- function(ss, llf, lpf,
                        psampler)
  {
    #inclomplete
    psampler(ss[which.min(ss$ll)],
             lpf)
  }
              

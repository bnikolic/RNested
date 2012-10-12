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
          if(lp >= s$lpr)
            {
              return (TRUE);
            }
          else
            {
              if ( exp(lp-s$lpr) > runif(1))
                {
                  print ("accepted ")
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

##' Make one step of the nested sampler
##'
##' .. content for \details{} ..
##' @title 
##' @param cs 
##' @param llf 
##' @param lpf 
##' @param psampler 
##' @return list (new current set, eliminated row)
##' @author bnikolic
nested.step <- function(cs,
                        llf, lpf,
                        psampler)
  {
    worsti <- which.min(cs$ll)
    worst <- cs[worsti,]
    cs[worsti,] <- psampler(worst, llf, lpf)
    return (list(cs, worst));
  }


nested.sample <- function(cs,
                          llf, lpf,
                          psampler,
                          cout=NA,
                          N=1)
  {
    for (i in 1:N)
      {
        r <- nested.step(cs, llf, lpf, psampler)
        cs <- r[[1]]

        if(is.na(cout))
          {
            cout <- rbind(r[[2]])
          }
        else
          {
            cout <- rbind(cout, r[[2]])
          }
      }
    return (list(cs, data.frame(cout)));    
  }
              

# (C) 2012 Bojan Nikolic <bojan@bnikolic.co.uk>
#
# GPL v2 License
#
# Output of results from nested sampling

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
                                    main=paste("Marginal probability distribution of parameter", i))
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


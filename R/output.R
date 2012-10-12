# (C) 2012 Bojan Nikolic <bojan@bnikolic.co.uk>
#
# GPL v2 License
#
# Output of results from nested sampling

library(plotrix)
library(fields)

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
              ff <- data.frame(x=r$p[,i],
                               y=r$p[,j])
              # This does not do an equivalent of a histogram!
              #contour(smooth.2d( w, x=ff))
            }
        }
  }



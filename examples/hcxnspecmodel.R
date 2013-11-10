library('devtools');
library('ggplot2');
load_all("/home/bnikolic/d/n/RNested/")



hcxnS <- function(p, J)
  {
    
    ccb <- 290.51832257;
    ccm <- 5.205;
    ccd <- 0.87478e-6;
    ccn <- 2.86e-3;
          
    kB <-  1.380e-16;
    hP <-  6.626e-27;
    n <-  10**p[[1]];
    T <- 10**p[[2]];
    
    JUp <- J+1;
    twoJone <- (2*JUp) + 1;
    twoJthree <- (2*J) + 3;

    jFactor <- n * twoJone;

    nu <- function(J)
      {
        JUp <- J+1;
        2*ccb*JUp - 4*ccd*JUp**3;
      }

    jNu <- nu(J);

    oneOverkT <-  1e6/(kB * T);

    Z <- function ()
      {
        z <- 0;
        for (j in 0:1000)
          {
            efact <- -hP * ccb * (j+1.)*j * oneOverkT;
            z <- z +((2.*j) + 1.)*exp(efact);
          }
        z
      }

    
    Aj <- jNu**3 * ccm**2 * JUp / twoJthree;

    nj <-  jFactor * Aj * ccn *
      exp( -hP * ccb * JUp * (JUp+1)  * oneOverkT) / Z();
    nj;
  }

mkLkl <- function()
  {
    dd <- read.csv("hc9ndata.csv", header=T)
    nested.GaussLkl(dd$J, dd$S, function(x,p) {hcxnS(p,x) }, dd$SErr);
  }


# optim( c(-9,0.5) , function(p) { -1*xx(p)}) check convergance with.
# ggplot(dd) + geom_errorbar(aes(L, S, ymin=S-SErr, ymax=S+SErr, colour='red')) + geom_point(aes(L, hcxnS(c(-9.333890,  0.659134), dd$L)))

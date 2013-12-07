library('devtools');
library('ggplot2');
library(densityvis)
load_all("/home/bnikolic/d/n/RNested/")

nuX <- function(X)
  {
    xx=read.csv("hcxconsts.csv", header=T, row.names=1 )
    ccb <- xx[X, "ccb"]
    ccd <- xx[X, "ccd"]
    nu <- function(J)
      {
        JUp <- J+1;
        2*ccb*JUp - 4*ccd*JUp**3;
      }
    nu
  }

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
        zi <- function(j)
          {
            efact <- -hP * ccb * (j+1.)*j * oneOverkT;
            ((2.*j) + 1.)*exp(efact);
          }
        sum(zi(0:1000))
      }

    
    Aj <- jNu**3 * ccm**2 * JUp / twoJthree;

    nj <-  jFactor * Aj * ccn *
      exp( -hP * ccb * JUp * (JUp+1)  * oneOverkT) / Z();
    nj;
  }

mkLkl <- function()
  {
    dd <- read.csv("hc9ndata.csv", header=T, comment.char="#")
    nested.GaussLkl(dd$J, dd$S, function(x,p) {hcxnS(p,x) }, dd$SErr);
  }

prepp <- function()
  {
    pp <- c( -12, -9, 0, 1)
    lpf <- boxp(pp)
    ss <- sset.box(array( pp , dim=c(2,length(pp)/2)),
                   100, 
                   mkLkl())
    list(lpf=lpf, ss=ss)
  }


# optim( c(-9,0.5) , function(p) { -1*xx(p)}) check convergance with.
# ggplot(dd) + geom_errorbar(aes(L, S, ymin=S-SErr, ymax=S+SErr, colour='red')) + geom_point(aes(L, hcxnS(c(-9.333890,  0.659134), dd$L)))

# ggplot(dd) + geom_errorbar(aes(L, S, ymin=S-SErr, ymax=S+SErr, colour='red')) + geom_point(aes(L, hcxnS(c(-9.333890,  0.659134), dd$L))) + facet_grid(. ~ L, scales="free_x") + theme(axis.text.x=element_blank() ,axis.ticks.x = element_blank() ,panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())

#ggplot(dd) + scale_y_continuous(limits=c(0, 1.25))+ geom_errorbar(aes(J, S, ymin=S-SErr, ymax=S+SErr, colour='red')) + geom_point(aes(J, hcxnS(c(-9.333890,  0.659134), dd$J))) + facet_grid(. ~ J, scales="free_x") +  theme_bw()+ theme(axis.text.x=element_blank() ,axis.ticks.x = element_blank() ,panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.title.x = element_text(face="bold", size=12),          axis.title.y = element_text(face="bold", size=12, angle=90), legend.title=element_blank())+     stat_abline(intercept=0, slope=0, linetype="dotted")


#p <- prepp()
#r <- nested.sample(p$ss, xx,  p$lpf,   mkCovarianceSampler(),                    N=1000)

# local_density_2d(r$cout$p[ww >0,1],  r$cout$p[ww >0,2], weight=(exp(r$cout$ll)*r$cout$w)[ww >0] )

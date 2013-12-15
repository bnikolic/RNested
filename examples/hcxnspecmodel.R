library('devtools');
library('ggplot2');
library('grid');
library(Cairo)
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

##' Turn an density function to dataframe suitable for tile plot
##'
##' Suitable for use with output from local_density_2d. (TODO: Move
##' this function to suitable file)
##' @title 
##' @param de density estimated function
##' @param n  Number of points on each axis to evalute
##' @return 
##' @author bnikolic
densitydf <- function(de, n)
  {
    ax <- function(lims) {seq(lims[[1]], lims[[2]], length.out=n)}
    gg <- expand.grid(x=ax(attributes(de)$xlim),
                      y=ax(attributes(de)$ylim))
    gg <- data.frame( x=gg$x, y=gg$y, z= de(gg$x, gg$y))
    gg
  }

# optim( c(-9,0.5) , function(p) { -1*xx(p)}) check convergance with.
# ggplot(dd) + geom_errorbar(aes(L, S, ymin=S-SErr, ymax=S+SErr, colour='red')) + geom_point(aes(L, hcxnS(c(-9.333890,  0.659134), dd$L)))

# ggplot(dd) + geom_errorbar(aes(L, S, ymin=S-SErr, ymax=S+SErr, colour='red')) + geom_point(aes(L, hcxnS(c(-9.333890,  0.659134), dd$L))) + facet_grid(. ~ L, scales="free_x") + theme(axis.text.x=element_blank() ,axis.ticks.x = element_blank() ,panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())

#ggplot(dd) + scale_y_continuous(limits=c(0, 1.25))+ geom_errorbar(aes(J, S, ymin=S-SErr, ymax=S+SErr, colour='red')) + geom_point(aes(J, hcxnS(c(-9.333890,  0.659134), dd$J))) + facet_grid(. ~ J, scales="free_x") +  theme_bw()+ theme(axis.text.x=element_blank() ,axis.ticks.x = element_blank() ,panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.title.x = element_text(face="bold", size=12),          axis.title.y = element_text(face="bold", size=12, angle=90), legend.title=element_blank())+     stat_abline(intercept=0, slope=0, linetype="dotted")


#p <- prepp()
#r <- nested.sample(p$ss, xx,  p$lpf,   mkCovarianceSampler(),                    N=1000)

#ww= exp(r$cout$ll)*r$cout$w
# local_density_2d(r$cout$p[ww >0,1],  r$cout$p[ww >0,2], weight=(exp(r$cout$ll)*r$cout$w)[ww >0] )


# See http://stackoverflow.com/questions/11546256/two-way-density-plot-combined-with-one-way-density-plot-with-selected-regions-in/11590446#11590446
# for how to combine plots

dplot <- function(cc)
  {
    ww  <- exp(cc$ll) * cc$w
    dfn <- local_density_2d(cc$p[ww >0, 1],
                            cc$p[ww >0, 2],
                            weight=ww[ww >0] )
    mm <- max(ww)
    xl = c(min(cc$p[ww > mm * 1e-10, 1]),
             max(cc$p[ww >mm * 1e-10, 1]))
    yl = c(min(cc$p[ww >mm * 1e-10, 2]),
           max(cc$p[ww >mm * 1e-10, 2]))

    th <- function() {
      theme_bw() + theme(
                     #axis.text.x=element_blank() ,
#            axis.ticks.x = element_blank() ,
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            axis.title.x = element_text(face="bold", size=12),
            axis.title.y = element_text(face="bold", size=12, angle=90),
            legend.title=element_blank() )
    }
    
    p1 <- ggplot_gtable(ggplot_build(ggplot(densitydf(dfn, 300),
                                            aes(x, y))+
                                     geom_tile(aes(fill=z))+ 
                                     geom_contour(aes(z=z))+
                                     coord_cartesian(xlim=xl, ylim=yl)+
                                     scale_fill_gradient(low = "white",
                                                         high = "steelblue")+
                                     th()+
                                     xlab(expression(T))+
                                     ylab(expression(rho))))

    gt1 <- gtable_add_cols(p1, unit(0.3, "null"), pos = -1)
    gt1 <- gtable_add_rows(gt1, unit(0.3, "null"), pos = 0)

    margp <- function(pii,
                      l,
                      flip)
      {
        d11 <- local_density_1d(cc$p[ww >0, pii],
                                weight=ww[ww >0] )

        pp <- ggplot(data.frame(x = l), aes(x)) +
          stat_function(fun = d11, n=500) +
            th() 
        
        if (flip)
          pp <- pp + coord_flip(l) + theme(axis.ticks.x = element_blank(),
                                                axis.text.x=element_blank())
        else
          pp <- pp + coord_cartesian(l) + theme(axis.ticks.y = element_blank(),
                                                axis.text.y=element_blank())

        ggplot_gtable(ggplot_build(pp))
      }

    p2 <- margp(1, xl, FALSE)

    gt1 <- gtable_add_grob(gt1, p2$grobs[[which(p2$layout$name == "panel")]],
                           1, 4, 1, 4)
    gt1 <- gtable_add_grob(gt1, p2$grobs[[which(p2$layout$name == "axis-l")]],
                           1, 3, 1, 3, clip = "off")

    p3 <- margp(2, yl, TRUE)
    gt1 <- gtable_add_grob(gt1,
                           p3$grobs[[which(p3$layout$name == "panel")]],
                           4, 6, 4, 6)
    gt1 <- gtable_add_grob(gt1,
                           p3$grobs[[which(p3$layout$name == "axis-b")]],
                           5, 6, 5, 6, clip = "off")
    

    #grid.newpage()
    grid.draw(gt1)
    
  }

# CairoPDF("test.pdf"); dplot(r$cout); dev.off()

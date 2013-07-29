## Functions from the paper "Analysis and synthesis of the variability
## of irradiance and PV power time series with the wavelet transform",
## Solar Energy, 85:1 (188-197), 2010
## Copyright (C) 2011, 2010 Oscar Perpiñán Lamigueiro
 
 ## This program is free software; you can redistribute it and/or
 ## modify it under the terms of the GNU General Public License as
 ## published by the Free Software Foundation; either version 2 of the
 ## License, or (at your option) any later version.
 
 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.
 
 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, write to the Free Software
 ## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 ## /

### The sinc function
sinc <- function(x){y=sin(pi*x)/(pi*x); y[is.na(y)] <- 1; y}

### Aggregate a time series
agg <- function(x, period=5){
  time <- time(x)
  av.time<-time-as.numeric(time)%%period
  av<-aggregate(x, av.time, mean)
}

### Synthesis of a irradiance time series with different methods
clear <- function(x, lambda){
  k=coredata(x$kt)
  Bo0=coredata(x$Bo0)

  ##With a wavelet bootstrap
  k.boot<-wavBootstrap(k)

  ##With a ARIMA simulation
  media.k=mean(k)
  rango.k<-diff(range(k))
  k=k-media.k                           
    
  k.ar<-ar(ts(k), na.action=na.exclude)
  ##Distribución exponencial con signo generado con una gaussiana
  rexp.sign<-function(n, A, lambda){A/lambda*sign(rnorm(n))*rexp(n,lambda)}

  k.exp<-arima.sim(n=length(k), model=list(ar=k.ar$ar),
                   rand=rexp.sign, A=1e-3, lambda=lambda) #UTILIZAR FISTDISTR de MASS
  rango.exp<-diff(range(k.exp))
  k.exp<-rango.k/rango.exp*(k.exp)
  k.exp=k.exp+media.k
  
  G0.exp<-k.exp*Bo0
  G0.boot<-k.boot*Bo0
 
  result <- zoo(data.frame(k.exp, k.boot, G0.exp, G0.boot), index(x))
}


### Read data from file and create a zoo object
data2zoo <- function(file, thresh=50){
  z <- read.zoo(file,
                colClasses=c('NULL', 'character', 'character',
                  rep('numeric', 18)),
                col.names=  c("x", "Fch","Hor","Irms1","Irms2","Irms3",
                  "Urms1","Urms2","Urms3",
                  "P1","P2","P3","Q1","Q2","Q3",
                  "T_amb","T_pnl_fij","T_pnl_mov",
                  "G0","Gef","Vel"),
                index=c(1, 2),
                FUN=function(fecha, hora){
                  as.POSIXct(paste(fecha, hora), format="%Y/%m/%d %H:%M:%S")
                },
                sep=',', skip=1, header=FALSE, fill=TRUE
                )
  ##Irradiance threshold
  z <- z[z$G0>thresh,]
  ##Active power (from a personal communication with Marcos Navarra)
  Pr=with(z,P1*cos(pi/6)+Q1*sin(pi/6))
  Ps=with(z,P2*cos(pi/6)+Q2*sin(pi/6))
  Pt=with(z,P3*cos(pi/6)+Q3*sin(pi/6))
  P=(Pr+Ps+Pt)*13200/400;
  z$P <- P

  z <- z[,c('G0', 'Gef', 'P', 'Vel')]
  z
}

### Horizon time plot of a wavelet decomposition
wavHorizon <- function(x, name='x', hour.out=TRUE, plotit=TRUE,...){
  wave<-wavMODWT(x) 
  rec<-wavMRD(wave)                     # reconstruction of the signal

  delta=deltat(x)
  tau=wave$scales
  lt=length(tau)
  scales=tau*delta
  h <- scales%/%3600
  aux <- scales-h*3600
  m <- aux%/%60
  s <- aux-m*60

  hms<-paste(ifelse(h>0, paste(h,'h',sep=''), ''),
             ifelse(m>0, paste(m,'m', sep=''), ''),
             paste(s,'s',sep=''),
             sep='')

  waveZ<-as.data.frame(wave$data)
  waveZ<-zoo(waveZ, index(x))
  waveZ$x<-diff(x)
  names(waveZ)<-c(hms, 'S', name)


  recDF <- data.frame(as.matrix(rec))
  recZ <- zoo(recDF, index(x))
  recZ$x <- x
  names(recZ)<-c(hms, 'smooth', name)
  N <- dim(recZ)[2]

  ## Conventional plot of the time series of the wavelet coefficients
  xy<-xyplot(waveZ, ..., layout=c(1,lt+2), strip=FALSE, strip.left=TRUE)

  out <- c(which(h>0), lt+1)
  if (hour.out) {
    waveZ<-waveZ[,-out]
    rows<-lt+2-length(out)
  } else {
    rows <- lt+2
  }

  ##Horizon plot of the time series of the wavelet coefficients
  hor<-horizonplot(waveZ, ..., layout=c(1, rows), origin=0,
                   col.regions=brewer.pal(n=11, name='RdBu'),
                   par.strip.text = list(cex = 0.35),
                   scales = list(y = list(relation = "free"))
                   )

  ## Conventional plot of the reconstruction of the time series
  xy.rec<-xyplot(recZ, ..., layout=c(1,lt+2), strip=FALSE, strip.left=TRUE)

  ## Horizon plot of the reconstruction of the time series
  hor.rec<-horizonplot(recZ[,-c(N-1, N)], ..., layout=c(1, lt), origin=0,
                       col.regions=brewer.pal(n=11, name='RdBu'),
                       par.strip.text = list(cex = 0.35),
                       bg='gray',
                       scales = list(y = list(relation = "free"))
                       )

  hor.rec <- c(hor.rec,
               xyplot(recZ[,c(N-1, N)],
                      superpose=TRUE, col=c('gray', 'black'),
                      strip=FALSE, strip.left=FALSE),
               layout=c(1,lt+1))
  
  if (plotit) {print(hor.rec)}
  
  result<-list(wave=waveZ, rec=recZ, xy=xy, hor=hor, xy.rec=xy.rec, hor.rec=hor.rec)
 
}

### Wavelet variance of a time series (only a wrapper of wmtsa::wavVar)
wavVar.uni <- function(x){
    aux <- wavVar(coredata(x))
    unbiased.data <- aux$block$unbiased
    period <- deltat(x)
    scales <- as.numeric(aux$scales)*period
    conf.low <- aux$confidence$n3$low
    conf.high <- aux$confidence$n3$high
    df <- data.frame(var=unbiased.data, scale=scales, low=conf.low, high=conf.high)
    rownames(df) <- NULL
    df
  }

### Wavelet variance of a multivariate time series (using wavVar.uni)
wavVar.multi <- function(x){
  result0 <- lapply(x, wavVar.uni)
  result <- do.call('rbind', result0)
  rownames(result) <- NULL
  result$ID <- rep(names(result0), each=dim(result)[1]/length(result0))
  result
}


### Plot of the wavelet variance

##Auxiliar panels
my.panel.bands <- function(x, y, upper, lower, fill, col,
                           subscripts, ..., font, fontface)
    {
      upper <- upper[subscripts]
      lower <- lower[subscripts]
      panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                    col = fill, border = FALSE,
                    ...)
    }
my.panel.lines <- function(x, y, upper, lower, subscripts, ...)
    {
      upper <- upper[subscripts]
      lower <- lower[subscripts]
      panel.xyplot(x, upper, ...)
      panel.xyplot(x, lower, ...)
    }

drop.levels <- function(dat){
  # Drop unused factor levels from all factors in a data.frame
  # Author: Kevin Wright.  Idea by Brian Ripley.
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}

##Main function
plotVar <- function(formula, data, groups=ID, subset=TRUE, ci='bands',
                    xlab='Physical Scale (s)', ylab='Wavelet Variance',
                    pos.lab=1, cex.lab=0.6)
{
  environment(formula) <- environment() ##Necesario para que xyplot acepte los argumentos
  sbt=eval(substitute(subset), data, environment(formula))
  data=drop.levels(data[sbt,])
  groups=eval(substitute(groups), data, environment(formula))
  if (is.null(groups)) {groups <- rep(factor(''), dim(data)[1])}
  ci.panel=switch(ci,
    bands="my.panel.bands",
    lines="my.panel.lines"
    )
  if (ci=='none') {
    data$high <- NA
    data$low <- NA
  }
  p <- xyplot(formula, data=data, groups=groups,
              layout=c(NA, 1),subscripts=TRUE,
              xlab=xlab, ylab=ylab, pch=c(20, 18, 17, 3:14), 
              scales=list(x=list(log=TRUE), y=list(log=TRUE)),
              strip=strip.custom(bg='gray', par.strip.text=list(cex=0.8)),
              xscale.components = xscale.components.log10ticks,
              yscale.components= yscale.components.log10ticks,
              fill = 'lightgray',
              upper = log10(data$high),
              lower = log10(data$low),
              cex.lab=cex.lab,#Se lo pasa a glayer
              pos.lab=pos.lab,#Se lo pasa a glayer
              panel = function(x, y, ...){
                if (ci!='none') {
                  panel.superpose(x, y, panel.groups = ci.panel, type='l', col='gray',...)
                }
                panel.xyplot(x, y, type='b', cex=0.6, lty=1,...)
              }
              )
  if (!is.null(groups))
    {p <- p+glayer(panel.text(x[1], y[1], group.value, cex=cex.lab, pos=pos.lab))}
    p
  }


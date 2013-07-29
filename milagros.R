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


library(latticeExtra)
library(zoo)
library(wmtsa)
library(solaR)

source('funciones.R')

## Customization of the lattice package
lattice.options(default.theme=standard.theme(color=FALSE))

custom.theme.gray <- function(n){
  x <- custom.theme(symbol=gray.colors(n, start=0, end=0.7))
  x$strip.background$col='gray'
  result <- x
  }

######################
###Location data   ###
######################

lat=42+15/60
lon=-(1+46/60)#NEGATIVA por estar al oeste del meridiano de Greenwich (meridiano 0)
Beta=45#Seguidor acimutal
#LH=0 porque el huso horario de la peninsula coincide con el meridiano 0
period=5

#######################
###Source data      ###
#######################
##Read the data
mil13102009<-data2zoo('Milagro_13_10_2009.csv')
mil14102009 <- data2zoo('Milagro_14_10_2009.csv')
mil15102009 <- data2zoo('Milagro_15_10_2009.csv')
mil16102009 <- data2zoo('Milagro_16_10_2009.csv')
## aggregation with a 5 seconds period
mil13102009.av <- agg(mil13102009, period=period)
mil14102009.av <- agg(mil14102009, period=period)
mil15102009.av <- agg(mil15102009, period=period)
mil16102009.av <- agg(mil16102009, period=period)

##Horizon plot 13/10/2009
## trellis.device(pdf,
##                file="ArticuloSimulador/Figuras/GefHorizon.pdf",
##                width=6,height=4.5,
##                family='Courier')

hor <- wavHorizon(mil13102009.av$Gef, name='Gef', hour.out=FALSE, plotit=FALSE)
print(hor$hor.rec)

## dev.off()

## Evolution of the irradiance during these four days
data1 <- scale(mil13102009.av[,c('Gef','P')], center=TRUE, scale=c(1000, 9.5e6))
data2 <- scale(mil14102009.av[,c('Gef','P')], center=TRUE, scale=c(1000, 9.5e6))
data3 <- scale(mil15102009.av[,c('Gef','P')], center=TRUE, scale=c(1000, 9.5e6))
data4 <- scale(mil16102009.av[,c('Gef','P')], center=TRUE, scale=c(1000, 9.5e6))

p1 <- xyplot(data1, superpose=TRUE, auto.key=list(corner=c(1,0), cex=0.6))
p2 <- xyplot(data2, superpose=TRUE)
p3 <- xyplot(data3, superpose=TRUE)
p4 <- xyplot(data4, superpose=TRUE)
p <- c(p1, p2, p3, p4, layout=c(4,1))

## trellis.device(pdf,
##                theme=custom.theme.gray(2),
##                file="ArticuloSimulador/Figuras/EvolucionRadiacionPotencia.pdf",
##                width=6,height=4.5,
##                family='Courier',
##                title='Evolución de la irradiación en Octubre')
update(p, strip=strip.custom(factor.levels=c('13th Oct',
                               '14th Oct',
                               '15th Oct',
                               '16th Oct')),
            scales=list(y=list(draw=FALSE), x=list(rot=15)),
            lwd=0.8, ylab='')
## dev.off()

########################
###Clearness index   ###
########################

## Sun calculations with solaR::calcSol
sol13102009 <- calcSol(lat=lat, BTi=index(mil13102009))
sol13102009 <- cbind(mil13102009,
                     as.zooI(sol13102009)[,c('w', 'cosThzS', 'AlS', 'AzS', 'Bo0')])
sol13102009$kt <- with(sol13102009, G0/Bo0)

##Synthesis of a irradiance time series
sol13.k<-clear(sol13102009, lambda=800)
sol13102009 <- cbind(sol13102009, sol13.k)

sol13102009.av <- agg(sol13102009[,c('G0', 'G0.boot')], period=5)

## trellis.device(pdf, 
##                theme=custom.theme.gray(2),
##                file="ArticuloSimulador/Figuras/G0Boot.pdf",
##                width=6,height=4.5,
##                family='Courier')
xyplot(sol13102009.av, superpose=TRUE,
       auto.key = list(x = .7, y = .9,
      text=c(expression(G(0)), expression(G(0)^B))),
       ylab='Irradiance (W/m²)')
## dev.off()


## trellis.device(pdf, 
##                theme=custom.theme.gray(2),
##                file="ArticuloSimulador/Figuras/G0BootDiff.pdf",
##                width=6,height=4.5,
##                family='Courier')
ecdfplot(~G0+G0.boot, 
    data=as.data.frame(coredata(diff(sol13102009))),
    xlab='Fluctuations of G(0) (W/m²)',
         xlim=c(-5, 5),
    auto.key = list(x = .7, y = .8,
      text=c(expression(G(0)), expression(G(0)^B))))
## dev.off()


#########################################################
###Simulation of effective irradiance time series     ###
#########################################################

mil13102009$detail <- wavMRDSum(coredata(mil13102009$Gef),
                                keep.smooth=FALSE, levels=1:9)
mil13102009$envol <- with(mil13102009, Gef-detail)
mil13102009$detailBoot <- wavBootstrap(coredata(mil13102009$detail))
mil13102009$GefBoot <- with(mil13102009, envol+detailBoot)

x <- agg(mil13102009[,c('Gef', 'GefBoot','envol', 'detail', 'detailBoot')], period=5)

## trellis.device(pdf, 
##                theme=custom.theme.gray(3),
##                file="ArticuloSimulador/Figuras/GefBoot.pdf",
##                width=6,height=4.5,
##                family='Courier')

p <- xyplot(x, superpose=TRUE, ylab='Irradiance (W/m²)',
             auto.key=list(x=0.8, y=0.5, cex=0.6,
               text=c(expression(G[ef]),
                 expression(G[ef]^B),
                 expression(G[ef]^T),
                 expression(G[ef]^D),
                 expression(G[ef]^{DB}))))
print(p)

## dev.off()

ecdfplot(~Gef+GefBoot, data=as.data.frame(diff(mil13102009)),
         auto.key = list(x = .7, y = .8, corner = c(0, 1)),
         xlab='Fluctuaction of Gef (W/m²)',
         xlim=c(-20, 20))

## Simulate fluctuactions with wavBootstrap
foo <- function(x){
  Gef <- x$Gef
  detail <- wavMRDSum(coredata(Gef),
                      keep.smooth=FALSE, levels=1:9)
  envol <- Gef-detail
  detailBoot <- wavBootstrap(coredata(detail))
  GefBoot <- envol+detailBoot
  difGefDF <- as.data.frame(diff(cbind(Gef, GefBoot))) # The result is only the fluctuations
  dia <- paste(unique(dom(index(x))), 'th', sep='')
  result <- cbind(difGefDF, dia)
}

lista <- list(mil13102009, mil14102009, mil15102009, mil16102009)
milDF <- lapply(lista, foo)
milDF <- do.call(rbind, milDF)

## trellis.device(pdf,
##                theme=custom.theme.gray(2),
##                file="ArticuloSimulador/Figuras/ecdfGefBoot.pdf",
##                width=6,height=4.5,
##                family='Courier')

ecdfplot(~Gef+GefBoot|dia, data=milDF, xlim=c(-10,10), layout=c(4,1),
                  xlab='Fluctuaction of irradiance (W/m²)',
                  auto.key = list(space='right', cex=0.7, 
                    text=c(expression(G[ef]), expression(G[ef]^B))))

## dev.off()

######################
###Wavelet variance###
######################

z131009 <-mil13102009.av[,c('Gef', 'P')]
z131009 <- scale(z131009, center=TRUE, scale=c(1000, 9.5e6))

varZ131009 <- wavVar.multi(z131009)

## trellis.device(pdf, 
##                theme=custom.theme(),
##                file="ArticuloSimulador/Figuras/VarPotencia.pdf",
##                width=6,height=4.5,
##                family='Courier')
p <- plotVar(var~scale, group=ID, data=varZ131009, cex.lab=0.7, pos.lab=1)
p
## dev.off()

##Evolution during several days
foo <- function(x){
  dat <- x[,c('Gef', 'P')]
  dat <- agg(dat, period)
  dat <- scale(dat, center=TRUE, scale=c(1000, 9.5e6))
  varZ <- wavVar.multi(dat)
  varZ$dia=paste(unique(dom(index(x))), 'th', sep='')
  varZ
}

dataList <- list(mil13102009, mil14102009, mil15102009, mil16102009)
varZ <- lapply(dataList, foo)
varZ <- do.call(rbind, varZ)
varZ$ID[varZ$ID=='Gef']='G'

## trellis.device(pdf,
##                file="ArticuloSimulador/Figuras/VarIrradianciaPotencia.pdf",
##                width=6,height=4.5,
##                family='Courier',
##                title='Varianza de la irradiancia y la potencia con intervalos de confianza')

p <- plotVar(var~scale|dia, group=ID, data=varZ, cex.lab=0.7, pos.lab=1)
print(p)
## dev.off()

##Filter
Gef <- subset(varZ, ID=='G', select=c('var', 'scale', 'dia'))
P <- subset(varZ, ID=='P', select=c('var', 'scale', 'dia'))
H <- data.frame(scale=Gef$scale)
H$var <- P$var/Gef$var
H$dia <- Gef$dia
H$high <- NA
H$low <- NA

## trellis.device(pdf,
##                file="ArticuloSimulador/Figuras/VarPlantaFiltro.pdf",
##                width=6,height=4.5,
##                family='Courier',
##                title='Filtro equivalente de una planta FV')

p <- plotVar(var~scale, groups=dia, data=H,
             ylab='Ratio between wavelet variances',
             ci='none', pos.lab=3, cex.lab=0.7)
print(p)
## dev.off()

###Filter and moving average
data2<-mil14102009
data2 <- data2[,c('Gef', 'P')]
data2 <- agg(data2, 5)
data2 <- scale(data2, center=TRUE, scale=c(1000, 9.5e6))
varZ2 <- wavVar.multi(data2)
varZ2$dia='14th Oct'
Gef <- subset(varZ2, ID=='Gef', select=c('var', 'scale', 'dia'))
P <- subset(varZ2, ID=='P', select=c('var', 'scale', 'dia'))
H <- data.frame(scale=Gef$scale)
H$var <- P$var/Gef$var
H$var <- H$var/max(H$var)
H$ID <- 'Real'

k <- 7:11
Tk <- lapply(k, function(x)rollmean(data2$Gef, x))
Tk <- do.call(cbind, Tk)
names(Tk) <- paste('T', 5*k, sep='')
xvar <- wavVar.multi(Tk)
Gefvar <- wavVar.uni(data2$Gef)
ref <- Gefvar$var
xvar$var <- xvar$var/ref
xvar$low <- NULL
xvar$high <- NULL

Ht <- rbind(xvar, H)

## trellis.device(pdf, 
##                file="ArticuloSimulador/Figuras/Filtro_MovingAverage.pdf",
##                width=6,height=4.5,
##                family='Courier',
##                title='Media movil y filtro de una planta FV')
p <- plotVar(var~scale, groups=ID, data=Ht, ci='none',
             ylab='Ratio between wavelet variances', pos.lab=2)
print(p)
## dev.off()

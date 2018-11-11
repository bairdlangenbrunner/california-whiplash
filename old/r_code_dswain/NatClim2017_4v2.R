dir = "/Users/dswain/Dropbox/Work-Home Synch/Academic_Projects/California Whiplash (NatClimChg 2017)/Data"
setwd(dir)

library(ncdf4)
library(Hmisc)
library(RColorBrewer)
library(gplots)
library(ppcor)
library(pracma)
library(plotrix)
library(graphics)
library(corrplot)
library(zoo)


#Do for 25year wet event
infile1 = nc_open("time_series_event_counts_ALL_LATLON_07262017.nc")
wet.25.ens = ncvar_get(infile1,"25yr_1season_wet_all_events_40_separate_ens")

#extract relevant gridpoint cluster
wet.25.i = wet.25.ens[9,c(17:18),,] #LA
#wet.25.i = wet.25.ens[13,c(14:15),,]  #SF
wet.25 = apply(wet.25.i,c(2,3),mean)
wet.25.r30 = apply(wet.25,1,rollsum,30)

#calculate relative change
wet.25.pic.raw = ncvar_get(infile1,"25yr_1season_wet_pic_count")
wet.25.pic = array(wet.25.pic.raw/(1798),dim(wet.25))
wet.25.pic.r30 = apply(wet.25.pic,1,rollsum,30)
wet.25.pct.chg = 100*((wet.25-wet.25.pic)/wet.25.pic)

wet.25.chg.r30 = apply(wet.25.pct.chg,1,rollmean,30)

#calculate percentiles for thresholds of interest
rcp.pctile.83.wet25 = rep(NA,151)
rcp.pctile.16.wet25 = rep(NA,151)
rcp.pctile.75.wet25 = rep(NA,151)
rcp.pctile.25.wet25 = rep(NA,151)

for (i in 1:151) {
rcp.pctile.83.wet25[i] = quantile(as.vector(wet.25.chg.r30[i,]),probs=.833333)
rcp.pctile.16.wet25[i] = quantile(as.vector(wet.25.chg.r30[i,]),probs=.166666)
rcp.pctile.75.wet25[i] = quantile(as.vector(wet.25.chg.r30[i,]),probs=.75)
rcp.pctile.25.wet25[i] = quantile(as.vector(wet.25.chg.r30[i,]),probs=.25)
}

#plot raw 30-year count timeseries for Los Angeles
cols.wetdry = brewer.pal(11,"BrBG")
cols.whiplash = brewer.pal(11,"PRGn")

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}


cols.wetdry.a = add.alpha(cols.wetdry,alpha=0.4)

#TIMESERIES PLOT V1
xyears=seq(1935,2085,10)
plot(y=as.vector(colMeans(wet.25.pct.chg)),x=(1935:2085),type="l",ylim=(c(-100,1000)),col="forestgreen",main="Seasonal wet extremes (25 yr RI)",xlab="Years",ylab="Percent change (%)",lwd=3.5,las=1,xaxt="n") #,xaxs="i"
grid()
#polygon(c(1935:2085,2085:1935),c(unlist(rcp.pctile.10[1:151]),rev(unlist(rcp.pctile.90[1:151]))),col=cols.trans[1],border=NA)
#polygon(c(1935:2085,2085:1935),c(unlist(rcp.pctile.25[1:151]),rev(unlist(rcp.pctile.75[1:151]))),col=cols.trans[1],border=NA)
boxplot(wet.25.pct.chg[151,], use.cols = TRUE,add=TRUE,at=c(2085),pars=list(boxwex=5.5),ylim=c(-100,1000),names=FALSE,axes=FALSE,outline=FALSE,col="forestgreen")

abline(h=0,lwd=2.5,lty=2,col="forestgreen")

axis(side=1,at=xyears)
grid()



##############################
wet.200.ens = ncvar_get(infile1,"200yr_40d_wet_all_events_40_separate_ens")

#extract relevant gridpoint cluster
wet.200.i = wet.200.ens[10,18,,]
wet.200 = apply(wet.200.i,c(2,3),mean)
wet.200.r30 = apply(wet.200,1,rollsum,30)

#calculate relative change
wet.200.pic.raw = ncvar_get(infile1,"200yr_40d_wet_pic_count")
wet.200.pic = array(wet.200.pic.raw/(1798),dim(wet.200.r30))
wet.200.pic.r30 = apply(wet.200.pic,1,rollsum,30)
wet.200.pct.chg = 100*((wet.200.r30-wet.200.pic)/wet.200.pic)

#calculate percentiles for thresholds of interest
rcp.pctile.83.wet200 = rep(NA,151)
rcp.pctile.16.wet200 = rep(NA,151)
rcp.pctile.75.wet200 = rep(NA,151)
rcp.pctile.25.wet200 = rep(NA,151)

for (i in 1:151) {
rcp.pctile.83.wet200[i] = quantile(as.vector(wet.200.pct.chg[i,]),probs=.83333)
rcp.pctile.16.wet200[i] = quantile(as.vector(wet.200.pct.chg[i,]),probs=.16666)
rcp.pctile.75.wet200[i] = quantile(as.vector(wet.200.pct.chg[i,]),probs=.75)
rcp.pctile.25.wet200[i] = quantile(as.vector(wet.200.pct.chg[i,]),probs=.25)
}

#TIMESERIES PLOT FOR 25YR WET
xyears=seq(1935,2085,10)
#plot(y=as.vector(rowMeans(wet.25.pct.chg)),x=(1935:2085),type="l",ylim=(c(-100,2500)),col=cols.wetdry[8],main="Seasonal wet extremes (25 yr RI), Los Angeles",xlab="Years",ylab="Percent change (%)",lwd=3.5,las=1,xaxt="n",,xaxs="i")
#lines(y=as.vector(rowMeans(wet.25.pct.chg)),x=(1935:2085),type="l",ylim=(c(-100,2500)),col=cols.wetdry[9],lwd=3.5,las=1,xaxt="n")
#grid()
#points(x=rep(2085,40),y=wet.25.pct.chg[151,],col=cols.wetdry.a[9])
#points(x=rep(2085,40),y=wet.200.pct.chg[151,],col=cols.wetdry.a[8])
#polygon(c(1935:2085,2085:1935),c(unlist(rcp.pctile.16.wet25[1:151]),rev(unlist(rcp.pctile.83.wet25[1:151]))),col=cols.wetdry.a[8],border=NA)
#polygon(c(1935:2085,2085:1935),c(unlist(rcp.pctile.16.wet200[1:151]),rev(unlist(rcp.pctile.75.wet25[1:151]))),col=cols.wetdry.a[9],border=NA)

plot(y=as.vector(rowMeans(wet.25.chg.r30)),x=(1935:2085),type="l",ylim=(c(-100,600)),col=cols.wetdry[10],main="Seasonal wet extremes (25 yr RI), Los Angeles",xlab="Years",ylab="Change in Frequency (%)",lwd=3.5,las=1,xaxt="n",xaxs="i")
polygon(c(1935:2085,2085:1935),c(unlist(rcp.pctile.16.wet25[1:151]),rev(unlist(rcp.pctile.83.wet25[1:151]))),col=cols.wetdry.a[10],border=NA)
lines(y=as.vector(rowMeans(wet.25.pct.chg)),x=(1935:2085),ylim=(c(-100,500)),col=cols.wetdry[10],lwd=3.5,las=1,xaxt="n",,xaxs="i")

axis(side=1,at=xyears)
grid()
abline(h=0,lwd=2.5,lty=2,col="black")


#######################
#DO FOR 100YR DRY
wet.200.ens = ncvar_get(infile1,"100yr_1season_dry_all_events_40_separate_ens")

#extract relevant gridpoint cluster
#wet.200.i = wet.25.ens[9,c(17:18),,] #LA
wet.200.i = wet.25.ens[13,c(14:15),,]  #SF
wet.200 = apply(wet.200.i,c(2,3),mean)
wet.200.r30 = apply(wet.200,1,rollmean,30)

#calculate relative change
wet.200.pic.raw = ncvar_get(infile1,"100yr_1season_pic_count")
wet.200.pic = array(wet.200.pic.raw/(1798),dim(wet.200.r30))
wet.200.pic.r30 = apply(wet.200.pic,1,rollmean,30)
wet.200.pct.chg = 100*((wet.200.r30-wet.200.pic)/wet.200.pic)

#calculate percentiles for thresholds of interest
rcp.pctile.83.wet200 = rep(NA,151)
rcp.pctile.16.wet200 = rep(NA,151)
rcp.pctile.75.wet200 = rep(NA,151)
rcp.pctile.25.wet200 = rep(NA,151)

for (i in 1:151) {
rcp.pctile.83.wet200[i] = quantile(as.vector(wet.200.pct.chg[i,]),probs=.83333)
rcp.pctile.16.wet200[i] = quantile(as.vector(wet.200.pct.chg[i,]),probs=.16666)
rcp.pctile.75.wet200[i] = quantile(as.vector(wet.200.pct.chg[i,]),probs=.75)
rcp.pctile.25.wet200[i] = quantile(as.vector(wet.200.pct.chg[i,]),probs=.25)
}

#TIMESERIES PLOT V2
xyears=seq(1935,2085,10)
#plot(y=as.vector(rowMeans(wet.200.pct.chg)),x=(1935:2085),type="l",ylim=(c(-100,2500)),col=cols.wetdry[8],main="Seasonal wet extremes (25 yr RI)",xlab="Years",ylab="Percent change (%)",lwd=3.5,las=1,xaxt="n",,xaxs="i")
#lines(y=as.vector(rowMeans(wet.25.pct.chg)),x=(1935:2085),type="l",ylim=(c(-100,2500)),col=cols.wetdry[9],lwd=3.5,las=1,xaxt="n")
#grid()
#polygon(c(1935:2085,2085:1935),c(unlist(rcp.pctile.16.wet200[1:151]),rev(unlist(rcp.pctile.83.wet200[1:151]))),col=cols.wetdry.a[8],border=NA)
#polygon(c(1935:2085,2085:1935),c(unlist(rcp.pctile.16.wet200[1:151]),rev(unlist(rcp.pctile.75.wet25[1:151]))),col=cols.wetdry.a[9],border=NA)

plot(y=as.vector(rowMeans(wet.200.pct.chg)),x=(1935:2085),type="l",ylim=(c(-100,2500)),col=cols.wetdry[2],main="Seasonal dry extremes (100 yr RI), San Francisco",xlab="Years",ylab="Change in Frequency (%)",lwd=3.5,las=1,xaxt="n",xaxs="i")
polygon(c(1935:2085,2085:1935),c(unlist(rcp.pctile.16.wet200[1:151]),rev(unlist(rcp.pctile.83.wet200[1:151]))),col=cols.wetdry.a[3],border=NA)
lines(y=as.vector(rowMeans(wet.200.pct.chg)),x=(1935:2085),ylim=(c(-100,2500)),col=cols.wetdry[2],lwd=3.5,las=1,xaxt="n",,xaxs="i")

axis(side=1,at=xyears)
grid()
abline(h=0,lwd=2.5,lty=2,col="black")



###################################
#CUMULATIVE OCCURENCE OF 200YR EVENT TIMESERES
wet.200.ens = ncvar_get(infile1,"200yr_40d_wet_all_events_40_separate_ens")

#extract relevant gridpoint cluster
wet.200.i = wet.25.ens[9,c(17:18),,97:180] #LA
#wet.200.i = wet.25.ens[13,c(14:15),,]  #SF
wet.200 = apply(wet.200.i,c(2,3),mean)
wet.200.r30 = apply(wet.200,1,cumsum)

#calculate relative change
wet.200.pic.raw = ncvar_get(infile1,"200yr_40d_wet_pic_count")
wet.200.pic = array(wet.200.pic.raw/(1798),dim(wet.200.r30))
wet.200.pic.r30 = apply(wet.200.pic,2,cumsum)
wet.200.pct.chg = 100*((wet.200.r30-wet.200.pic)/wet.200.pic)

#calculate percentiles for thresholds of interest
rcp.pctile.83.wet200 = rep(NA,84)
rcp.pctile.16.wet200 = rep(NA,84)
rcp.pctile.75.wet200 = rep(NA,84)
rcp.pctile.25.wet200 = rep(NA,84)

for (i in 1:84) {
rcp.pctile.83.wet200[i] = quantile(as.vector(wet.200.r30[i,]),probs=.83333)
rcp.pctile.16.wet200[i] = quantile(as.vector(wet.200.r30[i,]),probs=.16666)
rcp.pctile.75.wet200[i] = quantile(as.vector(wet.200.r30[i,]),probs=.75)
rcp.pctile.25.wet200[i] = quantile(as.vector(wet.200.r30[i,]),probs=.25)
}

xyears=seq(2020,2100,10)
yvals = seq(0,15,2.5)
plot(y=as.vector(rowMeans(wet.200.r30)),x=(2017:2100),type="l",ylim=(c(0,15)),col=cols.wetdry[2],main="Cumulative occurrence, 40-day wet events (200 yr RI), Los Angeles",xlab="Years",ylab="Number of events",lwd=3.5,las=1,xaxt="n",yaxt="n",xaxs="i")
polygon(c(2017:2100,2100:2017),c(unlist(rcp.pctile.16.wet200[1:84]),rev(unlist(rcp.pctile.83.wet200[1:84]))),col=cols.wetdry.a[3],border=NA)
lines(y=as.vector(rowMeans(wet.200.pic.r30)),x=(2017:2100),ylim=(c(0,15)),col="black",lwd=2.5,las=1,xaxt="n",,xaxs="i",lty=2)
lines(y=as.vector(rowMeans(wet.200.r30)),x=(2017:2100),ylim=(c(0,15)),col=cols.wetdry[2],lwd=3.5,las=1,xaxt="n",,xaxs="i")


axis(side=1,at=xyears)
axis(side=2,at=yvals,las=2)

grid()
#abline(h=0,lwd=2.5,lty=2,col="black")



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


#READ IN FILE DATA AND RELEVANT VARIABLES
#25year wet event
infile         = nc_open("time_series_event_counts_ALL_LATLON.nc")
wet.25.ens     = ncvar_get(infile,"25yr_1season_wet_all_events_40_separate_ens")
wet.25.pi.raw  = ncvar_get(infile,"25yr_1season_wet_pic_count")

#extract relevant gridpoint corresponding to Los Angeles from RCP ensemble (40x180)
wet.25 = wet.25.ens[9,17,,]

#convert PI to correct units (equivalent counts per ensemble member per year)
wet.25.pi.i = wet.25.pi.raw/1798

#convert PI to array of correct size for matrix division (40x180)
wet.25.pi = array(wet.25.pi.i,dim(wet.25))

#calculate relative change, RCP vs PI (40x180)
wet.25.pct.chg = 100*((wet.25-wet.25.pi)/wet.25.pi)

#plot ensemble mean without 30year smoothing
plot(colMeans(wet.25.pct.chg),type="l")

#calculate 30 year rolling mean (centered) (151x1)
wet.25.ensmean.roll30 = rollmean(colMeans(wet.25.pct.chg),30)

#plot ensemble mean without and with 30 year smoothing
plot(x=c(1921:2100),y=colMeans(wet.25.pct.chg),type="l",col="blue",main="Raw (blue) and 30y smoothed (red) Percent Change, Los Angeles (9,17 gridpoint)")
lines(y=wet.25.ensmean.roll30,x=c(1935:2085),type="l",col="red")

#test plot all smoothed ensemble members
wet.25.pct.chg.ens.rollmean30 = apply(wet.25.pct.chg, 1, rollmean, 30)
matplot(wet.25.pct.chg.ens.rollmean30,type="l",col="black",lty=1)



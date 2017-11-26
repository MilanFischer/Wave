###########################################################################################
#                                                                                         #
# WAVETEST Example R script for WAVELET, using NINO3 SST dataset                          #
#                                                                                         #
# See 'http://paos.colorado.edu/research/wavelets/'                                       #
# The original script written in Matlab in January 1998 by C. Torrence                    #
#                                                                                         #
# Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,                #
#   changed all 'log' to 'log2', changed logarithmic axis on GWS to a normal axis.        #
#                                                                                         #
#                                                                                         #
# Translated to R by Milan Fischer in October 2016 (fischer.milan[at]gmail[dot]com)       #
#                                                                                         #
###########################################################################################

	rm(list=ls())
	
	library(pracma)
	library(matlab)
	library(signal)
	library(graphics)

	#Functions to be used
	source('Wave_bases.R')
	source('Wavelet.R')
	source('Wave_signif.R')
	source('Chisquare_inv.R')
	source('Chisquare_solve.R')

	#Move one directory up
	setwd('..')

data<-read.table(paste(getwd(),'/Data/sst_nino3.dat',sep=''),header=FALSE)		# input SST time series
sst=data$V1

#------------------------------------------------------ Computation

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# 'http://paos.colorado.edu/research/wavelets/plot/'

variance=sd(sst)^2
sst=(sst-mean(sst))/sqrt(variance) 

n=length(sst)
dt=0.25
time=seq(from=1871,by=dt,
length.out=length(sst))							# construct time array
xlim=c(1870,2000)							# plotting range
pad=1									# pad the time series with zeroes (recommended)
dj=0.25									# this will do 4 sub-octaves per octave
s0=2*dt									# this says start at a scale of 6 months
j1=7/dj									# this says do 7 powers-of-two with dj sub-octaves each
lag1=0.72								# lag-1 autocorrelation for red noise background
mother='MORLET'

# Wavelet transform:
out=wavelet(sst,dt,pad,dj,s0,j1,mother)
wave=out$wave
period=out$period
scale=out$scale
coi=out$coi

power=(abs(wave))^2							# compute wavelet power spectrum 

# Significance levels: (variance=1 for the normalized SST)
out=wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother)
signif=out$signif
fft_theor=out$fft_theor

sig95=matrix(signif,length(signif),n)					# expand signif --> (J+1)x(N) array
sig95=power/sig95							# where ratio > 1, power is significant

# Global wavelet spectrum & significance levels:
global_ws=variance*(sum(t(power))/n)					# time-average over all times
dof=n-scale								# the -scale corrects for padding at edges
global_signif=wave_signif(variance,dt,scale,1,lag1,-1,dof,mother)$signif

# Scale-average between El Nino periods of 2--8 years
avg=which((scale>=2)&(scale<8))
Cdelta=0.776								# this is for the MORLET wavelet
scale_avg=matrix(scale,length(scale),n)					# expand scale --> (J+1)x(N) array
scale_avg=power/scale_avg						# [Eqn(24)]
scale_avg=variance*dj*dt/Cdelta*sum(scale_avg[avg,])			#[Eqn(24)]
scaleavg_signif=wave_signif(variance,dt=dt,scale=scale,sigtest=2,lag1=lag1,siglvl=0.95,dof=c(2,7.9),mother=mother,param=-1)$signif


# Print the list of available variables
ls()


#------------------------------------------------------ Plotting

tiff(paste(getwd(),'/Output/sst_nino3_',mother,'.tiff',sep=''),width=350,height=220,units='mm',res=200,type='cairo')

# or
# dev.new()

#--- Plot time series
suppressWarnings(par(fig=c(0.02,0.75,0.6,1),new=TRUE))
plot(time,sst,type='l',xlab=NA,ylab=NA,xaxs='i',yaxs='r',ylim=c(-4,4),las=1)
title('a) NINO3 Sea Surface Temperature (seasonal)',line=1,cex.main=0.9)
title(ylab='NINO3 SST (°C)',cex.lab=0.9)


#--- Contour plot wavelet power spectrum
par(fig=c(0.02,0.75,0.3,0.7),new=TRUE)
levels=c(0.0625,0.125,0.25,0.5,1,2,4,8,16) 
hmcols<-colorRampPalette(c('blue','green','yellow','orange','red','brown','black'))(length(levels))
contour(time,log2(period),t(log2(power)),drawlabels=FALSE,xaxs='i',yaxs='i',
levels=log2(levels),ylim=log2(c(70,0.5)),col=hmcols,lwd=2,
xlab=NA,ylab=NA,yaxt='n')
title('b) NINO3 SST Wavelet Power Spectrum',line=1,cex.main=0.9)
title(ylab='Period (years)',cex.lab=0.9)
lab<-c(1,2,4,8,16,32,64)
axis(2,at=c(0,1,2,3,4,5,6),labels=lab,las=1)
#*** or use 'filled.contour'

#--- add 95% significance contour
par(fig=c(0.02,0.75,0.3,0.7),new=TRUE)
contour(time,log2(period),t(sig95),drawlabels=FALSE,xaxs='i',yaxs='i',
levels=log2(levels),ylim=log2(c(70,0.5)),col='black',lwd=2.2,
xlab=NA,ylab=NA,xaxt='n',yaxt='n')
														
#--- add cone-of-influence, anything 'below' is dubious
lines(time,log2(coi),col='black',lwd=2,lty=2)


#--- Plot global wavelet spectrum
par(fig=c(0.7,1,0.3,0.7),new=TRUE)
plot(global_ws,log2(period),ylim=log2(c(70,0.5)),xlim=c(0,3),type='l',xlab=NA,ylab=NA,xaxs='r',yaxs='i',yaxt='n')
title('c) Global Wavelet Spectrum',line=1,cex.main=0.9)
axis(2,at=c(0,1,2,3,4,5,6),labels=lab,las=1)
lines(global_signif,log2(period),lty=2)
title(xlab=expression('Power '~(degree*C^{2})))


#--- Plot 2--8 yr scale-average time series
par(fig=c(0.02,0.75,0,0.4),new=TRUE)
plot(time,scale_avg,type='l',xlab=NA,ylab=NA,xaxs='i',yaxs='r',las=1)
abline(h=scaleavg_signif,lty=2)
title(main='d) 2-8 yr Scale-average Time Series',line=1,cex.main=0.9)
title(xlab='Time (year)',cex.lab=0.9)
title(ylab=expression('Avg variance '~(degree*C^{2})),line=2.5)

dev.off()
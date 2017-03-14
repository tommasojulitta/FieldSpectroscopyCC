\name{RadiometricCalibration}
\alias{RadiometricCalibration}
\title{RadiometricCalibration}
\description{Compute the VIS-NIR radiometric calibration coefficients based on outdoor measurements }
\usage{RadiometricCalibration(type, wl, radiance, DN, ReferenceCoeff)}
\arguments{
  \item{type}{value: 1 if the same white panel has been used; 2 if the spectrometer connected to the cosine receptor has to be cross calibrated with measurements collected over  a white panel}
  \item{wl}{numeric vector: wavelength vector}
  \item{radiance}{numeric vector: mean vector of n radiance spectra acquired}
  \item{DN}{numeric vector: mean vector of n spectra in Digital Number
acquired over the same reference standard during the same time interval 
(DN spectrum has to be already dark current subtracted and divided by integration time used) }
  \item{ReferenceCoeff}{numeric numeric vector:
vector fo the calibration coefficients of the refrence standard used. }
}
\details{Outdoor measurements are acquired simoultaneously over the same reference standard (e.g. white panel) 
or over different reference standard (e.g. white panel and cosine receptor) 
by a reference spectrometer and the spectrometer to be calibrated.
All the provided data in input need to be pre-processed in order to be homogenized in terms of wavelength
(see \code{\link{LinearResample}})}
\value{numeric vector: conversion factor from digital counts to radiance for each wavelength }

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

#
#
#for Indoor radiometric claibration
#
#

data("indoor_rad_cal_data")
integration_time<-250
#Remove DC fro spectra
sub_dat_m<-DCSubtraction(signal=indoor_rad_cal_data$data_lamp,DarkSignal = indoor_rad_cal_data$data_dc,type=1)
#Create matrix for radiometric calibration
DN_mat<-DNSpectralMatrixRadCal(spectra = sub_dat_m,IntegrationTime = integration_time)
#calculate mean of several spectra
DN_mean<-StatsOnSpectra(wl=indoor_rad_cal_data$DN_wl,spectra=DN_mat,fun='mean')
#linear resample at refrence radiance wavelength 
radiance_lamp_res<-LinearResample(indoor_rad_cal_data$lamp_radiance$Wavelength,indoor_rad_cal_data$lamp_radiance$rad,indoor_rad_cal_data$DN_wl)
#calculate calibration coefficients
rad_cal<-RadiometricCalibration(type=1,wl=indoor_rad_cal_data$DN_wl,radiance = radiance_lamp_res,DN = DN_mean)
#exclude regions of the spectrum, optically black pixels
range_to_exclude<-data.frame(indoor_rad_cal_data$DN_wl[1],indoor_rad_cal_data$DN_wl[30])
range_to_exclude<-rbind(range_to_exclude,c(indoor_rad_cal_data$DN_wl[length(indoor_rad_cal_data$DN_wl)-35],indoor_rad_cal_data$DN_wl[length(indoor_rad_cal_data$DN_wl)]))
#exclude noisy pixels
excluded_noisy_pixels<-ExcludeSpectralRegions(wl=indoor_rad_cal_data$DN_wl,spectrum = rad_cal,SpectralRegionToExclude = range_to_exclude)
rad_cal_coeff<-SplineSmoothGapfilling(wl=indoor_rad_cal_data$DN_wl,excluded_noisy_pixels) 
#plot results
x11()
plot(indoor_rad_cal_data$DN_wl,rad_cal,ylim=c(0,0.02),xlab="WL [nm]",ylab="Conversion coefficients")
points(indoor_rad_cal_data$DN_wl,excluded_noisy_pixels,pch=20,col="red")
lines(indoor_rad_cal_data$DN_wl,rad_cal_coeff,col="green",lwd=2)

#
#
#for Outdoor radiometric claibration
#
#

data("outdoor_rad_cal_data")
data("atmospheric_absorption_regions")
integration_time<-450
#Create matrix for radiometric calibration
DN_mat<-DNSpectralMatrixRadCal(spectra = outdoor_rad_cal_data$DN_matrix,IntegrationTime = integration_time)
#calculate mean of several spectra
radiance_mean<-StatsOnSpectra(wl=outdoor_rad_cal_data$radiance_wl,spectra=outdoor_rad_cal_data$radiance_matrix,fun='mean')
#linear resample at refrence radiance wavelength 
radiance_mean_res<-LinearResample(wl = outdoor_rad_cal_data$radiance_wl,spectrum = radiance_mean,wlRef = outdoor_rad_cal_data$DN_wl)
wp_coeff_res<-LinearResample(outdoor_rad_cal_data$wp_coef$V1,outdoor_rad_cal_data$wp_coef$V2,outdoor_rad_cal_data$DN_wl)
#calculate mean of several spectra
DN_mean<-StatsOnSpectra(wl=outdoor_rad_cal_data$DN_wl,spectra=outdoor_rad_cal_data$DN_matrix,fun='mean')
#calculate calibration coefficients
rad_cal<-RadiometricCalibration(type=1,wl=outdoor_rad_cal_data$DN_wl,radiance = radiance_mean_res,DN = DN_mean)
#exclude regions of the spectrum affected by atmospheric absorptions and noisy pixels
range_to_exclude<-data.frame(wl_start=c(outdoor_rad_cal_data$DN_wl[1],outdoor_rad_cal_data$DN_wl[length(outdoor_rad_cal_data$DN_wl)-35]),wl_end=c(outdoor_rad_cal_data$DN_wl[30],outdoor_rad_cal_data$DN_wl[length(outdoor_rad_cal_data$DN_wl)]))
atmospheric_absorption_regions<-rbind(atmospheric_absorption_regions,c(range_to_exclude))
exclude_atmospheric_absorption_features<-ExcludeSpectralRegions(wl=outdoor_rad_cal_data$DN_wl, spectrum = rad_cal,SpectralRegion = atmospheric_absorption_regions)
#smooth results
rad_cal_coeff<-SplineSmoothGapfilling(wl=outdoor_rad_cal_data$DN_wl,spectrum = exclude_atmospheric_absorption_features) 
#plot results
x11()
plot(outdoor_rad_cal_data$DN_wl,rad_cal,ylim=c(0,0.000025),xlab="WL [nm]",ylab="Conversion coefficients")
points(outdoor_rad_cal_data$DN_wl,exclude_atmospheric_absorption_features,pch=20,col="red")
lines(outdoor_rad_cal_data$DN_wl,rad_cal_coeff,col="green",lwd=2)
}

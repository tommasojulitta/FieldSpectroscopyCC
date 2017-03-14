#packageDir <- system.file(package='FieldSpectroscopyCC')	# only works if package has been installed
#dataDir <- file.path(packageDir,"data") 

.tmp.f <- function(){
	source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyCC\\R\\AdvancedFunctions.R")
	source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyCC\\R\\BaseFunctions.R")
	source("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyCC\\R\\IOFunctions.R")
	setwd("C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyCC\\data\\")
	dataDir <- "C:\\Users\\Tommi\\Documents\\R_workspace\\FieldSpectroscopyCC\\data\\"
	return(dataDir)
}
dataDir<-.tmp.f()
#---------------------------------------------------------------------------------------------------------------------------------------
########################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------------
##################################################TEST ON SNR DETERMINATION#############################################################
load(file.path(dataDir,"snr_data.Rdata"))
#perform statistics on spectra
specStats <- SummarizeSpectra(snr_data)
SNR<-SignalToNoiseRatio(specStats)

lamp_mean<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='mean')
lamp_sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='sd')
dc_mean<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_dc,fun='mean')
dc_sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_dc,fun='sd')
#calculate Signal to Noise Ratio
SNR<-SignalToNoiseRatio(average_light=lamp_mean,sd_light=lamp_sd,average_dark=dc_mean,sd_dark=dc_sd)
#plot results
x11();plot(snr_data$wl,SNR,type="l")
#---------------------------------------------------------------------------------------------------------------------------------------
##################################################TEST ON OUTDOOR RAD CAL###############################################################
#---------------------------------------------------------------------------------------------------------------------------------------
load(file.path(dataDir,"outdoor_rad_cal_data.Rdata"))
load(file.path(dataDir,"atmospheric_absorption_regions.Rdata"))
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
plot(outdoor_rad_cal_data$DN_wl,rad_cal,ylim=c(0,0.000025))
points(outdoor_rad_cal_data$DN_wl,exclude_atmospheric_absorption_features,pch=20,col="red")
lines(outdoor_rad_cal_data$DN_wl,rad_cal_coeff,col="green",lwd=2)
#---------------------------------------------------------------------------------------------------------------------------------------
##################################################TEST ON INDOOR RAD CAL###############################################################
#---------------------------------------------------------------------------------------------------------------------------------------
load("indoor_rad_cal_data.Rdata")
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
plot(indoor_rad_cal_data$DN_wl,rad_cal,ylim=c(0,0.02))#,xlim=c(650,700))
points(indoor_rad_cal_data$DN_wl,excluded_noisy_pixels,pch=20,col="red")
lines(indoor_rad_cal_data$DN_wl,rad_cal_coeff,col="green",lwd=2)
#---------------------------------------------------------------------------------------------------------------------------------------
##################################################TEST ON INDOOR WL CAL###############################################################
#---------------------------------------------------------------------------------------------------------------------------------------
load("indoor_wl_cal_data.Rdata")
lamp_spectrum_dc_sub<-DCSubtraction(signal=indoor_wl_cal_data$lamp_spectrum,DarkSignal = indoor_wl_cal_data$dc_spectrum,type=1)
region_to_analyze<-SelectSpectralRegion(wl = indoor_wl_cal_data$DN_wl,spectrum = lamp_spectrum_dc_sub,WlSelection = indoor_wl_cal_data$emission_lines$peak,buffer=3)
n_peaks_fit<-length(names(region_to_analyze))
wl_peaks<-as.numeric(names(region_to_analyze))

for(n in 1:n_peaks_fit)
{print(n)
  n_pixels<-data.frame(region_to_analyze[n])[,1]
  wl<-data.frame(region_to_analyze[n])[,2]
  DN<-data.frame(region_to_analyze[n])[,3]  
  wl_param<-GaussFit(n_pixels,wl,DN,plot=TRUE)
  if(n==1){wl_cal<-wl_param}else{
  wl_cal<-rbind(wl_cal,wl_param)}
}
#extarct coefficients for wl calibration
wl_coeff<-WlCal(pix_center = wl_cal$CenterPixel,wl_peaks = wl_peaks)
wl_cal<-predict(wl_coeff, data.frame(x=1:length(indoor_wl_cal_data$DN_wl)))
#plot results
x11()
plot(indoor_wl_cal_data$DN_wl,lamp_spectrum_dc_sub,xlab="wl",ylab="DN",type="l",cex.lab=1.5,cex=2)
lines(wl_cal,lamp_spectrum_dc_sub,col="red")
legend("topright",col=c("black","red"),lty=1,cex=1.2,legend=c("before calibration","after calibration"),box.col="white")
box()
#####################################################################################################################################################################

StatsOnSpectra<-function(
  ### Compute statistic across several spectra subset to a specified spectral range
  wl					##<< numeric vector: wavelength for each row in spectra 
  ,wlStart=wl[1]  		##<< numeric value or vector: first wavelength of the spectral range selection. Default value: first wavelength of wl vector
  ,wlEnd=wl[length(wl)]  ##<< numeric value or vector: last wavelength of the spectral range selection. Default value: last wavelength of wl vector
  ,spectra=spectra		##<< numeric matrix or data.frame: collection of several spectra acquired by several columns (Digital numbers, Radiance, Reflectance).
  ,fun='mean' 			##<< character or function: function to be applied on each row of the data.frame. Default function: mean
  ,margin=1  			##<< numeric value: a vector giving the subscripts which the function will be applied over. 1 indicates rows, 2 indicates columns
)
{
  
  if(length(wlStart)>1&length(wlEnd)>1)
  {
    ### In case vector of wlStart and wl end
    stat_spectrum<-numeric(length(wlStart))
    for( n in 1:length(wlStart)){
      spectral_subset<-which(wl>=wlStart[n]&wl<=wlEnd[n])
      spectra_subset<-data.frame(spectra[spectral_subset,])
      stat<-apply(spectra_subset,margin,fun,na.rm=TRUE)
      stat_spectrum<-as.numeric(stat_spectrum)
      stat_spectrum[n]<-stat[n]
    }
    ### loop on wl start and end
  }else{
    ### In case wlStart and wl end are single value
    spectral_subset<-which(wl>=wlStart&wl<=wlEnd)
    spectra_subset<-data.frame(spectra[spectral_subset,])
    stat_spectrum<-apply(spectra_subset,margin,fun,na.rm=TRUE)
    stat_spectrum<-as.numeric(stat_spectrum)
  }
  ##value<< numeric data.frame containing the computed statistic across the selected spectral range.
  
  return(stat_spectrum)
}


attr(StatsOnSpectra,"ex") <- function(){
  
    data("snr_data")
    #perform statistics on spectra
    #calculate mean of the  signal between 700 and 800 nm
    mean<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='mean',wlStart=700,wlEnd = 800)
    #calculate stadard deviation of the signal between 700 and 800 nm
    sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='sd',wlStart=700,wlEnd = 800)
}


#####################################################################################################################################################################

# checkSpecraFormat <- function(
#   data
# ){
#   ##TODO check whether wl, data_lamp, and data_dc conform 
#   return(TRUE)
# }
# attr(StatsOnSpectra,"ex") <- function(){
#   # This example appears in the generated help.
#   load(file.path(dataDir,"snr_data.Rdata"))
#   MeanLampSpectrum<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='mean')
# }

#####################################################################################################################################################################

SummarizeSpectra <- function(
  ### compute common summary statistics on light and dark field spectra 
  spectra	##<< list with entries wl (wavelength), data_lamp (matrix of several light spectra), and data_dc (matrix of several dark current)
  ,...	##<< further arguments to \code{\link{StatsOnSpectra}}, such as \code{wlStart} and \code{wlEnd} 	
){
  lamp_mean<-StatsOnSpectra(wl=spectra$wl,spectra=spectra$data_lamp,fun=mean,...)
  lamp_sd<-StatsOnSpectra(wl=spectra$wl,spectra=spectra$data_lamp,fun=sd,...)
  dc_mean<-StatsOnSpectra(wl=spectra$wl,spectra=spectra$data_dc,fun=mean,...)
  dc_sd<-StatsOnSpectra(wl=spectra$wl,spectra=spectra$data_dc,fun=sd,...)
  ##value<< data.frame with entries
  ans <- data.frame(
    meanLamp=lamp_mean			##<< mean light spectrum  
    ,sdLamp=lamp_sd			##<< standard deviation across light spectra
    ,nLamp=ncol(spectra$data_lamp)##<< number of summarized light spectra
    ,meanDC=dc_mean			##<< mean dark spectrum
    ,sdDC=dc_sd				##<< standard deviation across dark spectra
    ,nDar=ncol(spectra$data_dc)	##<< number of summarized dark spectra
  )
}


attr(SummarizeSpectra,"ex") <- function(){
  

    data("snr_data")
    #perform statistics on spectra
    specStats <- SummarizeSpectra(snr_data)
    #print header of specStats
    head(specStats)
}



#####################################################################################################################################################################

LinearResample<-function(
  ### Compute linear resample of a spectrum based on specific wavelength of a defined spectrometer
  wl 			##<< numeric vector: wavelength vector to be resamped
  ,spectrum 	##<< numeric vector: spectrum to be resampled (e.g. radiance, reflectance, digital numbers)
  ,wlRef 		##<< numeric vector: reference wavelength vector to be be used for resampling
){
  SpectrumResampled<-approx(wl,spectrum,wlRef);SpectrumResampled<-SpectrumResampled$y
  ### Compute linear resample
  ##value<< numeric vector containing the spectrum linearly resample according to the wavelngth vector wlRef.
  return(SpectrumResampled)
}

attr(LinearResample,"ex") <- function(){
  

    data("outdoor_rad_cal_data")
    data("atmospheric_absorption_regions")
    integration_time<-450
    #Create matrix for radiometric calibration
    DN_mat<-DNSpectralMatrixRadCal(spectra = outdoor_rad_cal_data$DN_matrix,IntegrationTime = integration_time)
    #calculate mean of several spectra
    radiance_mean<-StatsOnSpectra(wl=outdoor_rad_cal_data$radiance_wl,spectra=outdoor_rad_cal_data$radiance_matrix,fun='mean')
    #linear resample at refrence radiance wavelength 
    radiance_meanRes<-LinearResample(wl = outdoor_rad_cal_data$radiance_wl,spectrum = radiance_mean,wlRef = outdoor_rad_cal_data$DN_wl)
    #plot results
    x11()
    plot(outdoor_rad_cal_data$radiance_wl,radiance_mean,type="l",xlab="WL [nm]",ylab=expression("Radiance [W m"^-2* "sr"^-1* "nm"^-1*"]"),ylim=c(0,5))
    lines(outdoor_rad_cal_data$DN_wl,radiance_meanRes,col="blue")
    
}


#####################################################################################################################################################################

SplineSmoothGapfilling<-function(
  ### Apply a spline smoothing and gap filling in the region where NA are found
  wl  			##<< numeric vector: wavelength vector
  ,spectrum  	##<< numeric data.frame: First column wavelength vector, second column spectrum vector
  ,df = length(wl)/50 ##<< numeric value: the desired equivalent number of degrees of freedom (trace of the smoother matrix)
)
{
  data<-data.frame(wl,spectrum)
  temp<-na.omit(data)
  spline<-smooth.spline(temp$wl,temp$spectrum,df=df)
  spline_gapfilling<-predict(spline, data$wl);spline_gapfilling<-spline_gapfilling$y
  ##value<< numeric vector containing the smoothed spectrum.
  return(spline_gapfilling)
}


attr(SplineSmoothGapfilling,"ex") <- function(){
  

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
    x11()
    plot(outdoor_rad_cal_data$DN_wl,rad_cal,ylim=c(0,0.000025),xlab="WL [nm]",ylab="Conversion coefficients")
    points(outdoor_rad_cal_data$DN_wl,exclude_atmospheric_absorption_features,pch=20,col="red")
    lines(outdoor_rad_cal_data$DN_wl,rad_cal_coeff,col="green",lwd=2)
    
}

#####################################################################################################################################################################

PolSmoothGapfilling<-function(
  ### Apply a polynomial fitting and gap filling in the region where NA are found
  wl  ##<< numeric vector: wavelength vector
  ,spectrum  ##<< numeric data.frame: First column wavelength vector, second column spectrum vector
  ,deg = 4 ##<< numeric value: the degree of the polynomial
)
{
  data<-data.frame(wl,spectrum)
  temp<-na.omit(data)
  y<-temp[,2]
  x<-temp[,1]
  fit <- lm(y ~ poly(x, deg, raw=TRUE))
  predpol<-predict(fit,data.frame(x=wl))
  ##value<< numeric vector containing the smoothed spectrum.
  return(predpol)
}



attr(PolSmoothGapfilling,"ex") <- function(){
  

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
    rad_cal_coeff<-PolSmoothGapfilling(wl=outdoor_rad_cal_data$DN_wl,spectrum = exclude_atmospheric_absorption_features) 
    x11()
    plot(outdoor_rad_cal_data$DN_wl,rad_cal,ylim=c(0,0.000025),xlab="WL [nm]",ylab="Conversion coefficients")
    points(outdoor_rad_cal_data$DN_wl,exclude_atmospheric_absorption_features,pch=20,col="red")
    lines(outdoor_rad_cal_data$DN_wl,rad_cal_coeff,col="green",lwd=2)
    
}
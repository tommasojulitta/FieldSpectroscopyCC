SignalToNoiseRatio<-function(
  ### Compute the Signal to Noise Ratio of a spectrometer
  specStats		##<< alternative parameterization data.frame with columns meanLight, sdLight, meanDark, sdDark
  ,avgLamp=specStats$meanLamp	##<< numeric vector: mean spectrum of n scans of stable lamp
  ,sdLamp=specStats$sdLamp  		##<< numeric vector: sd of the spectra acquired with a stable lamp
  ,avgDC=specStats$meanDC  	##<< numeric vector: mean spectrum of dark current measurements acquired with the same integration time as the measurements over the lamp
  ,sdDC=specStats$sdDC 			##<< numeric vector: sd of the dark current spectra acquired with with the same integration time as the measurements over the lamp
  )
{
  signal<-avgLamp-avgDC
  # The signal substracting the dark current
  noise<-sqrt((sdLamp^2)+(sdDC^2))
  # The noise, based on standard deviation of lamp signal and dark current
  spectralSNR<-signal/noise
  ##value<< numeric vector: ratio of signal and noise for each wavelength
  return(spectralSNR)
}


attr(SignalToNoiseRatio,"ex") <- function(){
  

    data("snr_data")
    #perform statistics on spectra
    lamp_mean<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='mean')
    #calculate mean of the lamp signal
    lamp_sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='sd')
    #calculate standard deviation of the lamp signal
    dc_mean<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_dc,fun='mean')
    #calculate mean of the dark current signal
    dc_sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_dc,fun='sd')
    #calculate standard deviation of the dark current signal
    SNR<-SignalToNoiseRatio(avgLamp =lamp_mean,sdLamp = lamp_sd,avgDC = dc_mean,sdDC = dc_sd)
    #plot results
    x11();plot(snr_data$wl,SNR,type="l",ylab="SNR",xlab="WL [nm]")
}


#####################################################################################################################################################################

NoiseEquivalentRadiance<-function(
  ### Compute the Signal to Noise Ratio of a spectrometer
  sdLamp=specStats$sdLamp  		##<< numeric vector: sd of the spectra acquired with a stable lamp
  ,sdDC=specStats$sdDC##<< numeric vector: sd of the dark current spectra acquired with with the same integration time as the measurements over the lamp
  ,RadCalCoeff  ##<< numeric vector: wavelength dependent vector of coefficient for calibration
  ,IntegrationTime  ##<< numeric value: integration time used for the acqisition of the spectra.
)
{
  noise<-sqrt((sdLamp^2)+(sdDC^2))/IntegrationTime
  # The noise, based on standard deviation of lamp signal and dark current
  spectralNER<-noise*RadCalCoeff
  ##value<< numeric vector: noise equivalent delta radiance for each wavelength
  return(spectralNER)
}

attr(NoiseEquivalentRadiance,"ex") <- function(){
  

    data("snr_data")
    data("rad_cal")
    
    #perform statistics on spectra
    #calculate standard deviation of the lamp signal
    lamp_sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='sd')
    #calculate standard deviation of the dark spectrum
    dc_sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_dc,fun='sd')
    #calculate noise equivalent delta radiance
    NER<-NoiseEquivalentRadiance(sdLamp = lamp_sd,sdDC = dc_sd,IntegrationTime =  4000,RadCalCoeff = rad_cal$cal)
    x11();plot(rad_cal$wl,NER*1000,type="l",xlab="WL [nm]",ylab=expression("Radiance [mW m"^-2* "sr"^-1* "nm"^-1*"]"),ylim=c(0,5))
}


#####################################################################################################################################################################

RadiometricCalibration<-function(
  ### Compute the VIS-NIR radiometric calibration coefficients based on outdoor measurements 
  type  	##<< value: 1 if the same white panel has been used; 2 if the spectrometer connected to the cosine receptor has to be cross calibrated with measurements collected over  a white panel 
  ,wl		##<< numeric vector: wavelength vector
  ,radiance	##<< numeric vector: mean vector of n radiance spectra acquired 
  ,DN		##<< numeric vector: mean vector of n spectra in Digital Number 
		## acquired over the same reference standard during the same time interval 
		## (DN spectrum has to be already dark current subtracted and divided by integration time used) 
  ,ReferenceCoeff  ##<< numeric numeric vector: 
		## vector fo the calibration coefficients of the refrence standard used. 
){
  ##details<<
  ## Outdoor measurements are acquired simoultaneously over the same reference standard (e.g. white panel) 
  ## or over different reference standard (e.g. white panel and cosine receptor) 
  ## by a reference spectrometer and the spectrometer to be calibrated.
  ## All the provided data in input need to be pre-processed in order to be homogenized in terms of wavelength
  ## (see \code{\link{LinearResample}})
  if (type==1) gain<-radiance/DN
  if (type==2){upwelling_radiance<-radiance/ReferenceCoeff;gain<-upwelling_radiance/DN}
  ##value<< numeric vector: conversion factor from digital counts to radiance for each wavelength 
  return(gain)
}


attr(RadiometricCalibration,"ex") <- function(){
  

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


#####################################################################################################################################################################

ExcludeSpectralRegions<-function(
  ### Give NA values to specific spectral regions of the spectrum
  wl##<< numeric vector: wavelength vector
  ,spectrum  ##<< numeric vector: vector reporting a specific value at each wavelenght. 
  ,SpectralRegion  ##<< numeric data.frame: data.frame containing the absorption atmospheric regions. First column refer to the first wl of the region, second column to the last wl of the region. Rows represent the number of regions to be considered
)
{
for(n in 1:dim(SpectralRegion)[1])
{
  RangeToExclude<-which(wl>= SpectralRegion[n,1] & wl<= SpectralRegion[n,2])
  if(length(RangeToExclude)>0) {spectrum[RangeToExclude]<-NA}
}
  spectrum[which(spectrum==Inf)]<-NA
  spectrum[which(spectrum==-Inf)]<-NA
  ##value<< numeric vector: vector of values (same units as input) with NA values in corrispondance of the selected spectral range
  return(spectrum)
}


attr(ExcludeSpectralRegions,"ex") <- function(){
  

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
}


#####################################################################################################################################################################

DNSpectralMatrixRadCal<-function(
  ### Create a numeric matrix where the spectra (in Digital Number) are divided by the integration time of the measurement 
  spectra  				##<< numeric data.frame or matrix: n columns correspond to n spectra acquired during the cross calibration measurements
  ,IntegrationTime  	##<< numeric value or vector: if all the measurement collected during the field cross calibration have been acquired using the same integration time one value is needed. 
  	## If the integration time changed within the measurements a vector of intergration time can be provided (length of vector must to be equal to the number of columns of the spectra data frame) 
)
{
  spectra<-as.matrix(spectra)
  if (length(IntegrationTime)==1){
	  DN_matrix<-spectra/IntegrationTime  
  } else {
    if(length(IntegrationTime) != dim(spectra)[2]) warning(
				"length of integration time did not correspond to the number of columns in the provided spectra")
     DN_matrix<-t(apply(spectra, 1, function(x) x/IntegrationTime))
  }
  ##value<< numeric matrix of the same size as input spectra, with each entry devided by the corresponding integration time
  return(DN_matrix)
}


attr(DNSpectralMatrixRadCal,"ex") <- function(){
  
    data("outdoor_rad_cal_data")
    data("atmospheric_absorption_regions")
    integration_time<-450
    #Create matrix for radiometric calibration
    DN_mat<-DNSpectralMatrixRadCal(spectra = outdoor_rad_cal_data$DN_matrix,IntegrationTime = integration_time)
    
}



#####################################################################################################################################################################

DCSubtraction<-function(
  ### Remove the dark current signal from the signal itself. Preferably remove the dark spectrum from the signal spectrum. In case no dark spectrum is acquired it use the average value of the optically black pixels
  signal  		##<< numeric data.frame or vector: in case of data frame the number of colums refer to the number of spectra acquired
  ,DarkSignal	##<< numeric data.frame or vector: in case of data frame the number of colums refer to the number of dark spectra acquired
  ,type = 1		##<< numeric value: 1 in case dark spectra have been acquired; 2 in case optically black pixels have to be used
  ,BlackPixels = 1:4##<<numeric vector: vector containing the number of the optically black pixels
){
  if (type == 1) 
  {
    if (class(signal)== "data.frame" & class(DarkSignal)== "data.frame") {SignalDcSubtracted<-signal-DarkSignal}
    if (class(signal)== "matrix" & class(DarkSignal)== "matrix") {SignalDcSubtracted<-signal-DarkSignal}
    if (class(signal)== "numeric" & class(DarkSignal)== "numeric") {SignalDcSubtracted<-signal-DarkSignal}
    if (class(signal)== "integer" & class(DarkSignal)== "integer") {SignalDcSubtracted<-signal-DarkSignal}
    if (class(signal)!=class(DarkSignal))
      {
      warning("class of signal and class of dark current differ. Na returned")
      SignalDcSubtracted<-NA
      }
    return(SignalDcSubtracted)
  }
  if (type == 2)
  {
    if (class(signal)== "numeric") 
    {
      DcObPixel<-mean(signal[BlackPixels],na.rm=TRUE)
      SignalDcSubtracted<-as.vector(signal-DcObPixel)
    }
    
    if(class(signal)== "data.frame")
    {
      n_spectra<-dim(signal)[2]
      for(n in 1:n_spectra)
      {
        DcObPixel<-mean(signal[BlackPixels,n],na.rm=TRUE)
        signal_dc_sub<-signal[,n]-DcObPixel
        if(n==1) {SignalDcSubtracted<-signal_dc_sub}else{
          SignalDcSubtracted<-cbind(SignalDcSubtracted,signal_dc_sub) 
        }
      }
    }
    ##value<< numeric data.frame or vector: spectra with durk current removed
    return(SignalDcSubtracted)
  }
}


attr(DCSubtraction,"ex") <- function(){
  
    data("snr_data")
    data("rad_cal")
    
    dataDCsubtracted<-DCSubtraction(signal = snr_data$data_lamp,DarkSignal = snr_data$data_dc)
    x11();plot(rad_cal$wl,snr_data$data_lamp[,1],type="l",ylim=c(0,16000),ylab="Digital Counts",xlab="WL [nm]")
    lines(rad_cal$wl,snr_data$data_dc[,1],col="red");lines(rad_cal$wl,dataDCsubtracted[,1],col="blue")
    legend("topleft",col=c("black","red","blue"),lty=1,cex=1.2,legend=c("Signal","dark current","Signal DC subtracted"),box.col="white")
    box()
}



#####################################################################################################################################################################

SelectSpectralRegion<-function(
  ### Select spectral region around defined wavelenths with a specified buffer
  wl  ##<< numeric vector: wavelength vector
  ,spectrum ##<< numeric vector: spectrum to be processed
  ,WlSelection  ##<< numeric vector or value: vector or value of defined wavelength to be selected
  ,buffer  ##<< numeric value: width (in nanometer) of the buffer around the defined wavelength
)
{
  WlSelection<-WlSelection[which(WlSelection>=wl[1]&WlSelection<=wl[length(wl)])]
  NPixels<-1:length(wl)
  n_wl_sel<-length(WlSelection)
  SelectedRegions <- list()
  for(n in 1:n_wl_sel)
  {
    SingleRegionSelection<-which(wl>= WlSelection[n]-(buffer/2) & wl<=WlSelection[n]+(buffer/2))
    name <- as.character(WlSelection[n])
    region<-data.frame(NPixels[SingleRegionSelection],wl[SingleRegionSelection],spectrum[SingleRegionSelection]);names(region)<-c("n_pixel","wl","DN")
    SelectedRegions[[name]] <- region
  }
  return(SelectedRegions)
  ##value<< numeric list containing for each selected range a data.frame with: pixel numbers,wavelengths and digital count values.
}


attr(SelectSpectralRegion,"ex") <- function(){
  
    data("indoor_wl_cal_data")
    lamp_spectrum_dc_sub<-DCSubtraction(signal=indoor_wl_cal_data$lamp_spectrum,DarkSignal = indoor_wl_cal_data$dc_spectrum,type=1)
    region_to_analyze<-SelectSpectralRegion(wl = indoor_wl_cal_data$DN_wl,spectrum = lamp_spectrum_dc_sub,WlSelection = indoor_wl_cal_data$emission_lines$peak,buffer=3)
    
}


#####################################################################################################################################################################

GaussFit<-function(
  ### fit a gaussian function to points and extract the corresponidng center, in terms of wl and pixel position, and Full Width at Half Maximum
  NPixels ##<< numeric vector: number of pixel vector to be fitted
  ,wl ##<< numeric vector: wavelength vector
  ,DN  ##<< numeric vector: DN vector to be fitted
  ,plot=FALSE ##<< a logical value indicating whether the output is plotted or not
)
{
  if(length(DN)<5) 
  {
    warning("too few points to perform a fitting. NAs returned")
    CenterPixel<-NA;center_wl<-NA;fwhm<-NA
  }else{
    
    #control on double peak
    dpcontrol<-diff(DN)
    min_der<-which(dpcontrol==min(dpcontrol))
    max_before_peak_pos<-which(dpcontrol==max(dpcontrol[1:min_der]))
    max_before_peak<-dpcontrol[max_before_peak_pos]
    max_after_peak_pos<-which(dpcontrol==max(dpcontrol[min_der:length(dpcontrol)]))
    max_after_peak<-dpcontrol[max_after_peak_pos]
    if((max_before_peak-max_after_peak)<0.1*max_before_peak)
    {
      warning("Double peak found in the specified range. NA values are returned. Please revise the buffer size.")
      CenterPixel<-NA;center_wl<-NA;fwhm<-NA
    }else{
      
      #subset around the peak
      dpcontrol<-abs(dpcontrol)
      dpcontrol<-(dpcontrol-min(dpcontrol))/(max(dpcontrol)-min(dpcontrol))
      gret_than_th<-which(dpcontrol>0.1)
      if(gret_than_th[1]>4 & (length(dpcontrol)- gret_than_th[length(gret_than_th)])>4)
      {
        Selected_range<-c((gret_than_th[1]-4):(gret_than_th[length(gret_than_th)]+4))
        DN<-DN[Selected_range];NPixels<-NPixels[Selected_range];wl<-wl[Selected_range]
      }
      # Estimate some starting values.
      mu <- wl[which.max(DN)]; sigma <- (max(wl)-min(wl))/4; k <-max(DN)-min(DN)
      tab <- data.frame(x=wl, r=DN)
      res <- nls( r ~ k*exp(-1/2*(x-mu)^2/sigma^2), start=c(mu=mu,sigma=sigma,k=k) , data = tab)
      center_wl<-round(summary(res)$parameters["mu", 1],2)
      fwhm<-round((2*sqrt(2*log(2)))*summary(res)$parameters["sigma", 1],2)
      v <- summary(res)$parameters[,"Estimate"]
      
      if(fwhm>(wl[length(wl)]-wl[1])) warning("FWHM estiate exceeds spectral the considered spectral range. Possible error.")
      
      if(plot==TRUE)
      {
        x11()
        plot(DN~wl,ylim=c(0,max(DN)+0.3*max(DN)),xlab="wl",ylab="DN",type="o",lty=2,lwd=0.5,pch=20)
        plot(function(x) v[3]*exp(-1/2*(x-v[1])^2/v[2]^2),col=2,add=T,xlim=range(tab$x) )
        legend("topright",pch=NA,cex=1.2,legend=c(paste("CENTER: ",center_wl,sep=""),paste("FWHM: ",fwhm,sep="")),box.col="white")
        box()
      }
      
      mu <- NPixels[which.max(DN)]; sigma <- (max(NPixels)-min(NPixels))/4; k <-max(DN)-min(DN)
      tab <- data.frame(x=NPixels, r=DN)
      res <- nls( r ~ k*exp(-1/2*(x-mu)^2/sigma^2), start=c(mu=mu,sigma=sigma,k=k) , data = tab)
      CenterPixel<-round(summary(res)$parameters["mu", 1],2)
    }
  }
  out_fit<-data.frame(CenterPixel,center_wl,fwhm)
  
  ##value<< numeric data.frame containing pixel number of the peak, the corresponing wavelength and the full width at half maximum.
  return(out_fit)
}



attr(GaussFit,"ex") <- function(){
  

    data("indoor_wl_cal_data")
    #extract dc subracted spectra
    lamp_spectrum_dc_sub<-DCSubtraction(signal=indoor_wl_cal_data$lamp_spectrum,DarkSignal = indoor_wl_cal_data$dc_spectrum,type=1)
    #Select spectral region to analyse
    region_to_analyze<-SelectSpectralRegion(wl = indoor_wl_cal_data$DN_wl,spectrum = lamp_spectrum_dc_sub,WlSelection = indoor_wl_cal_data$emission_lines$peak,buffer=1)
    #define data.frame containing the expected results
    n_peaks_fit<-length(names(region_to_analyze))
    wl_peaks<-as.numeric(names(region_to_analyze))
    #loop on spectral region to analyze
    for(n in 1:n_peaks_fit)
    { print(n)
      n_pixels<-data.frame(region_to_analyze[n])[,1]
      wl<-data.frame(region_to_analyze[n])[,2]
      DN<-data.frame(region_to_analyze[n])[,3]  
      wl_param<-GaussFit(n_pixels,wl,DN,plot=TRUE)
      if(n==1){wl_cal<-wl_param}else{
        wl_cal<-rbind(wl_cal,wl_param)}
    }
    
}




####################################################################################################################################################################

WlCal<-function(
  ### fit a 3 order polynomial to pix value and wl for wavelength calibration
  pix_center ##<< numeric vector: center pixel corresponding to emission line peak
  ,wl_peaks ##<< numeric vector: wavelength vector of emission lines lamp
)
{
  data_lm<-data.frame(pix_center,wl_peaks);data_lm<-na.omit(data_lm)
  x<-as.array(data_lm[,1])
  y<-as.array(data_lm[,2])
  model <- lm(y ~ x + I(x^2) + I(x^3))
  ##value<< an object of class "lm" containing the coefficients needed for the wavelength calibration.
  
  return(model)
}


attr(WlCal,"ex") <- function(){
  

    data("indoor_wl_cal_data")
    #extract dc subracted spectra
    lamp_spectrum_dc_sub<-DCSubtraction(signal=indoor_wl_cal_data$lamp_spectrum,DarkSignal = indoor_wl_cal_data$dc_spectrum,type=1)
    #Select spectral region to analyse
    region_to_analyze<-SelectSpectralRegion(wl = indoor_wl_cal_data$DN_wl,spectrum = lamp_spectrum_dc_sub,WlSelection = indoor_wl_cal_data$emission_lines$peak,buffer=1)
    #define data.frame containing the expected results
    n_peaks_fit<-length(names(region_to_analyze))
    wl_peaks<-as.numeric(names(region_to_analyze))
    #loop on spectral region to analyze
    for(n in 1:n_peaks_fit)
    { print(n)
      n_pixels<-data.frame(region_to_analyze[n])[,1]
      wl<-data.frame(region_to_analyze[n])[,2]
      DN<-data.frame(region_to_analyze[n])[,3]  
      wl_param<-GaussFit(n_pixels,wl,DN,plot=TRUE)
      if(n==1){wl_cal<-wl_param}else{
        wl_cal<-rbind(wl_cal,wl_param)}
    }
    
    #extarct coefficients for wl calibration
    wl_coeff<-WlCal(pix_center = wl_cal$CenterPixel,wl_peaks = wl_peaks)
}


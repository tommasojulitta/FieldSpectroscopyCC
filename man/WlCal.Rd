\name{WlCal}
\alias{WlCal}
\title{WlCal}
\description{fit a 3 order polynomial to pix value and wl for wavelength calibration}
\usage{WlCal(pix_center, wl_peaks)}
\arguments{
  \item{pix_center}{numeric vector: center pixel corresponding to emission line peak}
  \item{wl_peaks}{numeric vector: wavelength vector of emission lines lamp}
}

\value{an object of class "lm" containing the coefficients needed for the wavelength calibration.}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

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

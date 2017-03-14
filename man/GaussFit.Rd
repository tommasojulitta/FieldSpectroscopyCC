\name{GaussFit}
\alias{GaussFit}
\title{GaussFit}
\description{fit a gaussian function to points and extract the corresponidng center, in terms of wl and pixel position, and Full Width at Half Maximum}
\usage{GaussFit(NPixels, wl, DN, plot = FALSE)}
\arguments{
  \item{NPixels}{numeric vector: number of pixel vector to be fitted}
  \item{wl}{numeric vector: wavelength vector}
  \item{DN}{numeric vector: DN vector to be fitted}
  \item{plot}{a logical value indicating whether the output is plotted or not}
}

\value{numeric data.frame containing pixel number of the peak, the corresponing wavelength and the full width at half maximum.}

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

}

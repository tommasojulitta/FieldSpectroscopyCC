\name{NoiseEquivalentRadiance}
\alias{NoiseEquivalentRadiance}
\title{NoiseEquivalentRadiance}
\description{Compute the Signal to Noise Ratio of a spectrometer}
\usage{NoiseEquivalentRadiance(sdLamp = specStats$sdLamp, sdDC = specStats$sdDC, 
    RadCalCoeff, IntegrationTime)}
\arguments{
  \item{sdLamp}{numeric vector: sd of the spectra acquired with a stable lamp}
  \item{sdDC}{numeric vector: sd of the dark current spectra acquired with with the same integration time as the measurements over the lamp}
  \item{RadCalCoeff}{numeric vector: wavelength dependent vector of coefficient for calibration}
  \item{IntegrationTime}{numeric value: integration time used for the acqisition of the spectra.}
}

\value{numeric vector: noise equivalent delta radiance for each wavelength}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

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

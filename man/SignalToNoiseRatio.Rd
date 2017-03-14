\name{SignalToNoiseRatio}
\alias{SignalToNoiseRatio}
\title{SignalToNoiseRatio}
\description{Compute the Signal to Noise Ratio of a spectrometer}
\usage{SignalToNoiseRatio(specStats, avgLamp = specStats$meanLamp, sdLamp = specStats$sdLamp, 
    avgDC = specStats$meanDC, sdDC = specStats$sdDC)}
\arguments{
  \item{specStats}{alternative parameterization data.frame with columns meanLight, sdLight, meanDark, sdDark}
  \item{avgLamp}{numeric vector: mean spectrum of n scans of stable lamp}
  \item{sdLamp}{numeric vector: sd of the spectra acquired with a stable lamp}
  \item{avgDC}{numeric vector: mean spectrum of dark current measurements acquired with the same integration time as the measurements over the lamp}
  \item{sdDC}{numeric vector: sd of the dark current spectra acquired with with the same integration time as the measurements over the lamp}
}

\value{numeric vector: ratio of signal and noise for each wavelength}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

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

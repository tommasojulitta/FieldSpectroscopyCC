\name{StatsOnSpectra}
\alias{StatsOnSpectra}
\title{StatsOnSpectra}
\description{Compute statistic across several spectra subset to a specified spectral range}
\usage{StatsOnSpectra(wl, wlStart = wl[1], wlEnd = wl[length(wl)], spectra = spectra, 
    fun = "mean", margin = 1)}
\arguments{
  \item{wl}{numeric vector: wavelength for each row in spectra}
  \item{wlStart}{numeric value or vector: first wavelength of the spectral range selection. Default value: first wavelength of wl vector}
  \item{wlEnd}{numeric value or vector: last wavelength of the spectral range selection. Default value: last wavelength of wl vector}
  \item{spectra}{numeric matrix or data.frame: collection of several spectra acquired by several columns (Digital numbers, Radiance, Reflectance).}
  \item{fun}{character or function: function to be applied on each row of the data.frame. Default function: mean}
  \item{margin}{numeric value: a vector giving the subscripts which the function will be applied over. 1 indicates rows, 2 indicates columns}
}

\value{numeric data.frame containing the computed statistic across the selected spectral range.}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  
data("snr_data")
#perform statistics on spectra
#calculate mean of the  signal between 700 and 800 nm
mean<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='mean',wlStart=700,wlEnd = 800)
#calculate stadard deviation of the signal between 700 and 800 nm
sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='sd',wlStart=700,wlEnd = 800)
}

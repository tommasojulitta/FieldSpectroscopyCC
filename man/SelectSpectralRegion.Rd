\name{SelectSpectralRegion}
\alias{SelectSpectralRegion}
\title{SelectSpectralRegion}
\description{Select spectral region around defined wavelenths with a specified buffer}
\usage{SelectSpectralRegion(wl, spectrum, WlSelection, buffer)}
\arguments{
  \item{wl}{numeric vector: wavelength vector}
  \item{spectrum}{numeric vector: spectrum to be processed}
  \item{WlSelection}{numeric vector or value: vector or value of defined wavelength to be selected}
  \item{buffer}{numeric value: width (in nanometer) of the buffer around the defined wavelength}
}

\value{numeric list containing for each selected range a data.frame with: pixel numbers,wavelengths and digital count values.}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  
data("indoor_wl_cal_data")
lamp_spectrum_dc_sub<-DCSubtraction(signal=indoor_wl_cal_data$lamp_spectrum,DarkSignal = indoor_wl_cal_data$dc_spectrum,type=1)
region_to_analyze<-SelectSpectralRegion(wl = indoor_wl_cal_data$DN_wl,spectrum = lamp_spectrum_dc_sub,WlSelection = indoor_wl_cal_data$emission_lines$peak,buffer=3)

}

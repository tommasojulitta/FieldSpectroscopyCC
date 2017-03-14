\name{LinearResample}
\alias{LinearResample}
\title{LinearResample}
\description{Compute linear resample of a spectrum based on specific wavelength of a defined spectrometer}
\usage{LinearResample(wl, spectrum, wlRef)}
\arguments{
  \item{wl}{numeric vector: wavelength vector to be resamped}
  \item{spectrum}{numeric vector: spectrum to be resampled (e.g. radiance, reflectance, digital numbers)}
  \item{wlRef}{numeric vector: reference wavelength vector to be be used for resampling}
}

\value{numeric vector containing the spectrum linearly resample according to the wavelngth vector wlRef.}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

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

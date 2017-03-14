\name{PolSmoothGapfilling}
\alias{PolSmoothGapfilling}
\title{PolSmoothGapfilling}
\description{Apply a polynomial fitting and gap filling in the region where NA are found}
\usage{PolSmoothGapfilling(wl, spectrum, deg = 4)}
\arguments{
  \item{wl}{numeric vector: wavelength vector}
  \item{spectrum}{numeric data.frame: First column wavelength vector, second column spectrum vector}
  \item{deg}{numeric value: the degree of the polynomial}
}

\value{numeric vector containing the smoothed spectrum.}

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

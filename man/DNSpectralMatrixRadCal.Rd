\name{DNSpectralMatrixRadCal}
\alias{DNSpectralMatrixRadCal}
\title{DNSpectralMatrixRadCal}
\description{Create a numeric matrix where the spectra (in Digital Number) are divided by the integration time of the measurement }
\usage{DNSpectralMatrixRadCal(spectra, IntegrationTime)}
\arguments{
  \item{spectra}{numeric data.frame or matrix: n columns correspond to n spectra acquired during the cross calibration measurements}
  \item{IntegrationTime}{numeric value or vector: if all the measurement collected during the field cross calibration have been acquired using the same integration time one value is needed.
If the integration time changed within the measurements a vector of intergration time can be provided (length of vector must to be equal to the number of columns of the spectra data frame) }
}

\value{numeric matrix of the same size as input spectra, with each entry devided by the corresponding integration time}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  
data("outdoor_rad_cal_data")
data("atmospheric_absorption_regions")
integration_time<-450
#Create matrix for radiometric calibration
DN_mat<-DNSpectralMatrixRadCal(spectra = outdoor_rad_cal_data$DN_matrix,IntegrationTime = integration_time)

}

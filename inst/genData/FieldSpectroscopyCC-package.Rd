\name{FieldSpectroscopyCC-package}
\alias{FieldSpectroscopyCC-package}
\alias{FieldSpectroscopyCC}
\title{Spectrometer Characterization and Calibratraion}
\description{
Characterization of a spectrometer, 
radiometric calibration (indoor, outdoor),spectral calibration (indoor).
}
\details{

The basic functions are: 
\itemize{
\item{ Compute statistic on spectra data.frame over a specified spectral range: \code{\link{StatsOnSpectra}}  }
\item{ Compute common summary statistics on light and dark field spectra: \code{\link{SummarizeSpectra}}  }
\item{ Compute linear resample of a spectrum based on specific wavelength of a defined spectrometer: \code{\link{LinearResample}}  }
\item{ Apply a spline smoothing and gap filling in the region where NA are found: \code{\link{SplineSmoothGapfilling}}  }
\item{ Apply a polynomial fitting and gap filling in the region where NA are found: \code{\link{PolSmoothGapfilling}}  }
}%itemize

There are functions helping wiht input/output (TODO):
\itemize{
\item{ Create matrix of spectra from single measurements: \code{\link{CreateDataMatrix}}  }
\item{ Load a vector form an asci files: \code{\link{DataLoadVector}}  }
\item{ Read the header of an asci file: \code{\link{DataOpenHeader}}  }
\item{ Create .Rdata file containing calibration coefficients for spectral and radiometric calibration.: \code{\link{SaveCalFile}}  }

}%itemize

Further functionality (TODO add other functions maybe grouped by topic): 
\itemize{
\item{ Compute the Signal to Noise Ratio of a spectrometer: \code{\link{SignalToNoiseRatio}}  }
\item{ Compute the Noise equivalent in Radiance of a spectrometer: \code{\link{NoiseEquivalentRadiance}}  }
\item{ Compute the VIS-NIR radiometric calibration coefficients based on outdoor measurements: \code{\link{RadiometricCalibration}}  }
\item{ Give NA values to specific spectral regions of the spectrum: \code{\link{ExcludeSpectralRegions}}  }
\item{ Create a numeric matrix where the spectra in Digital Number are divided by the integration time of the measurement: \code{\link{DNSpectralMatrixRadCal}}  }
\item{ Remove the dark current signal from the signal itself. Preferably remove the dark spectrum from the signal spectrum. In case no dark spectrum is acquired it use the average value of the optically black pixels, if any: \code{\link{DCSubtraction}}  }
\item{ Select spectral region around defined wavelenths with a specified buffer: \code{\link{SelectSpectralRegion}}  }
\item{ Fit a gaussian function to points and extract the corresponidng center, in terms of wl and pixel position, and Full Width at Half Maximum: \code{\link{GaussFit}}  }
\item{ Fit a 3 order polynomial to pix value and wl for wavelength calibration: \code{\link{WlCal}}  }

}%itemize

}%details

\author{Thomas Wutzler}

\keyword{ package }



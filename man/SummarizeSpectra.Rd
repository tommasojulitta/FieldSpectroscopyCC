\name{SummarizeSpectra}
\alias{SummarizeSpectra}
\title{SummarizeSpectra}
\description{compute common summary statistics on light and dark field spectra }
\usage{SummarizeSpectra(spectra, ...)}
\arguments{
  \item{spectra}{list with entries wl (wavelength), data_lamp (matrix of several light spectra), and data_dc (matrix of several dark current)}
  \item{\dots}{further arguments to \code{\link{StatsOnSpectra}}, such as \code{wlStart} and \code{wlEnd}}
}

\value{data.frame with entries
\item{meanLamp}{mean light spectrum}
\item{sdLamp}{standard deviation across light spectra}
\item{nLamp}{number of summarized light spectra}
\item{meanDC}{mean dark spectrum}
\item{sdDC}{standard deviation across dark spectra}
\item{nDar}{number of summarized dark spectra
}}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

data("snr_data")
#perform statistics on spectra
specStats <- SummarizeSpectra(snr_data)
#print header of specStats
head(specStats)
}

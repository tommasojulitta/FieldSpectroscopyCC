\name{DCSubtraction}
\alias{DCSubtraction}
\title{DCSubtraction}
\description{Remove the dark current signal from the signal itself. Preferably remove the dark spectrum from the signal spectrum. In case no dark spectrum is acquired it use the average value of the optically black pixels}
\usage{DCSubtraction(signal, DarkSignal, type = 1, BlackPixels = 1:4)}
\arguments{
  \item{signal}{numeric data.frame or vector: in case of data frame the number of colums refer to the number of spectra acquired}
  \item{DarkSignal}{numeric data.frame or vector: in case of data frame the number of colums refer to the number of dark spectra acquired}
  \item{type}{numeric value: 1 in case dark spectra have been acquired; 2 in case optically black pixels have to be used}
  \item{BlackPixels}{numeric vector: vector containing the number of the optically black pixels}
}

\value{numeric data.frame or vector: spectra with durk current removed}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  
data("snr_data")
data("rad_cal")

dataDCsubtracted<-DCSubtraction(signal = snr_data$data_lamp,DarkSignal = snr_data$data_dc)
x11();plot(rad_cal$wl,snr_data$data_lamp[,1],type="l",ylim=c(0,16000),ylab="Digital Counts",xlab="WL [nm]")
lines(rad_cal$wl,snr_data$data_dc[,1],col="red");lines(rad_cal$wl,dataDCsubtracted[,1],col="blue")
legend("topleft",col=c("black","red","blue"),lty=1,cex=1.2,legend=c("Signal","dark current","Signal DC subtracted"),box.col="white")
box()
}

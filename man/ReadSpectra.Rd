\name{ReadSpectra}
\alias{ReadSpectra}
\title{ReadSpectra}
\description{Read Spectra from ASCII files and list them }
\usage{ReadSpectra(filenames, WLname = NULL, WLpos = NULL, Ename = NULL, 
    E2name = NULL, dcEname = NULL, Lname = NULL, dcLname = NULL, 
    Epos = NULL, E2pos = NULL, dcEpos = NULL, Lpos = NULL, dcLpos = NULL, 
    ...)}
\arguments{
  \item{filenames}{character value or vector: names of the file(s) to be opened}
  \item{WLname}{character value: string of the name in the ASCII file of the wavelength vector, if any}
  \item{WLpos}{numeric value: value of the position (number of column) in the ASCII file of the wavelength vector, if any}
  \item{Ename}{character value: string of the name in the ASCII file of the solar radiance vector, if any}
  \item{E2name}{character value: string of the name in the ASCII file of the second solar radiance vector (in case of sandwicth acquisition), if any}
  \item{dcEname}{character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of solar irradiance, if any}
  \item{Lname}{character value: string of the name in the ASCII file of the reflected radiance vector, if any}
  \item{dcLname}{character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of reflected radiance vector, if any}
  \item{Epos}{numeric value: value of the position (number of column) in the ASCII file of the solar radiance vector, if any}
  \item{E2pos}{numeric value: value of the position (number of column) in the ASCII file of the second solar radiance vector (in case of sandwicth acquisition), if any}
  \item{dcEpos}{numeric value: value of the position (number of column) in the ASCII file of dark current vector acquired with same integration time of solar irradiance, if any}
  \item{Lpos}{numeric value: value of the position (number of column) in the ASCII file of the the reflected radiance vector, if any}
  \item{dcLpos}{numeric value: value of the position (number of column) in the ASCII file of dark current vector acquired with same integration time of reflected irradiance, if any}
  \item{\dots}{Further arguments to be passed to read.table.}
}

\value{list containing all the info specified in the opened file(s).}

\author{Tommaso Julitta, Mirco Migliavacca, Thomas Wutzler}




\examples{
  

path<-"..."
#list files 
list_files<-list.files(path = path,pattern = "h4",full.names = TRUE)
#open an list them in a single file
dat<-ReadSpectra(filenames = list_files,WLpos = 1,Epos = 2,Lpos = 3,header =TRUE, sep=";")
}



ReadSpectra<-function(
  ### Read Spectra from ASCII files and list them 
  filenames ##<< character value or vector: names of the file(s) to be opened
  ,WLname =NULL  ##<< character value: string of the name in the ASCII file of the wavelength vector, if any
  ,WLpos =NULL  ##<< numeric value: value of the position (number of column) in the ASCII file of the wavelength vector, if any
  ,Ename =NULL  ##<< character value: string of the name in the ASCII file of the solar radiance vector, if any
  ,E2name =NULL  ##<< character value: string of the name in the ASCII file of the second solar radiance vector (in case of sandwicth acquisition), if any
  ,dcEname =NULL  ##<< character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of solar irradiance, if any
  ,Lname =NULL  ##<< character value: string of the name in the ASCII file of the reflected radiance vector, if any
  ,dcLname =NULL  ##<< character value: string of the name in the ASCII file of the dark current vector acquired with same integration time of reflected radiance vector, if any
  ,Epos =NULL  ##<< numeric value: value of the position (number of column) in the ASCII file of the solar radiance vector, if any
  ,E2pos =NULL  ##<< numeric value: value of the position (number of column) in the ASCII file of the second solar radiance vector (in case of sandwicth acquisition), if any
  ,dcEpos =NULL  ##<< numeric value: value of the position (number of column) in the ASCII file of dark current vector acquired with same integration time of solar irradiance, if any
  ,Lpos =NULL  ##<< numeric value: value of the position (number of column) in the ASCII file of the the reflected radiance vector, if any
  ,dcLpos =NULL  ##<< numeric value: value of the position (number of column) in the ASCII file of dark current vector acquired with same integration time of reflected irradiance, if any
  ,...  ##<< Further arguments to be passed to read.table.
)
{
  meas_data<-lapply(filenames, function(x) read.table(..., x))
  ### Open the ASCII file
  Data<-list()
  ### Create an empty list
  if(!is.null(WLname))  {Wl<-lapply(meas_data, "[[", WLame);Wl<-matrix(unlist(Wl),ncol=length(sapply(Wl,dim)));Wl<-Wl[,1];Data[["Wl"]]<- Wl}
  if(!is.null(WLpos))  {Wl<-lapply(meas_data, "[[", WLpos);Wl<-matrix(unlist(Wl),ncol=length(sapply(Wl,dim)));Wl<-Wl[,1];Data[["Wl"]]<- Wl}
  if(!is.null(Ename))  {E<-lapply(meas_data, "[[", Ename);E<-matrix(unlist(E),ncol=length(sapply(E,dim)));Data[["E"]]<- E}
  if(!is.null(Epos))  {E<-lapply(meas_data, "[[", Epos);E<-matrix(unlist(E),ncol=length(sapply(E,dim)));Data[["E"]]<- E}
  if(!is.null(E2name))  {E2<-lapply(meas_data, "[[", E2name);E2<-matrix(unlist(E2),ncol=length(sapply(E2,dim)));Data[["E2"]]<- E2}
  if(!is.null(E2pos))  {E2<-lapply(meas_data, "[[", E2pos);E2<-matrix(unlist(E2),ncol=length(sapply(E2,dim)));Data[["E2"]]<- E2}
  if(!is.null(dcEname))  {dcE<-lapply(meas_data, "[[", dcEname);dcE<-matrix(unlist(dcE),ncol=length(sapply(dcE,dim)));Data[[ "dcE" ]] <- dcE}
  if(!is.null(dcEpos))  {dcE<-lapply(meas_data, "[[", dcEpos);dcE<-matrix(unlist(dcE),ncol=length(sapply(dcE,dim)));Data[[ "dcE" ]] <- dcE}
  if(!is.null(Lname))  {L<-lapply(meas_data, "[[", Lname);L<-matrix(unlist(L),ncol=length(sapply(L,dim)));Data[[ "L" ]] <- L}
  if(!is.null(Lpos))  {L<-lapply(meas_data, "[[", Lpos);L<-matrix(unlist(L),ncol=length(sapply(L,dim)));Data[[ "L" ]] <- L}
  if(!is.null(dcLname))  {dcL<-lapply(meas_data, "[[", dcLname);dcL<-matrix(unlist(dcL),ncol=length(sapply(dcL,dim)));Data[[ "dcL" ]] <- dcL}
  if(!is.null(dcLpos))  {dcL<-lapply(meas_data, "[[", dcLpos);dcL<-matrix(unlist(dcL),ncol=length(sapply(dcL,dim)));Data[[ "dcL" ]] <- dcL}
  ### Fill the list with data contained in the ASCII file, attributing specific names (E,wl, L, dc etc...) to the list  
  ##value<< list containing all the info specified in the opened file(s).
  return(Data)
}  

attr(ReadSpectra,"ex") <- function(){
  

    path<-"..."
    #list files 
    list_files<-list.files(path = path,pattern = "h4",full.names = TRUE)
    #open an list them in a single file
    dat<-ReadSpectra(filenames = list_files,WLpos = 1,Epos = 2,Lpos = 3,header =TRUE, sep=";")
}

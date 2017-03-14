#require(testthat)
context("testBaseFun")

packageDir <- system.file(package='FieldSpectroscopyCC')	# only works if package has been installed
load(file.path(packageDir,"data","snr_data.RData"))

test_that("StatsOnSpectra",{
		lamp_mean<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='mean')			
		lamp_mean<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun=mean)	# proving function directly is usually faster			
		expect_that( length(lamp_mean), equals(length(snr_data$wl)) )
		lamp_sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_lamp,fun='sd')
		dc_mean<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_dc,fun='mean')
		dc_sd<-StatsOnSpectra(wl=snr_data$wl,spectra=snr_data$data_dc,fun='sd')
	})

test_that("StatsOnSpectra on subSpectra",{
			#TODO
		})


test_that("SignalToNoiseRatio",{
		specStats <- SummarizeSpectra(snr_data)
		SNR<-SignalToNoiseRatio(specStats)
		expect_that( length(SNR), equals(nrow(specStats)))
	})



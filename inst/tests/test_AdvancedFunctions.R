#require(testthat)
context("testBaseFun")

packageDir <- system.file(package='FieldSpectroscopyCC')	# only works if package has been installed
dataDir <- file.path(packageDir,"data") 

load(file.path(dataDir,"outdoor_rad_cal_data.Rdata"))
load(file.path(dataDir,"atmospheric_absorption_regions.Rdata"))


test_that("DN_spectral_matrix_for_rad_cal single integration time",{
			integration_time<-450
			DN_mat<-DN_spectral_matrix_for_rad_cal(outdoor_rad_cal_data$DN_matrix,integration_time=integration_time)
			expect_equal( DN_mat, as.matrix(outdoor_rad_cal_data$DN_matrix)/integration_time )
	})

test_that("DN_spectral_matrix_for_rad_cal vector of  integration time",{
			#TODO
		})
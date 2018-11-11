import numpy
import scipy.stats
import sys
from netCDF4 import Dataset
import netCDF4
import numpy.linalg
import datetime
import time
import matplotlib
#matplotlib.use('Qt4Agg')
matplotlib.use('TkAgg')
import matplotlib.pyplot as mp



##########################################################################################
##########################################################################################
##########################################################################################

season_names = ['djf','mam','jja','son','ondjfm']
#season='djf'
#season='mam'
#season='jja'
#season='son'
#season='annual'
season='ndjfm'

ncfile = Dataset('/ninod/baird/cmip5/concat_nc_files/pr_regrid_2.5x2.5/pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', 'r', format='NETCDF4')
pr_data_orig = ncfile.variables['pr'][:]
pr_lat = ncfile.variables['lat'][:]
pr_lon = ncfile.variables['lon'][:]

global_nlat, global_nlon = pr_data_orig.shape[1:3]
global_lat_vals = pr_lat[:]
global_lon_vals = pr_lon[:]

##########################################################################################
##########################################################################################
##########################################################################################

model_names = numpy.array(( \
'ACCESS1-0', \
'ACCESS1-3', \
'bcc-csm1-1', \
'bcc-csm1-1-m', \
'BNU-ESM', \
'CanESM2', \
'CCSM4', \
'CESM1-BGC', \
'CESM1-CAM5', \
'CMCC-CESM', \
'CMCC-CM', \
'CMCC-CMS', \
'CNRM-CM5', \
'CSIRO-Mk3-6-0', \
'EC-EARTH', \
'FGOALS-g2', \
'GFDL-CM3', \
'GFDL-ESM2G', \
'GFDL-ESM2M', \
'GISS-E2-H', \
'GISS-E2-R', \
'HadGEM2-AO', \
'HadGEM2-CC', \
'HadGEM2-ES', \
'inmcm4', \
'IPSL-CM5A-LR', \
'IPSL-CM5A-MR', \
'IPSL-CM5B-LR', \
'MIROC5', \
'MIROC-ESM-CHEM', \
'MIROC-ESM', \
'MPI-ESM-LR', \
'MPI-ESM-MR', \
'MRI-CGCM3', \
'NorESM1-ME', \
'NorESM1-M' ))

file_root = '/ninod/baird/cmip5/concat_nc_files/pr_regrid_2.5x2.5/'

pr_hist = numpy.array(( \
'pr_Amon_ACCESS1-0_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_ACCESS1-3_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_bcc-csm1-1_historical_r1i1p1_185001-201212_2.5x2.5regrid.nc', \
'pr_Amon_bcc-csm1-1-m_historical_r1i1p1_185001-201212_2.5x2.5regrid.nc', \
'pr_Amon_BNU-ESM_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_CanESM2_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_CCSM4_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_CESM1-BGC_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_CESM1-CAM5_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_CMCC-CESM_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_CMCC-CM_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_CMCC-CMS_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_CNRM-CM5_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_EC-EARTH_historical_r8i1p1_185001-201212_2.5x2.5regrid.nc', \
'pr_Amon_FGOALS-g2_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_GFDL-CM3_historical_r1i1p1_186001-200512_2.5x2.5regrid.nc', \
'pr_Amon_GFDL-ESM2G_historical_r1i1p1_186101-200512_2.5x2.5regrid.nc', \
'pr_Amon_GFDL-ESM2M_historical_r1i1p1_186101-200512_2.5x2.5regrid.nc', \
'pr_Amon_GISS-E2-H_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_GISS-E2-R_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_HadGEM2-AO_historical_r1i1p1_186001-200512_2.5x2.5regrid.nc', \
'pr_Amon_HadGEM2-CC_historical_r1i1p1_185912-200511_2.5x2.5regrid.nc', \
'pr_Amon_HadGEM2-ES_historical_r1i1p1_185912-200511_2.5x2.5regrid.nc', \
'pr_Amon_inmcm4_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_IPSL-CM5A-MR_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_IPSL-CM5B-LR_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_MIROC5_historical_r1i1p1_185001-201212_2.5x2.5regrid.nc', \
'pr_Amon_MIROC-ESM-CHEM_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_MIROC-ESM_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_MPI-ESM-MR_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_MRI-CGCM3_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_NorESM1-ME_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc', \
'pr_Amon_NorESM1-M_historical_r1i1p1_185001-200512_2.5x2.5regrid.nc' ))

pr_rcp85 = numpy.array(( \
'pr_Amon_ACCESS1-0_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_ACCESS1-3_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_bcc-csm1-1_rcp85_r1i1p1_200601-209912_2.5x2.5regrid.nc', \
'pr_Amon_bcc-csm1-1-m_rcp85_r1i1p1_200601-209912_2.5x2.5regrid.nc', \
'pr_Amon_BNU-ESM_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_CanESM2_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_CCSM4_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_CESM1-BGC_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_CESM1-CAM5_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_CMCC-CESM_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_CMCC-CM_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_CMCC-CMS_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_CNRM-CM5_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_CSIRO-Mk3-6-0_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_EC-EARTH_rcp85_r8i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_FGOALS-g2_rcp85_r1i1p1_200601-210112_2.5x2.5regrid.nc', \
'pr_Amon_GFDL-CM3_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_GFDL-ESM2G_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_GFDL-ESM2M_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_GISS-E2-H_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_GISS-E2-R_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_HadGEM2-AO_rcp85_r1i1p1_200601-209912_2.5x2.5regrid.nc', \
'pr_Amon_HadGEM2-CC_rcp85_r1i1p1_200512-210012_2.5x2.5regrid.nc', \
'pr_Amon_HadGEM2-ES_rcp85_r1i1p1_200512-209912_2.5x2.5regrid.nc', \
'pr_Amon_inmcm4_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_IPSL-CM5A-LR_rcp85_r1i1p1_200601-230012_2.5x2.5regrid.nc', \
'pr_Amon_IPSL-CM5A-MR_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_IPSL-CM5B-LR_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_MIROC5_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_MIROC-ESM-CHEM_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_MIROC-ESM_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_MPI-ESM-LR_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_MPI-ESM-MR_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_MRI-CGCM3_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_NorESM1-ME_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc', \
'pr_Amon_NorESM1-M_rcp85_r1i1p1_200601-210012_2.5x2.5regrid.nc' ))

##########################################################################################
##########################################################################################
##########################################################################################

# (YYYY,MM,DD)
hist_start = netCDF4.datetime(1970,1,1)
hist_end = netCDF4.datetime(2001,3,1)

rcp_start = netCDF4.datetime(2069,11,1)
rcp_end = netCDF4.datetime(2099,3,1)

pr_histmonths_list = []
pr_rcpmonths_list = []

pr_histclim_list = []
pr_rcpclim_list = []
pr_anoms_list = []

pr_histstdevs_list = []
pr_rcpstdevs_list = []

pr_pvals_list = []
	
for i in range(len(model_names)):
#for i in [19,20]:

	print("opening model", model_names[i])

	# OPEN HISTORICAL PERIOD PR DATA
	modelname = model_names[i]
	ncfile = Dataset(file_root+pr_hist[i], 'r', format='NETCDF4')
	pr_hist_data = ncfile.variables['pr'][:,:,:]*86400.
	time_variable = ncfile.variables['time']
	print(time_variable.units)
	date_start = netCDF4.date2num(hist_start, time_variable.units, time_variable.calendar)
	date_end = netCDF4.date2num(hist_end, time_variable.units, time_variable.calendar)
	model_time = time_variable[:]
	
	time_variable_converted = netCDF4.num2date(time_variable[:], time_variable.units, time_variable.calendar)
		
	#timespan = model_time[numpy.where((model_time[:]>=date_start)&(model_time<date_end))]
	#pr_hist_data = pr_hist_data[numpy.where((model_time[:]>=date_start)&(model_time<date_end))[0], :,:]
	#ncfile.close()
	
	if season=='ndjfm':		
		time_indices = numpy.array([(t.month in [11,12,1,2,3])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
		pr_hist_data_seas = pr_hist_data[time_indices,:,:]
	
	elif season=='annual':
		pr_hist_data_seas = pr_hist_data[:,:,:]
	pr_histmonths_list.append(pr_hist_data_seas[:,:,:])

	# OPEN RCP PERIOD PR DATA
	#print "opening model", model_names[i]
	modelname = model_names[i]
	ncfile = Dataset(file_root+pr_rcp85[i], 'r', format='NETCDF4')
	pr_rcp_data = ncfile.variables['pr'][:,:,:]*86400.
	time_variable = ncfile.variables['time']
	print(time_variable.units)
	date_start = netCDF4.date2num(rcp_start, time_variable.units, time_variable.calendar)
	date_end = netCDF4.date2num(rcp_end, time_variable.units, time_variable.calendar)
	model_time = time_variable[:]

	time_variable_converted = netCDF4.num2date(time_variable[:], time_variable.units, time_variable.calendar)
	
	if season=='ndjfm':
		time_indices = numpy.array([(t.month in [11,12,1,2,3])&(t.year in range(rcp_start.year, rcp_end.year+1)) for t in time_variable_converted])
		pr_rcp_data_seas = pr_rcp_data[time_indices,:,:]

	elif season=='annual':
		pr_rcp_data_seas = pr_rcp_data[:,:,:]
		
	pr_rcpmonths_list.append(pr_rcp_data_seas[:,:,:])
	print(pr_hist_data_seas.shape)
	print(pr_rcp_data_seas.shape)


print("calculating climatologies")
for i in range(len(model_names)):
	pr_histclim_list.append(numpy.mean(pr_histmonths_list[i], axis=0)[:,:])
	pr_rcpclim_list.append(numpy.mean(pr_rcpmonths_list[i], axis=0)[:,:])
	pr_histstdevs_list.append(numpy.std(pr_histmonths_list[i], axis=0, ddof=1)[:,:])
	pr_rcpstdevs_list.append(numpy.std(pr_rcpmonths_list[i], axis=0, ddof=1)[:,:])
	pr_anoms_list.append(pr_rcpclim_list[-1] - pr_histclim_list[-1])
	pr_pvals_list.append(scipy.stats.ttest_ind(pr_rcpclim_list[i], pr_histmonths_list[i], axis=0))

print("calculating pvals for ENTIRE collection of models")
# must collect all data first...
all_model_data_hist = numpy.zeros((0, global_nlat, global_nlon))
all_model_data_rcp = numpy.zeros((0, global_nlat, global_nlon))
for i in range(len(model_names)):
	all_model_data_hist = numpy.vstack((all_model_data_hist, pr_histmonths_list[i]))
	all_model_data_rcp = numpy.vstack((all_model_data_rcp, pr_rcpmonths_list[i]))

hist_mmem_clim = numpy.mean(numpy.array(pr_histclim_list), axis=0)
rcp_mmem_clim = numpy.mean(numpy.array(pr_rcpclim_list), axis=0)
rcp_minus_hist = rcp_mmem_clim - hist_mmem_clim

all_model_data_hist_clim = numpy.array((pr_histclim_list))
all_model_data_rcp_clim = numpy.array((pr_rcpclim_list))
agreement_positive = numpy.sum(all_model_data_rcp_clim-all_model_data_hist_clim>0, axis=0)

pvals_array = scipy.stats.ttest_ind(all_model_data_hist, all_model_data_rcp, axis=0, equal_var=False)[1]

#cf=mp.contourf(pr_pvals_array)
#mp.colorbar(cf)
#mp.show()
#exit()




# save hist clim
filename = '/home/baird/CMIP5_climatologies_anomalies_agreement.nc'
ncfile = netCDF4.Dataset(filename, 'w', format='NETCDF4')

lat_dim = ncfile.createDimension('lat', pr_lat.size)
lon_dim = ncfile.createDimension('lon', pr_lon.size)

lat_var = ncfile.createVariable('lat', 'f4', ('lat',))
lon_var = ncfile.createVariable('lon', 'f4', ('lon',))
lat_var[:] = pr_lat
lon_var[:] = pr_lon
lat_var.units = 'degrees_north'
lon_var.units = 'degrees_east'

hist_clim_var = ncfile.createVariable('hist_clim', 'f4', ('lat','lon'))
hist_clim_var[:] = hist_mmem_clim
hist_clim_var.units = 'mm day-1'
hist_clim_var.description = 'Historical climatology for NDJFM during 1970-2000, calculated from monthly data'

rcp_clim_var = ncfile.createVariable('rcp_clim', 'f4', ('lat','lon'))
rcp_clim_var[:] = rcp_mmem_clim
rcp_clim_var.units = 'mm day-1'
rcp_clim_var.description = 'RCP8.5 climatology for NDJFM during 2070-2100, calculated from monthly data'

rcp_minus_hist_var = ncfile.createVariable('rcp_minus_hist', 'f4', ('lat','lon'))
rcp_minus_hist_var[:] = rcp_minus_hist
rcp_minus_hist_var.units = 'mm day-1'
rcp_minus_hist_var.description = 'RCP8.5 minus historical climatologies (2070-2100 minus 1970-2000) for NDJFM, calculated from monthly data'

ttest_pvals_var = ncfile.createVariable('ttest_pvals', 'f4', ('lat','lon'))
ttest_pvals_var[:] = pvals_array
ttest_pvals_var.units = 'unitless'
ttest_pvals_var.description = 'P-values from a t-test on the difference of the mean of independent samples (scipy.stats.ttest_ind()) with unequal variances'

agreement_var = ncfile.createVariable('agreement', 'f4', ('lat','lon'))
agreement_var[:] = agreement_positive
agreement_var.units = 'number of ensemble members (out of 40)'
agreement_var.description = 'Number of ensemble members that agrees on a positive end-of-century change'

ncfile.history = 'Created ' + time.ctime(time.time())
ncfile.close()
print(filename, "saved")
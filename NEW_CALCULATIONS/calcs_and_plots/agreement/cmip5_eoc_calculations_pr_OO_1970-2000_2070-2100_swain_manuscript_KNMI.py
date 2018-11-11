import numpy
import scipy
import scipy.stats
import netCDF4
import datetime
import matplotlib
import time
matplotlib.use('TkAgg')
import matplotlib.pyplot as mp

# pr_Amon_ens_rcp85_028.nc

hist_start = datetime.datetime(1970,1,1)
hist_end = datetime.datetime(1999,12,31)

rcp_start = datetime.datetime(2070,1,1)
rcp_end = datetime.datetime(2099,12,31)

working_dir = '/ninod/cmip5_baird/knmi_explorer/'
ensemble_indices = ['{:03}'.format(i) for i in range(0,78)]

ncfile = netCDF4.Dataset(working_dir+'pr_Amon_ens_rcp85_'+ensemble_indices[0]+'.nc', 'r', 'NetCDF4')
lat = ncfile.variables['lat'][:]
lon = ncfile.variables['lon'][:]

ensemble_data = []
ensemble_names = []
ensemble_time_var = []
ensemble_time = []
for i in range(78):
	ncfile = netCDF4.Dataset(working_dir+'pr_Amon_ens_rcp85_'+ensemble_indices[i]+'.nc', 'r', 'NetCDF4')
	ensemble_names.append(ncfile.model_id)
	print(ensemble_names[-1])
	ensemble_data.append(ncfile.variables['pr'][:]*86400.)
	ensemble_time_var.append(ncfile.variables['time'])
	ensemble_time.append(ncfile.variables['time'][:])

ensemble_names = numpy.array((ensemble_names))
ensemble_names_unique = numpy.unique(ensemble_names)

# pull out seasonal data
ensemble_hist_seasons = []
ensemble_rcp_seasons = []
for i in range(ensemble_names.size):
	time_variable_converted = netCDF4.num2date(ensemble_time[i], ensemble_time_var[i].units, ensemble_time_var[i].calendar)

	hist_time_indices = numpy.array([(t.month in [11,12,1,2,3])&(t.year in range(hist_start.year, hist_end.year+1)) for t in time_variable_converted])
	
	rcp_time_indices = numpy.array([(t.month in [11,12,1,2,3])&(t.year in range(rcp_start.year, rcp_end.year+1)) for t in time_variable_converted])
	
	ensemble_hist_seasons.append( ensemble_data[i][hist_time_indices,:,:] )
	ensemble_rcp_seasons.append( ensemble_data[i][rcp_time_indices,:,:] )

# calculate climatologies
ensemble_hist_clim_seasons = [i.mean(axis=0) for i in ensemble_hist_seasons]
ensemble_rcp_clim_seasons = [i.mean(axis=0) for i in ensemble_rcp_seasons]

"""
# save climatologies

for i in range(ensemble_names.size):
	filename = '/ninod/cmip5_baird/knmi_explorer/climatology_calculations/' + ensemble_indices[i] + '_' + ensemble_names[i] + '_pr_1970-2000_and_2070-2100_climatology_NDJFM.nc'
	ncfile = netCDF4.Dataset(filename, 'w', format='NETCDF4')
	ncfile.model_id = ensemble_names[i]
	
	lat_dim = ncfile.createDimension('lat', lat.size)
	lon_dim = ncfile.createDimension('lon', lon.size)
	lat_var = ncfile.createVariable('lat', 'f4', ('lat',))
	lon_var = ncfile.createVariable('lon', 'f4', ('lon',))
	lat_var[:] = lat
	lon_var[:] = lon

	pr_var = ncfile.createVariable('pr_hist', 'f4', ('lat','lon'))
	pr_var[:] = ensemble_hist_clim_seasons[i]
	lat_var.units = 'degrees_north'
	lon_var.units = 'degrees_east'
	pr_var.units = 'mm/day'
	pr_var.description = str(hist_start.year)+'-'+str(hist_end.year)
		
	pr_var = ncfile.createVariable('pr_rcp', 'f4', ('lat','lon'))
	pr_var[:] = ensemble_rcp_clim_seasons[i]
	lat_var.units = 'degrees_north'
	lon_var.units = 'degrees_east'
	pr_var.units = 'mm/day'
	pr_var.description = str(rcp_start.year)+'-'+str(rcp_end.year)
	
	pr_var = ncfile.createVariable('pr_anom', 'f4', ('lat','lon'))
	pr_var[:] = ensemble_rcp_clim_seasons[i] - ensemble_hist_clim_seasons[i]
	lat_var.units = 'degrees_north'
	lon_var.units = 'degrees_east'
	pr_var.units = 'mm/day'
	pr_var.description = str(rcp_start.year)+'-'+str(rcp_end.year)+' minus '+str(hist_start.year)+'-'+str(hist_end.year)

	ncfile.history = 'Created ' + time.ctime(time.time())
	ncfile.close()
	print(filename, "saved")
"""
# 
# 
# # take mean of realizations
# ensemble_fields_unique = numpy.zeros((ensemble_names_unique.size, lat.size, lon.size))
# hist_clim_unique = numpy.zeros((ensemble_names_unique.size, lat.size, lon.size))
# rcp_clim_unique = numpy.zeros((ensemble_names_unique.size, lat.size, lon.size))
# #for i in range(model_names_unique.size):
# n_members = numpy.array(([numpy.sum(model_names==model_names_unique[i]) for i in range(model_names_unique.size)]))
# print(n_members)
# 
# for i in range(n_members.size):
#     start_index = numpy.sum(n_members[0:i])
#     end_index = numpy.sum(n_members[0:i+1])
#     #print(start_index, end_index)
#     ensemble_fields_unique[i,:,:] = numpy.mean(ensemble_fields[start_index:end_index,:,:], axis=0)
# 


# now take average of the ensembles that matter






ensemble_hist_seasons = numpy.array((ensemble_hist_seasons))
ensemble_rcp_seasons = numpy.array((ensemble_rcp_seasons))

ensemble_hist_clim_seasons = numpy.array((ensemble_hist_clim_seasons))
ensemble_rcp_clim_seasons = numpy.array((ensemble_rcp_clim_seasons))
anomaly_clim_seasons = ensemble_rcp_clim_seasons - ensemble_hist_clim_seasons

hist_clim_unique = numpy.zeros((ensemble_names_unique.size, lat.size, lon.size))
rcp_clim_unique = numpy.zeros((ensemble_names_unique.size, lat.size, lon.size))
anomaly_clim_unique = numpy.zeros((ensemble_names_unique.size, lat.size, lon.size))

n_members = numpy.array(([numpy.sum(ensemble_names==ensemble_names_unique[i]) for i in range(ensemble_names_unique.size)]))


print(ensemble_hist_seasons.shape, ensemble_rcp_seasons.shape) # 78x150x72x144

ensemble_hist_seasons_unique = numpy.zeros((ensemble_names_unique.size, ensemble_hist_seasons.shape[1], ensemble_hist_seasons.shape[2], ensemble_hist_seasons.shape[3])) #35x150x72x144

ensemble_rcp_seasons_unique = numpy.zeros((ensemble_names_unique.size, ensemble_rcp_seasons.shape[1], ensemble_rcp_seasons.shape[2], ensemble_rcp_seasons.shape[3])) #35x150x72x144

print(ensemble_rcp_seasons_unique.shape)

# take ensemble means of full fields
for i in range(n_members.size):
	start_index = numpy.sum(n_members[0:i])
	end_index = numpy.sum(n_members[0:i+1])
	ensemble_hist_seasons_unique[i,:,:,:] = numpy.mean(ensemble_hist_seasons[start_index:end_index,:,:,:], axis=0)
	ensemble_rcp_seasons_unique[i,:,:,:] = numpy.mean(ensemble_rcp_seasons[start_index:end_index,:,:,:], axis=0)

# take means of climatological fields
for i in range(n_members.size):
	start_index = numpy.sum(n_members[0:i])
	end_index = numpy.sum(n_members[0:i+1])
	hist_clim_unique[i,:,:] = numpy.mean(ensemble_hist_clim_seasons[start_index:end_index,:,:], axis=0)
	rcp_clim_unique[i,:,:] = numpy.mean(ensemble_rcp_clim_seasons[start_index:end_index,:,:], axis=0)
	anomaly_clim_unique[i,:,:] = numpy.mean(anomaly_clim_seasons[start_index:end_index,:,:], axis=0)

# now take means of unique fields and do statistical significance test
hist_mmem_clim = numpy.mean(hist_clim_unique, axis=0)
rcp_mmem_clim = numpy.mean(rcp_clim_unique, axis=0)
rcp_minus_hist = rcp_mmem_clim - hist_mmem_clim
agreement_positive = numpy.sum(anomaly_clim_unique>0, axis=0)
pvals_array = scipy.stats.ttest_ind(ensemble_hist_seasons_unique.reshape((-1,72,144)), ensemble_rcp_seasons_unique.reshape((-1,72,144)), axis=0, equal_var=False)[1]

# NOW MAKE MULTIMODEL ENSEMBLE MEAN
# save hist clim

filename = '/home/baird/CMIP5_climatologies_anomalies_agreement_KNMI.nc'
ncfile = netCDF4.Dataset(filename, 'w', format='NETCDF4')

ncfile.description = '35 separate models, consisting of 78 total runs, were used in these calculations.  The information is contained in model_id and ensemble_count below.'
ncfile.model_id = ', '.join([i for i in ensemble_names_unique])
ncfile.ensemble_count = ', '.join([str(i) for i in n_members])

lat_dim = ncfile.createDimension('lat', lat.size)
lon_dim = ncfile.createDimension('lon', lon.size)

lat_var = ncfile.createVariable('lat', 'f4', ('lat',))
lon_var = ncfile.createVariable('lon', 'f4', ('lon',))
lat_var[:] = lat
lon_var[:] = lon
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
import numpy
import scipy
import netCDF4
import matplotlib.pyplot as mp
import matplotlib.ticker
import matplotlib.colors
import scipy.stats
import pandas
import itertools
import datetime
import os

mp.rcParams.update({'mathtext.default': 'regular'})
#get_ipython().magic('matplotlib inline')

# ======================================================================

PRECT_nlat = 26
PRECT_nlon = 25

latlon_indices = list(itertools.product(range(PRECT_nlat), range(PRECT_nlon)))
region = 'whole_domain'
window=1

PRECT_lat = numpy.load('/Users/baird/Dropbox/_analysis/attribution_2017/NEW_CALCULATIONS/npy_files/PRECT_lat.npy')
PRECT_lon = numpy.load('/Users/baird/Dropbox/_analysis/attribution_2017/NEW_CALCULATIONS/npy_files/PRECT_lon.npy')

# ======================================================================

#working_dir = '/Users/baird/Dropbox/_analysis/attribution_2017/NEW_CALCULATIONS/npy_files/'

#latlon_indices = numpy.load(working_dir + 'ncal_latlon_indices_array.npy'); region='ncal'
#latlon_indices = numpy.load(working_dir + 'ccal_latlon_indices_array.npy'); region='ccal'
#latlon_indices = numpy.load(working_dir + 'scal_latlon_indices_array.npy'); region='scal'

threshold = 0.1

year_start_pic = 402 #time_subsets[chunk,0]
year_end_pic = 2200 #time_subsets[chunk,1]

year_start = 1920
year_end = 2100

lo_perc = 20
hi_perc = 80

#year_start_whiplash = 1920
#year_end_whiplash = 1950

year_start_list = numpy.arange(1920,2100)
#year_end_list = numpy.arange(1950,2101)

whiplash_ratios_all = numpy.zeros(( year_start_list.size, len(latlon_indices) ))
whiplash_counts_pic_all = numpy.zeros(( year_start_list.size, len(latlon_indices) ))
whiplash_counts_rcp_all = numpy.zeros(( year_start_list.size, len(latlon_indices) ))

# create list of names of members '001','002','003', ...
ensemble_members = numpy.hstack((numpy.arange(1,36), numpy.arange(101,106)))
ensemble_names = ['{:03d}'.format(i) for i in ensemble_members]

working_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/calculated_npy_files/whole_domain/'

# create season strings for PIC
years_pic = numpy.arange(year_start_pic, year_end_pic+1, 1).astype(numpy.int)
half_years_pic = numpy.arange(year_start_pic+0.75, year_end_pic, 1)

season_strings_pic = [str(years_pic[i])+'-'+str(years_pic[i+1]) for i in range(years_pic.size-1)]
member_strings_pic = ['{:03d}'.format(i) for i in range(1,36)]
n_seasons_pic=year_end_pic-year_start_pic

# ========== open pic
ncfile = netCDF4.Dataset('/Users/baird/Dropbox/_analysis/attribution_2017/NEW_CALCULATIONS/calcs_and_plots/create_ncfiles_seasonal/'+\
                       'seasonal_totals_pic.nc', 'r', 'NETCDF4')
seasonal_totals_pic = ncfile.variables['seasonal_total'][:]
seasonal_totals_pic = seasonal_totals_pic.reshape((1798,-1))

# ========== open hist+rcp
ncfile = netCDF4.Dataset('/Users/baird/Dropbox/_analysis/attribution_2017/NEW_CALCULATIONS/calcs_and_plots/create_ncfiles_seasonal/'+\
                       'seasonal_totals_hist_rcp.nc', 'r', 'NETCDF4')
seasonal_totals_hist_rcp = ncfile.variables['seasonal_total'][:]
seasonal_totals_hist_rcp = seasonal_totals_hist_rcp.reshape((180,40,-1))

whiplash_hist_rcp_markers = numpy.zeros((180,40,len(latlon_indices)), dtype=numpy.int)
whiplash_pic_markers = numpy.zeros((1798,len(latlon_indices)), dtype=numpy.int)


for latlon_idx in range(len(latlon_indices)):

	if latlon_idx%10==0:
		print('latlon_idx', latlon_idx)

	# ======================================================================
	# Open preindustrial control info (dict_pic)
	seasonal_total_pic = seasonal_totals_pic[:,latlon_idx]
	#print(seasonal_total_pic)
	
	# ======================================================================
	# ======================================================================
	# ======================================================================
	# whiplash calculation setup
	# get top and bottom 10th percentiles in PIC
	# then see how often it transitions from 1 or 2 seasons with that

	pic_hi = numpy.percentile(seasonal_total_pic, hi_perc)

	# take 3 year running mean for pic_lo
	seasonal_total_pic_window_mean = numpy.array(pandas.Series(seasonal_total_pic).rolling(window=window).mean())
	pic_lo = numpy.nanpercentile(seasonal_total_pic_window_mean, lo_perc)

	#print(seasonal_totals_hist_rcp.shape)
	#print(whiplash_hist_rcp_markers.shape)

	#for ens_idx in range(40):
	#	whiplash_hist_rcp_markers[:,:,latlon_idx][ seasonal_totals_hist_rcp[:,ens_idx,latlon_idx]<pic_lo ] = -1
	whiplash_hist_rcp_markers[:,:,latlon_idx][ seasonal_totals_hist_rcp[:,:,latlon_idx]<pic_lo ] = -1
	whiplash_hist_rcp_markers[:,:,latlon_idx][ seasonal_totals_hist_rcp[:,:,latlon_idx]>pic_hi ] = 1


# calculate percent change
# if the pic starts at 20th percentile across entire domain, then find the 20th percentile across all ensemble members
# same with 80th percentile

ratio_wet_dry = numpy.zeros((180, len(latlon_indices) ))
wet_total = numpy.zeros((180, len(latlon_indices) ))
dry_total = numpy.zeros((180, len(latlon_indices) ))

for latlon_idx in range(len(latlon_indices)):
	wet_total[:,latlon_idx] = numpy.sum(whiplash_hist_rcp_markers[:,:,latlon_idx]==1, axis=1)
	dry_total[:,latlon_idx] = numpy.sum(whiplash_hist_rcp_markers[:,:,latlon_idx]==-1, axis=1)

wet_total_fraction = wet_total/40.
dry_total_fraction = dry_total/40.

wet_total_fraction = wet_total_fraction.reshape((180, PRECT_nlat, PRECT_nlon))
dry_total_fraction = dry_total_fraction.reshape((180, PRECT_nlat, PRECT_nlon))

print(wet_total_fraction)
print(wet_total_fraction.shape)
print()
print(dry_total_fraction)









# ========= save netcdf files of -1,1 whiplash values for hist+RCP8.5
year_list = numpy.arange(1921,2101)
time_datetime = [datetime.datetime(i,1,15) for i in year_list]
time_nc = netCDF4.date2num(time_datetime, units='days since 1920-01-01', calendar='standard')





# save hist clim
filename = 'hist_rcp_20_low_80_high.nc'

if os.path.exists(filename):
    print('file exists')
    os.remove(filename)

ncfile = netCDF4.Dataset(filename, 'w', format='NETCDF4')

time_dim = ncfile.createDimension('time', None)
time_var = ncfile.createVariable('time', 'f4', ('time',))
time_var[:] = time_nc
time_var.units = 'days since 1920-01-01'

lat_dim = ncfile.createDimension('lat', PRECT_nlat)
lat_var = ncfile.createVariable('lat', 'f4', ('lat',))
lat_var[:] = PRECT_lat
lat_var.units = 'degrees North'

lon_dim = ncfile.createDimension('lon', PRECT_nlon)
lon_var = ncfile.createVariable('lon', 'f4', ('lon',))
lon_var[:] = PRECT_lon
lon_var.units = 'degrees East'

ens_dim = ncfile.createDimension('ensemble', len(ensemble_names))
ens_var = ncfile.createVariable('ensemble', 'f4', ('ensemble',))
ens_var[:] = ensemble_names
ens_var.units = 'NCAR LENS ensemble member'

data_var = ncfile.createVariable('seasonal_total', 'f4', ('time','ensemble','lat','lon',))
data_var.units = 'a 20/80 whiplash event is when a gridpoint has a -1 followed by a +1'
data_var[:] = whiplash_hist_rcp_markers

ncfile.close()

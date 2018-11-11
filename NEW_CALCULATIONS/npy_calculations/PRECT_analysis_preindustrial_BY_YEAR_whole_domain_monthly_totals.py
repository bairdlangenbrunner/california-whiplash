import numpy
import scipy
import netCDF4
import pandas
import os
#import matplotlib.pyplot as mp
#import matplotlib.ticker

#mp.rcParams.update({'mathtext.default': 'regular'})
#get_ipython().magic('matplotlib inline')

################################################################################
################################################################################
################################################################################

time_subsets = numpy.array((
[402,499], \
))

################################################################################
################################################################################
################################################################################
file_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/B1850C5CN_CA_REGION/'

file = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.04020101-22001231_CA_REGION.nc'
#file_list = numpy.array(( \
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.04020101-04991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.05000101-05991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.06000101-06991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.07000101-07991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.08000101-08991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.09000101-09991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.10000101-10991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.11000101-11991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.12000101-12991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.13000101-13991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.14000101-14991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.15000101-15991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.16000101-16991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.17000101-17991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.18000101-18991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.19000101-19991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.20000101-20991231_SUBSET.nc
#b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.21000101-22001231_SUBSET.nc ))

################################################################################
################################################################################
################################################################################
#d_oct = 31
d_nov = 30
d_dec = 31
d_jan = 31
d_feb = 28
d_mar = 31
days_per_season = d_nov+d_dec+d_jan+d_feb+d_mar
#print(days_per_season)

################################################################################
################################################################################
################################################################################
#i = 0

#ncfile = netCDF4.Dataset(file_dir + file)
#PRECT_lat = ncfile.variables['lat'][:]
#PRECT_lon = ncfile.variables['lon'][:]
#PRECT_time_var = ncfile.variables['time']

#PRECT_time_dates = netCDF4.num2date(PRECT_time_var[:], PRECT_time_var.units, PRECT_time_var.calendar)
#time_indices_NDJFM = numpy.array([t.month in [11,12,1,2,3] for t in PRECT_time_dates], dtype=bool)
#PRECT_time_dates_NDJFM = PRECT_time_dates[time_indices_NDJFM]

#PRECT_NDJFM_CA = ncfile.variables['PRECT'][time_indices_NDJFM,(LA_lat_idx-2):(LA_lat_idx+2), (LA_lon_idx-2):(LA_lon_idx+2)]*86400.*1000
#PRECT_NDJFM_CA = ncfile.variables['PRECT'][time_indices_NDJFM,LA_lat_idx, LA_lon_idx]*86400.*1000

################################################################################
################################################################################
################################################################################
working_dir = '/Users/baird/Dropbox/_analysis/attribution_2017/NEW_CALCULATIONS/npy_files/'
save_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/calculated_npy_files/whole_domain_MONTHLY_TOTALS/'

threshold=0.1
region='whole_domain'

################################################################################
################################################################################
################################################################################
year_start = 402 #time_subsets[chunk,0]
year_end = 2200 #time_subsets[chunk,1]

# create season strings
years = numpy.arange(year_start, year_end+1, 1).astype(numpy.int)
season_strings = numpy.empty(years.size-1, dtype=numpy.str)

season_strings = [str(years[i])+'-'+str(years[i+1]) for i in range(years.size-1)]

################################################################################
################################################################################
################################################################################
# open up example file
ncfile = netCDF4.Dataset(file_dir + file)
PRECT_lat = ncfile.variables['lat'][:]
PRECT_lon = ncfile.variables['lon'][:]
PRECT_time_var = ncfile.variables['time']
PRECT_time_dates = netCDF4.num2date(PRECT_time_var[:], PRECT_time_var.units, PRECT_time_var.calendar)
PRECT_total = ncfile.variables['PRECT'][:, :, :]*86400.*1000

n_seasons = year_end-year_start
years = numpy.arange(year_start, year_end+1)
PRECT_monthly_data = numpy.zeros((PRECT_lat.size, PRECT_lon.size, year_end-year_start+1, 12))

for i in [0]:
	save_file = 'member_005_years_'+'{:04d}'.format(year_start)+'-'+'{:04d}'.format(year_end)+'_threshold_'+str(threshold)+'mmday_'+region+'.npy'
	if os.path.isfile(save_dir+save_file):
		print(save_file + ' exists')
		pass#continue

	for yr_idx in range(year_end-year_start+1):
		print(yr_idx+1)
		time_indices_year = numpy.array([(t.year==years[yr_idx]) for t in PRECT_time_dates], dtype=bool)
		time_dates_year = PRECT_time_dates[time_indices_year]
		PRECT_region = PRECT_total[time_indices_year, :, :]
		PRECT_region[PRECT_region<threshold] = 0
		
		for lat_idx in range(PRECT_lat.size):
			for lon_idx in range(PRECT_lon.size):
				for month_idx in range(12):
					time_indices_month = numpy.array([(t.month==month_idx+1) for t in time_dates_year], dtype=bool)
					PRECT_monthly_data[lat_idx,lon_idx,yr_idx,month_idx] = numpy.sum(PRECT_region[time_indices_month,lat_idx,lon_idx])
					
	numpy.save(save_dir + save_file, PRECT_monthly_data)
	print('Saved '+save_file)
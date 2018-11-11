import numpy
import scipy
import netCDF4
import multiprocessing
#import matplotlib.pyplot as mp
#import matplotlib.ticker

#mp.rcParams.update({'mathtext.default': 'regular'})
#get_ipython().magic('matplotlib inline')

################################################################################
################################################################################
################################################################################
def geo_idx(dd, dd_array):
   """
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
    """
   geo_idx = (numpy.abs(dd_array - dd)).argmin()
   return geo_idx

LA_lat = 34.0522
LA_lon = 118.2437 # deg west
LA_lon = 180. + (180-LA_lon)

Oroville_dam_lat = 39.5380
Oroville_dam_lon = 121.4831 # deg west
Oroville_dam_lon = 360 - Oroville_dam_lon

################################################################################
################################################################################
################################################################################
file_dir = '/ninod/NCAR_LENS/daily/PRECT/hist_and_rcp85/'

file_list = numpy.array(( 'hist_and_rcp85.001.PRECT.18500101-21001231.nc', \
'hist_and_rcp85.002.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.003.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.004.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.005.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.006.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.007.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.008.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.009.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.010.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.011.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.012.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.013.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.014.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.015.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.016.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.017.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.018.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.019.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.020.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.021.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.022.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.023.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.024.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.025.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.026.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.027.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.028.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.029.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.030.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.031.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.032.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.033.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.034.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.035.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.101.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.102.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.103.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.104.PRECT.19200101-21001231.nc', \
'hist_and_rcp85.105.PRECT.19200101-21001231.nc' ))

# ### Pull out proper lat/lon (doing this for one grid point right now)

################################################################################
################################################################################
################################################################################
d_oct = 31
d_nov = 30
d_dec = 31
d_jan = 31
d_feb = 28
d_mar = 31
days_per_season = d_oct+d_nov+d_dec+d_jan+d_feb+d_mar

################################################################################
################################################################################
################################################################################

i = 0

ncfile = netCDF4.Dataset(file_dir + file_list[i])
PRECT_lat = ncfile.variables['lat'][:]
PRECT_lon = ncfile.variables['lon'][:]
PRECT_time_var = ncfile.variables['time']

PRECT_time_dates = netCDF4.num2date(PRECT_time_var[:], PRECT_time_var.units, PRECT_time_var.calendar)
time_indices_ONDJFM = numpy.array([t.month in [10,11,12,1,2,3] for t in PRECT_time_dates], dtype=bool)
PRECT_time_dates_ONDJFM = PRECT_time_dates[time_indices_ONDJFM]

LA_lat_idx = geo_idx(LA_lat, PRECT_lat)
LA_lon_idx = geo_idx(LA_lon, PRECT_lon)

#PRECT_ONDJFM_CA = ncfile.variables['PRECT'][time_indices_ONDJFM,(LA_lat_idx-2):(LA_lat_idx+2), (LA_lon_idx-2):(LA_lon_idx+2)]*86400.*1000
PRECT_ONDJFM_CA = ncfile.variables['PRECT'][time_indices_ONDJFM,LA_lat_idx, LA_lon_idx]*86400.*1000


################################################################################
################################################################################
################################################################################

# California lat/lon combos
#ca_latlon_array = numpy.load('ca_latlon_array.npy')

# CALCULATE GLOBAL VALUES WHERE THESE EXIST
#global_lat_indices = [numpy.where(numpy.isclose(ca_latlon_array[i,0], PRECT_lat))[0][0] for i in range(ca_latlon_array.shape[0])]
#global_lon_indices = [numpy.where(numpy.isclose(ca_latlon_array[i,1], PRECT_lon))[0][0] for i in range(ca_latlon_array.shape[0])]

#ca_latlon_indices_array = numpy.column_stack((global_lat_indices,global_lon_indices))
#numpy.save('ca_latlon_indices_array.npy', full_ca_latlon_array)
#PRECT_ONDJFM_CA = PRECT_ONDJFM[:, (LA_lat_idx-2):(LA_lat_idx+2), (LA_lon_idx-2):(LA_lon_idx+2)]

working_dir = '/ninod/baird/cmip5/cmip5_calculations/attribution_2017/storm_counting/'
ca_latlon_indices_array = numpy.load(working_dir + 'ca_latlon_indices_array.npy')

threshold=0.5

# for each season, do counts, store them in a list

################################################################################
################################################################################
################################################################################
year_start = 1920 #time_subsets[chunk,0]
year_end = 2100 #time_subsets[chunk,1]

# create season strings
years = numpy.arange(year_start, year_end+1, 1).astype(numpy.int)
season_strings = numpy.empty(years.size-1, dtype=numpy.str)

season_strings = [str(years[i])+'-'+str(years[i+1]) for i in range(years.size-1)]

################################################################################
################################################################################
################################################################################
def main():
	pool = multiprocessing.Pool()
	input = zip(ca_latlon_indices_array[:,0], ca_latlon_indices_array[:,1])
	pool.map(storm_count, input)

n_seasons = year_end-year_start
save_dir = '/ninod/baird/cmip5/cmip5_calculations/attribution_2017/storm_counting_parallel/'


def storm_count(args):
	lat_idx, lon_idx = args
	print(lat_idx, lon_idx)
	#for latlon in range(ca_latlon_indices_array.shape[0]):
	#lat_idx = ca_latlon_indices_array[latlon,0]
	#lon_idx = ca_latlon_indices_array[latlon,1]
	save_dict = {}
	for i in range(file_list.size):
		#print(file_list[i])
		ensemble_member = file_list[i].split('.')[1]
		ncfile = netCDF4.Dataset(file_dir + file_list[i])
		PRECT_lat = ncfile.variables['lat'][:]
		PRECT_lon = ncfile.variables['lon'][:]
		PRECT_time_var = ncfile.variables['time']
		PRECT_time_dates = netCDF4.num2date(PRECT_time_var[:], PRECT_time_var.units, PRECT_time_var.calendar)
		time_indices_ONDJFM = numpy.array([(t.month in [10,11,12,1,2,3])&(t.year in range(year_start,year_end+1)) for t in PRECT_time_dates], dtype=bool)
		PRECT_time_dates_ONDJFM = PRECT_time_dates[time_indices_ONDJFM]
		PRECT_ONDJFM_region = ncfile.variables['PRECT'][time_indices_ONDJFM, lat_idx, lon_idx]*86400.*1000
		ncfile.close()
		# get index of first October
		# and index of last March
		t=0
		while t < PRECT_time_dates_ONDJFM.size:
			if PRECT_time_dates_ONDJFM[t].month!=10:
				t+=1
			else:
				first_October_index = t
				break
		t=0
		while i < PRECT_time_dates_ONDJFM.size:
			if PRECT_time_dates_ONDJFM[::-1][t].month!=3:
				t+=1
			else:
				last_March_index = PRECT_time_dates_ONDJFM.size-t-1
				break
		PRECT_time_dates_ONDJFM_fullseasons = PRECT_time_dates_ONDJFM[first_October_index:last_March_index+1]
		PRECT_ONDJFM_region_fullseasons = PRECT_ONDJFM_region[first_October_index:last_March_index+1]
		#years_array = numpy.array(([t.year for t in PRECT_time_dates_ONDJFM_fullseasons]), dtype=numpy.int)
		storm_count=0 # per season
		storm_length=0 # for each storm
		storm_magnitude=0.0 # for each storm
		for s in range(n_seasons):
			storm_count=0
			storm_magnitude_list = []
			storm_length_list = []
			seas_PRECT = PRECT_ONDJFM_region_fullseasons[s*days_per_season:(s*days_per_season+days_per_season)]
			for d in range(days_per_season):
				if seas_PRECT[d]>=threshold:
					storm_magnitude+=seas_PRECT[d]
					storm_length+=1 # days
				elif (seas_PRECT[d]<threshold)&(storm_magnitude>0.0):
					storm_magnitude_list.append(storm_magnitude)
					storm_length_list.append(storm_length)
					storm_magnitude=0.0
					storm_length=0
					storm_count+=1
			seasonal_total = numpy.sum(seas_PRECT[seas_PRECT>threshold])
			precipitation_days = seas_PRECT[seas_PRECT>threshold]
			save_dict[season_strings[s]] = {'storm_count' : storm_count}
			save_dict[season_strings[s]]['storm_magnitude_list'] = storm_magnitude_list
			save_dict[season_strings[s]]['storm_length_list'] = storm_length_list
			save_dict[season_strings[s]]['seasonal_total'] = seasonal_total
			save_dict[season_strings[s]]['precipitation_days'] = precipitation_days

		save_file = 'member_'+ensemble_member+'_latidx_'+str(lat_idx)+'_lonidx_'+str(lon_idx)+'_years_'+str(year_start)+'-'+str(year_end)+'_threshold_'+str(threshold)+'mmday_PARALLEL.npy'
		#numpy.save(save_dir + save_file, save_dict)
		#print('Saved '+save_file)

main()

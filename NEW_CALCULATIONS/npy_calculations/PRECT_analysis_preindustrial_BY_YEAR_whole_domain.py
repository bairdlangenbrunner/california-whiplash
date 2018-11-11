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
save_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/calculated_npy_files/whole_domain/'

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
time_indices_NDJFM = numpy.array([(t.month in [11,12,1,2,3])&(t.year in range(year_start,year_end+1)) for t in PRECT_time_dates], dtype=bool)
PRECT_time_dates_NDJFM = PRECT_time_dates[time_indices_NDJFM]


n_seasons = year_end-year_start

for lat_idx in range(PRECT_lat.size):

	for lon_idx in range(PRECT_lon.size):

		save_dict = {}

		save_file = 'member_005_latidx_'+'{:02d}'.format(lat_idx)+'_lonidx_'+'{:02d}'.format(lon_idx)+'_years_'+'{:04d}'.format(year_start)+'-'+'{:04d}'.format(year_end)+'_threshold_'+str(threshold)+'mmday_'+region+'.npy'
		if os.path.isfile(save_dir+save_file):
			print(save_file + ' exists')
			pass#continue

		PRECT_NDJFM_region = ncfile.variables['PRECT'][time_indices_NDJFM, lat_idx, lon_idx]*86400.*1000

		# get index of first November
		# and index of last March
		t=0
		while t < PRECT_time_dates_NDJFM.size:
			if PRECT_time_dates_NDJFM[t].month!=11:
				t+=1
			else:
				first_November_index = t
				break
		t=0
		while t < PRECT_time_dates_NDJFM.size:
			if PRECT_time_dates_NDJFM[::-1][t].month!=3:
				t+=1
			else:
				last_March_index = PRECT_time_dates_NDJFM.size-t-1
				break

		PRECT_time_dates_NDJFM_fullseasons = PRECT_time_dates_NDJFM[first_November_index:last_March_index+1]
		PRECT_NDJFM_region_fullseasons = PRECT_NDJFM_region[first_November_index:last_March_index+1]

		storm_count=0 # per season
		storm_length=0 # for each storm
		storm_magnitude=0.0 # for each storm

		# loop through seasons
		# count storms
	
		for s in range(n_seasons):
			storm_count=0
			storm_magnitude_list = []
			storm_length_list = []
			seas_PRECT = PRECT_NDJFM_region_fullseasons[s*days_per_season:(s*days_per_season+days_per_season)]
			
			# weird fact:  the 402-499 data are missing Dec. 31st, so need to fix that by looping over days_per_season-1
			for d in range(seas_PRECT.size):
				if seas_PRECT[d]>=threshold:
					storm_magnitude+=seas_PRECT[d]
					storm_length+=1 # days
				elif (seas_PRECT[d]<threshold)&(storm_magnitude>0.0):
					storm_magnitude_list.append(storm_magnitude)
					storm_length_list.append(storm_length)
					storm_magnitude=0.0
					storm_length=0
					storm_count+=1
			
			# consecutive dry days (use 1mm version)
			cdd_count=0
			cdd_list=[]
			cdd_threshold=1 #mm/day
			# weird fact:  the 402-499 data are missing Dec. 31st, so need to fix that by looping over days_per_season-1
			for d in range(seas_PRECT.size):
				if seas_PRECT[d]<cdd_threshold:
					cdd_count+=1 # days
				elif (seas_PRECT[d]>cdd_threshold)&(cdd_count>0):
					cdd_list.append(cdd_count)
					cdd_count=0

			seas_PRECT_with_zeros = seas_PRECT.copy()
			seas_PRECT_with_zeros[seas_PRECT_with_zeros<threshold] = 0
			dry_days = numpy.sum(seas_PRECT<threshold)
			running_40d_sum = list(pandas.Series(seas_PRECT_with_zeros).rolling(window=40).sum())
			seasonal_total = numpy.sum(seas_PRECT[seas_PRECT>=threshold])
			#precipitation_days = seas_PRECT[seas_PRECT>=threshold]
		
			# calculate monthly sums
			sum_nov = sum(seas_PRECT_with_zeros[0 : d_nov])
			sum_dec = sum(seas_PRECT_with_zeros[d_nov : d_nov+d_dec])
			sum_jan = sum(seas_PRECT_with_zeros[d_nov+d_dec : d_nov+d_dec+d_jan])
			sum_feb = sum(seas_PRECT_with_zeros[d_nov+d_dec+d_jan : d_nov+d_dec+d_jan+d_feb])
			sum_mar = sum(seas_PRECT_with_zeros[d_nov+d_dec+d_jan+d_feb : d_nov+d_dec+d_jan+d_feb+d_mar])
			
			# calculate monthly dry days
			sum_nov_dd = sum(seas_PRECT[0 : d_nov]<threshold)
			sum_dec_dd = sum(seas_PRECT[d_nov : d_nov+d_dec]<threshold)
			sum_jan_dd = sum(seas_PRECT[d_nov+d_dec : d_nov+d_dec+d_jan]<threshold)
			sum_feb_dd = sum(seas_PRECT[d_nov+d_dec+d_jan : d_nov+d_dec+d_jan+d_feb]<threshold)
			sum_mar_dd = sum(seas_PRECT[d_nov+d_dec+d_jan+d_feb : d_nov+d_dec+d_jan+d_feb+d_mar]<threshold)
		
			save_dict[season_strings[s]] = {'storm_count' : storm_count}
			save_dict[season_strings[s]]['storm_magnitude_list'] = storm_magnitude_list
			save_dict[season_strings[s]]['storm_length_list'] = storm_length_list
			save_dict[season_strings[s]]['seasonal_total'] = seasonal_total
			#save_dict[season_strings[s]]['precipitation_days'] = precipitation_days
			save_dict[season_strings[s]]['running_40d_sum'] = running_40d_sum
			save_dict[season_strings[s]]['monthly_totals'] = [sum_nov, sum_dec, sum_jan, sum_feb, sum_mar]
			save_dict[season_strings[s]]['monthly_dry_days'] = [sum_nov_dd, sum_dec_dd, sum_jan_dd, sum_feb_dd, sum_mar_dd]
			save_dict[season_strings[s]]['seasonal_dry_days'] = dry_days
			if len(cdd_list)==0:
				save_dict[season_strings[s]]['cdd_max'] = 0
			else:
				save_dict[season_strings[s]]['cdd_max'] = numpy.max(cdd_list)

		save_file = 'member_005_latidx_'+'{:02d}'.format(lat_idx)+'_lonidx_'+'{:02d}'.format(lon_idx)+'_years_'+'{:04d}'.format(year_start)+'-'+'{:04d}'.format(year_end)+'_threshold_'+str(threshold)+'mmday_'+region+'.npy'	
		numpy.save(save_dir + save_file, save_dict)
		print('Saved '+save_file)
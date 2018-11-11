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
file_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/hist_and_rcp85_CA_REGION/'

file_list = numpy.array(( \
'hist_and_rcp85.001.PRECT.18500101-21001231_CA_REGION.nc', \
'hist_and_rcp85.002.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.003.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.004.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.005.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.006.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.007.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.008.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.009.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.010.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.011.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.012.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.013.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.014.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.015.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.016.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.017.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.018.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.019.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.020.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.021.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.022.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.023.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.024.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.025.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.026.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.027.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.028.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.029.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.030.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.031.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.032.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.033.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.034.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.035.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.101.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.102.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.103.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.104.PRECT.19200101-21001231_CA_REGION.nc', \
'hist_and_rcp85.105.PRECT.19200101-21001231_CA_REGION.nc' ))

# ### Pull out proper lat/lon (doing this for one grid point right now)

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

################################################################################
################################################################################
################################################################################
working_dir = '/Users/baird/Dropbox/_analysis/attribution_2017/NEW_CALCULATIONS/npy_files/'
save_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/calculated_npy_files/whole_domain/'

#ca_latlon_regional_indices_array = numpy.load(working_dir + 'ncal_latlon_indices_array.npy'); region='ncal'
#ca_latlon_regional_indices_array = numpy.load(working_dir + 'ccal_latlon_indices_array.npy'); region='ccal'
#ca_latlon_regional_indices_array = numpy.load(working_dir + 'scal_latlon_indices_array.npy'); region='scal'

threshold=0.1
region='whole_domain'
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
# open up example file
ncfile = netCDF4.Dataset(file_dir + file_list[1])
PRECT_lat = ncfile.variables['lat'][:]
PRECT_lon = ncfile.variables['lon'][:]
PRECT_time_var = ncfile.variables['time']
PRECT_time_dates = netCDF4.num2date(PRECT_time_var[:], PRECT_time_var.units, PRECT_time_var.calendar)
print(PRECT_time_dates.size)
time_indices_NDJFM = numpy.array([(t.month in [11,12,1,2,3])&(t.year in range(year_start,year_end+1)) for t in PRECT_time_dates], dtype=bool)
PRECT_time_dates_NDJFM = PRECT_time_dates[time_indices_NDJFM]
ncfile.close()

n_seasons = year_end-year_start

print(time_indices_NDJFM.shape)

for lat_idx in range(PRECT_lat.size):

	for lon_idx in range(PRECT_lon.size):

		save_dict = {}

		for i in range(30,40):#range(40):#file_list.size):
			print(file_list[i])
		
			ensemble_member = file_list[i].split('.')[1]
		
			save_file = 'member_'+ensemble_member+'_latidx_'+'{:02d}'.format(lat_idx)+'_lonidx_'+'{:02d}'.format(lon_idx)+'_years_'+str(year_start)+'-'+str(year_end)+'_threshold_'+str(threshold)+'mmday_'+region+'.npy'
			if os.path.isfile(save_dir+save_file):
				print(save_file + ' exists')
				pass#continue

			ncfile = netCDF4.Dataset(file_dir + file_list[i])
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

			#years_array = numpy.array(([t.year for t in PRECT_time_dates_NDJFM_fullseasons]), dtype=numpy.int)

			storm_count=0 # per season
			storm_length=0 # for each storm
			storm_magnitude=0.0 # for each storm

			for s in range(n_seasons):
				storm_count=0
				storm_magnitude_list = []
				storm_length_list = []
				seas_PRECT = PRECT_NDJFM_region_fullseasons[s*days_per_season:(s*days_per_season+days_per_season)]
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


			save_file = 'member_'+ensemble_member+'_latidx_'+'{:02d}'.format(lat_idx)+'_lonidx_'+'{:02d}'.format(lon_idx)+'_years_'+str(year_start)+'-'+str(year_end)+'_threshold_'+str(threshold)+'mmday_'+region+'.npy'
			
			numpy.save(save_dir + save_file, save_dict)
			
			print('Saved '+save_file)
			
			
			
			
			
			
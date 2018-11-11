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

	whiplash_count_pic = 0
	whiplash_hi_seasons = []
	whiplash_lo_seasons = []
	lo_count = 0

	s=0
	# check and see if the 3 year running mean is less than the necessary low
	# if so, consider it a low count
	# and then ask if the NEXT season is above the hi perc
	while s<seasonal_total_pic.size:
		if (seasonal_total_pic_window_mean[s]<pic_lo):
			lo_count+=1
			s+=1
		elif (seasonal_total_pic[s]>pic_hi)&(lo_count>0):
			whiplash_count_pic+=1
			whiplash_hi_seasons.append(s)
			whiplash_lo_seasons.append(s-1)
			lo_count=0
			whiplash_pic_markers[s,latlon_idx]+=1
			whiplash_pic_markers[s-1,latlon_idx]-=1
			s+=window
		else:
			lo_count=0
			s+=1

	#flatten this list
	whiplash_lo_seasons = [item for item in whiplash_lo_seasons]

	#whiplash_ratio = []
	
	#for yr_idx in range(180):
	for ens_idx in range(40):
		seasonal_total_hist_rcp = seasonal_totals_hist_rcp[:,ens_idx,latlon_idx]

		seasonal_total_separate_ensembles_window_mean = numpy.array(pandas.Series(seasonal_total_hist_rcp).rolling(window=window).mean() )

		# get top and bottom 10th percentiles in PIC
		# then see how often it transitions from 1 or 2 seasons with that
		whiplash_count = 0
		whiplash_hi_seasons = []
		whiplash_lo_seasons = []
		lo_count = 0
		whiplash_lo_seasons_ens = []
		whiplash_hi_seasons_ens = []

		s=0
		while s<180:
			if (seasonal_total_separate_ensembles_window_mean[s]<pic_lo):
				lo_count += 1
				s+=1
			elif (seasonal_total_separate_ensembles_window_mean[s]>pic_hi)&(lo_count>0):
				whiplash_count += 1
				whiplash_hi_seasons.append(s)
				whiplash_lo_seasons.append(s-1)
				lo_count = 0
				whiplash_hist_rcp_markers[s,ens_idx,latlon_idx]+=1
				whiplash_hist_rcp_markers[s-1,ens_idx,latlon_idx]-=1
				s+=window
			else:
				s+=1
				lo_count=0

			#whiplash_lo_seasons = [item for item in whiplash_lo_seasons]
			#whiplash_lo_seasons_ens.append(whiplash_lo_seasons)
			#whiplash_hi_seasons_ens.append(whiplash_hi_seasons)
	
			#whiplash_hi_seasons = []
			#whiplash_lo_seasons = []

		## times per century for PREINDUSTRIAL whiplash event
		#pic_freq = (whiplash_count_pic/(1798))*100
		## times per century for RCP8.5-like warming
		#rcp_freq = (whiplash_count/(len(ensemble_names)*(180)))*100
		## store this value
#
#		whiplash_ratio = rcp_freq/pic_freq
#		whiplash_ratios_all[yr_idx, latlon_idx] = whiplash_ratio
#		
#		whiplash_counts_pic_all[yr_idx, latlon_idx] = whiplash_count_pic
#		whiplash_counts_rcp_all[yr_idx, latlon_idx] = whiplash_count



# whiplash_ratios = numpy.zeros(( 180, 40, len(latlon_indices) ))
# for latlon_idx in range(len(latlon_indices)):
# 	whiplash_count_pic = numpy.sum( whiplash_pic_markers[:,latlon_idx]==1 )
# 	
# 	for ens_idx in range(40):
# 		whiplash_count_hist_rcp = numpy.sum( whiplash_hist_rcp_markers[:, :, latlon_idx]==1 )
# 		print(whiplash_count_hist_rcp.shape)
# for each grid box
	# for each ens member
		# calculate the ratios
	# then do same across all ensemble members as an average










# ========= save netcdf files of -1,1 whiplash values for PIC
year_list = numpy.arange(1,1798+1)
time_datetime = [datetime.datetime(i,1,15) for i in year_list]
time_nc = netCDF4.date2num(time_datetime, units='days since 0001-01-01', calendar='standard')

# save hist clim
filename = 'seasonal_whiplash_PIC_low_to_high.nc'

if os.path.exists(filename):
    print('file exists')
    os.remove(filename)

ncfile = netCDF4.Dataset(filename, 'w', format='NETCDF4')

time_dim = ncfile.createDimension('time', None)
time_var = ncfile.createVariable('time', 'f4', ('time',))
time_var[:] = time_nc
time_var.units = 'days since 0001-01-01'

lat_dim = ncfile.createDimension('lat', PRECT_nlat)
lat_var = ncfile.createVariable('lat', 'f4', ('lat',))
lat_var[:] = PRECT_lat
lat_var.units = 'degrees North'

lon_dim = ncfile.createDimension('lon', PRECT_nlon)
lon_var = ncfile.createVariable('lon', 'f4', ('lon',))
lon_var[:] = PRECT_lon
lon_var.units = 'degrees East'

data_var = ncfile.createVariable('whiplash_events', 'f4', ('time','lat','lon',))
data_var.units = 'a 20/80 whiplash event is when a gridpoint has a -1 followed by a +1'
data_var[:] = whiplash_pic_markers

ncfile.close()





# ========= save netcdf files of -1,1 whiplash values for hist+RCP8.5
year_list = numpy.arange(1921,2101)
time_datetime = [datetime.datetime(i,1,15) for i in year_list]
time_nc = netCDF4.date2num(time_datetime, units='days since 1920-01-01', calendar='standard')

# save hist clim
filename = 'seasonal_whiplash_hist_rcp_low_to_high.nc'

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


# 
# print(whiplash_counts_rcp_all)
# exit()
# print('Now saving all files!')
# # changed on 6/21 to save whiplash counts
# for yr_idx in range(year_start_list.size):
# 
# 	whiplash_label = str(year_start_list[yr_idx])+'-'+str(year_start_list[yr_idx]+1)
# 	
# 	whiplash_ratios_yearly = numpy.column_stack((whiplash_ratios_all[yr_idx,:], whiplash_counts_pic_all[yr_idx,:], whiplash_counts_rcp_all[yr_idx, :]))
# 	
# 	whiplash_ratios_all_df = pandas.DataFrame(whiplash_ratios_yearly, columns=[whiplash_label, whiplash_label, whiplash_label])
# 	whiplash_ratios_all_df.to_csv('/Users/baird/Dropbox/_analysis/attribution_2017/NEW_CALCULATIONS/calcs_and_plots/whiplash/csv_files_ALL_TIME/store_all_whiplash_ratios_dataframe_' + str(hi_perc)+str(lo_perc) + '_'+region+'_'+str(year_start_list[yr_idx]) + '-'+str(year_end_list[yr_idx])+'_windowsize_'+str(window)+'.csv')
# 


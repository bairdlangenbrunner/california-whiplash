import numpy
import scipy
import netCDF4
import matplotlib.pyplot as mp
import matplotlib.ticker
import matplotlib.colors
import scipy.stats
import pandas
import itertools

mp.rcParams.update({'mathtext.default': 'regular'})
#get_ipython().magic('matplotlib inline')

# ======================================================================

PRECT_nlat = 26
PRECT_nlon = 25

latlon_indices = list(itertools.product(range(PRECT_nlat), range(PRECT_nlon)))
region = 'whole_domain'
window=1

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

year_start_list = numpy.arange(1920,2071)
year_end_list = numpy.arange(1950,2101)

whiplash_ratios_all = numpy.zeros(( year_start_list.size, len(latlon_indices) ))
whiplash_counts_pic_all = numpy.zeros(( year_start_list.size, len(latlon_indices) ))
whiplash_counts_rcp_all = numpy.zeros(( year_start_list.size, len(latlon_indices) ))

# create list of names of members '001','002','003', ...
ensemble_members = numpy.hstack((numpy.arange(1,36), numpy.arange(101,106)))
ensemble_names = ['{:03d}'.format(i) for i in ensemble_members]

working_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/calculated_npy_files/whole_domain/'

# import monthly pic data
working_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/calculated_npy_files/whole_domain_MONTHLY_TOTALS/'
monthly_data_pic = numpy.load(working_dir+'member_005_years_0402-2200_threshold_0.1mmday_whole_domain.npy')
#print(monthly_data_pic.shape)
# (26, 25, 1799, 12)

monthly_data_pic_NDJFM = numpy.zeros((26,25,1798,5))

monthly_data_pic_NDJFM[:,:,:,0] = monthly_data_pic[:,:,:1798,10]
monthly_data_pic_NDJFM[:,:,:,1] = monthly_data_pic[:,:,:1798,11]
monthly_data_pic_NDJFM[:,:,:,2] = monthly_data_pic[:,:,1:,0]
monthly_data_pic_NDJFM[:,:,:,3] = monthly_data_pic[:,:,1:,1]
monthly_data_pic_NDJFM[:,:,:,4] = monthly_data_pic[:,:,1:,2]

# import monthly RCP data
monthly_data_hist_rcp = []
for ens_name in ensemble_names:
	monthly_data_hist_rcp.append(numpy.load(working_dir+'PRECT_monthly_data_member_' +ens_name+ '_years_1920-2100_threshold_0.1mmday_whole_domain.npy'))
monthly_data_hist_rcp = numpy.array((monthly_data_hist_rcp))
#print(monthly_data_hist_rcp.shape)
# (40, 26, 25, 181, 12)

# store just NDJFM
monthly_data_hist_rcp_NDJFM = numpy.zeros((40,26,25,180,12))
monthly_data_hist_rcp_NDJFM[:,:,:,:,0] = monthly_data_hist_rcp[:,:,:,:180,10]
monthly_data_hist_rcp_NDJFM[:,:,:,:,1] = monthly_data_hist_rcp[:,:,:,:180,11]
monthly_data_hist_rcp_NDJFM[:,:,:,:,2] = monthly_data_hist_rcp[:,:,:,1:,0]
monthly_data_hist_rcp_NDJFM[:,:,:,:,3] = monthly_data_hist_rcp[:,:,:,1:,1]
monthly_data_hist_rcp_NDJFM[:,:,:,:,4] = monthly_data_hist_rcp[:,:,:,1:,2]

# calculate 80th, 20th percentiles for PIC for each month
hi_perc = 80
lo_perc = 20
pic_hi_perc_NDJFM = numpy.percentile(monthly_data_pic_NDJFM, hi_perc, axis=2) # 26x25x5
pic_lo_perc_NDJFM = numpy.percentile(monthly_data_pic_NDJFM, lo_perc, axis=2) # 26x25x5

# whiplash events for pic runs
whiplash_lotohi_count_pic = numpy.zeros((26,25))
whiplash_hitolo_count_pic = numpy.zeros((26,25))
for lat_idx in range(26):
	for lon_idx in range(25):
		for s in range(180):
			for m in range(5-1):
				# lo to hi
				m_val = monthly_data_pic_NDJFM[lat_idx,lon_idx,s,m]
				m_plus_one_val = monthly_data_pic_NDJFM[lat_idx,lon_idx,s,m+1]
				if ( (m_val<pic_lo_perc_NDJFM[lat_idx,lon_idx,m]) and (m_plus_one_val>pic_hi_perc_NDJFM[lat_idx,lon_idx,m+1]) ):
					whiplash_lotohi_count_pic[lat_idx,lon_idx]+=1
				elif ( (m_val>pic_hi_perc_NDJFM[lat_idx,lon_idx,m]) and (m_plus_one_val<pic_lo_perc_NDJFM[lat_idx,lon_idx,m+1]) ):
					whiplash_hitolo_count_pic[lat_idx,lon_idx]+=1

# whiplash events for RCP runs
whiplash_lotohi_count_hist_rcp = numpy.zeros((180,40,26,25))
whiplash_hitolo_count_hist_rcp = numpy.zeros((180,40,26,25))
for ens_idx in range(40):
	print(ens_idx)
	for lat_idx in range(26):
		for lon_idx in range(25):
			for s in range(180):
				for m in range(5-1):
					# lo to hi
					m_val = monthly_data_hist_rcp_NDJFM[ens_idx,lat_idx,lon_idx,s,m]
					m_plus_one_val = monthly_data_hist_rcp_NDJFM[ens_idx,lat_idx,lon_idx,s,m+1]
					if ( (m_val<pic_lo_perc_NDJFM[lat_idx,lon_idx,m]) and (m_plus_one_val>pic_hi_perc_NDJFM[lat_idx,lon_idx,m+1]) ):
						whiplash_lotohi_count_hist_rcp[s,ens_idx,lat_idx,lon_idx]+=1
					elif ( (m_val>pic_hi_perc_NDJFM[lat_idx,lon_idx,m]) and (m_plus_one_val<pic_lo_perc_NDJFM[lat_idx,lon_idx,m+1]) ):
						whiplash_hitolo_count_hist_rcp[s,ens_idx,lat_idx,lon_idx]+=1

#print(whiplash_lotohi_count_hist_rcp)
whiplash_total_count_hist_rcp = whiplash_lotohi_count_hist_rcp + whiplash_hitolo_count_hist_rcp
whiplash_total_count_hist_rcp = numpy.sum(whiplash_total_count_hist_rcp, axis=0)
whiplash_total_count_hist_rcp = numpy.sum(whiplash_total_count_hist_rcp, axis=0)
mp.contourf(whiplash_total_count_hist_rcp)
mp.show()
					

exit()
for latlon_idx in range(len(latlon_indices)):

	if latlon_idx%10==0:
		print('latlon_idx', latlon_idx)

	# ======================================================================
	# Open preindustrial control info (dict_pic)
	# create season strings for PIC
	years_pic = numpy.arange(year_start_pic, year_end_pic+1, 1).astype(numpy.int)
	half_years_pic = numpy.arange(year_start_pic+0.75, year_end_pic, 1)

	season_strings_pic = [str(years_pic[i])+'-'+str(years_pic[i+1]) for i in range(years_pic.size-1)]
	member_strings_pic = ['{:03d}'.format(i) for i in range(1,36)]
	n_seasons_pic=year_end_pic-year_start_pic

	filename = 'member_005_latidx_'+'{:02d}'.format(latlon_indices[latlon_idx][0])+'_lonidx_'+'{:02d}'.format(latlon_indices[latlon_idx][1])+'_years_'+'{:04d}'.format(year_start_pic)+'-'+'{:04d}'.format(year_end_pic)+'_threshold_'+str(threshold)+'mmday_'+region+'.npy'

	dict_pic = numpy.load(working_dir + filename).item()

	seasonal_total_pic = numpy.array([dict_pic[s]['seasonal_total'] for s in season_strings_pic])

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
	# if so
	while s<seasonal_total_pic.size:
		if (seasonal_total_pic_window_mean[s]<pic_lo):
			lo_count+=1
			s+=1
		elif (seasonal_total_pic[s]>pic_hi)&(lo_count>0):
			whiplash_count_pic+=1
			whiplash_hi_seasons.append(s)
			whiplash_lo_seasons.append(s-1)
			lo_count=0
			s+=window-1
		else:
			lo_count=0
			s+=1

	#flatten this list
	whiplash_lo_seasons = [item for item in whiplash_lo_seasons]

	whiplash_ratio = []
	
	for yr_idx in range(year_start_list.size):
		if latlon_idx==0:
			print('year index is', yr_idx)
	
		year_start_whiplash=year_start_list[yr_idx]
		year_end_whiplash=year_end_list[yr_idx]

		# create season strings for whiplash
		years = numpy.arange(year_start_whiplash, year_end_whiplash+1, 1).astype(numpy.int)
		season_strings = [str(years[i])+'-'+str(years[i+1]) for i in range(years.size-1)]
		member_strings_rcp = ['{:03d}'.format(i) for i in range(1,36)]

		# ======================================================================
		# open all rcp and hist data; place in dict_list_hist_rcp
		
		dict_list_hist_rcp = []
		for i in range(len(ensemble_names)):
			ensemble_member=ensemble_names[i]
			filename = 'member_'+ensemble_member+'_latidx_'+'{:02d}'.format(latlon_indices[latlon_idx][0])+'_lonidx_'+'{:02d}'.format(latlon_indices[latlon_idx][1])+'_years_'+'{:04d}'.format(year_start)+'-'+'{:04d}'.format(year_end)+'_threshold_'+str(threshold)+'mmday_'+region+'.npy'

			dict_list_hist_rcp.append(numpy.load(working_dir + filename).item())

		# ======================================================================
		# Import seasonal accumulation data

		seasonal_total_list_rcp = []
		for i in range(len(ensemble_names)):
			seasonal_total = [dict_list_hist_rcp[i][s]['seasonal_total'] for s in season_strings]
			#seasonal_total = numpy.array([item for sublist in seasonal_total for item in sublist])  
			seasonal_total_list_rcp.append(seasonal_total)
		seasonal_total_rcp = numpy.array([item for sublist in seasonal_total_list_rcp for item in sublist])

		# ======================================================================
		# ======================================================================
		# ======================================================================
		# now count the number of events that happen in the RCP data set
		# for each season, count number of events gte the value in the pic

		# import all rcp data NOT BUNCHED TOGETHER but separately as each ensemble member
		seasonal_total_separate_ensembles = []

		# ======================================================================

		for i in range(len(ensemble_names)):
			seasonal_total = [dict_list_hist_rcp[i][s]['seasonal_total'] for s in season_strings]
			#seasonal_total = [item for sublist in seasonal_total for item in sublist]
			seasonal_total_separate_ensembles.append(numpy.array(seasonal_total))

		# now calculate windowed running mean for each ensemble member
		seasonal_total_separate_ensembles_window_mean = []
		for data in seasonal_total_separate_ensembles:
			seasonal_total_separate_ensembles_window_mean.append( numpy.array(pandas.Series(data).rolling(window=window).mean() ))

		# get top and bottom 10th percentiles in PIC
		# then see how often it transitions from 1 or 2 seasons with that

		whiplash_count = 0
		whiplash_hi_seasons = []
		whiplash_lo_seasons = []
		lo_count = 0
		whiplash_lo_seasons_ens = []
		whiplash_hi_seasons_ens = []

		for i in range(len(ensemble_names)):
			s=0
			while s<len(season_strings):
				if (seasonal_total_separate_ensembles_window_mean[i][s]<pic_lo):
					lo_count += 1
					s+=1
				elif (seasonal_total_separate_ensembles[i][s]>pic_hi)&(lo_count>0):
					whiplash_count += 1
					whiplash_hi_seasons.append(s)
					whiplash_lo_seasons.append(s-1)
					lo_count = 0
					s+=window-1
				else:
					s+=1
					lo_count=0

			whiplash_lo_seasons = [item for item in whiplash_lo_seasons]
			whiplash_lo_seasons_ens.append(whiplash_lo_seasons)
			whiplash_hi_seasons_ens.append(whiplash_hi_seasons)
	
			whiplash_hi_seasons = []
			whiplash_lo_seasons = []

		# times per century for PREINDUSTRIAL whiplash event
		pic_freq = (whiplash_count_pic/(1798-window+1))*100
		# times per century for RCP8.5-like warming
		rcp_freq = (whiplash_count/(len(ensemble_names)*(year_end_whiplash-year_start_whiplash-window+1)))*100
		# store this value

		whiplash_ratio = rcp_freq/pic_freq
		whiplash_ratios_all[yr_idx, latlon_idx] = whiplash_ratio
		
		whiplash_counts_pic_all[yr_idx, latlon_idx] = whiplash_count_pic
		whiplash_counts_rcp_all[yr_idx, latlon_idx] = whiplash_count

print('Now saving all files!')
# changed on 6/21 to save whiplash counts
for yr_idx in range(year_start_list.size):

	whiplash_label = str(year_start_list[yr_idx])+'-'+str(year_end_list[yr_idx])
	
	whiplash_ratios_yearly = numpy.column_stack((whiplash_ratios_all[yr_idx,:], whiplash_counts_pic_all[yr_idx,:], whiplash_counts_rcp_all[yr_idx, :]))
	
	whiplash_ratios_all_df = pandas.DataFrame(whiplash_ratios_yearly, columns=[whiplash_label, whiplash_label, whiplash_label])
	whiplash_ratios_all_df.to_csv('/Users/baird/Dropbox/_analysis/attribution_2017/NEW_CALCULATIONS/calcs_and_plots/whiplash/csv_files_ALL_TIME/store_all_whiplash_ratios_dataframe_' + str(hi_perc)+str(lo_perc) + '_'+region+'_'+str(year_start_list[yr_idx]) + '-'+str(year_end_list[yr_idx])+'_windowsize_'+str(window)+'.csv')
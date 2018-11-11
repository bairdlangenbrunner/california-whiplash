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

lo_perc = 10
hi_perc = 90

year_start_whiplash = 1920
year_end_whiplash = 1950

year_start_list = [1920,1950,1980,2010,2040,2070]
year_end_list = [1950,1980,2010,2040,2070,2100]

whiplash_ratios_all = numpy.zeros((len(latlon_indices)))
#whiplash_labels_all = [str(year_start_list[i])+'-'+str(year_end_list[i]) for i in range(6)]

for yr_idx in [5]:#range(6):
	print('year index is', yr_idx)
	
	whiplash_label = str(year_start_list[yr_idx])+'-'+str(year_end_list[yr_idx])
	
	year_start_whiplash=year_start_list[yr_idx]
	year_end_whiplash=year_end_list[yr_idx]

	whiplash_ratio = []
	# ======================================================================
	working_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/calculated_npy_files/'

	for latlon_idx in range(len(latlon_indices)):
		#print('==================================================')
		print(latlon_idx)

		#print('opening preindustrial control data')

		# Open preindustrial control info
		# create season strings for PIC
		years_pic = numpy.arange(year_start_pic, year_end_pic+1, 1).astype(numpy.int)
		half_years_pic = numpy.arange(year_start_pic+0.75, year_end_pic, 1)

		season_strings_pic = [str(years_pic[i])+'-'+str(years_pic[i+1]) for i in range(years_pic.size-1)]
		member_strings_pic = ['{:03d}'.format(i) for i in range(1,36)]
		n_seasons_pic=year_end_pic-year_start_pic

		filename = 'member_005_latidx_'+'{:02d}'.format(latlon_indices[latlon_idx][0])+'_lonidx_'+'{:02d}'.format(latlon_indices[latlon_idx][1])+'_years_'+'{:04d}'.format(year_start_pic)+'-'+'{:04d}'.format(year_end_pic)+'_threshold_'+str(threshold)+'mmday_'+region+'.npy'
		working_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/calculated_npy_files/whole_domain/'

		dict_pic = numpy.load(working_dir + filename).item()

		# ======================================================================
		# create list of season strings
		year_start_rcp = 1920
		year_end_rcp = 2100
		# create season strings
		years = numpy.arange(year_start_rcp, year_end_rcp+1, 1).astype(numpy.int)
		half_years_rcp = numpy.arange(year_start_rcp+0.75, year_end_rcp, 1)
		season_strings = [str(years[i])+'-'+str(years[i+1]) for i in range(years.size-1)]
		member_strings_rcp = ['{:03d}'.format(i) for i in range(1,36)]
		n_seasons_rcp=year_end_rcp-year_start_rcp

		# create list of names of members '001','002','003', ...
		ensemble_members = numpy.hstack((numpy.arange(1,36), numpy.arange(101,106)))
		ensemble_names = ['{:03d}'.format(i) for i in ensemble_members]

		# ======================================================================
		# open all rcp and hist data; place in dict_list_hist_rcp
		working_dir = '/Users/baird/Dropbox/_data_original/NCAR_LENS/daily/PRECT/calculated_npy_files/whole_domain/'
		#print('opening hist and rcp data')

		dict_list_hist_rcp = []
		for i in range(len(ensemble_names)):
			ensemble_member=ensemble_names[i]
			filename = 'member_'+ensemble_member+'_latidx_'+'{:02d}'.format(latlon_indices[latlon_idx][0])+'_lonidx_'+'{:02d}'.format(latlon_indices[latlon_idx][1])+'_years_'+'{:04d}'.format(year_start)+'-'+'{:04d}'.format(year_end)+'_threshold_'+str(threshold)+'mmday_'+region+'.npy'

			dict_list_hist_rcp.append(numpy.load(working_dir + filename).item())

		# ======================================================================
		# Import seasonal accumulation data

		#print('importing all seasonal accumulation data')
		seasonal_total_list_hist = []
		for i in range(len(ensemble_names)):
			seasonal_total = [dict_list_hist_rcp[i][s]['seasonal_total'] for s in season_strings] 
			seasonal_total_list_hist.append(seasonal_total)
		seasonal_total_hist = numpy.array([item for sublist in seasonal_total_list_hist for item in sublist])

		seasonal_total_list_rcp = []
		for i in range(len(ensemble_names)):
			seasonal_total = [dict_list_hist_rcp[i][s]['seasonal_total'] for s in season_strings]
			#seasonal_total = numpy.array([item for sublist in seasonal_total for item in sublist])  
			seasonal_total_list_rcp.append(seasonal_total)
		seasonal_total_rcp = numpy.array([item for sublist in seasonal_total_list_rcp for item in sublist])

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
		#print(whiplash_lo_seasons)

		#flatten this list
		#whiplash_lo_seasons = [item for sublist in whiplash_lo_seasons for item in sublist]
		whiplash_lo_seasons = [item for item in whiplash_lo_seasons]

		# ======================================================================
		# ======================================================================
		# ======================================================================
		# now count the number of events that happen in the RCP data set
		# for each season, count number of events gte the value in the pic

		# ======================================================================
		# ======================================================================
		# ======================================================================
		#print('whiplash calculations for hist and rcp')

		# import all rcp data NOT BUNCHED TOGETHER but separately as each ensemble member
		seasonal_total_separate_ensembles = []

		# calculate index
		year_start_index = 180 - (2100 - year_start_whiplash)
		year_end_index = 180 - (2100 - year_end_whiplash)

		# ======================================================================

		for i in range(len(ensemble_names)):
			seasonal_total = [dict_list_hist_rcp[i][s]['seasonal_total'] for s in season_strings[year_start_index:year_end_index]]
			#seasonal_total = [item for sublist in seasonal_total for item in sublist]
			seasonal_total_separate_ensembles.append(numpy.array(seasonal_total))

		# now calculate 3 yr running mean for each ensemble member
		seasonal_total_separate_ensembles_window_mean = []
		for data in seasonal_total_separate_ensembles:
			seasonal_total_separate_ensembles_window_mean.append( numpy.array(pandas.Series(data).rolling(window=window).mean() ))


		# get top and bottom 10th percentiles in PIC
		# then see how often it transitions from 1 or 2 seasons with that

		whiplash_count = 0
		whiplash_hi_seasons = []
		whiplash_lo_seasons = []
		lo_count = 0
		hi_count = 0
		whiplash_lo_seasons_ens = []
		whiplash_hi_seasons_ens = []

		for i in range(len(ensemble_names)):
			s=0
			while s<len(season_strings[year_start_index:year_end_index])-1:
				if ((seasonal_total_separate_ensembles[i][s]<pic_lo)&(seasonal_total_separate_ensembles[i][s+1]>pic_hi)):
					whiplash_count+=1
					s+=1
				elif ((seasonal_total_separate_ensembles[i][s]>pic_hi)&(seasonal_total_separate_ensembles[i][s+1]<pic_lo)):
					whiplash_count+=1
					s+=1
				else:
					s+=1
					lo_count=0

			whiplash_lo_seasons = [item for item in whiplash_lo_seasons]
			whiplash_lo_seasons_ens.append(whiplash_lo_seasons)
			whiplash_hi_seasons_ens.append(whiplash_hi_seasons)
	
			whiplash_hi_seasons = []
			whiplash_lo_seasons = []

		# ======================================================================
		# ======================================================================
		# ======================================================================
		# now count the number of events that happen in the RCP data set
		# for each season, count number of events gte the value in the pic
		#fig.savefig('./figs/all_hist_rcp_seasonal_accumulations_threshold_WHIPLASH_'+str(threshold)+'_mmday_loperc_'+'{:.0f}'.format(lo_perc)+'_hiperc_'+'{:.0f}'.format(hi_perc)+'_dryseasons_'+str(n_lo_events)+'_'+location+'.pdf', transparent=True, bbox_inches='tight')

		# times per century for PREINDUSTRIAL whiplash event
		pic_freq = (whiplash_count_pic/(1798-window+1))*100
		# times per century for RCP8.5-like warming
		rcp_freq = (whiplash_count/(len(ensemble_names)*(year_end_whiplash-year_start_whiplash-window+1)))*100
		# store this value
		whiplash_ratio.append(rcp_freq/pic_freq)

	print('========== '+region+' '+str(year_start_whiplash)+'-'+str(year_end_whiplash)+' ==========')
	for i in whiplash_ratio:  print(i)

	#print('========== whiplash counts PIC ==========')
	#print(whiplash_count_pic)
	#print('========== whiplash counts hist+rcp ==========')
	#print(whiplash_count)

	whiplash_ratios_all[:] = whiplash_ratio
	#whiplash_ratios_all[0:5] = whiplash_ratio

	whiplash_ratios_all_df = pandas.DataFrame(whiplash_ratios_all, columns=[whiplash_label])
	whiplash_ratios_all_df.to_csv('store_all_whiplash_ratios_dataframe_'+str(hi_perc)+str(lo_perc)+'_'+region+'_'+str(year_start_list[yr_idx])+'-'+str(year_end_list[yr_idx])+'_reversible.csv')
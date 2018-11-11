import numpy
import scipy
import netCDF4
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

time_subsets = numpy.array((
[402,499], \
))

################################################################################
################################################################################
################################################################################
file_dir = '/ninod/NCAR_LENS/daily/PRECT/B1850C5CN/'

file = 'b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.04020101-22001231_SUBSET.nc'
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
d_oct = 31
d_nov = 30
d_dec = 31
d_jan = 31
d_feb = 28
d_mar = 31
days_per_season = d_oct+d_nov+d_dec+d_jan+d_feb+d_mar
#print(days_per_season)

################################################################################
################################################################################
################################################################################
i = 0

ncfile = netCDF4.Dataset(file_dir + file)
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
working_dir = '/ninod/baird/cmip5/cmip5_calculations/attribution_2017/storm_counting/'
ca_latlon_regional_indices_array = numpy.load(working_dir + 'ca_latlon_regional_indices_array.npy')

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
threshold= 1.0

n_seasons = year_end-year_start

for latlon in range(ca_latlon_regional_indices_array.shape[0]):
    
    lat_idx = ca_latlon_regional_indices_array[latlon,0]
    lon_idx = ca_latlon_regional_indices_array[latlon,1]
    
    save_dict = {}

    print(file)

    ncfile = netCDF4.Dataset(file_dir + file)
    PRECT_lat = ncfile.variables['lat'][:]
    PRECT_lon = ncfile.variables['lon'][:]
    PRECT_time_var = ncfile.variables['time']

    PRECT_time_dates = netCDF4.num2date(PRECT_time_var[:], PRECT_time_var.units, PRECT_time_var.calendar)
    time_indices_ONDJFM = numpy.array([(t.month in [10,11,12,1,2,3])&(t.year in range(year_start,year_end+1)) for t in PRECT_time_dates], dtype=bool)
    PRECT_time_dates_ONDJFM = PRECT_time_dates[time_indices_ONDJFM]

    LA_lat_idx = geo_idx(LA_lat, PRECT_lat)
    LA_lon_idx = geo_idx(LA_lon, PRECT_lon)
    print(LA_lat_idx, LA_lon_idx, PRECT_lat, PRECT_lon)
    PRECT_ONDJFM_region = ncfile.variables['PRECT'][time_indices_ONDJFM, lat_idx, lon_idx]*86400.*1000

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

    storm_count=0 # per season
    storm_length=0 # for each storm
    storm_magnitude=0.0 # for each storm

    for s in range(n_seasons):
        storm_count=0
        storm_magnitude_list = []
        storm_length_list = []
        seas_PRECT = PRECT_ONDJFM_region_fullseasons[s*days_per_season:(s*days_per_season+days_per_season)]
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
        seasonal_total = numpy.sum(seas_PRECT[seas_PRECT>=threshold])
        precipitation_days = seas_PRECT[seas_PRECT>=threshold]
        save_dict[season_strings[s]] = {'storm_count' : storm_count}
        save_dict[season_strings[s]]['storm_magnitude_list'] = storm_magnitude_list
        save_dict[season_strings[s]]['storm_length_list'] = storm_length_list
        save_dict[season_strings[s]]['seasonal_total'] = seasonal_total
        save_dict[season_strings[s]]['precipitation_days'] = precipitation_days
    
    save_file = 'member_005_latidx_'+'{:02d}'.format(lat_idx)+'_lonidx_'+'{:02d}'.format(lon_idx)+'_years_'+'{:04d}'.format(year_start)+'-'+'{:04d}'.format(year_end)+'_threshold_'+str(threshold)+'mmday.npy'
    numpy.save(working_dir + save_file, save_dict)
    print('Saved '+save_file)
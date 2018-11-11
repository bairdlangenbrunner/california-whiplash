Gridded 5km GHCN-Daily Temperature and  Precipitation Dataset (nClimGrid/CLIMGRID)

The CLIMGRID dataset tracks four climate variables derived from the GHCN-D v3.12 dataset:  maximum temperature (TMAX), minimum temperature (TMIN), average temperature (TAVG) and precipitation (PRCP).  Each file provides monthly values with a lat/lon centerpoint for a nominally 5 x 5 km lat/lon grid of the Continental United States.  Data are available from Jan 1895-present.  

Each element is in a separate ASCII text file under the file name:  YYYYMM.element_name.pnt
Where, 
YYYYMM are the year and month of the data within the file
element_name are TMAX, TMIN, TAVE, or PRCP
	example record from data file: "48.9792  -99.6458   114.10"


Within each file, the data elements are organized by latitude, longitude (center of each grid square), and element value.  These are formatted as real numbers separated by one "blank" space. The latitude element has a width of 8 characters (up to 4 characters can be decimals); the longitude element has a width of 9 characters (up to 4 decimal places); the element value has a width of 8 characters (up to 2 decimal places) (in FORTRAN - F8.4,1X,F9.4,1X,F8.2).  Temperature values are reported in Celsius and precipitation in Millimeters. A description of data adjustments can be found in "Improved Historical Temperature and Precipitation Time Series for U.S. Climate Divisions" (referenced below).

These monthly files are immediately available after monthly processing on the 4th of each month and are considered "preliminary."  Files use the daily version of GHCN-D produced the day before the create date of the CLIMGRID file (for example, a CLIMGRID file created on 2014-09-04 uses the GHCN-D file created 2014-09-03). The "final" version will be released annually and will include late and reprocessed data. The most recently produced file will be available to users as soon as it is produced.

For more information on the gridding methods used in the production of this dataset, refer to Part 3 of Vose et al. 

Further documentation can be found:
1)	Transitioning from the traditional divisional dataset to the Global Historical Climatology Network-Daily gridded divisional dataset, Fenimore, et. El. 2011:  ftp://ftp.ncdc.noaa.gov/pub/data/cmb/GrDD-Transition.pdf
2)	Vose, Russel S., Scott Applequist, Mike Squires, Imke Durre, Matthew J. Menne, Claude N. Williams, Jr., Chris Fenimore, Karin Gleason, and Derek Arndt (NOAA/NESIDES National Climatic Data Center).  "Improved Historical Temperature and Precipitation Time Series for U.S. Climate Divisions." Journal of Applied Meteorology and Climatology. 2014: ftp://ftp.ncdc.noaa.gov/pub/data/cmb/vose-et-al.pdf

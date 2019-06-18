# Exercise 5
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

input_dir  = Path("C:/Users/Alex/Documents/GitHub/exercise-5-aschupp/DS_II_Data")
output_dir = Path("C:/Users/Alex/Documents/GitHub/exercise-5-aschupp/DS_II_output")

# 1. Go to http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles
#    and download the 0.25 deg. file for daily mean temperature.
#    Save the file into the data directory but don't commit it to github!!! [2P]


# 2. Read the file using xarray. Get to know your data. What's in the file?
#    Calculate monthly means for the reference periode 1981-2010 for Europe (Extent: Lon_min:-13, Lon_max: 25, Lat_min: 30, Lat_max: 72). [2P]

data = input_dir / "tg_ens_mean_0.25deg_reg_v19.0e.nc"
exc = xr.open_dataset(data)
exc_full = exc.sel(latitude = slice(30,72), longitude = slice(-13,25), time = slice("1981-01-01","2010-12-31"))
exc_slice_mean = exc_full.groupby("time.month").mean("time")


# 3. Calculate monthly anomalies from the reference period for the year 2018 (use the same extent as in #2).
#    Make a quick plot of the anomalies for the region. [2P]

exc_2018 = exc.sel(latitude = slice(30,72), longitude = slice(-13,25), time = slice("2018","2018"))
per_2018 = exc_2018.groupby("time.month").mean("time")
anom_2018 = per_2018 - exc_slice_mean 
anom_2018["tg"].plot()


# 4. Calculate the mean anomaly for the year 2018 for Europe (over all pixels of the extent from #2) 
#    Compare this overall mean anomaly to the anomaly of the pixel which contains Marburg. 
#    Is the anomaly of Marburg lower or higher than the one for Europe? [2P] 

exc_mean_full =  exc.sel(latitude = slice(30,72), longitude = slice(-13,25), time = slice("1981-01-01","2010-12-31")).mean("time")
exc_mean_2018 =  exc.sel(latitude = slice(30,72), longitude = slice(-13,25), time = slice("2018","2018")).mean("time")
anom_full_2018 = exc_mean_2018 - exc_mean_full


# location and comparison

anom_mr = anom_full_2018.sel(latitude = 50.80, longitude = 8.77, method = "nearest")

if anom_mr > anom_full_2018:
    print("The anomaly of Marburg is higher than the one for Europe.")
else:
    print("The anomaly of Marburg is lower than the one for Europe.")


# 5. Write the monthly anomalies from task 3 to a netcdf file with name "europe_anom_2018.nc" to the solution directory.
#    Write the monthly anomalies for Marburg to a csv file with name "marburg_anom_2018.csv" to the solution directory. [2P]

anom_2018.to_netcdf(output_dir / "europe_anom_2018.nc")
marburg_anom_mon = anom_2018.sel(latitude = 50.80, longitude = 8.77, method = "nearest")
mr_df = marburg_anom_mon.to_dataframe()
mr_df.to_csv(output_dir / "marburg_anom_2018.csv")
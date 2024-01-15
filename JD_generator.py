import os
from astropy.io import fits
from astropy.time import Time
import pandas as pd


ut_time_string, ut_time =[],[]
header_data,obs_data, obs_data_JD = [],[],[]

PARENT_DIR = os.getcwd()

##### Read the list of fits files in the directory ######
fits_list_path = os.path.join(PARENT_DIR, 'fits_list')
with open(fits_list_path,'r') as file1:
    filelist = [line.rstrip() for line in file1]

for f in range(len(filelist)):
        header_data.append(fits.getheader(filelist[f]))
        obs_data.append(header_data[f]["DATE-OBS"])

        ut_time_string.append(obs_data[f].replace('T',' '))                   ## Replace the T string between the date and time if present else comment out the next statement##

        ut_time.append(Time(ut_time_string[f], format='iso', scale='utc'))    ## Convert UT time string to Astropy Time object ##
        
        obs_data_JD.append(ut_time[f].jd)                                     ## Get Julian Date (JD) from the Astropy Time object ##

#### Save the JD dates in a desirable file in the same directory ####
df = pd.DataFrame({"DATE-OBS-JD": obs_data_JD})
csv_output_path = os.path.join(PARENT_DIR, 'date_obs_JD')
df.to_csv(csv_output_path, index=False)

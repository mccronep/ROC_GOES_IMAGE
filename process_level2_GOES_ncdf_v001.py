# Load GOES imagery and plot
#
#!/home/pmccrone/anaconda3/bin/python
# -*- coding: utf-8 -*-
#==============================================================
#
#==============================================================
#
#==-ROC/FRB PYTHON PROGRAM DEFINITION-==========================================
#
#/home/pmccrone/python/src/CSVMaker
#
# NAME:
# :::::::::::::::::::::::::::::::::::::::::::::::
# process_level2_GOES_ncdf_v001.py
# :::::::::::::::::::::::::::::::::::::::::::::::
#
#  PROGRAM OVERVIEW:
#       (0) The PYTHON CODE reads GOES information from a NETCDF file. 
#       (1) The information is used to decode GOES data for further post analysis.
#
#--------------------------------------------------------------------------------------------------
# PARAMETER TABLE:
#--------------------------------------------------------------------------------------------------
#
# I/O           NAME                               TYPE            FUNCTION
#--------------------------------------------------------------------------------------------------
#  I            NCDF Level I                       input           INPUT           DATA FROM GOES
#  O            Formatted table of data            output          formatted information
#_________________________________________________________________________________________________
#=================================================================================================
#
#=================================================================================================
#-
#
# Programmer: Mr. Paul McCrone     12 Dec 2025
# Modification  :  BELOW
#========================================================================================
#  Version 1.0   , Dated 2025-Dec-12
#                  Initial Build.
#========================================================================================
#  NOTE: THIS PROGRAM ASSUMES THE USE OF Python version 3.8.8+ for RHEL.
#---------------------------------------------------------------
#!/home/pmccrone/anaconda3/bin/pip install s3fs xarray metpy cartopy matplotlib h5netcdf
INSTALLERROR=0
try:
    import s3fs
    import xarray as xr
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    #
    import cartopy.feature as cfeature
    from datetime import datetime
    import numpy as np
    import metpy.plots as mpp # metpy is used for some plotting utilities
    from io import BytesIO
    import sys
    from datetime import datetime,date
    import boto3
    import os
    from botocore import UNSIGNED
    from netCDF4 import Dataset
    from botocore.config import Config
    

    INSTALLERROR=1
except:
    print("ERROR: Cound not install libraries!")
    INSTALLERROR=0

if INSTALLERROR==0:
     
    print("Errors occurred when starting the program. Make sure all modules were installed.")
    print("\n")
    sys.exit()
else:
    print("All python libraries installed without an issue.")

# Variables I will use in all functions:
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

goes_data_direct='/import/frb_archive/pmccrone/goes_data'

DADASH='-----------------------------------------------------'
dadash='-----------------------------------------------------'
dadashes='-----------------------------------------------------'
#
DAEQUALS='==--==--==--==--==--==--==--==--==--==--==--==--==--'
#
DADASHES='----------------------------------------------------'
PRTERR="--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--"
PRTOK='--OK--OK--OK--OK--OK--OK--OK--OK--OK--OK--OK--OK--OK--'


def printbn():
    #
    print('\n')
    #END OF Function

def printok():
    #
    print(PRTOK)
    #END OF Function

def printerr():
    #
    print(PRTERR)
    #END OF Function

def printds():
    #
    print('--------------------------------------------')
    #END

#
# Function to convert year month and date into ordinal.
#

def convert_ordinal_date(year, month, day):
    """
    Converts a year, month, and day into an ordinal date (1-366).
    Accounts for leap years automatically.
    """
    try:
        # Create a date object
        d = date(year, month, day)
        
        # %j is the format code for day of the year (001-366)
        ordinal_str = d.strftime('%j')
        
        return ordinal_str
    except ValueError as e:
        return f"Error: {e}"

#
#
#
"""
Function get_valid_user_date
Ask the user to input a year, then month, date, and time in minutes and seconds. 
Put this in a loop that ensures the user does not ask for a month and date that 
doesn't exist (September 31) and that the selected date is a current or past date 
before January 2016

To ensure the user provides a valid date and that the date falls within your specific 
range (after the beginning of time but before January 1, 2016), we can use a while 
True loop combined with a try/except block.

The datetime module will automatically handle the "September 31st" logic (raising 
a ValueError), and we can use a simple comparison for the date limit.

"""
def get_valid_user_date():
    # The upper limit constraint
    limit_date = datetime(2026, 1, 1)

    while True:
        try:
            print("\n--- Enter Date Details ---")
            year = int(input("Year: "))
            month = int(input("Month (1-12): "))
            day = int(input("Day: "))
            hour = int(input("Hour (0-23): "))
            minute = int(input("Minute (0-59): "))
            second = int(input("Second (0-59): "))

            # This line will fail if the date is impossible (e.g., Feb 30)
            user_dt = datetime(year, month, day, hour, minute, second)

            # Check if the date is before January 1st, 2016
            if user_dt >= limit_date:
                print("Error: The date must be BEFORE January 1, 2016. Please try again.")
                continue
            
            # If we reach here, the date is valid and within range
            return user_dt

        except ValueError as e:
            # This catches invalid numbers (like Month 13) or non-integers
            print(f"Invalid input: {e}. Please check your numbers and try again.")

#
# Function to plot GOES imagery
#

def plot_goes_data(file_path, channel_id):
    """
    Loads a GOES Level 2 NetCDF file, selects a specific channel, 
    applies the map projection, and displays the image.

    Args:
        file_path (str): The path to the NetCDF file (can be a local path or an S3 path).
        channel_id (int or str): The channel number or name (e.g., 2, 13, "CMI").
    """
    
    
    # Handle S3 paths using s3fs
    if file_path.startswith("s3://"):
        fs = s3fs.S3FileSystem(anon=True)
        file_obj = fs.open(file_path, mode='rb')
    else:
        file_obj = file_path


    print('Trying to read GOES data.')
        
    try:
        # Open the dataset with xarray
        ds = xr.open_dataset(file_obj, engine='h5netcdf')
        
        # Determine the data variable name
        if isinstance(channel_id, int):
            # For specific channel files (e.g., L1b products)
            data_var = [var for var in ds.data_vars if f"C{channel_id:02d}" in var or "Rad" in var]
            if not data_var:
                print(f"Could not find data variable for channel {channel_id}.")
                return
            data_var = data_var[0]
        else:
            # For MCMIP files (all channels in one)
            data_var = 'CMI' 

        data = ds[data_var].squeeze()
        
        # Get the projection information from the file using netCDF4 library
        # xarray can sometimes be finicky with the projection info attributes
        nc_file = Dataset(file_path if not file_path.startswith("s3://") else file_obj.full_name)
        proj_info = nc_file.variables['goes_imager_projection']

        printds()
        print("Able to ingest image file.") 

        # 2. Extract values specifically, ensuring they are floats
        # Note: some versions of xarray/netCDF4 require .item() to get the raw number

        longg_origin = proj_info.longitude_of_projection_origin
        sat_height = proj_info.perspective_point_height
        sweep = proj_info.sweep_angle_axis

        # 3. Construct the projection object
        geos = ccrs.Geostationary(
            central_longitude=longg_origin,
            satellite_height=sat_height,
            sweep_axis=sweep
        )

        print("Projection object created successfully!")

        # Plot the data using matplotlib and cartopy
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, projection=geos)
        
        # Add coastlines and borders
        ax.coastlines(resolution='50m', color='black', linewidth=1)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5) # Add country borders
        ax.add_feature(cfeature.STATES, linewidth=0.5) # Add US states
        ax.add_feature(cfeature.COASTLINE)           # Add ocean coastlines
#        ax.add_feature(cfeature.LAND, facecolor='lightgray') # Fill land color
#        ax.add_feature(cfeature.OCEAN)               # Fill ocean color
        ax.add_feature(cfeature.LAKES, alpha=0.5)    # Add lakes with transparency

        # Display the image
        #img = ax.imshow(data.values, origin='upper', extent=(data.x.min(), data.x.max(), data.y.min(), data.y.max()), 
        #                  transform=geos, interpolation='none', cmap='viridis') # Use 'viridis' or appropriate cmap
        img = ax.imshow(data.values, origin='upper', extent=(data.x.min(), data.x.max(), data.y.min(), data.y.max()), 
                          transform=geos, interpolation='none', cmap='Greys_r') # Use 'viridis' or appropriate cmap



        plt.colorbar(img, ax=ax, orientation='horizontal', pad=.05, label=f'{data.long_name} ({data.units})')
        plt.title(f'GOES Satellite Data - Channel {channel_id}')
        plt.show()

    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        if file_path.startswith("s3://"):
            file_obj.close()


def main(exitvalue):
    
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    PROGRAM_NAME="process_level2_GOES_ncdf_v001"
    #
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    exitvalue=0

    printds()
    print("Running... "+str(PROGRAM_NAME))

    # Create an anonymous s3 filesystem object
    fs = s3fs.S3FileSystem(anon=True)

    # Define the bucket and product path
    # Example: GOES-16 ABI-L2-CMIP Full Disk for a specific date and hour
    # You'll need to update the date and time to find available files
    year = '2024'
    day_of_year = '300' # Example: October 26, 2024 # Ordinal date
    hour = '18' # Example: 6 PM UTC

    print("Please enter the following details to create a timestamp:")

    # 1. Invoke the function
    selected_date = get_valid_user_date()

    # 2. Store the components into your specific variables
    syear  = selected_date.year
    smonth = selected_date.month
    sdate  = selected_date.day
    shour  = selected_date.hour
    smin   = selected_date.minute

    # 3. Verification output
    print("\n--- Variables Stored ---")
    print(f"Year (syear): {syear}")
    print(f"Month (smonth): {smonth}")
    print(f"Date (sdate): {sdate}")
    print(f"Hour (shour): {shour}")
    print(f"Minute (smin): {smin}")

    day_of_year = convert_ordinal_date(syear, smonth, sdate)    


    # Satellite Selection
    print("\nAvailable Satellites: 16, 17, 18, 19")
    sat_choice = input("Select GOES Satellite number: ").strip()
    while sat_choice not in ['16', '17', '18', '19']:
        sat_choice = input("Invalid choice. Please enter 16, 17, 18, or 19: ").strip()

    # 3. Channel Selection
    print("\n--- ABI Channel Options ---")
    print("01-02: Visible (Blue, Red)")
    print("Channel 02: Visual (Red - 0.64 µm)")
    print("03-06: Near-Infrared")
    print("Channel 07: Shortwave Window (3.9 µm - Fire/Fog detection)")
    print("Channel 13: Clean Infrared (10.3 µm - Cloud Top Temp)")
    print("07-16: Infrared")
    schan = input("Select Channel (01-16): ").strip().zfill(2) # Ensures '1' becomes '01'

    printds()
    print('Channel selected:'+str(schan))    

    # 4. Construct the S3 path
    # Note: The 'C' in the path often corresponds to the channel (e.g., C01, C02)
    path = f'noaa-goes{sat_choice}/ABI-L2-CMIPF/{syear}/{day_of_year}/{shour}/'

    # List the files in the directory
    file_list = fs.ls(path)

    # Filter for a specific file (e.g., the first one, or latest one)
    # The file name includes scan start time and other metadata
    if file_list:
        goes_file = file_list[-1] # Get the latest file
        print(f"Using file: {goes_file}")
    else:
        print("No files found for the specified time. Please update the year, day_of_year, or hour.")
        # Exit or handle error as needed
        exitvalue=99
        return exitvalue

    # 1. Setup Local Directory
    download_dir = goes_data_direct
    if not os.path.exists(download_dir):
        printerr()
        print("Problem, the goes data directory doesnt exist. Please resolve this")        
        printerr()
        exitvalue=99
        return exitvalue
        #os.makedirs(download_dir)

    # 2. Initialize S3 Client (using UNSIGNED for public access)
    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

    # 3. Define Bucket and Prefix
    # Using variables from your previous steps (sat_choice, syear, day_of_year, shour, schan)
    bucket_name = f'noaa-goes{sat_choice}'
    prefix = f'ABI-L2-CMIPF/{syear}/{day_of_year}/{shour}/'

    print(f"\nSearching S3 bucket: {bucket_name}...")

    # 4. List and Filter Files
    response = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)

    if 'Contents' in response:
        # Filter files that match the specific channel (e.g., C02, C07, C13)
        matching_files = [obj['Key'] for obj in response['Contents'] if f"C{schan}" in obj['Key']]
    
        if matching_files:
            # For this example, we take the first matching file for that hour
            s3_file_key = matching_files[0]
            filename = os.path.basename(s3_file_key)
            local_path = os.path.join(download_dir, filename)
             
            
            print(f"Downloading: {filename}...")
            s3.download_file(bucket_name, s3_file_key, local_path)
            print(f"Finished! File saved to: {local_path}")
            goes_file=local_path
        
            # --- YOUR PROCESSING CODE STARTS HERE ---
            # Now you can use variables like 'local_path' to open the file with NetCDF4 or xarray
            # example: ds = xarray.open_dataset(local_path)
       
        else:
            print(f"No files found for Channel {schan} in that hour.")
            exitvalue=99
            return exitvalue
    else:
        print("No files found at that S3 path. Check your date/satellite selection.")
        exitvalue=99
        return exitvalue

    

    # Plot data on map
    print("Starting GOES image process.")
    plot_goes_data(goes_file, schan)

    return exitvalue
#end main function

state=0

execute_value= main(state)

print("State of program = "+str(state))

#
# END 





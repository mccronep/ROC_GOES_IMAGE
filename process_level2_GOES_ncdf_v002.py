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
#/home/pmccrone/python/src/GOES
#
# NAME:
# :::::::::::::::::::::::::::::::::::::::::::::::
# process_level2_GOES_ncdf_v002.py
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
#  Version 1.0   , Dated 2025-Dec-12                   Initial Build.
#  Version 2.0   , Dated 2025-Dec-23rd                 AI version Build.

#========================================================================================
#  NOTE: THIS PROGRAM ASSUMES THE USE OF Python version 3.8.8+ for RHEL.
#---------------------------------------------------------------
#!/home/pmccrone/anaconda3/bin/pip install s3fs xarray metpy cartopy matplotlib h5netcdf

# IMPORTANT VERSIONING NOTES:
# - Version 001 (process_level2_GOES_ncdf_v001.py, a.k.a. "001") is FROZEN.
#   It is kept as a reference for original behavior and must not be modified.
# - This file is version 002 ("002"). It is allowed to be refactored and
#   extended, but should initially preserve the core scientific behavior of 001.
#
# HIGH-LEVEL DESIGN FOR 002:
# - 002 preserves the notion of downloading GOES Level 2 NetCDF data from
#   the NOAA GOES S3 bucket (e.g. noaa-goes16).
# - When 002 fetches files from S3, it ALWAYS saves them to local system disk
#   for later re-use. We do not want to redownload the same files every time.
# - The actual "processing" step in 002 operates on LOCAL NetCDF files on disk,
#   not directly on remote S3 objects.
# - In later versions (e.g. v003 and beyond), we will optionally allow the
#   script to:
#       1) Prefer already-downloaded local files (acting as a local cache).
#       2) Only go back to S3 when needed files are not present locally.
# - To prepare for that, 002 is structured so that:
#       - There is a clear "acquire phase" that talks to S3 and writes files
#         to disk (always, in 002).
#       - There is a clear "process phase" that only cares about local files.
#
# PROJECT STRUCTURE INTENT:
# - This script is the main pipeline / entry point for Level 2 processing.
# 
#

INSTALLERROR=0
try:
    import s3fs
    import xarray as xr
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from datetime import datetime, timedelta,date
    import numpy as np
    import metpy.plots as mpp # metpy is used for some plotting utilities
    from io import BytesIO
    import sys
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

WARNING_INIT_ERROR=1

BIGERRORLIST=[]

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

goes_data_direct='/import/frb_archive/pmccrone/goes_data'

DADASHES='-----------------------------------------------------'
DADASH='-----------------------------------------------------'
PRTERR="--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--"
PRTOK='--OK--OK--OK--OK--OK--OK--OK--OK--OK--OK--OK--OK--OK--'
PRTWAR="--WARNING--WARNING--WARNING--WARNING--WARNING--WARNING--WARNING--"
#
DAEQUALS='==--==--==--==--==--==--==--==--==--==--==--==--==--'


#-----------------------------------------------------------------------------
# FUNCTION: # printing funtions.
#-----------------------------------------------------------------------------
#============================================================================================
# functions
#============================================================================================
#
# The printxx functions are just simple functions to make it easy to format prints.
#

def printbn():
    """
    This function is used for making printed output readable
    """
    #
    print('\n')
    #END OF Function

def printok():
    """
    This function is used for making printed output readable
    """
    #
    print(PRTOK)
    #END OF Function

def printerr():
    """
    This function is used for making printed output readable
    """
    #
    print(PRTERR)
    #END OF Function

def printwarn():
    """
    This function is used for making printed output readable
    """
    #
    print(PRTWAR)
    #END OF Function

def printds():
    """
    This function is used for making printed output readable
    """
    #
    print('--------------------------------------------')
    #END
#
def printeq():
    """
    This function is used for making printed output readable
    """
    #
    print('============================================')
    #END

#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----
#######  Begin Function Print_Current_Time
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----
#
def Print_Current_Time(now):
    #-----
    ###import datetime
    #-----
    now = datetime.now()
    #-----
    print ()
    print( "Current date and time using str method of datetime object:")
    print( str(now))
    #-----
    print( " \n")
    print( "Current date and time using instance attributes:")
    print( "Current year: %d" % now.year)
    print( "Current month: %d" % now.month)
    print( "Current day: %d" % now.day)
    print( "Current hour: %d" % now.hour)
    print( "Current minute: %d" % now.minute)
    print( "Current second: %d" % now.second)
    print( "Current microsecond: %d" % now.microsecond)
    #-----
    print( " \n")
    print( "Current date and time using strftime:")
    #print now.strftime("%Y-%m-%d %H:%M")
    print( now.strftime("%Y-%m-%d...%H:%M"))
    #-----
    print( " \n")
    print( "Current date and time using isoformat:")
    print( now.isoformat())
    return( now.strftime("%Y-%m-%d...%H:%M"))
    #return now
    #
    #-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----
    #### END OF Print_Current_Time FUNCTION
    #-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----

#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----
#######  Begin Function small_print_current_time
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----
#
def small_print_current_time(now):
    """
    Displays current date and time.
    """
    #-----
    ###import datetime
    #-----
    now = datetime.now()
    #-----
    print ()
    print( "Current date and time using str method of datetime object:")
    print( str(now))
    #-----
    print( " \n")
    print( "Current date and time using strftime:")
    print( now.strftime("%Y-%m-%d...%H:%M"))
    print( " \n")
    print( "Current date and time using isoformat:")
    print( now.isoformat())
    return  now.strftime("%Y-%m-%d...%H:%M")
    #-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----
    #### END OF small_print_current_time FUNCTION
    #-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----


#

#.......................................................................
# Function to convert year month and date into ordinal.
#

#-----------------------------------------------------------------------------
# FUNCTION: #
#-----------------------------------------------------------------------------#
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
#-----------------------------------------------------------------------------
# FUNCTION: #
#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
# FUNCTION: compute_lat_lon_from_goes
#-----------------------------------------------------------------------------
def compute_lat_lon_from_goes(ds):
    """
    Compute latitude and longitude arrays for a GOES-R ABI fixed grid product.

    Uses:
      - x, y scan angle coordinates (radians)
      - goes_imager_projection attributes:
          * semi_major_axis
          * semi_minor_axis
          * perspective_point_height
          * longitude_of_projection_origin

    Returns
    -------
    lat : xarray.DataArray  (dims: y, x)
    lon : xarray.DataArray  (dims: y, x)
    """

    import numpy as np
    import xarray as xr

    if "goes_imager_projection" not in ds.variables:
        raise ValueError("Dataset does not contain 'goes_imager_projection' variable.")

    proj = ds["goes_imager_projection"]

    # Projection parameters
    r_eq = proj.semi_major_axis
    r_pol = proj.semi_minor_axis
    h_sat = proj.perspective_point_height + r_eq
    lon0_rad = np.deg2rad(proj.longitude_of_projection_origin)

    # 1D fixed grid scan angle coordinates (radians)
    if "x" not in ds.coords or "y" not in ds.coords:
        raise ValueError("Dataset does not contain 'x' and 'y' scan angle coordinates.")

    x_1d = ds["x"].values
    y_1d = ds["y"].values

    # 2D grids
    X, Y = np.meshgrid(x_1d, y_1d)

    # GOES-R fixed grid to latitude/longitude (NOAA standard)
    a = (np.sin(X) ** 2 +
         np.cos(X) ** 2 *
         (np.cos(Y) ** 2 + (r_eq ** 2 / r_pol ** 2) * np.sin(Y) ** 2))
    b = -2.0 * h_sat * np.cos(X) * np.cos(Y)
    c = h_sat ** 2 - r_eq ** 2

    # Discriminant can go slightly negative from roundoff; clip at 0 to
    # avoid sqrt of negative values.
    disc = b ** 2 - 4.0 * a * c
    disc = np.maximum(disc, 0.0)

    r_s = (-b - np.sqrt(disc)) / (2.0 * a)

    s_x = r_s * np.cos(X) * np.cos(Y)
    s_y = -r_s * np.sin(X)
    s_z = r_s * np.cos(X) * np.sin(Y)

    lon = lon0_rad + np.arctan2(s_y, s_x)
    lat = np.arctan((r_eq ** 2 / r_pol ** 2) *
                    (s_z / np.sqrt((h_sat - s_x) ** 2 + s_y ** 2)))

    lon_deg = np.rad2deg(lon)
    lat_deg = np.rad2deg(lat)

    # Wrap into DataArrays aligned to y/x dims
    lat_da = xr.DataArray(
        lat_deg,
        dims=("y", "x"),
        coords={"y": ds["y"], "x": ds["x"]},
        name="lat",
    )
    lon_da = xr.DataArray(
        lon_deg,
        dims=("y", "x"),
        coords={"y": ds["y"], "x": ds["x"]},
        name="lon",
    )

    return lat_da, lon_da

#
#-----------------------------------------------------------------------------
# FUNCTION: # Function to plot GOES imagery
#-----------------------------------------------------------------------------
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
                          transform=geos, interpolation='none', cmap='gray_r') # Use 'viridis' or appropriate cmap



        plt.colorbar(img, ax=ax, orientation='horizontal', pad=.05, label=f'{data.long_name} ({data.units})')
        plt.title(f'GOES Satellite Data - Channel {channel_id}')
        plt.show()

    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        if file_path.startswith("s3://"):
            file_obj.close()


#-----------------------------------------------------------------------------
# FUNCTION:plot_goes_conus_region
#-----------------------------------------------------------------------------
def plot_goes_conus_region(file_path, channel_id):
    """
    Plot GOES data for an approximate CONUS region:
        Latitude:  24N to 47N
        Longitude: -125W to -66W

    This version:
      - Computes lat/lon from GOES projection if not present.
      - Slices by index (no where(drop=True)).
      - Flips y so that north is up when using origin="lower".
    """

    import xarray as xr
    import numpy as np
    import matplotlib.pyplot as plt

    ds = xr.open_dataset(file_path)

    if "CMI" in ds.data_vars:
        data_var_name = "CMI"
    else:
        ds.close()
        raise ValueError("Could not find expected data variable (e.g., 'CMI') in dataset.")

    data = ds[data_var_name]

    if "lat" in ds.variables and "lon" in ds.variables:
        lat = ds["lat"]
        lon = ds["lon"]
    else:
        lat, lon = compute_lat_lon_from_goes(ds)

    lat_min = 24.0
    lat_max = 47.0
    lon_min = -125.0
    lon_max = -66.0

    region_mask = (
        (lat >= lat_min) & (lat <= lat_max) &
        (lon >= lon_min) & (lon <= lon_max)
    )

    mask_vals = region_mask.values
    valid_y, valid_x = np.where(mask_vals)

    if valid_y.size == 0 or valid_x.size == 0:
        ds.close()
        print("No data points found inside CONUS bounds for this file.")
        return

    y_min = int(valid_y.min())
    y_max = int(valid_y.max())
    x_min = int(valid_x.min())
    x_max = int(valid_x.max())

    if "y" in data.dims and "x" in data.dims:
        data_conus = data.isel(
            y=slice(y_min, y_max + 1),
            x=slice(x_min, x_max + 1)
        )
    else:
        ds.close()
        raise ValueError("Expected 'y' and 'x' in data dimensions, got " + str(data.dims))

    # Flip north-south so that north is up in the image
    data_conus_plot = data_conus.isel(y=slice(None, None, -1))

    plt.figure(figsize=(10, 6))
    plt.imshow(
        data_conus_plot,
        origin="lower",
        cmap="gray_r"
    )
    plt.title("GOES Channel " + str(channel_id) + " - CONUS Region")
    plt.colorbar(label=data_var_name)
    plt.tight_layout()
    plt.show()

    ds.close()


#-----------------------------------------------------------------------------
# FUNCTION: plot_tropical_atlantic_region
#-----------------------------------------------------------------------------
def plot_tropical_atlantic_region(file_path, channel_id):
    """
    Plot GOES data for the Tropical Atlantic region:
        Latitude:   0N  to 30N
        Longitude: -90W to 20E
    """

    import xarray as xr
    import numpy as np
    import matplotlib.pyplot as plt

    ds = xr.open_dataset(file_path)

    # Main data variable (adjust if needed)
    if "CMI" in ds.data_vars:
        data_var_name = "CMI"
    else:
        ds.close()
        raise ValueError("Could not find expected data variable (e.g., 'CMI') in data set for Tropical Atlantic plotting.")

    data = ds[data_var_name]

    # Lat/lon, either from file or computed
    if "lat" in ds.variables and "lon" in ds.variables:
        lat = ds["lat"]
        lon = ds["lon"]
    else:
        lat, lon = compute_lat_lon_from_goes(ds)

    # Tropical Atlantic bounds
    lat_min = 0.0
    lat_max = 30.0
    lon_min = -90.0
    lon_max = 20.0

    region_mask = (
        (lat >= lat_min) & (lat <= lat_max) &
        (lon >= lon_min) & (lon <= lon_max)
    )

    mask_vals = region_mask.values
    valid_y, valid_x = np.where(mask_vals)

    if valid_y.size == 0 or valid_x.size == 0:
        ds.close()
        print("No data points found inside Tropical Atlantic bounds.")
        return

    y_min = int(valid_y.min())
    y_max = int(valid_y.max())
    x_min = int(valid_x.min())
    x_max = int(valid_x.max())

    if "y" in data.dims and "x" in data.dims:
        data_tatl = data.isel(y=slice(y_min, y_max + 1),
                              x=slice(x_min, x_max + 1))
    else:
        ds.close()
        raise ValueError("Expected 'y' and 'x' in data dimensions, got " + str(data.dims))




    plt.figure(figsize=(10, 6))
 
    # after you have data_tatl from isel(...)
    data_tatl_plot = data_tatl.isel(y=slice(None, None, -1))

    plt.imshow(
        data_tatl_plot,
        origin="lower",
        cmap="gray_r"
    )

    plt.title("GOES Channel " + str(channel_id) + " - Tropical Atlantic Region")
    plt.colorbar(label=data_var_name)
    plt.tight_layout()
    plt.show()

    ds.close()

#-----------------------------------------------------------------------------
# FUNCTION: get_region_selection
#-----------------------------------------------------------------------------
def get_region_selection():
    """
    Ask the user (via keyboard) which region they want to view:
        1) Full Disk
        2) CONUS
        3) Tropical Atlantic
        4) Specific State

    If "Specific State" is selected, a second menu is presented where the
    user can choose a particular U.S. state.

    Returns
    -------
    region_type : str
        One of: "FULL_DISK", "CONUS", "TROPICAL_ATLANTIC", "STATE"
    region_value : str or None
        For:
          - "STATE": the name of the state selected (e.g. "Colorado")
          - others: None
    """

    while True:
        print("")
        print("==============================================")
        print(" Select Region To View")
        print("==============================================")
        print("  1) Full Disk")
        print("  2) CONUS")
        print("  3) Tropical Atlantic")
        print("  4) Specific State")
        print("")

        choice = input("Enter choice (1-4): ").strip()

        if choice == "1":
            return "FULL_DISK", None
        elif choice == "2":
            return "CONUS", None
        elif choice == "3":
            return "TROPICAL_ATLANTIC", None
        elif choice == "4":
            # Go to separate state menu
            region_type, region_value = get_state_selection()
            return region_type, region_value
        else:
            print("")
            print(" Invalid selection. Please enter a number from 1 to 4.")
            print("")


#-----------------------------------------------------------------------------
# FUNCTION: get_state_selection
#-----------------------------------------------------------------------------
def get_state_selection():
    """
    Present a menu of U.S. states and ask the user to select one.

    This function is called when the user chooses "Specific State" from
    get_region_selection(). The list of states can be expanded or modified
    as needed.

    Returns
    -------
    region_type : str
        Always "STATE" for this function.
    state_name : str
        The name of the selected state (e.g. "Colorado").
    """

    # You can expand or customize this list as needed.
    states = [
        "Alabama", "Alaska", "Arizona", "Arkansas", "California",
        "Colorado", "Connecticut", "Delaware", "Florida", "Georgia",
        "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa",
        "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland",
        "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri",
        "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey",
        "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio",
        "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina",
        "South Dakota", "Tennessee", "Texas", "Utah", "Vermont",
        "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming"
    ]

    while True:
        print("")
        print("==============================================")
        print(" Select Specific State")
        print("==============================================")

        # Show states in a simple numbered list
        for idx, state in enumerate(states, start=1):
            print("  " + str(idx).rjust(2) + ") " + state)
        print("")

        choice = input("Enter state number (1-" + str(len(states)) + "): ").strip()

        # Validate numeric choice and range
        if not choice.isdigit():
            print("")
            print(" Invalid input. Please enter a number.")
            print("")
            continue

        idx = int(choice)
        if idx < 1 or idx > len(states):
            print("")
            print(" Invalid selection. Please enter a number between 1 and " + str(len(states)) + ".")
            print("")
            continue

        # Valid selection
        state_name = states[idx - 1]
        return "STATE", state_name

#-----------------------------------------------------------------------------
# FUNCTION: latlonstate
#-----------------------------------------------------------------------------
def latlonstate(state_name):
    """
    Given the name of a U.S. state, return an approximate latitude and
    longitude corresponding to the center of that state's land area.

    This is intended to work with the state names returned by the
    get_state_selection() function.

    Parameters
    ----------
    state_name : str
        Full state name, e.g. "Colorado", "Florida", "New York".

    Returns
    -------
    (lat, lon) : tuple of float
        Approximate center latitude and longitude in degrees.

    Notes
    -----
    - These coordinates are approximate geographic centers, not
      population-weighted centers.
    - If the state is not recognized, this function raises a ValueError.
    """

    # Dictionary of approximate geographic centers for U.S. states.
    # Latitudes and longitudes are in decimal degrees.
    state_centers = {
        "Alabama":        (32.806671, -86.791130),
        "Alaska":         (61.370716, -152.404419),
        "Arizona":        (33.729759, -111.431221),
        "Arkansas":       (34.969704, -92.373123),
        "California":     (36.116203, -119.681564),
        "Colorado":       (39.059811, -105.311104),
        "Connecticut":    (41.597782, -72.755371),
        "Delaware":       (39.318523, -75.507141),
        "Florida":        (27.766279, -81.686783),
        "Georgia":        (33.040619, -83.643074),
        "Hawaii":         (21.094318, -157.498337),
        "Idaho":          (44.240459, -114.478828),
        "Illinois":       (40.349457, -88.986137),
        "Indiana":        (39.849426, -86.258278),
        "Iowa":           (42.011539, -93.210526),
        "Kansas":         (38.526600, -96.726486),
        "Kentucky":       (37.668140, -84.670067),
        "Louisiana":      (31.169546, -91.867805),
        "Maine":          (44.693947, -69.381927),
        "Maryland":       (39.063946, -76.802101),
        "Massachusetts":  (42.230171, -71.530106),
        "Michigan":       (43.326618, -84.536095),
        "Minnesota":      (45.694454, -93.900192),
        "Mississippi":    (32.741646, -89.678696),
        "Missouri":       (38.456085, -92.288368),
        "Montana":        (46.921925, -110.454353),
        "Nebraska":       (41.125370, -98.268082),
        "Nevada":         (38.313515, -117.055374),
        "New Hampshire":  (43.452492, -71.563896),
        "New Jersey":     (40.298904, -74.521011),
        "New Mexico":     (34.840515, -106.248482),
        "New York":       (42.165726, -74.948051),
        "North Carolina": (35.630066, -79.806419),
        "North Dakota":   (47.528912, -99.784012),
        "Ohio":           (40.388783, -82.764915),
        "Oklahoma":       (35.565342, -96.928917),
        "Oregon":         (44.572021, -122.070938),
        "Pennsylvania":   (40.590752, -77.209755),
        "Rhode Island":   (41.680893, -71.511780),
        "South Carolina": (33.856892, -80.945007),
        "South Dakota":   (44.299782, -99.438828),
        "Tennessee":      (35.747845, -86.692345),
        "Texas":          (31.054487, -97.563461),
        "Utah":           (40.150032, -111.862434),
        "Vermont":        (44.045876, -72.710686),
        "Virginia":       (37.769337, -78.169968),
        "Washington":     (47.400902, -121.490494),
        "West Virginia":  (38.491226, -80.954453),
        "Wisconsin":      (44.268543, -89.616508),
        "Wyoming":        (42.755966, -107.302490),
    }

    if state_name not in state_centers:
        raise ValueError("Unknown state name for latlonstate: " + str(state_name))

    return state_centers[state_name]

#-----------------------------------------------------------------------------
# FUNCTION: plot_goes_state_scale
#-----------------------------------------------------------------------------
def plot_goes_state_scale(file_path, channel_id, center_lat, center_lon):
    """
    Plot GOES Level 2 data for a circular region (radius ~620 km) centered
    on a given latitude/longitude (e.g. state center).
    """

    import xarray as xr
    import numpy as np
    import matplotlib.pyplot as plt

    ds = xr.open_dataset(file_path)

    # Main data variable (adjust name if needed)
    if "CMI" in ds.data_vars:
        data_var_name = "CMI"
    else:
        ds.close()
        raise ValueError("Could not find expected data variable (e.g., 'CMI') in data set.")

    data = ds[data_var_name]

    # Lat/lon, either from file or computed
    if "lat" in ds.variables and "lon" in ds.variables:
        lat = ds["lat"]
        lon = ds["lon"]
    else:
        lat, lon = compute_lat_lon_from_goes(ds)

    # Great-circle distance from center
    R_earth_km = 6371.0
    radius_km = 620.0

    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    lat0_rad = np.deg2rad(center_lat)
    lon0_rad = np.deg2rad(center_lon)

    dlat = lat_rad - lat0_rad
    dlon = lon_rad - lon0_rad

    a_gc = (np.sin(dlat / 2.0) ** 2 +
            np.cos(lat0_rad) * np.cos(lat_rad) * np.sin(dlon / 2.0) ** 2)
    c_gc = 2.0 * np.arctan2(np.sqrt(a_gc), np.sqrt(1.0 - a_gc))
    distance_km = R_earth_km * c_gc

    region_mask = distance_km <= radius_km

    # Find bounding box indices of points inside the radius
    mask_vals = region_mask.values
    valid_y, valid_x = np.where(mask_vals)

    if valid_y.size == 0 or valid_x.size == 0:
        ds.close()
        print("No data points found within 620 km of the chosen center.")
        return

    y_min = int(valid_y.min())
    y_max = int(valid_y.max())
    x_min = int(valid_x.min())
    x_max = int(valid_x.max())

    # Subset data via slices (no drop=True, avoids 0-length dims)
    if "y" in data.dims and "x" in data.dims:
        data_state = data.isel(y=slice(y_min, y_max + 1),
                               x=slice(x_min, x_max + 1))
    else:
        ds.close()
        raise ValueError("Expected 'y' and 'x' in data dimensions, got " + str(data.dims))

    plt.figure(figsize=(8, 6))
    data_state_plot = data_state.isel(y=slice(None, None, -1))

    plt.imshow(
        data_state_plot,
        origin="lower",
        cmap="gray_r"
     )
    plt.title(
        "GOES Channel " + str(channel_id) +
        " - ~620 km Radius around (" +
        str(round(center_lat, 2)) + ", " + str(round(center_lon, 2)) + ")"
    )
    plt.colorbar(label=data_var_name)
    plt.tight_layout()
    plt.show()

    ds.close()

#-----------------------------------------------------------------------------
# FUNCTION: main
#-----------------------------------------------------------------------------

def main(exitvalue):
    
    now=0
    Print_Current_Time(now)

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    PROGRAM_NAME="process_level2_GOES_ncdf_v002.py"
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
    # Select the image region to display

    region_type, region_value = get_region_selection()

    print("Region type selected:", region_type)
    print("Region value:", region_value)

    if region_type == "FULL_DISK":
        # set full disk extent / projection
        plot_goes_data(goes_file, schan)
        exitvalue=1
    elif region_type == "CONUS":
        # set CONUS extent
        plot_goes_conus_region(goes_file, schan)
        exitvalue=1
    elif region_type == "TROPICAL_ATLANTIC":
        # set Tropical Atlantic extent
        plot_tropical_atlantic_region(goes_file, schan)
        exitvalue=1
    elif region_type == "STATE":
        # region_value is the state name, e.g. "Colorado"
        # use that to look up a bounding box, shapefile, etc.
        state_name=get_state_selection()
        center_lat, center_lon=latlonstate(state_name)
        plot_goes_state_scale(goes_file, schan, center_lat, center_lon)
        exitvalue=1
    else:
        # Error
        printerr()
        print("Invalid region selection. Program ends.")
        printerr()
        exitvalue=99
        sys.exit()

    return exitvalue
#end main function

state=0

execute_value= main(state)

print("State of program = "+str(state))
now=0
small_print_current_time(now)
#
# END 





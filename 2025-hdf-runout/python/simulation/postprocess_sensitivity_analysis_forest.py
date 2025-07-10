######################################################################
# Copyright (C) 2025 BFH
#
# Script for post processing the runout simulations for the 
# sensitivity analyis for historical landslides in the forest.
#
# Author: Christoph Schaller, BFH-HAFL, June 2025
######################################################################

#
# Imports
#

from osgeo import gdal
from osgeo.gdalconst import *
from osgeo.gdalnumeric import *
from osgeo import gdal_array

gdal.UseExceptions()

import numpy as np
import pandas as pd
import glob
import os
import sys 

# Adding the path to the utilities.py to the path and import
utility_path = "E:/GIS_Projekte/Paper_4/code/utilities"
# utility_path = os.environ["GHK_UTILITY_PATH"]
sys.path.append(utility_path)
import utilities

# Default entry point
if __name__ == "__main__":

    #
    # Configs
    #

    input_base_path = "E:/GIS_Projekte/Paper_4/data/sensitivity_forest"

    #
    # Extract results from simulations
    #

    # Iterate the folders for all landslides
    slide_folders = glob.glob(os.path.join(input_base_path,"slide_*"))
    for slide_folder_path in slide_folders:
        # Skip files
        if not os.path.isdir(slide_folder_path):
            continue

        # Skip if extracted stats already exist
        sim_stat_path = os.path.join(slide_folder_path,"sim_stats.csv")
        if os.path.exists(sim_stat_path):
            continue

        slide_folder_name = os.path.basename(slide_folder_path)

        # Skip log folder (if present)
        if slide_folder_name == "log":
            continue

        # Extract slide ID from folder name
        slide_id = int(slide_folder_name.split("_")[-1])

        # Process all simulations
        sim_records = []
        sim_folders = glob.glob(os.path.join(slide_folder_path,"sim_*"))
        print(len(sim_folders))
        i = 0
        for sim_folder_path in sim_folders:
            # Skip files
            if not os.path.isdir(sim_folder_path):
                continue

            sim_folder_name = os.path.basename(sim_folder_path)

            # Skip log folder (if present)
            if sim_folder_name=="log":
                continue

            # Extract combination index from folder name
            sim_id = int(sim_folder_name.split("_")[-1])

            # Read the needed output rasters
            runout_path = os.path.join(sim_folder_path,"reach_probability.tif")
            velocity_path = os.path.join(sim_folder_path,"V_mean.tif")
            fh_path = os.path.join(sim_folder_path,"Fh_mean.tif")

            runout_values = gdal_array.LoadFile(runout_path)
            v_values = gdal_array.LoadFile(velocity_path)
            fh_values = gdal_array.LoadFile(fh_path)



            # Determine the global mean/max values for velocity, and flow height 
            if np.sum(v_values)>0:
                v_max = np.max(v_values[v_values>0])
                v_mean = np.mean(v_values[v_values>0])
            else:
                v_max = 0
                v_mean = 0

            if np.sum(fh_values[fh_values>0])>0:
                fh_max = np.max(fh_values[fh_values>0])
                fh_mean = np.mean(fh_values[fh_values>0])
            else: 
                fh_max = 0
                fh_mean = 0
    
            # Determine dimensions of the runout
            x_min = 100000
            x_max = -1
            y_min = 100000
            y_max = -1
            max_spread = 0
            n_runout = np.sum(runout_values>0)

            for x in range(np.shape(runout_values)[1]):
                if np.sum(runout_values[:,x])>0:
                    x_min = min(x_min,x)
                    x_max = max(x_max,x)
                    
                    y_min_curr = 10000
                    y_max_curr = -1
                    for y in range(np.shape(runout_values)[0]):
                        if runout_values[y,x]>0:
                            y_min_curr = min(y_min_curr,y)
                            y_max_curr = max(y_max_curr,y)

                    y_spread_curr = y_max_curr-y_min_curr
                    max_spread = max(max_spread,y_spread_curr)

                    y_min = min(y_min,y_min_curr)
                    y_max = max(y_max,y_max_curr)

            rec = {"slide_id":slide_id, "step":sim_id, 
                "x_min":x_min,"x_max":x_max,"y_min":y_min,"y_max":y_max,"w":x_max-x_min,"h":y_max-y_min,"max_spread":max_spread,"n_runout":n_runout,"v_max":v_max,"v_mean":v_mean,"fh_max":fh_max,"fh_mean":fh_mean}
            sim_records.append(rec)

        # Skip saving results if no results were collected
        if len(sim_records)==0:
            continue        

        # Create a CSV with the results
        sim_stat_path = os.path.join(slide_folder_path,"sim_stats.csv")
        if os.path.exists(sim_stat_path):
            os.remove(sim_stat_path)

        sim_stats_df = pd.DataFrame(sim_records)
        sim_stats_df.to_csv(sim_stat_path,sep=",",index=False)


    #
    # Collect and merge CSV files with the results per landslide
    #

    # Get list of folders
    subfolders = [ f.path for f in os.scandir(input_base_path) if f.is_dir() ]

    # Collect results from individual simulations
    records = []
    for slide_folder_path in subfolders:
        print(slide_folder_path)
        # Skip log folder (if present)
        if slide_folder_path.endswith("log"):
            continue

        sim_stats_in_path = os.path.join(slide_folder_path,"sim_stats.csv")

        df_sim = pd.read_csv(sim_stats_in_path,sep=',')
        if len(df_sim)==0:
            continue
        records.append(df_sim)

    # Merge the individual results and save as CSV
    df_merged = pd.concat(records)
    sim_stats_out_path = os.path.join(input_base_path,"sim_stats_merged.csv")
    df_merged = df_merged.sort_values(by="step")
    df_merged.to_csv(sim_stats_out_path,sep=",",index=False)

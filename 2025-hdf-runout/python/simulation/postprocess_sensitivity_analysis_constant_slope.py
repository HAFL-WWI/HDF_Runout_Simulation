######################################################################
# Copyright (C) 2025 BFH
#
# Script for post processing the runout simulations for the 
# sensitivity analyis on a constant gradient slope.
#
# Author: Christoph Schaller, BFH-HAFL, June 2025
######################################################################

#
# Imports
#

import os
import pandas as pd

# Default entry point
if __name__ == "__main__":

    #
    # Configs
    #

    input_path = "E:/GIS_Projekte/Paper_4/data/sensitivity_constant_slope/"
    
    # Get list of folders
    subfolders = [ f.path for f in os.scandir(input_path) if f.is_dir() ]

    # Collect results from individual simulations
    records = []
    for sim_path in subfolders:
        print(sim_path)

        # Skip log folder (if present)
        if sim_path.endswith("log"):
            continue

        sim_stats_in_path = os.path.join(sim_path,"sim_stats.csv")

        df_sim = pd.read_csv(sim_stats_in_path,sep=',')
        if len(df_sim)==0:
            continue
        records.append(df_sim)

    # Merge the individual results and save as CSV
    df_merged = pd.concat(records)
    sim_stats_out_path = os.path.join(input_path,"sim_stats_merged.csv")
    df_merged = df_merged.sort_values(by="step")
    df_merged.to_csv(sim_stats_out_path,sep=",",index=False)

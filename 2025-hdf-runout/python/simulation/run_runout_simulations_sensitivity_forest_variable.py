######################################################################
# Copyright (C) 2025 BFH
#
# Script for running the runout simulations for the sensitivity 
# analyis for historic landslides in forests using variable friction.
#
# Author: Christoph Schaller, BFH-HAFL, June 2025
######################################################################

#
# Imports
#

from osgeo import gdal
from osgeo.gdalconst import *
from osgeo.gdalnumeric import *
gdal.UseExceptions()

import numpy as np
import pandas as pd
import geopandas as gpd
import os
import sys 

from datetime import datetime
import time

# Adding the path to the utilities.py to the path and import
utility_path = "E:/GIS_Projekte/GHK/dataprocessing_be/utilities"
# utility_path = os.environ["GHK_UTILITY_PATH"]
sys.path.append(utility_path)
import utilities

import glob

from datetime import datetime, date, time, timedelta
import time

from multiprocessing import Process, Pool, Queue, JoinableQueue, Manager, current_process, freeze_support
from queue import Empty


#
# Parallel processing functions
#

# Worker function that will handle the processing of the individual records
def worker(q, r, work_function, cfg):
    
    while True:
        #Consume work as long as there are items in the queue
        try:
            job_record = q.get()
            if job_record == None:
                q.task_done()
                print("Queue End")
                break

            res = work_function(job_record, cfg)
            if res:
                # print("queuing result")
                r.put([res])

            q.task_done()
        except Empty:
            # print("Queue empty")
            break
    #No more work available
    print("Exit:",current_process())
    return

# Process records in parallel fashion 
def process_records(records, process_function, cfg, num_processes = 1):

    with Manager() as manager:
        # Create queues
        job_queue = JoinableQueue()
        result_queue = manager.Queue()

        if len(records) < num_processes:
            num_processes = len(records)

        #Insert records into queue
        for r in records:
            job_queue.put(r)


        #Create and start worker processes
        processes = [] 
        for i in range(num_processes):
            job_queue.put(None)
            proc = Process(target=worker, args=(job_queue,result_queue,process_function,cfg,))
            processes.append(proc)
            print("Start: ",proc)
            proc.start()

        job_queue.join()

        print("Processing finished")
        print("Queue empty:",result_queue.empty())

        res = []
        if not result_queue.empty():
            for r in records:
                item = result_queue.get()
                res.extend(item)
        return res

# Process records in a linear fashion -> good for debugging the script
def process_records_linear(records, process_function, cfg,  num_processes = 1):
    # Create queues
    job_queue = JoinableQueue()
    result_queue = Queue()

    #Insert records into queue
    for r in records:
        job_queue.put(r)

    job_queue.put(None)

    print("Start:")
    worker(job_queue,result_queue,process_function,cfg)

    print("Processing finished")
    print("Queue empty:",result_queue.empty())
    res = []
    if not result_queue.empty():
        for r in records:
            item = result_queue.get()
            res.extend(item)

    return res

#
# Business logic
#

# Setup the records describing the work units for parallel processing
def prepare_tiles (cfg, num_processes = 1):
    output_path = cfg["output_path"]
    
    # Read the parameter values
    param_in_path = os.path.join(output_path,"sensitivity_param_values_with_forest_variable.txt")
    param_values = np.loadtxt(param_in_path) 

    # For each folder/historic landslide
    folders = glob.glob(os.path.join(output_path,"slide_*"))
    print(len(folders))
    tile_records = []
    for slide_folder_path in folders:
        slide_folder_name = os.path.basename(slide_folder_path)
        slide_id = int(slide_folder_name.replace("slide_",""))

        # Read only start area polygon
        slide_gpkg = os.path.join(slide_folder_path,"slide.gpkg")
        slide_layer = "slide"
        slide_df = gpd.GeoDataFrame.from_file(slide_gpkg,sql="SELECT * FROM {0} WHERE slide_id={1} AND area_type=1".format(slide_layer,slide_id),engine="pyogrio",)

        slide_df["in_forest"] = slide_df["in_forest"].apply(lambda x: x if not np.isnan(x) else 0)

        # Create a record for each parameter combination
        for i in range(len(param_values)):
            [thickness, density, height, type, mu, deposit] = param_values[i]

            tile_record = {
                "slide_id": slide_id,
                "slide_name": slide_folder_name,
                "slide_thickness": thickness,
                "in_forest": 0,
                "output_path": output_path,
                "slide_path": slide_folder_path,
                "slide_type": type,
                "max_flow_height": height,
                "moving_mass_density": density,
                "mu": mu,
                "deposit": deposit,
                "step":i,
                "output_name":"sim_{0}".format(i),

                "dem_5m_path": os.path.join(slide_folder_path,"dem_5_filled.tif"),
                "slide_mask_path": os.path.join(slide_folder_path,"slide_mask.tif"),
                "start_mask_path": os.path.join(slide_folder_path,"start_mask.tif"),

                "tree_deposit_path": os.path.join(slide_folder_path,"tree_deposit.tif"),
                "forest_friction_path": os.path.join(slide_folder_path,"forest_dn.tif"),
            }
            tile_records.append(tile_record)

    return tile_records

# Worker function to setup and run simulations
def process_simulation(sim_record, cfg):

    # Get parameter values for sensitivity from record
    moving_mass_density = int(round(sim_record["moving_mass_density"]))
    slide_type = int(sim_record["slide_type"])
    slide_thickness = sim_record["slide_thickness"]
    max_flow_height = sim_record["max_flow_height"]
    friction = sim_record["mu"]
    deposit = sim_record["deposit"]

    # Get parameters for simulation from config
    num_simulations = cfg["num_simulations"]
    pressure_model = cfg["pressure_model"]
    slideforce_path = cfg["slideforce_path"]

    slide_name = sim_record["slide_name"]    
    slide_path = sim_record["slide_path"]
    in_forest = sim_record["in_forest"]
    dem_path = sim_record["dem_5m_path"]
    slide_mask_path = sim_record["slide_mask_path"]
    start_mask_path = sim_record["start_mask_path"]
    tree_deposit_path = sim_record["tree_deposit_path"]
    forest_friction_path = sim_record["forest_friction_path"]

    # Create output directory
    output_path = os.path.join(slide_path,sim_record["output_name"])
    utilities.ensure_dir(output_path)

    # Test if simulation was already run and skip if so
    reach_probability_path = os.path.join(output_path,"reach_probability.tif")
    if os.path.exists(reach_probability_path):
        return slide_path

    # Available options in SlideForce:
    #   -h [ --help ]              produce help message
    #   -o [ --outputdir ] arg     output directory where generated files will be put
    #   -m [ --dem ] arg           input digital elevation model
    #   -a [ --start ] arg         input raster for start cells
    #   -q [ --depthraster ] arg   input raster for depth of start cells
    #   -r [ --friction ] arg      input raster for friction modification of cells
    #   -f [ --forest ] arg        input raster for forest in cells
    #   -d [ --depth ] arg         slide depth for all start cells (default is 1.0)
    #   -e [ --density ] arg       moving mass density (kg/m3) (default is 1600)
    #   -g [ --height ] arg        maximum flow height (default is 0.5)
    #   -t [ --slidetype ] arg     slide type (0=Loose material/high fluidity,
    #                              1=Loose material/normal fluidity, 2=Adhesive
    #                              material/high fluidity, 3=Adhesive material/normal
    #                              fluidity) (default is 0)
    #   -p [ --pressuremodel ] arg pressure model (0=Scheidegger, 1=VKG/Fine-grained,
    #                              2=VKG/Blocky, 3=Hungr) (default is 0)
    #   -s [ --sim ] arg           number of simulations (default is 10)

    if deposit==0: # Simulation without tree deposit
        if friction==0: # Simulation without additional friction
            slideforce_cmd = "\"{0}/slideForceCmd_64_mu_dynamic.exe\" --outputdir {1} --dem {2} --start {3} --density {4} --slidetype {5} --depth {6} --sim {7} --pressuremodel {8} --height {9}"
            sf_cmd = slideforce_cmd.format(slideforce_path, output_path, dem_path, start_mask_path, moving_mass_density, slide_type, slide_thickness, num_simulations, pressure_model, max_flow_height)
        else: # Simulation with additional friction
            slideforce_cmd = "\"{0}/slideForceCmd_64_mu_dynamic.exe\" --outputdir {1} --dem {2} --start {3} --density {4} --slidetype {5} --depth {6} --sim {7} --pressuremodel {8} --height {9} --friction {10}"
            sf_cmd = slideforce_cmd.format(slideforce_path, output_path, dem_path, start_mask_path, moving_mass_density, slide_type, slide_thickness, num_simulations, pressure_model, max_flow_height,  forest_friction_path)
    else: # Simulation with tree deposit
        if friction==0: # Simulation without additional friction
            slideforce_cmd = "\"{0}/slideForceCmd_64_mu_dynamic.exe\" --outputdir {1} --dem {2} --start {3} --density {4} --slidetype {5} --depth {6} --sim {7} --pressuremodel {8} --height {9} --forest {10}"
            sf_cmd = slideforce_cmd.format(slideforce_path, output_path, dem_path, start_mask_path, moving_mass_density, slide_type, slide_thickness, num_simulations, pressure_model, max_flow_height, tree_deposit_path)
        else: # Simulation with additional friction
            slideforce_cmd = "\"{0}/slideForceCmd_64_mu_dynamic.exe\" --outputdir {1} --dem {2} --start {3} --density {4} --slidetype {5} --depth {6} --sim {7} --pressuremodel {8} --height {9} --forest {10} --friction {11}"
            sf_cmd = slideforce_cmd.format(slideforce_path, output_path, dem_path, start_mask_path, moving_mass_density, slide_type, slide_thickness, num_simulations, pressure_model, max_flow_height, tree_deposit_path,  forest_friction_path)


    print(sf_cmd)
    os.system(sf_cmd)

    # Delete unneeded outputs to save space
    for filepath in glob.iglob(os.path.join(output_path,"*.tif")):
        if "reach_probability" in filepath or "V_mean" in filepath or "Fh_mean" in filepath:
            continue
        else:
            os.remove(filepath)

    return slide_path


# Entry function that initiates and coordinates the process
def process_sf_simulations(cfg):
    start_time = time.time()

    print("Start Setup",datetime.now())
    
    # Preparing records for parallel processing describing the hisorical slides
    output_path = cfg["output_path"]
    utilities.ensure_dir(output_path)
    slide_records = prepare_tiles(cfg)

    # Running simulations
    num_processes = cfg["num_processes"]

    print("Start Simulating", datetime.now())
    # process_records_linear(slide_records, process_simulation, cfg, num_processes = num_processes)
    process_records(slide_records, process_simulation, cfg, num_processes = num_processes)

    print("TOTAL PROCESSING TIME: %s (h:min:sec)" % str(timedelta(seconds=(time.time() - start_time))))


# Default entry point
if __name__ == "__main__":
    freeze_support()

    cfg = {
        "output_path": "E:/GIS_Projekte/Paper_4/data/sensitivity_forest/", # path where the simulations will be performed
        "slideforce_path": "E:/GIS_Projekte/Paper_4/slideforce", # path to SlideForce installation  
        "num_processes": 30, # number of parallel simulations
        "moving_mass_density": 1600, # moving mass density for simulations
        "num_simulations": 100, # number of repeated simulations
        "pressure_model": 0, # pressure model for simulations (Scheidegger)
    }

    process_sf_simulations(cfg)

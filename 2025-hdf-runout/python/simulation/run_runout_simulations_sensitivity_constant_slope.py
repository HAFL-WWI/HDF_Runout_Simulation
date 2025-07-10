######################################################################
# Copyright (C) 2025 BFH
#
# Script for running the runout simulations for the sensitivity 
# analyis on a constant gradient slope using variable friction.
#
# Author: Christoph Schaller, BFH-HAFL, June 2025
######################################################################

#
# Imports
#

from osgeo import gdal
from osgeo import gdal_array
from osgeo.gdalconst import *
from osgeo.gdalnumeric import *
gdal.UseExceptions()

import numpy as np
import pandas as pd
import geopandas as gpd
import os
import sys 

import math

from datetime import datetime
import time

# Adding the path to the utilities.py to the path and import
utility_path = "E:/GIS_Projekte/GHK/dataprocessing_be/utilities"
# utility_path = os.environ["GHK_UTILITY_PATH"]
sys.path.append(utility_path)
import utilities

import glob

from shapely import geometry
from shapely.geometry import shape, Point, mapping

from fiona import collection

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
    param_in_path = os.path.join(output_path,"sensitivity_param_values_infinite_slope_variable.txt")
    print(param_in_path)
    param_values = np.loadtxt(param_in_path) 

    folders = glob.glob(os.path.join(output_path,"slide_*"))
    print(len(folders))
    tile_records = []

    # Create a record for each parameter combination
    for i in range(len(param_values)):
        [thickness, density, height, type, friction, deposit, start_size,slope,dbh,stemdensity] = param_values[i]

        output_name = "sim_{0}".format(i)
        sim_path  = os.path.join(output_path,output_name)

        tile_record = {
            "slide_thickness": thickness,
            "output_path": output_path,
            "sim_path": sim_path,
            "step":i,
            "output_name":output_name,

            "slide_type": int(type),
            "max_flow_height": height,
            "moving_mass_density": density,
            "friction": friction,
            "deposit": int(deposit),
            "start_size": int(start_size),
            "slope": slope,
            "dbh": dbh,
            "stemdensity": stemdensity,
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
    friction = sim_record["friction"]
    deposit = sim_record["deposit"]

    start_size = sim_record["start_size"]
    slope = sim_record["slope"]
    dbh = sim_record["dbh"]
    stemdensity = sim_record["stemdensity"]
    step = sim_record["step"]
    

    # Get parameters for simulation from config
    num_simulations = cfg["num_simulations"]
    pressure_model = cfg["pressure_model"]
    slideforce_path = cfg["slideforce_path"]

    output_path = cfg["output_path"]
    sim_path = sim_record["sim_path"]
    step = sim_record["step"]
    output_name = sim_record["output_name"]

    # Create output directory
    output_base_path = sim_path
    utilities.ensure_dir(output_base_path)

    #
    # Generate simulation inputs
    #

    raster_resolution = cfg["resolution"]
    crs = "EPSG:{0}".format(cfg["epsg"])
    tree_seed = 42

    # Dimensions of the slope
    z_max = 1000
    x_min = 0
    x_max = 1500
    y_min = 0
    y_max = 205

    # Generate points representing the pixels of the DEM
    points = []
    for x in range(0,x_max//raster_resolution):
        if x>0:
            # Decrease height based on slope parameter
            z -= raster_resolution*math.tan(math.radians(slope))
        else:
            z = z_max
        for y in range(0,y_max//raster_resolution):
            point = [x,y,x*raster_resolution+raster_resolution/2,y*raster_resolution+raster_resolution/2,z]
            points.append(point)

    # Save points to a geopackage and rasterize to generate the DEM
    points_df = pd.DataFrame(points, columns=["col","row","x_coord","y_coord","z"])
    points_gdf = gpd.GeoDataFrame(points_df, geometry=gpd.points_from_xy(points_df["x_coord"], points_df["y_coord"]), crs=crs)
    points_gpkg_path = os.path.join(output_base_path,"points.gpkg")

    points_gdf.to_file(points_gpkg_path, driver="GPKG")

    dem_path = os.path.join(output_base_path,"dem.tif")
    with gdal.Rasterize(destNameOrDestDS=dem_path, srcDS=points_gpkg_path, outputSRS=crs, format="GTiff", outputType=gdalconst.GDT_Float32, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], attribute="z", initValues=0) as ds:
        ds=None

    # Generate start area
    slide_radii = {0:5.64,1:7.98,2:13.02,3:16.93}
    all_touched = False
    x_center = 27.5
    y_center = 102.5

    # Get start area matching the parameter value
    radius = slide_radii[start_size]

    # Generate a shapefile ith a circular start area located at the top of the slope
    start_area_shp_name = "start_area.shp"
    start_area_shp_path = os.path.join(output_base_path,start_area_shp_name)

    schema = { 'geometry': 'Polygon', 'properties': { 'name': 'str' } }
    with collection(start_area_shp_path, "w", "ESRI Shapefile", schema) as output:
        center_point = Point(x_center, y_center)
        output.write({
                'properties': {
                    'name': "{0} {1}".format(x_center, y_center)
                },
                'geometry': mapping(shape(center_point).buffer(radius))
            })

    # Generate start cell raster by rasterizing the shapefile
    start_area_rast_name = "start_area.tif"
    start_mask_path = os.path.join(output_base_path,start_area_rast_name)

    with gdal.Rasterize(destNameOrDestDS=start_mask_path, srcDS=start_area_shp_path, outputSRS=crs, format="GTiff", burnValues = 1, initValues=-9999, noData =-9999, outputType=gdalconst.GDT_Int16, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], allTouched=all_touched) as ds:
        ds=None

    utilities.delete_shapefile(start_area_shp_path)

    # Generate forest deposit raster
    # Randomly sample points to match the stem density
    plot_bbox_polygon = geometry.box(minx=x_min,maxx=x_max,miny=y_min,maxy=y_max)
    print(plot_bbox_polygon)
    area_ha = (x_max-x_min)*(y_max-y_min)/10000
    n_points = int(area_ha*stemdensity)

    plot_df = pd.DataFrame({"id":[1],"dbh_mean":[dbh],"dbhstdev":[7.5],"npoints":[n_points]})
    plot_gdf = gpd.GeoDataFrame(plot_df, geometry=[plot_bbox_polygon], crs=crs)
    sample_trees = plot_gdf["geometry"].sample_points(n_points,method="uniform", rng=tree_seed )
    
    # Apply deposit capacity formula to points/trees
    sample_trees_df = sample_trees.explode(index_parts=True).reset_index()
    sample_trees_df["deposit"] = dbh*1.13277/100-0.18478
    sample_trees_df["deposit"] = sample_trees_df["deposit"].apply(lambda x: max(0,x))

    # Calculate the corner coordinate of the raster cell containing the tree 
    sample_trees_df["x"] = sample_trees_df.geometry.x
    sample_trees_df["y"] = sample_trees_df.geometry.y
    sample_trees_df["x_corner"] = sample_trees_df["x"].apply(lambda x: int(math.floor(x/raster_resolution)*raster_resolution))
    sample_trees_df["y_corner"] = sample_trees_df["y"].apply(lambda x: int(math.floor(x/raster_resolution)*raster_resolution))

    # Agregate deposition capacity per cell and save to geopackage
    sample_tree_grouped_df = sample_trees_df.groupby(["x_corner","y_corner"], as_index =False)["deposit"].sum()
    sample_tree_grouped_df.reset_index(inplace=True)
    geom = gpd.points_from_xy(sample_tree_grouped_df["x_corner"]+2.5, sample_tree_grouped_df["y_corner"]+2.5)
    sample_tree_grouped_gdf = gpd.GeoDataFrame(sample_tree_grouped_df,geometry=geom)

    # Rasterize the deposition capacity
    sample_treefile_gpk = os.path.join(output_base_path,"sample_treefile.gpkg")
    tree_grouped_layer = "tree_grouped_layer"
    sample_tree_grouped_gdf.to_file(sample_treefile_gpk, driver="GPKG", layer=tree_grouped_layer,crs=crs)

    tree_deposit_path = os.path.join(output_base_path,"sample_tree_deposit.tif")
    with gdal.Rasterize(destNameOrDestDS=tree_deposit_path, srcDS=sample_treefile_gpk, outputSRS=crs, layers=tree_grouped_layer, format="GTiff", outputType=gdalconst.GDT_Float32, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], attribute="deposit", initValues=0) as ds:
        ds=None
    utilities.unset_nodata(tree_deposit_path)

    # Generate friction input raster
    # This assumes forest on the entire slope with a mean stem density and dbh
    stem_density_sample = stemdensity/10000
    dbh_mean_sample = dbh/100 
    # Product of the mean stem count and DBH within the perimeter 
    dn_sample = stem_density_sample*dbh_mean_sample

    # Generate points representing the pixels in the raster
    for x in range(0,x_max//raster_resolution):
        for y in range(0,y_max//raster_resolution):
            point = [x,y,x*raster_resolution+raster_resolution/2,y*raster_resolution+raster_resolution/2,dn_sample]
            points.append(point)


    # Save points to a geopackage and rasterize to generate the friction input raster
    points_df = pd.DataFrame(points, columns=["col","row","x_coord","y_coord","dn_sample"])
    points_gdf = gpd.GeoDataFrame(points_df, geometry=gpd.points_from_xy(points_df["x_coord"], points_df["y_coord"]), crs=crs)
    points_gpkg_path = os.path.join(output_base_path,"points_dn.gpkg")

    points_gdf.to_file(points_gpkg_path, driver="GPKG")

    dn_sample_path = os.path.join(output_base_path,"forest_dn.tif")
    with gdal.Rasterize(destNameOrDestDS=dn_sample_path, srcDS=points_gpkg_path, outputSRS=crs, format="GTiff", outputType=gdalconst.GDT_Float32, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], attribute="dn_sample", initValues=0) as ds:
        ds=None

    utilities.delete_gpkg(points_gpkg_path)

    utilities.unset_nodata(dn_sample_path)


    #
    # Run simulation
    #

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
            sf_cmd = slideforce_cmd.format(slideforce_path, output_base_path, dem_path, start_mask_path, moving_mass_density, slide_type, slide_thickness, num_simulations, pressure_model, max_flow_height)
        else: # Simulation with additional friction
            slideforce_cmd = "\"{0}/slideForceCmd_64_mu_dynamic.exe\" --outputdir {1} --dem {2} --start {3} --density {4} --slidetype {5} --depth {6} --sim {7} --pressuremodel {8} --height {9} --friction {10}"
            sf_cmd = slideforce_cmd.format(slideforce_path, output_base_path, dem_path, start_mask_path, moving_mass_density, slide_type, slide_thickness, num_simulations, pressure_model, max_flow_height,  dn_sample_path)
    else:  # Simulation with tree deposit
        if friction==0: # Simulation without additional friction
            slideforce_cmd = "\"{0}/slideForceCmd_64_mu_dynamic.exe\" --outputdir {1} --dem {2} --start {3} --density {4} --slidetype {5} --depth {6} --sim {7} --pressuremodel {8} --height {9} --forest {10}"
            sf_cmd = slideforce_cmd.format(slideforce_path, output_base_path, dem_path, start_mask_path, moving_mass_density, slide_type, slide_thickness, num_simulations, pressure_model, max_flow_height, tree_deposit_path)
        else: # Simulation with additional friction
            slideforce_cmd = "\"{0}/slideForceCmd_64_mu_dynamic.exe\" --outputdir {1} --dem {2} --start {3} --density {4} --slidetype {5} --depth {6} --sim {7} --pressuremodel {8} --height {9} --forest {10} --friction {11}"
            sf_cmd = slideforce_cmd.format(slideforce_path, output_base_path, dem_path, start_mask_path, moving_mass_density, slide_type, slide_thickness, num_simulations, pressure_model, max_flow_height, tree_deposit_path,  dn_sample_path)

    print(sf_cmd)
    os.system(sf_cmd)


    #
    # Post-process simulation
    #

    # Read the needed output rasters
    runout_path = os.path.join(output_base_path,"reach_probability.tif")
    velocity_path = os.path.join(output_base_path,"V_mean.tif")
    fh_path = os.path.join(output_base_path,"Fh_mean.tif")
    p_path = os.path.join(output_base_path,"P_mean.tif")
    
    runout_values = gdal_array.LoadFile(runout_path)
    v_values = gdal_array.LoadFile(velocity_path)
    fh_values = gdal_array.LoadFile(fh_path)
    p_values = gdal_array.LoadFile(p_path)


    # Determine the global mean/max values for velocity, pressure, and flow height 
    n_runout = np.sum(runout_values>0)
    if np.sum(v_values)>0:
        v_max = np.max(v_values[v_values>0])
        v_mean = np.mean(v_values[v_values>0])
    else:
        v_max = 0
        v_mean = 0

    if np.sum(p_values[p_values>0])>0:
        p_max = np.max(p_values[p_values>0])
        p_mean = np.mean(p_values[p_values>0])
    else:
        p_max = 0
        p_mean = 0

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

    # Create a CSV with the results
    rec = {"step":step,"thickness":slide_thickness,"density":moving_mass_density,"height":max_flow_height,"type":slide_type,"friction":friction,"deposit":deposit,"size":start_size,"slope":slope,"dbh":dbh,"stemdensity":stemdensity,
           "x_min":x_min,"x_max":x_max,"y_min":y_min,"y_max":y_max,"w":x_max-x_min,"h":y_max-y_min,"max_spread":max_spread,"n_runout":n_runout,"v_max":v_max,"v_mean":v_mean,"p_max":p_max,"p_mean":p_mean,"fh_max":fh_max,"fh_mean":fh_mean}
    
    sim_stats_df = pd.DataFrame([rec])
    sim_stats_df.to_csv(os.path.join(output_base_path,"sim_stats.csv"),sep=",",index=False)

    #
    # Cleanup 
    #
    for filepath in glob.iglob(os.path.join(output_base_path,"*.tif")):
        if "reach_probability" in filepath:
            continue
        else:
            os.remove(filepath)

    for filepath in glob.iglob(os.path.join(output_base_path,"*.gpkg")):
        os.remove(filepath)

    return sim_path


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
        "output_path": "E:/GIS_Projekte/Paper_4/data/sensitivity_constant_slope/", # path where the simulations will be performed
        "slideforce_path": "E:/GIS_Projekte/Paper_4/slideforce", # path to SlideForce installation  
        "num_processes": 30, # number of parallel simulations
        "moving_mass_density": 1600, # moving mass density for simulations
        "num_simulations": 100, # number of repeated simulations
        "pressure_model": 0, # pressure model for simulations (Scheidegger)
        "resolution": 5, # resolution of the generated rasters
        "epsg": 2056, # EPSG of coordinate system for the generated data; LV95  
    }

    process_sf_simulations(cfg)

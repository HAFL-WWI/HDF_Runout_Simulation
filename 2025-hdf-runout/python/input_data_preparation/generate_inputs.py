######################################################################
# Copyright (C) 2025 BFH
#
# Script for preparing all necessary inputs for the SlideForce
# simulations of the historical landslides.
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
import glob
import os
import sys 

import math

# Adding the path to the utilities.py to the path and import
utility_path = "E:/GIS_Projekte/Paper_4/code/utilities"
# utility_path = os.environ["GHK_UTILITY_PATH"]
sys.path.append(utility_path)
import utilities

#
# Configs and paths
#

saga_path = "C:/OSGeo4W64/apps/saga-ltr"

# Folder containing all current swissALTI3D DEM tiles at 0.5m resulution 
dem_input_folder = "P:/HAFL/9 Share/PermanenteOrdner/Geodaten/Nationale_Daten/Nationale_Daten_SWISSTOPO/swissALTI3D/DTM_05m/tiles"

slide_input_gpkg = "P:/HAFL/7 WWI/74b FF GNG/742b Aktuell/2022_2025_Diss_ChristophSchaller/Paper_4/data/Dataset_Runout_Merged.gpkg"
slide_input_layer = "inputs_merged"
slide_layer = "slide"

fintch_db_paths = {
    "BE":["F:/fint-ch/AP10__Validierung_der_EBD_Analysen/FINT_v2_BE.gpkg","processed_trees"],
    "GR":["F:/fint-ch/AP10__Validierung_der_EBD_Analysen/FINT_v2_GR.gpkg","processed_trees"],
    "SG":["F:/fint-ch/AP10__Validierung_der_EBD_Analysen/FINT_v2_SG.gpkg","sg_processed_tree"],
    "LU":["F:/fint-ch/AP10__Validierung_der_EBD_Analysen/FINT_v2_LU.gpkg","lu_processed_tree"],
    "TI":["F:/fint-ch/AP10__Validierung_der_EBD_Analysen/FINT_v2_TI.gpkg","ti_processed_tree"],
    "FL":["F:/fint-ch/AP10__Validierung_der_EBD_Analysen/FINT_v2_FL.gpkg","FL_processed_tree"],
}

tlm_gpkg = "E:/GIS_Projekte/Geodaten/swisstlm3d_2024-03_2056_5728/SWISSTLM3D_2024_LV95_LN02.gpkg"

raster_resolution = 5

forest_friction_factor = 0.02

slide_buffer = 500

# Path to generate slide inputs in
output_base_path = "E:/GIS_Projekte/Paper_4/data/simulation_inputs"

# # Uncomment to generate inputs for sensitivity analysis
# output_base_path = "E:/GIS_Projekte/Paper_4/data/sensitivity_forest"


utilities.ensure_dir(output_base_path)

#
# Generating the simulation inputs per slide
#     

# Read the reference data
slide_df = gpd.read_file(slide_input_gpkg,layer=slide_input_layer,engine="pyogrio")

for slide_id in slide_df["slide_id"].unique():
    # # Uncomment to generate only inputs for forest plots for sensitivity analysis
    # if slide_df[(slide_df["slide_id"]==slide_id)&(slide_df["area_type"]==1)].iloc[0]["in_forest"]==0:
    #     continue
    
    # Prepare output folder for the current slide
    slide_folder_name = "slide_{0}".format(slide_id)
    slide_folder_path = os.path.join(output_base_path,slide_folder_name)
    utilities.ensure_dir(slide_folder_path)

    # Create Geopackage with just the current slide
    runout_gpkg = os.path.join(slide_folder_path,"slide.gpkg")
    with gdal.VectorTranslate(destNameOrDestDS=runout_gpkg, srcDS=slide_input_gpkg, format="GPKG", layerName=slide_layer, SQLStatement="select * FROM {0} WHERE slide_id={1}".format(slide_input_layer,slide_id), SQLDialect="SQLite", accessMode="overwrite") as ds :
        ds = None


    #
    # DEM preparation
    #

    # Determine extent for DEM
    extent_gpkg = utilities.get_gpkg_extent(runout_gpkg,slide_layer)
    print(extent_gpkg)

    x_min = raster_resolution*(math.floor(extent_gpkg[0]/raster_resolution)) - slide_buffer
    x_max = raster_resolution*(math.ceil(extent_gpkg[1]/raster_resolution))  + slide_buffer
    y_min = raster_resolution*(math.floor(extent_gpkg[2]/raster_resolution)) - slide_buffer
    y_max = raster_resolution*(math.ceil(extent_gpkg[3]/raster_resolution))  + slide_buffer

    # Determine input tiles
    x_km_min = math.floor(x_min/1000)
    y_km_min = math.floor(y_min/1000)
    x_km_max = math.floor(x_max/1000)
    y_km_max = math.floor(y_max/1000)

    print(x_km_min,y_km_min,x_km_max,y_km_max)

    dem_tiles = []
    x_km_cur = x_km_min
    while x_km_cur<=x_km_max:
        y_km_cur = y_km_min
        while y_km_cur<=y_km_max:
            tile_pattern = "swissalti3d_*_{0}-{1}_0.5_2056_*.tif".format(x_km_cur,y_km_cur)
            tiles = glob.glob(os.path.join(dem_input_folder,tile_pattern))
            dem_tiles.append(tiles[-1])
            y_km_cur+=1
        x_km_cur+=1

    # Merge DEM input tiles and resample to 5m
    dem_0_5_vrt_path = os.path.join(slide_folder_path,"dem_0_5.vrt")
    gdal.BuildVRT(dem_0_5_vrt_path,srcDSOrSrcDSTab=dem_tiles)
    
    dem_5_path = os.path.join(slide_folder_path,"dem_5.tif")

    with gdal.Warp(destNameOrDestDS=dem_5_path,srcDSOrSrcDSTab=dem_0_5_vrt_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}, xRes=raster_resolution, yRes=raster_resolution, dstSRS="EPSG:2056", resampleAlg=gdal.GRA_Average, outputBounds=[x_min,y_min,x_max,y_max]) as ds:
        ds=None

    # Sink fill the DEM using Wang & Liu
    dem_5_sgrd_path = dem_5_path.replace(".tif",".sdat") 
    with gdal.Translate(destName=dem_5_sgrd_path, srcDS=dem_5_path, format="SAGA") as ds:
        ds=None

    fill_sinks_cmd = "{0}/saga_cmd ta_preprocessor 4 -ELEV {1} -FILLED {2}" #Wang & Liu

    dem_5_sgrd_filled_path = dem_5_sgrd_path.replace(".sdat","_filled.sdat") 

    dem_fill_cmd = fill_sinks_cmd.format(saga_path, dem_5_sgrd_path.replace(".sdat",".sgrd"), dem_5_sgrd_filled_path.replace(".sdat",".sgrd"))
    os.system(dem_fill_cmd)

    dem_5_filled_path = dem_5_sgrd_filled_path.replace(".sdat",".tif")
    with gdal.Translate(destName=dem_5_filled_path, srcDS=dem_5_sgrd_filled_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}) as ds:
        ds=None

    # Clip start area DEM 0.5m
    dem_start_05_path = os.path.join(slide_folder_path,"dem_start_05.tif")
    with gdal.Warp(destNameOrDestDS=dem_start_05_path,srcDSOrSrcDSTab=dem_0_5_vrt_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"},  dstSRS="EPSG:2056", resampleAlg=gdal.GRA_Average, cutlineDSName=runout_gpkg, cutlineSQL="SELECT * from {0} WHERE area_type=1".format(slide_layer)) as ds:
        ds=None

    
    #
    # Preparing slide masks
    #

    # Generate slide mask with start area and runout
    slide_mask_path = os.path.join(slide_folder_path,"slide_mask.tif")
    with gdal.Rasterize(destNameOrDestDS=slide_mask_path, srcDS=runout_gpkg, layers=slide_layer, format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0, noData=None, allTouched=True) as ds:
        ds=None
    utilities.unset_nodata(slide_mask_path)


    # Generate failure area mask (start cells)
    start_mask_path = os.path.join(slide_folder_path,"start_mask.tif")
    with gdal.Rasterize(destNameOrDestDS=start_mask_path, srcDS=runout_gpkg, SQLStatement="SELECT * from {0} WHERE area_type=1".format(slide_layer), format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0, noData=None, allTouched=True) as ds:
        ds=None
    utilities.unset_nodata(start_mask_path)

    #
    # Generate the tree deposit raster
    #

    treefile_gpk = os.path.join(slide_folder_path,"treefile.gpkg")
    fintch_db = fintch_db_paths[slide_df[slide_df["slide_id"]==slide_id]["canton"].unique()[0]]

    # Generate a geopackage containing the trees within the perimeter and calculate the deposition capacity per tree
    treefile_gpk_tmp = os.path.join(slide_folder_path, "treefile_tmp.gpkg")
    tree_layer = "tree_layer"
    tree_grouped_layer = "tree_grouped_layer"
    with gdal.VectorTranslate(destNameOrDestDS=treefile_gpk_tmp, srcDS=fintch_db[0], format="GPKG", layerName=tree_layer, layers=fintch_db[1], accessMode="overwrite", spatFilter=[x_min,y_min,x_max,y_max]) as ds : #FINT-CH
        gdal.VectorInfo(ds=ds, SQLStatement="ALTER TABLE {0} ADD COLUMN deposit FLOAT;".format(tree_layer))
        gdal.VectorInfo(ds=ds, SQLStatement="UPDATE {0} SET deposit=CASE WHEN (1.13277*bhd/100-0.18478)>0 THEN (1.13277*bhd/100-0.18478) ELSE 0 END".format(tree_layer) , SQLDialect="SQLite")
        ds = None

    # Read trees into a geodataframe
    tree_df = gpd.read_file(treefile_gpk_tmp,layer=tree_layer,engine="pyogrio")

    # Calculate the corner coordinate of the raster cell containing the tree 
    tree_df["x_corner"] = tree_df["x"].apply(lambda x: int(math.floor(x/raster_resolution)*raster_resolution))
    tree_df["y_corner"] = tree_df["y"].apply(lambda x: int(math.floor(x/raster_resolution)*raster_resolution))

    # Agregate deposition capacity per cell and save to geopackage
    tree_grouped_df = tree_df.groupby(["x_corner","y_corner"], as_index =False)["deposit"].sum()
    tree_grouped_df.reset_index(inplace=True)
    geometry = gpd.points_from_xy(tree_grouped_df["x_corner"]+2.5, tree_grouped_df["y_corner"]+2.5)
    tree_grouped_gdf = gpd.GeoDataFrame(tree_grouped_df,geometry=geometry)
    tree_grouped_gdf.to_file(treefile_gpk_tmp, driver="GPKG", layer=tree_grouped_layer,crs=tree_df.crs)

    # Rasterize the deposition capacity
    tree_deposit_path = os.path.join(slide_folder_path,"tree_deposit.tif")
    with gdal.Rasterize(destNameOrDestDS=tree_deposit_path, srcDS=treefile_gpk_tmp, layers=tree_grouped_layer, format="GTiff", outputType=gdalconst.GDT_Float32, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"No"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], attribute="deposit", initValues=0) as ds:
        ds=None
    utilities.unset_nodata(tree_deposit_path)

    # utilities.delete_gpkg(treefile_gpk_tmp)

    #
    # Generate swissTLM3D forest mask
    #
    forest_gpkg = os.path.join(slide_folder_path,"tlm_forest.gpkg")
    tlm_groundcover_layer = "tlm_bb_bodenbedeckung"
    forest_layer = "forest"
    with gdal.VectorTranslate(destNameOrDestDS=forest_gpkg, srcDS=tlm_gpkg, format="GPKG", layerName=forest_layer , SQLStatement="select * from {0} WHERE (OBJEKTART = 'Wald' OR  OBJEKTART = 'Wald offen')".format(tlm_groundcover_layer),clipSrc=[x_min,y_min,x_max,y_max], spatFilter=[x_min,y_min,x_max,y_max]) as ds :
        ds = None

    forest_mask_path = os.path.join(slide_folder_path,"forest_tlm_mask.tif")
    with gdal.Rasterize(destNameOrDestDS=forest_mask_path, srcDS=forest_gpkg, layers=forest_layer, format="GTiff", outputType=gdalconst.GDT_Byte, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"No"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=1, initValues=0) as ds:
        ds=None
    utilities.unset_nodata(forest_mask_path)

    #
    # Generate forest friction raster with fixed mu
    #
    friction_path = os.path.join(slide_folder_path,"forest_friction.tif")
    with gdal.Rasterize(destNameOrDestDS=friction_path, srcDS=forest_gpkg, layers=forest_layer, format="GTiff", outputType=gdalconst.GDT_Float32, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=forest_friction_factor, initValues=0, noData=None, allTouched=True) as ds:
        ds=None
    utilities.unset_nodata(friction_path)

    #
    # Generate forest deposit raster from mean DBH and stem density based on field survey
    #
    slide = slide_df[slide_df["slide_id"]==slide_id].iloc[0]

    # Only generate if surveyed DBH is present
    if slide["dbh_mean_cm"]>0:
        forest_df = gpd.read_file(forest_gpkg,layer=forest_layer,engine="pyogrio")
        forest_areas = forest_df.geometry.area/10000
        stem_density = forest_areas.apply(lambda x: int(x*slide["stemdensity_per_ha"]))

        # Randomly sample points within the forest mask to match the surveyed stem density
        sample_trees = forest_df.geometry.sample_points(size=stem_density, method="uniform", rng=42)  
        sample_trees_df = sample_trees.explode(index_parts=True).reset_index()
        sample_trees_df["deposit"] = slide["dbh_mean_cm"]*1.13277/100-0.18478
        sample_trees_df["deposit"] = sample_trees_df["deposit"].apply(lambda x: max(0,x))

        sample_trees_df["x"] = sample_trees_df.geometry.x
        sample_trees_df["y"] = sample_trees_df.geometry.y
        sample_trees_df["x_corner"] = sample_trees_df["x"].apply(lambda x: int(math.floor(x/raster_resolution)*raster_resolution))
        sample_trees_df["y_corner"] = sample_trees_df["y"].apply(lambda x: int(math.floor(x/raster_resolution)*raster_resolution))

        sample_tree_grouped_df = sample_trees_df.groupby(["x_corner","y_corner"], as_index =False)["deposit"].sum()
        sample_tree_grouped_df.reset_index(inplace=True)
        geometry = gpd.points_from_xy(sample_tree_grouped_df["x_corner"]+2.5, sample_tree_grouped_df["y_corner"]+2.5)
        sample_tree_grouped_gdf = gpd.GeoDataFrame(sample_tree_grouped_df,geometry=geometry)

        sample_treefile_gpk = os.path.join(slide_folder_path,"sample_treefile.gpkg")
        sample_tree_grouped_gdf.to_file(sample_treefile_gpk, driver="GPKG", layer=tree_grouped_layer,crs=forest_df.crs)

        sample_tree_deposit_path = os.path.join(slide_folder_path,"sample_tree_deposit.tif")
        with gdal.Rasterize(destNameOrDestDS=sample_tree_deposit_path, srcDS=sample_treefile_gpk, layers=tree_grouped_layer, format="GTiff", outputType=gdalconst.GDT_Float32, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"No"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], attribute="deposit", initValues=0) as ds:
            ds=None
        utilities.unset_nodata(sample_tree_deposit_path)

        utilities.delete_gpkg(sample_treefile_gpk)


    #
    # Generate forest friction raster with finput for variable mu
    #

    # Determine the mean stem density and dbh form the detected trees 
    stem_density_absolute = len(tree_df)
    forest_df = gpd.read_file(forest_gpkg,layer=forest_layer,engine="pyogrio")
    forest_areas = forest_df.geometry.area
    area_ha = forest_areas.sum()
    stem_density = stem_density_absolute/(area_ha)

    dbh_mean = tree_df["bhd"].mean()/100 

    # Product of the mean stem density and DBH within the perimeter 
    dn = stem_density*dbh_mean

    dn_path = os.path.join(slide_folder_path,"forest_dn.tif")
    with gdal.Rasterize(destNameOrDestDS=dn_path, srcDS=forest_gpkg, layers=forest_layer, format="GTiff", outputType=gdalconst.GDT_Float32, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=dn, initValues=0, noData=None, allTouched=True) as ds:
        ds=None
    utilities.unset_nodata(dn_path)


    #
    # Generate forest friction raster with finput for variable mu for sampled trees
    #

    # Only generate if surveyed value is present
    if slide["dbh_mean_cm"]>0:   
        stem_density_sample = slide["stemdensity_per_ha"]//10000
        dbh_mean_sample = slide["dbh_mean_cm"]/100 
        dn_sample = stem_density_sample*dbh_mean_sample
        dn_sample_path = os.path.join(slide_folder_path,"sample_forest_dn.tif")
        with gdal.Rasterize(destNameOrDestDS=dn_sample_path, srcDS=forest_gpkg, layers=forest_layer, format="GTiff", outputType=gdalconst.GDT_Float32, creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"YES"}, xRes=raster_resolution, yRes=raster_resolution, outputBounds=[x_min,y_min,x_max,y_max], burnValues=dn_sample, initValues=0, noData=None, allTouched=True) as ds:
            ds=None
        utilities.unset_nodata(dn_path)


    #
    # Cleanup of unneeded files
    #
    utilities.delete_raster(dem_0_5_vrt_path)
    utilities.delete_raster(dem_5_sgrd_path)
    utilities.delete_raster(dem_5_sgrd_filled_path)



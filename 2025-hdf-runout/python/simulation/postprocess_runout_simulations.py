######################################################################
# Copyright (C) 2025 BFH
#
# Script for post processing the runout simulations for the 
# historical landslides.
#
# Author: Christoph Schaller, BFH-HAFL, June 2025
######################################################################

#
# Imports
#

from osgeo import gdal
from osgeo.gdalconst import *
from osgeo.gdalnumeric import *
from osgeo_utils import gdal_calc, gdal_polygonize
from osgeo import gdal_array

gdal.UseExceptions()

import rasterio as rio

import numpy as np
import pandas as pd
import geopandas as gpd
import glob
import os
import sys 

import math

from sklearn.metrics import roc_auc_score
from scipy import ndimage

# Adding the path to the utilities.py to the path and import
utility_path = "E:/GIS_Projekte/Paper_4/code/utilities"
# utility_path = os.environ["GHK_UTILITY_PATH"]
sys.path.append(utility_path)
import utilities

def process_result(input_base_path, runout_threshold, runout_threshold_comparator, raster_resolution, forest_labels, slide_type_labels):
    # Lists to store results
    sim_records = []
    slide_poly_records = []
    sim_poly_records = []
    mask_poly_records = []

    # For all landslides
    slide_folders = glob.glob(os.path.join(input_base_path,"slide_*"))
    print(len(slide_folders))
    for slide_folder_path in slide_folders:
        # Skip files
        if not os.path.isdir(slide_folder_path):
            continue

        # Extract slide ID from folder name
        slide_folder_name = os.path.basename(slide_folder_path)
        print(slide_folder_path,slide_folder_name)
        slide_id = int(slide_folder_name.replace("slide_",""))

        print(slide_folder_name)

        # Get subfolder for the individual simulations
        sim_folders = glob.glob(os.path.join(slide_folder_path,"sim_*"))
        print(len(sim_folders))

        # Determine extent for historic slide
        slide_gpkg = os.path.join(slide_folder_path,"slide.gpkg")
        slide_layer = "slide"

        extent_slide_gpkg = utilities.get_gpkg_extent(slide_gpkg,slide_layer)
        print(extent_slide_gpkg)

        x_min = raster_resolution*(math.floor(extent_slide_gpkg[0]/raster_resolution)) #- raster_resolution
        x_max = raster_resolution*(math.ceil(extent_slide_gpkg[1]/raster_resolution))  #+ raster_resolution
        y_min = raster_resolution*(math.floor(extent_slide_gpkg[2]/raster_resolution)) #- raster_resolution
        y_max = raster_resolution*(math.ceil(extent_slide_gpkg[3]/raster_resolution))  #+ raster_resolution

        # Read slide mask as base for UNION mask
        slide_mask_path = os.path.join(slide_folder_path,"slide_mask.tif")
        mask_values = gdal_array.LoadFile(slide_mask_path)

        # Polygonize slide mask 
        utilities.reset_nodata(slide_mask_path,0) # NoData=0 -> only pixels with value 1 are polygonized
        slide_mask_gpkg = slide_mask_path.replace(".tif",".gpkg")
        slide_mask_layer = "slide_mask"
        utilities.delete_file_by_path(slide_mask_gpkg)
        gdal_polygonize.gdal_polygonize(src_filename=slide_mask_path, band_number=1, dst_filename=slide_mask_gpkg, dst_layername=slide_mask_layer, driver_name="GPKG", overwrite=True, connectedness8=True)
        utilities.unset_nodata(slide_mask_path) # Remove NoData value again
        slide_mask_df = gpd.read_file(slide_mask_gpkg,layer=slide_mask_layer,engine="pyogrio")
        slide_mask_df["slide_id"] = slide_id
        slide_poly_records.append(slide_mask_df)

        # Determine Union of all runout extents and historic slide extent
        for sim_folder_path in sim_folders:
            # Determine extent for simulation runout
            runout_path = os.path.join(sim_folder_path,"reach_probability.tif")

            # Binarize runout with specified threshold
            runout_binary_path = runout_path.replace(".tif","_binary_{0}.tif".format(runout_threshold))
            with gdal_calc.Calc(calc="(A {0} {1})*1".format(runout_threshold_comparator,runout_threshold),A=runout_path, outfile=runout_binary_path, creation_options = ["COMPRESS=DEFLATE", "BIGTIFF=False"], NoDataValue=0, overwrite=True) as ds:
                ds=None

            # Read binary mask and add add to UNION mask
            runout_binary_values = gdal_array.LoadFile(runout_binary_path)
            mask_values = mask_values+runout_binary_values

            # Polygonize runout mask
            runout_binary_gpkg = runout_binary_path.replace(".tif",".gpkg")
            runout_layer = "runout"
            utilities.delete_file_by_path(runout_binary_gpkg)
            gdal_polygonize.gdal_polygonize(src_filename=runout_binary_path, band_number=1, dst_filename=runout_binary_gpkg, dst_layername=runout_layer, driver_name="GPKG", overwrite=True, connectedness8=True)
            utilities.unset_nodata(runout_binary_path)

            # Get extent from runout mask
            extent_runout_gpkg = utilities.get_gpkg_extent(runout_binary_gpkg,runout_layer)
            print(extent_runout_gpkg)

            # Round to resolution
            x_min_runout = raster_resolution*(math.floor(extent_runout_gpkg[0]/raster_resolution)) #- raster_resolution
            x_max_runout = raster_resolution*(math.ceil(extent_runout_gpkg[1]/raster_resolution))  #+ raster_resolution
            y_min_runout = raster_resolution*(math.floor(extent_runout_gpkg[2]/raster_resolution)) #- raster_resolution
            y_max_runout = raster_resolution*(math.ceil(extent_runout_gpkg[3]/raster_resolution))  #+ raster_resolution

            # "UNION" of extent with other masks
            x_min = min(x_min,x_min_runout)
            x_max = max(x_max,x_max_runout)
            y_min = min(y_min,y_min_runout)
            y_max = max(y_max,y_max_runout)

        # Add a buffer of 1 pixel to extent
        x_min = x_min - raster_resolution
        x_max = x_max + raster_resolution
        y_min = y_min - raster_resolution
        y_max = y_max + raster_resolution
        print([x_min,x_max,y_min,y_max])

        # Buffer union mask by 1 cell
        se = ndimage.generate_binary_structure(2, 2)
        mask_values = (mask_values>0).astype(int) # Binarize UNION mask
        mask_values = ndimage.binary_dilation(mask_values , structure=se).astype(int)

        # Save UNION mask
        union_mask_path = os.path.join(slide_folder_path, "slide_union_mask.tif")
        with rio.open(slide_mask_path) as src:
            with rio.open(
                union_mask_path,
                'w',
                driver='GTiff',
                height=src.shape[0],
                width=src.shape[1],
                count=1,
                dtype=mask_values.dtype,
                compress="lzw",
                crs=src.crs,
                transform=src.transform,
            ) as dst:
                dst.write(mask_values, 1)

        # Polygonize UNION mask
        utilities.reset_nodata(union_mask_path,0)
        union_mask_gpkg = union_mask_path.replace(".tif",".gpkg")
        union_mask_layer = "union_mask"
        utilities.delete_file_by_path(union_mask_gpkg)
        gdal_polygonize.gdal_polygonize(src_filename=union_mask_path, band_number=1, dst_filename=union_mask_gpkg, dst_layername=union_mask_layer, driver_name="GPKG", overwrite=True, connectedness8=True)
        utilities.unset_nodata(union_mask_path)
        union_mask_df = gpd.read_file(union_mask_gpkg,layer=union_mask_layer,engine="pyogrio")
        union_mask_df["slide_id"] = slide_id
        mask_poly_records.append(union_mask_df)

        # Process each simulation to extract performance metrics
        for sim_folder_path in sim_folders:
            # Extract the slide type and forest influence type from folder name
            sim_folder_name = os.path.basename(sim_folder_path)
            tokens = sim_folder_name.split("_")
            slide_type = int(tokens[2])
            forest = int(tokens[4])

            print(sim_folder_path)


            # Clip the rasters to the extent of the UNION mask
            runout_path = os.path.join(sim_folder_path,"reach_probability.tif")
            runout_binary_path = runout_path.replace(".tif","_binary_{0}.tif".format(runout_threshold))

            slide_mask_clipped_path = os.path.join(sim_folder_path,"slide_mask_clipped.tif")
            with gdal.Warp(destNameOrDestDS=slide_mask_clipped_path,srcDSOrSrcDSTab=slide_mask_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}, dstSRS="EPSG:2056", outputBounds=[x_min,y_min,x_max,y_max]) as ds:
                ds=None

            union_mask_clipped_path = os.path.join(sim_folder_path,"slide_union_mask_clipped.tif")
            with gdal.Warp(destNameOrDestDS=union_mask_clipped_path,srcDSOrSrcDSTab=union_mask_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}, dstSRS="EPSG:2056", outputBounds=[x_min,y_min,x_max,y_max]) as ds:
                ds=None

            runout_binary_clipped_path = os.path.join(sim_folder_path,"runout_binary_clipped.tif")
            with gdal.Warp(destNameOrDestDS=runout_binary_clipped_path,srcDSOrSrcDSTab=runout_binary_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}, dstSRS="EPSG:2056", outputBounds=[x_min,y_min,x_max,y_max]) as ds:
                ds=None

            runout_clipped_path = os.path.join(sim_folder_path,"runout_clipped.tif")
            with gdal.Warp(destNameOrDestDS=runout_clipped_path,srcDSOrSrcDSTab=runout_path, format="GTiff", creationOptions={"COMPRESS":"DEFLATE", "BIGTIFF":"NO"}, dstSRS="EPSG:2056", outputBounds=[x_min,y_min,x_max,y_max]) as ds:
                ds=None

            # Read the clipped rasters
            slide_mask_values = gdal_array.LoadFile(slide_mask_clipped_path).flatten()
            union_mask_values = gdal_array.LoadFile(union_mask_clipped_path).flatten()
            runout_values = gdal_array.LoadFile(runout_clipped_path).flatten()
            runout_values = runout_values/100.0 # Convert percent to decimal 
            runout_binary_values = gdal_array.LoadFile(runout_binary_clipped_path).flatten()

            # Determine pixel counts
            n = len(slide_mask_values) # Pixels in extent
            n_runout = np.sum(runout_binary_values) # Pixels in runout
            n_slide = np.sum(slide_mask_values) # Pixels in observed slide
            n_mask = np.sum(union_mask_values) # Pixels in UNION mask

            # Calculate confusion matrix values within UNION mask
            tp = np.sum((slide_mask_values[union_mask_values==1]==1)&(runout_binary_values[union_mask_values==1]==1))
            tn = np.sum((slide_mask_values[union_mask_values==1]==0)&(runout_binary_values[union_mask_values==1]==0))
            fp = np.sum((slide_mask_values[union_mask_values==1]==0)&(runout_binary_values[union_mask_values==1]==1))
            fn = np.sum((slide_mask_values[union_mask_values==1]==1)&(runout_binary_values[union_mask_values==1]==0))

            # Calculate AUC
            auc = roc_auc_score(slide_mask_values,runout_binary_values)
            auc_prob = roc_auc_score(slide_mask_values,runout_values)

            # Read runout polygon as basis for result record
            runout_binary_gpkg = runout_binary_path.replace(".tif",".gpkg")
            runout_layer = "runout"
            runout_mask_df = gpd.read_file(runout_binary_gpkg,layer=runout_layer,engine="pyogrio")

            # Add attributes to record
            runout_mask_df = runout_mask_df.dissolve("DN")
            runout_mask_df["sim_path"] = sim_folder_path
            runout_mask_df["runout_path"] = runout_binary_gpkg
            runout_mask_df["slide_id"] = slide_id
            runout_mask_df["slide_type"] = slide_type
            runout_mask_df["forest"] = forest
            runout_mask_df["n"] = n
            runout_mask_df["n_slide"] = n_slide
            runout_mask_df["n_runout"] = n_runout
            runout_mask_df["n_mask"] = n_mask
            runout_mask_df["tp"] = tp
            runout_mask_df["tn"] = tn
            runout_mask_df["fp"] = fp
            runout_mask_df["fn"] = fn
            runout_mask_df["auc"] = auc
            runout_mask_df["auc_prob"] = auc_prob
            runout_mask_df["slide_type"] = slide_type
            runout_mask_df["slide_type"] = slide_type
            runout_mask_df["slide_type"] = slide_type
            runout_mask_df["slide_type"] = slide_type
            runout_mask_df["slide_type"] = slide_type

            runout_mask_df["slide_type_label"] = runout_mask_df["slide_type"].apply(lambda x: slide_type_labels[x])
            runout_mask_df["forest_label"] = runout_mask_df["forest"].apply(lambda x: forest_labels[x])

            # Calculate and add performance metrics
            runout_mask_df["ppv"] = runout_mask_df["tp"] / (runout_mask_df["tp"]+runout_mask_df["fp"])
            runout_mask_df["npv"] = runout_mask_df["tn"] / (runout_mask_df["tn"]+runout_mask_df["fn"])
            runout_mask_df["tpr"] = runout_mask_df["tp"] / (runout_mask_df["tp"]+runout_mask_df["fn"])
            runout_mask_df["fpr"] = runout_mask_df["fp"] / (runout_mask_df["fp"]+runout_mask_df["tn"])
            runout_mask_df["csi"] = runout_mask_df["tp"] / (runout_mask_df["tp"]+runout_mask_df["fn"]+runout_mask_df["fp"])
            runout_mask_df["tnr"] = runout_mask_df["tn"] / (runout_mask_df["tn"]+runout_mask_df["fp"])
            runout_mask_df["acc"] = (runout_mask_df["tp"] + runout_mask_df["tn"]) / (runout_mask_df["tp"]+runout_mask_df["fp"]+runout_mask_df["fn"]+runout_mask_df["tn"])
            runout_mask_df["for"] = runout_mask_df["fn"] / (runout_mask_df["tn"]+runout_mask_df["fn"])
            runout_mask_df["fdr"] = runout_mask_df["fp"] / (runout_mask_df["tp"]+runout_mask_df["fp"])
            runout_mask_df["rj"] = np.sqrt((1-runout_mask_df["tpr"])**2+runout_mask_df["fpr"]**2)

            sim_poly_records.append(runout_mask_df)

            # Create record for CSV 
            sim_records.append({"slide_id":slide_id,"slide_type":slide_type,"forest":forest,"n":n,"n_slide":n_slide,"n_runout":n_runout
                                ,"tp":tp,"tn":tn,"fp":fp,"fn":fn, "auc":auc, "auc_prob":auc_prob
                            })

    # Save CSV with results per slide and simulation
    sim_stat_path = os.path.join(input_base_path,"sim_stats_{0}.csv".format(runout_threshold))
    utilities.delete_file_by_path(sim_stat_path)

    sim_stats_df = pd.DataFrame(sim_records)
    sim_stats_df.to_csv(sim_stat_path,sep=",",index=False)

    # Save Geopackage with runout polygons containing results as attributes
    sim_stat_gpkg = os.path.join(input_base_path,"sim_stats_{0}.gpkg".format(runout_threshold))
    sim_stat_layer = "sim_stats"
    utilities.delete_file_by_path(sim_stat_gpkg)

    sim_stats_gdf = rdf = gpd.GeoDataFrame(pd.concat(sim_poly_records, ignore_index=True), crs=sim_poly_records[0].crs)
    sim_stats_gdf.to_file(sim_stat_gpkg, layer=sim_stat_layer, driver="GPKG")

    # Save Geopackage with slide mask polygons
    slide_mask_gpkg = os.path.join(input_base_path,"slide_mask_{0}.gpkg".format(runout_threshold))
    slide_mask_layer = "slide_mask"
    utilities.delete_file_by_path(slide_mask_gpkg)

    slide_mask_gdf = rdf = gpd.GeoDataFrame(pd.concat(slide_poly_records, ignore_index=True), crs=slide_poly_records[0].crs)
    slide_mask_gdf.to_file(slide_mask_gpkg, layer=slide_mask_layer, driver="GPKG")

    # Save Geopackage with UNION mask polygons
    union_mask_gpkg = os.path.join(input_base_path,"union_mask_{0}.gpkg".format(runout_threshold))
    union_mask_layer = "union_mask"
    utilities.delete_file_by_path(union_mask_gpkg)

    union_mask_gdf = rdf = gpd.GeoDataFrame(pd.concat(mask_poly_records, ignore_index=True), crs=mask_poly_records[0].crs)
    union_mask_gdf.to_file(union_mask_gpkg, layer=union_mask_layer, driver="GPKG")


# Default entry point
if __name__ == "__main__":

    #
    # Configs
    #

    input_base_path = "E:/GIS_Projekte/Paper_4/data/simulation_inputs"
    raster_resolution = 5

    forest_labels = {
        -1:"Without forest (open field)",
        0:"Without forest (all)",
        1:"With deposit",
        2:"With deposit, friction fix",
        3:"With deposit, friction variable",
        4:"With friction fix",
        5:"With friction variable",
    }

    slide_type_labels = {
        0:"Loose material/high water content",
        1:"Loose material/normal water content", 
        2:"Adhesive material/high water content",
        3:"Adhesive material/normal water content"
    }

    thresholds = [
        # 100th percentile of reach probability
        {"runout_threshold": 0, "runout_threshold_comparator": ">"},
        # 90th percentile of reach probability
        {"runout_threshold": 10, "runout_threshold_comparator": ">="},
        # 75th percentile of reach probability
        {"runout_threshold": 25, "runout_threshold_comparator": ">="}
    ]

    #
    # Processing
    #

    # Process for all defined thresholds
    for threshold in thresholds:
        runout_threshold = threshold["runout_threshold"]
        runout_threshold_comparator = runout_threshold_comparator["runout_threshold"]
        process_result(input_base_path, runout_threshold, runout_threshold_comparator, raster_resolution, forest_labels, slide_type_labels)





 
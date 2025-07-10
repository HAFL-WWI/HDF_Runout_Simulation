# Data
Data for the analysis
- __sensitivity_constant_slope__: Subfolder containing the parameter combinations and result data for the sensitivity analysis on a constant slope. 
- __sensitivity_forest__: Subfolder containing the parameter combinations and result data for the sensitivity analysis for the historical landslides in the forest.
- __simulation_inputs__: Subfolder containing the result data for the simulations.
- __Dataset_Runout_Merged.gpkg__: Geopackage containing the reference data with the failure and runout areas for the historical landslides. The data is stored in the "inputs_merged" layer:

|**Attribute**  | **Datatype** | **Description**| 
|:---|:---|:---|
|fid   | integer | Feature ID |
|slide_id   | integer | ID of the individual landslides |
|area_type   | integer | Type of polygon: 1 = failure area, 2 = runout area |
|STORMENR   | text | StorMe inventory number (where present) |
|year   | integer | Year the landslide occurred |
|canton   | text | Abbreviation of the canton where the landslide is located |
|in_forest   | integer | Location of landslide occurrence: 0 = open field, 1 = in forest |
|volume_ori   | double | Displaced volume of the landslide in m3 according to the source inventory |
|failure_area_ori   | double | Failure area the landslide in m2 according to the source inventory |
|thickness_ori   | double | Failure thickness the landslide in m according to the source inventory |
|stemdensity_per_ha   | double | Stem density in trees per hectare based on estimates in the field |
|dbh_mean_cm   | double | Mean diameter at breast height in cm based on estimates in the field |
|thickness_m   | double | Failure thickness in m. Where a volume was present in the source, it is calculated from the inventory volume and the area of the failure area polygon. Otherwise the thickness from the source is copied. Only set on failure area. |
|area_m2   | double | Area of the polygon in m2. |



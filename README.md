Name: Sara Cornejo
GIS 540 Geospatial Programming, Spring 2025 – Final Project

I.
Title: Identifying urban hot spots via land surface temperature

II.
Abstract:
Urban growth and climate change have significantly impacted surface temperatures, especially in metropolitan areas where built infrastructure is more consolidated, and population density is higher, factors that increase the effect on public health and environmental sustainability. There are several ways to assess the different vulnerabilities of cities, and one of the most important is analyzing urban heat islands, a phenomenon caused by replacing natural areas with heat-absorbing materials such as concrete, pavement, and buildings. The San Salvador Metropolitan Area, located in the tropical country of El Salvador, is a growing city that experiences two main seasons: a wet season from May to September and a dry season from October to April. The dry season is characterized by high temperatures, and in recent years, the frequency of heat waves, especially from March to May, has increased significantly. This project aims to automate the processing of publicly available Landsat 8 imagery to calculate land surface temperature (LST) in the San Salvador Metropolitan Area and analyze temperature trends over time. The maps produced will visualize these trends and can later be analyzed in relation to urban features. The results will help identify areas most affected by rising temperatures and serve as a support tool for urban planning, particularly for allocating green infrastructure, parks, and other public spaces in the most vulnerable locations to reduce urban heat.

III.
Input Data:
•Dataset1: Landsat 8 Satellite Images
▪Data format: .tar file when compressed: .TIF files for individual bands after extraction
▪Attributes used: Band 10 (Thermal Infrared), Band 4 (Red), Band 5 (Near Infrared), acquisition date and constants (from metadata).
▪This dataset will be used to calculate land surface temperature and analyze temperature variations across different years.
▪Data source: https://earthexplorer.usgs.gov/
Figure 1. Data set 1, bands needed for calculations.
Figure 2. Metadata text file provided when downloading from EarthExplorer. It contains the constants required for the calculations of each scene.
•Dataset2: Metropolitan Area Boundary
▪Data format: Polygon shapefile.
▪Attributes used: Geometry (area and perimeter); no specific attribute fields used.
▪This shapefile is used to clip the raster data to the San Salvador Metropolitan Area, focusing the analysis on the urban extent.
▪Data source: OPAMSS (Urban Planning Office of the San Salvador Metropolitan Area)
Figure 3. Dataset 2, Metropolitan Area of San Salvador Boundary

IV.
Data Products:
1.Clipped LST Rasters
  oFormat: .tif
  oDescription: LST rasters clipped to the Metropolitan Area shapefile per Landsat scene processed
  oPurpose: Shows the calculated land surface temperature for each date focusing only inside your area of interest.
2.Average Clipped LST Raster
  oFormat: .tif
  oDescription: mean temperature calculated from all available scenes processed.
  oPurpose: Useful for understanding typical heat patterns across multiple years.
3.Normalized Difference Vegetation Index (NDVI) Rasters
  oFormat: .tif
  oDescription: NDVI maps for each date processed.
  oPurpose: Indicates vegetation health and density; can support further analyses like correlating vegetation cover with surface temperatures.
4.Map Exports
  oHigh-quality exported JPEG maps with symbology applied (LST_Map.jpg or user-defined name).
5.HTML Report
  oSimple webpage summarizing:
    ▪Scene dates processed.
    ▪Minimum, maximum, and average LST values per date and the average for all the scene dates.
    ▪Shows the jpg exported map for visualization.
    
V.
Workflow
GET input_path (user can provide single .tar file or folder of .tar files) GET shapefile for clipping LST rasters (user input or default) GET output filename (user input or default "LST_Map.jpg")
FUNC date_from_filename
SPLIT part from file name
IF file name has more than 3 parts:
GET the raw date
GET year
GET month
GET day
ELSE return unknown date
SET a dictionary called constants
FUNC read_metadata
WITH open metadata file to read
FOR each line in the text file
IF line starts with RADIANCE_MULT_BAND_10 =
SPLIT the words in the line and get the last final number in the line
SAVE that as ‘ml’ in the dictionary with the specific number
IF line starts with RADIANCE_ADD_BAND_10 =
SPLIT the words in the line and get the last final number in the line
SAVE that as ‘al’ in the dictionary with the specific number
IF line starts with K1_CONSTANT_BAND_10 =
SPLIT the words in the line and get the last final number in the line
SAVE that as ‘k1’ in the dictionary with the specific number
IF line starts with K2_CONSTANT_BAND_10 =
SPLIT the words in the line and get the last final number in the line
SAVE that as ‘k2’ in the dictionary with the specific number
ENDFOR
RETURN constants
ENDFUNC
FUNC calculate_toa_radiance
GET constant ml
GET constant al
SET value of OI = 0.29
TOA = (BAND 10 * ml) + al – oi
SAVE TOA
RETURN TOA
ENDFUNC
FUNC calculate_bt
BT = k2 / (Ln (k2/TOA +1))
SAVE BT
RETURN BT
ENDFUNC
FUNC calculate_ndvi
NDVI = (BAND5 – BAND 4) / (BAND5 + BAND 4)
SAVE NDVI
RETURN NDVI
ENDFUNC
FUNC calculate_prop_veg
GET min value of NDVI
GET max value of NDVI
PV = NDVI – min_val / (max_value – min value) ^2
SAVE PV
RETURN PV
ENDFUNC
FUNC calculate_emissivity
emissivity = 0.004 * PV + 0.986
SAVE emissivity
RETURN emissivity
ENDFUNC
FUNC calculate_lst
LST_k = BT / (1+(10.8 * BT/14388) * ln (emissitvity)
LST_c = LST – 273.15
SAVE LST_c
RETURN LST_c
FUNC
FUNC extract_tar
WITH tarfile extract the files
ENDFUNC
FOR each metadata .txt file in Data/: CALL FUNC read_metadata READ radiometric constants (ml, al, k1, k2) ENDFUNC
SET paths for Band 10, Band 4, Band 5
CALL FUNC calculate_toa
CALL FUNC calculate_bt
CALL FUNC calculate_ndvi
CALL FUNC calculate_prop_veg
CALL FUNC calculate_emissivity CALL FUNC calculate_lst
COLLECT all clipped LST rasters
IF clipped rasters exist: SUM all clipped rasters DIVIDE by number of rasters to get average LST SAVE average raster as "lst_clipped_avg.tif" ENDIF
CREATE MAP OUTPUT
Set APRX to current
addData from path to map layer
Set symbology field of layer
Choose symbology color
INITIATE empty list
FOR each file in Output/ folder:
IF file name starts with "lst_clipped_" AND ends with ".tif" AND is NOT "lst_clipped_avg.tif":
GET minimum, maximum, and mean temperature values
ADD these statistics in a dictionary
ENDIF
ENDFOR
CHECK if "lst_clipped_avg.tif" exists:
GET minimum, maximum, and mean temperature values
ADD these statistics in a dictionary labeled "Average"
ENDIF
CREATE HTML start content
FOR each statistics in the dictionary
CREATE a table row with Scene, Min, Max, and Mean temperature values
ENDFOR
CREATE HTML end content
ADD jpg image inside the HTML
COMBINE start + end into a full HTML page
SAVE HTML file

VI.
Keywords:
Land surface temperature, minimum temperature, maximum temperature, mean temperature, batch raster analysis.

VII.
How to run:
It is necessary to install and import the arcpy package to run the script, the script tool was developed with ArcGis Pro and the script takes the following parameters:
Parameter label Data type Default value (all parameters must have defaults)
Select a single .tar file or a folder with multiple .tar files
Folder
C:\GIS540_Project\Data\compressed
Area of interest
Feature Class
C:\GIS540_Project\Data\perimetro_AMSS.shp
Output file name
Long
“LST_Map.jpg”
Figure 4. Screenshot of the script tool name LST Calculation

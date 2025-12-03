# LST.py
#
# Author: Sara Cornejo 26.04.2025
# Unity ID: sgcornej

# Purpose: Automate the extraction, processing, and analysis of Land Surface Temperature (LST)
# from Landsat imagery downloaded from EarthExplorer including NDVI calculation and
# clipping to an area of interest.

# Procedure Summary:
#           --Accept one or more .tar files containing Landsat imagery, or default to a predefined compressed
#           directory if no input is provided.
#           --Extract the Landsat data, identify relevant bands (Band 10, Band 4, Band 5) and metadata.
#           --Calculate TOA Radiance, Brightness Temperature, NDVI, proportion of vegetation,
#           emissivity, and final LST in Celsius.
#           --Clip the output LST to a user defined shapefile or default (area of interest).
#           --Compute the average LST raster across all dates.
#
#
# Main Steps:
# Step 1: Get a .tar file or a folder containing .tar files and also a get a shapefile
# defining the area of interest for clipping.If any of those not provided, use defaults.
#
# Step 2: Check if input is a .tar file or directory. Extract the .tar files into a working folder.
#
# Step 3: Read the metadata .txt file to get radiometric constants needed for LST and BT calculations.
#
# Step 4: LST Processing Pipeline: TOA, Brightness Temperature, NDVI, PV, emissivity and LST in celsius.
#
# Step 5: Clip Output to Area of Interest
#
# Step 6: calculate the average LST if multiple years are processed for the clipped area, saved as lst_clipped_avg.tif
#
# Step 7: Mapping
#
# Step 8: Create an HTML report.
#
# Software Requirements:
#       - ArcGIS Pro arcpy package must be installed.
#       - ArcGIS Pro version 3.3.1
# 
# Example input:
#
#

#import needed modules
import arcpy, os, sys, tarfile


# ---- Defining functions -----


# Setting the functions. 
def date_from_filename(metadata_path):
    """Extracts the acquisition date (YYYY_MM_DD) from the Landsat metadata filename.
    Example output: '2021_02_20' """
    filename = os.path.basename(metadata_path)
    parts = filename.split('_')
    if len(parts) > 3:
        raw_date = parts[3]  
        year = raw_date[0:4]
        month = raw_date[4:6]
        day = raw_date[6:8]
        return f"{year}_{month}_{day}"
    else:
        return "unknown_date"

# Making a dictionary to get the values from a .txt file which is the metadata. 
def read_metadata(metadata_path):
    """Reads the textfile metadata from Lansdat and extracts ml, al, K1, K2 constants required for calculations ahead"""

    constants = {}
    with open(metadata_path, 'r') as file:
        for line in file:
            if line.startswith("    RADIANCE_MULT_BAND_10 ="):
                words = line.split()
                constants['ml'] = float(words[-1])
            if line.startswith("    RADIANCE_ADD_BAND_10 ="):
                words = line.split()
                constants['al'] = float(words[-1])
            if line.startswith("    K1_CONSTANT_BAND_10 ="):
                words = line.split()
                constants['k1'] = float(words[-1])
            if line.startswith("    K2_CONSTANT_BAND_10 ="):
                words = line.split()
                constants['k2'] = float(words[-1])
    return constants

# Calculating Top of Atmospheric (TOA) spectral radiance
def calculate_toa_radiance(band10, ml, al, output_path):
    """Calculates TOA Radiance from thermal band"""
    ml = constants['ml']
    al = constants['al']
    oi = 0.29
    toa = (arcpy.sa.Float(band10) * ml) + al - oi
    toa.save(output_path)
    return arcpy.Raster(output_path)

# Calculating top of atmosphere Brightness-temperature
def calculate_bt(toa_raster, k1, k2, output_path):
    """Calculates Brightness Temperature."""
    bt = k2 / arcpy.sa.Ln((k1/arcpy.sa.Float(toa_raster))+1)
    bt.save(output_path)
    return arcpy.Raster(output_path)

# Calculating NDVI
def calculate_ndvi(BAND5, BAND4, output_path):
    """Calculates NDVI from NIR (band5) and Red (band4) bands."""
    nir_BAND5 = arcpy.Raster(BAND5)
    R_BAND4 = arcpy.Raster(BAND4)
    ndvi = (nir_BAND5-R_BAND4)/(nir_BAND5+R_BAND4)
    ndvi.save(output_path)
    return arcpy.Raster(output_path)
    
# Calculating the propotion of vegetation
def calculate_prop_veg(ndvi_raster, output_path):
    """Calculates proportion of vegetation (PV) from NDVI."""
    min_value = float(arcpy.GetRasterProperties_management(ndvi, "MINIMUM").getOutput(0))
    max_value = float(arcpy.GetRasterProperties_management(ndvi, "MAXIMUM").getOutput(0))
    pv = ((ndvi_raster - min_value) / (max_value - min_value)) ** 2
    pv.save(output_path)
    return arcpy.Raster(output_path)

#Calcutating the emissivity
def calculate_emissivity(pv_raster, output_path):
    """Calculates land surface emissivity."""
    emissivity = 0.004 * pv_raster + 0.986
    emissivity.save(output_path)
    return arcpy.Raster(output_path)

# Calculating LST
def calculate_lst(bt_raster, emissivity_raster, output_path):
    """Calculates Land Surface Temperature in Celsius."""
    lst_k = bt_raster / (1 + (10.8 * bt_raster / 14388) * arcpy.sa.Ln(emissivity))
    lst_c = lst_k - 273.15
    lst_c.save(output_path)
    return arcpy.Raster(output_path)


# Note, tarfile module must be imported. Thanks for code from U13-Forward for tar extraction to an specific path at
# https://stackoverflow.com/questions/59875978/tarfile-open-does-not-extract-into-the-right-directory-path/59876047
def extract_tar(tar_path, extract_to):
    """Extract a .tar file to the given path."""
    print(f"Extracting: {tar_path}")
    with tarfile.open(tar_path) as tar:
        tar.extractall(path=extract_to)
    print(f"Extracted to: {extract_to}")

# ---- Setting the workspace and directories -----

# Setting geoprocessing environment. 
arcpy.CheckOutExtension('Spatial')
arcpy.env.overwriteOutput = True

# Setting the directory based on the script's relative location. 
scriptPath = sys.argv[0]
codePath = os.path.dirname(scriptPath)
baseDir = os.path.dirname(codePath) + "/"
dataDir = baseDir + "Data/"

# Prepare the output workspace.
outputDir = baseDir + "Output/"

# Specify the template project.
projectPath = codePath + "/LST_calculation/LST_calculation.aprx"
layout_file = codePath + "/LST_calculation/LST_Layout.pagx"
# ---- Getting users arguments -----

# Step 1: Getting users arguments and setting default paths 
try:
    input_path = sys.argv[1]  # user can provide a single .tar file or a folder with multiple .tar files
except IndexError:
    input_path = dataDir + "compressed/"  # default 

try:
    clip_shp = sys.argv[2]
    output_filename = sys.argv[3]
    
    if output_filename.lower().endswith(".jpg"):
        output_map = outputDir + output_filename
    else:
        output_map = outputDir + output_filename + ".jpg"
    
except IndexError:
    clip_shp = dataDir + 'perimetro_AMSS.shp' # default
    output_map = outputDir + "LST_Map.jpg"

# ----  Starting data processing -----

# Step 2: Extracting .tar files from a single file of a folder with multiple .tar files
if os.path.isfile(input_path) and input_path.lower().endswith(".tar"):
    extract_tar(input_path, dataDir)

elif os.path.isdir(input_path):
    for file in os.listdir(input_path):
        if file.lower().endswith(".tar"):
            extract_tar(os.path.join(input_path, file), dataDir)
else:
    print(f"Invalid input: {input_path}")
    arcpy.AddMessage(f"Invalid input: {input_path}")


# Step 3. Read the metadata .txt file to get radiometric constants and identifying 
# .tif files (BAND 10, BAND 4  & BAND 5) needed for LST and BT calculations.
for filename in os.listdir(dataDir):
    if filename.endswith("_MTL.txt"):
        metadata_path = os.path.join(dataDir, filename)
        # Get the date from the metadata filename
        scene_date = date_from_filename(metadata_path)
        base = filename.split("_MTL")[0]
        
        band10_path = os.path.join(dataDir, base + "_B10.TIF")
        BAND4 = os.path.join(dataDir, base + "_B4.TIF")
        BAND5 = os.path.join(dataDir, base + "_B5.TIF")

        print(f"Processing scene: {scene_date}")
        """Prints the .txt in arcGIS"""
        arcpy.AddMessage(f'Processing scene: {scene_date}')
    
        # Read metadata constants
        constants = read_metadata(metadata_path)

        # Creating a date-based output filenames, first creating a list of the steps, then creating
        # a dictionary to hold the output paths and loop through and build each path 

        step_names = ["toa", "bt", "ndvi", "pv", "emissivity", "lst", "lst_classified"]
        output_paths = {}
        for step in step_names:
            output_paths[step] = outputDir + f"{step}_{scene_date}.tif"

# Step 4: LST Processing Pipeline

        # Calling the functions and pass the dynamic output paths
        toa = calculate_toa_radiance(band10_path, constants['ml'], constants['al'], output_paths['toa'])
        bt = calculate_bt(toa, constants['k1'], constants['k2'], output_paths['bt'])
        ndvi = calculate_ndvi(BAND5, BAND4, output_paths['ndvi'])
        pv = calculate_prop_veg(ndvi, output_paths['pv'])
        emissivity = calculate_emissivity(pv, output_paths['emissivity'])
        lst = calculate_lst(bt, emissivity, output_paths['lst'])

# Step 5: Clip Output to Area of Interest
        # Clip the LST rasters using the San Salvador metropolitan area perimeter shapefile as default
        # if not provided by the user
        clipped_lst_path = outputDir + f"lst_clipped_{scene_date}.tif"
        arcpy.management.Clip(
        in_raster=output_paths['lst'],
        rectangle="#",
        out_raster=clipped_lst_path,
        in_template_dataset=clip_shp,
        nodata_value="#",
        clipping_geometry="ClippingGeometry",
        maintain_clipping_extent="NO_MAINTAIN_EXTENT")       

        print(f"Finished processing {scene_date}")
        arcpy.AddMessage(f'Processing scene: Finished processing {scene_date}')
    
# Step 6: Average LST
# Average LST Calculation 
clipped_lst_rasters = []
for file in os.listdir(outputDir):
    if file.startswith("lst_clipped_") and file.endswith(".tif"):
        raster_path = os.path.join(outputDir, file)
        clipped_lst_rasters.append(raster_path)

if len(clipped_lst_rasters) == 0:
    print("No clipped LST rasters found.")
    arcpy.AddMessage("No clipped LST rasters found.")
else:
    print(f"Found {len(clipped_lst_rasters)} rasters. Calculating average...")
    arcpy.AddMessage(f"Found {len(clipped_lst_rasters)} rasters. Calculating average...")

    # Convert list to Raster objects and sum them
    raster_sum = arcpy.Raster(clipped_lst_rasters[0])
    i = 1
    while i < len(clipped_lst_rasters):
        raster_i = arcpy.Raster(clipped_lst_rasters[i])
        raster_sum = raster_sum + raster_i
        i = i + 1
    average_raster = raster_sum / len(clipped_lst_rasters)

    avg_output_path = outputDir + "lst_clipped_avg.tif"
    average_raster.save(avg_output_path)

    print("Average LST raster saved to:", avg_output_path)
    arcpy.AddMessage(f"Average LST raster saved to: {avg_output_path}")

# ----  Mapping -----
# Step 7: Mapping
try:
    # When working inside ArcGIS Pro,
    # use the project name "CURRENT"
    aprx = arcpy.mp.ArcGISProject("CURRENT")

except OSError:
    # When working outside ArcGIS Pro,
    # use the project full path file name.
    aprx = arcpy.mp.ArcGISProject(projectPath)


myMap = aprx.listMaps()[0]
layout = aprx.listLayouts("LST_layout")[0]

color_ramp = aprx.listColorRamps("Slope")[0]
avg_raster_path = outputDir + "lst_clipped_avg.tif"
if os.path.exists(avg_raster_path):
    layer = myMap.addDataFromPath(avg_raster_path)
    sym = layer.symbology
    sym.colorizer.classificationMethod = 'NaturalBreaks'
    sym.colorizer.breakCount = 7
    sym.colorizer.colorRamp = color_ramp
    layer.symbology = sym

    layout.exportToJPEG(output_map, resolution=300)
    print(f"Exported layout to: {output_map}")
    arcpy.AddMessage(f"Exported layout to: {output_map}")
   
else:
    print("Raster not found at expected path.")
    arcpy.AddMessage("Raster not found at expected path.")
    
# Step 8: Create an HTML report.

raster_stats = []

for file in os.listdir(outputDir):
    if file.startswith("lst_clipped_") and file.endswith(".tif") and file != "lst_clipped_avg.tif":
        raster_path = os.path.join(outputDir, file)
        raster = arcpy.Raster(raster_path)
        stats = {
            'scene': file.replace("lst_clipped_", "").replace(".tif", ""),
            'min': raster.minimum,
            'max': raster.maximum,
            'mean': raster.mean
        }
        raster_stats.append(stats)

avg_raster_path = os.path.join(outputDir, "lst_clipped_avg.tif")
if os.path.exists(avg_raster_path):
    avg_raster = arcpy.Raster(avg_raster_path)
    avg_stats = {
        'scene': "Average",
        'min': avg_raster.minimum,
        'max': avg_raster.maximum,
        'mean': avg_raster.mean
    }
    raster_stats.append(avg_stats)



html_start =f"""<html><head><title>Land Surface Temperature Report</title></head>
<body>
<center>
<h1>Land Surface Temperature Report</h1>
<table border='1'><tr><th>Scene</th><th>Min (°C)</th><th>Max (°C)</th><th>Mean (°C)</th></tr>"""

html_rows = []
for stat in raster_stats:
    row = f"""<tr><td>{stat['scene']}</td><td>{stat['min']:.2f}</td><td>{stat['max']:.2f}</td><td>{stat['mean']:.2f}</td></tr>"""
    html_rows.append(row)
    
html_end = f"""</table>
<br>
<img src="LST_map.jpg" alt="LST Map" width="600">
</center>
</body></html>"""

full_html = html_start + "".join(html_rows) + html_end

out_path = os.path.join(outputDir, 'lst_report.html')
outf = open(out_path, "w")
outf.write(full_html)
outf.close()
print(f"HTML report generated at:{out_path}")
arcpy.AddMessage(f"HTML report generated at: {out_path}")

        



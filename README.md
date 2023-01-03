# BowenIndex

This python tool enables calculation of Bowen Index and several other characteristics of evapotraspiration (albedo, land surface temperature, net radiation etc.)

## Requirements

* Python (tested on versions from 3.7 to 3.10)
* Python packages
    * Native: csv, json, math, os, sys, pprint, 
    * Installed: numpy, fiona, rasterio, tifffile, osgeo, sympy
 
## Input Data
* Landsat Images (tested on Landsat 8 and Landsat 9 images)
* Polygon mask in GeoPackage (same coordinate system as for the Landsat images)
* CSV with meterological data in structure (sensing date in YYYYMMDD)

## Output data

The script will generate a series of folders: *clipped_bands*, *albedo*, *bowenindex*, *flux*, *lst*, *netBudget*, *vegHeigth*, *vegIndices* containing respective variable for each day. To be sure, the name of variable is written at the end of original name. 

## Input Structure

The tool offers automatic search in folders. For the script to work properly you need to unzip the downloaded images and place all the unzipped folder into one folder with original name. Keep the JSON metadata. 

Make sure that the dates in the CSV file are in the same format as in the name of your downloaded images (**meaning they are in YYYYMMDD format**).

## Note

As this is an on-going project, the full optimlization has yet to be done. 

Research of this project has been supported by grant DSGC-2021-0083 *Modelling of spatiotemporal variability of selected bioclimatic factors of the landscape using open data in GIS* within framework of the project  *Improving schematics of Doctoral student grant competition and their pilot implementation CZ.02.2.69/0.0/0.0/19_073/0016713* ensured by Palacky University in Olomouc, Czech Republic.

import os
import fiona
import rasterio
import rasterio.mask
from pprint import pprint
import numpy as np
from sympy import Le
import tifffile as tf
import csv
from osgeo import gdal, osr
import math
import json



######################################################################
#   GENERAL FUNCTIONS
######################################################################

def getfilepath(input_folder, suffix):
    pathListFolder = []
    for root, dirs, files in os.walk(input_folder, topdown=False):
        for name in files:
            if name.endswith(suffix):
                if name not in pathListFolder:
                    pathListFolder.append(os.path.join(root, name))
                #'C:\\Users\\Tereza\\Documents\\PhD\\zkousky\\2_vyvoj_pro_FOSS\\snimky\\LC08_L1TP_190025_20170528_20200903_02_T1\\LC08_L1TP_190025_20170528_20200903_02_T1_B5.TIF'
    return pathListFolder

######################################################################

def createmeteodict(csv_file):
    with open(csv_file, mode = 'r') as infile:
        csv_list = [[val.strip() for val in r.split(",")] for r in infile.readlines()] #[['date', 'avg_temp', 'wind_sp', 'relHum', 'max_temp', 'min_temp'],
                                                                                       #['20220518', '14.25', '1.25', '65.52', '19.8', '8.7'],
                                                                                       #['20220612', '22.60', '0.83', '70.21', '28.4', '13.3']]

        (_, *header), *data = csv_list                                              #(('date', 'avg_temp', 'wind_sp', 'relHum', 'max_temp', 'min_temp'),
                                                                                    #['20220518', '14.25', '1.25', '65.52', '19.8', '8.7'],
                                                                                    #['20220612', '22.60', '0.83', '70.21', '28.4', '13.3'])
        csv_dict = {}
        for row in data:
            key, *values = row                                                      # ('20220612', '22.60', '0.83', '70.21', '28.4', '13.3')
            csv_dict[key] = {key : float(value) for key, value in zip(header, values)}     #{'date': [avg_temp': '22.60','max_temp': '28.4','min_temp': '13.3','relHum': '70.21','wind_sp': '0.83']}
    return  csv_dict

######################################################################

def load_json(jsonFile):
    with open(jsonFile, 'r') as j:
        data = json.load(j)
    return data

######################################################################

def clipimage(maskPath, inputBand, outImgPath):
    with fiona.open(maskPath, "r") as gpkg:
        shapes = [feature["geometry"] for feature in gpkg]

    with rasterio.open(inputBand) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta
    
    out_meta.update({"driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform})

    with rasterio.open(outImgPath, "w", **out_meta) as dest:
        dest.write(out_image)
    return 

######################################################################

def savetif(img, outputPath, epsg = 32633):
    new_dataset = rasterio.open(outputPath, "w", 
        driver = "GTiff",
        height = img.shape[0],
        width = img.shape[1],
        count = 1,
        nodata = -9999,
        dtype = img.dtype,
        crs = epsg,
        # upper left, pixel size
        transform = rasterio.transform.from_origin(656265, 5503485, 30, 30))
    new_dataset.write(img, 1)
    new_dataset.close()
    
    return

######################################################################
#   ATMOSHEPRIC FUNCTIONS
######################################################################
def e0(T):
    """
    # Partial Water Vapour Pressure [kPa]
    # Tetens Formula via https://www.omnicalculator.com/chemistry/vapour-pressure-of-water

    T = Air Temperature [˚C]
    """
    e0 = 0.6108 * math.exp(((17.27 * T) / (T + 237.3)))
    return e0

######################################################################

def es(Tamax, Tamin):
    """
    # Saturated Vapour Pressure [kPa]
    # via https://www.fao.org/3/x0490e/x0490e07.htm#atmospheric%20pressure%20(p)
    
    Tamax = Maximum Air Temperature [˚C]
    Tamin = Minimal Air Temperature [˚C]
    """

    e0_max = e0(Tamax)
    e0_min = e0(Tamin)

    es = (e0_max + e0_min) / 2  
    return es

######################################################################

def atmPress(Z):
    """
    # Atmospheric Pressure [kPa]
    # via https://designbuilder.co.uk/helpv3.4/Content/Calculation_of_Air_Density.htm

    Z = Elevation above sea level [m]
    """
    
    P = (101325 * ((1.0 - (Z * 0.0000225577)) ** 5.2559)) / 1000
    return P

######################################################################


def	atmemis(es, Ta):
	#(Brutsaert 1982).
    """
    # Emmisivity of the Atmosphere
    # via Brutsart, 1982 (https://link.springer.com/book/10.1007/978-94-017-1497-6)

    Ta = Air Temperature [˚C]
    es = Saturated vapour pressure [kPa]
    """
    
    atmEmis = 1.24 * (es * 10.0/(Ta + 273.15)) ** (1.0/7.0)

    return atmEmis
    
######################################################################

def densityair(Z, Ta, R = 278.05):
    """
    Density of Air [kg/m3]
    # via https://designbuilder.co.uk/helpv3.4/Content/Calculation_of_Air_Density.htm
    
    R  = Gas Constant (287.05 J/kg-K)
    Ta = Air Temperature in K (T [˚C] + 273.15)
    P = Standard Pressure
    """
   
    P = (101325 * ((1.0 - (Z * 0.0000225577)) ** 5.2559)) 
    air_density = P / (R * (Ta + 273.15))
    
    return air_density

######################################################################
def psychroCons(P):
    """
    # Psychrometric constant [kPa/˚C]
    # via https://www.fao.org/3/x0490e/x0490e07.htm#psychrometric%20constant%20(g)

    P = Atmospheric Pressure [kPa]
    """
    y =  0.00065 * P
    return y
######################################################################

def tpw(RU, P):
    """
    # Total Precipitable Water [kg/m2]
    # via https://www.scielo.br/j/rbeaa/a/sX6cJjNXWMfHQ5p4h33B8Zz/?lang=en&format=pdf

    RU = Relative Humidity [%]
    P = Atmospheric Pressure [kPa]  
    """
    
    tpw = 0.14 * RU * P + 2.1
    return tpw

######################################################################
# https://www.scielo.br/j/rbeaa/a/sX6cJjNXWMfHQ5p4h33B8Zz/?lang=en&format=pdf  
def toc(P, Z, tpw, Kt = 1.0): #Kt from 1.0 for clear air, 0.5 for extremely polluted air
    """
    # Atmospheric transmittance in the solar radiation domain
    # via https://www.scielo.br/j/rbeaa/a/sX6cJjNXWMfHQ5p4h33B8Zz/?lang=en&format=pdf
    # Unitless

    P = Atmospheric Pressure[kPa]
    Z = Elevation above sea level [m]
    tpw = Total Precipitable Water [kg/m2]
    Kt = Air Turbidity coefficient (Kt= 1.0 for clear air and Kt= 0.5 for extremely turbid or polluted air)
    """

    Toc = 0.35 + 0.627 * math.exp((-(0.0014 * P) / (Kt * math.cos(Z))) -0.075 * (tpw / math.cos(Z) ** 0.4))
    return Toc


def ndvi(nir, red, outputPath):
    """
    # Normalized Differential Vegetation Index
    # NDVI = (NIR - RED) / (NIR + RED)
    # Unitless
    # Range from -1 to +1
    """
    nir = np.array(tf.imread(nir))
    red = np.array(tf.imread(red))
    zero_except = np.seterr(all = "ignore")
    NDVI = (nir - red) / (nir + red)
    NDVI[NDVI > 1] = 0.99
    NDVI[NDVI < -1] = -0.99


    savetif(NDVI, outputPath)

    return NDVI

    ######################################################################

def bt(K1, K2, RadAddBand, RadMultBand, band):
    """
    # Top of Atmosphere Brightness Temperature [K]
    # via Landsat 8 Data Users Handbook (https://www.usgs.gov/media/files/landsat-8-data-users-handbook)

    K1 = Band-specific thermal conversion constant from the metadata (K1_CONSTANT_BAND_x, where x is the thermal band number)
    K2 = Band-specific thermal conversion constant from the metadata (K2_CONSTANT_BAND_x, where x is the thermal band number)
    RadAddBand = Radiation Additive Term from the Metadata (RADIANCE_ADD_BAND_X)
    RadMultBand =  Radiation Multiplicative Term from the Metadata (RADIANCE_MULT_BAND_X)
    band = Landsat Thermal Band
    TOARad = Top of Atmosphere Radiation
    """
    
    band = np.array(tf.imread(band))
    TOARad = RadMultBand * band + RadAddBand

    BT = (K2 / np.log(K1/TOARad + 1)) - 273.15
    BT[BT == -273.15] = 0
    return BT

#vegetation fraction (proportion of vegetation)
def pv(NDVI):
    """
    # Fraction of Vegetation
    # via https://giscrack.com/how-to-calculate-land-surface-temperature-with-landsat-8-images/

    NDVI = Normal Differential Vegetation Index
    """
    zero_except = np.seterr(all = "ignore")
    #PV = np.sqrt((NDVI - np.min(NDVI))/ (np.max(NDVI) - np.min(NDVI)))
    PV = np.divide(np.power(NDVI, 2), 0.3) 
    
    PV[PV > 1] = 0.99
    PV[PV < 0] = 0.1

    return PV

def emis(red, ndvi, Pv): #surface emis
    """
    # Surface Emmisivity
    # via https://giscrack.com/how-to-calculate-land-surface-temperature-with-landsat-8-images/

    red = Landsat Red Band
    ndvi = Normal Differential Vegetation Index
    Pv = Fraction of Vegetation
    """
    red = np.array(tf.imread(red))
    
    E = (0.004 * Pv) + 0.986
    E = np.where(ndvi < 0.2, 1 - red, E) 
    E = np.where(ndvi > 0.5, 0.99, E)
    E[E > 1] = 0.99
    E[E < 0.8] = 0.8
    return E


def LST(BT, emis, band, outputPath):
    """
    # Land Surface Temperature [˚C]
    # via https://giscrack.com/how-to-calculate-land-surface-temperature-with-landsat-8-images/

    BT - Brightness temperature [K]
    emis - emissivity [-]
    band - Landsat band [-]
    """
    band = np.array(tf.imread(band))
    
    LST = (BT / (1 + ((0.0015 * BT)/1.4488) * np.log(emis)))
    
    savetif(LST, outputPath)
    return LST

def rb_band(AddRef, MultRef, band, SE, outputPath):
   
    """
    # Band Reflectance 
    # via https://www.scielo.br/j/rbeaa/a/sX6cJjNXWMfHQ5p4h33B8Zz/?lang=en&format=pdf

    AddRef = Reflectance Additive Term from the Metadata (REFLECTANCE_ADD_BAND_X)
    MultRef =  Reflectance Multiplicative Term from the Metadata (REFLECTANCE_MULT_BAND_X)
    band = Landsat band
    SE = Sun Elevation Angle [˚]
    """
    band = np.array(tf.imread(band))
    
    toa_unitless = MultRef * band + AddRef
    rb = (toa_unitless / math.cos(SE))

    savetif(rb, outputPath)
    return rb


def kb(rb, Z, AddRad, MultRad, band, dES, outputPath):
    """
    # Solar Constant [W/m2]
    # via https://www.scielo.br/j/rbeaa/a/sX6cJjNXWMfHQ5p4h33B8Zz/?lang=en&format=pdf

    rb = Band reflectance
    Z = Elevation above sea level [m]
    AddRad = Radiance Additive Term from the Metadata (RADIANCE_ADD_BAND_X)
    MultRad = Radiance Multiplicative Term from the Metadata (RADIANCE_MULT_BAND_X)
    band = Landsat band
    dES = Earth-Sun Dinstance from Metadata ('EARTH_SUN_DISTANCE') [AU]

    Lb = Radiance of each pixel [W/m2]
    dr = correction of the eccentricity of the terrestrial orbit

    """
    
    band = np.array(tf.imread(band))
    
    Lb = AddRad + MultRad * band
    dr = (1/dES) ** 2 
    kb = math.pi * Lb / rb * math.cos(Z * dr)

    savetif(kb, outputPath) 
    return kb

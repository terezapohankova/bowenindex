import csv
import json
import math
import os
from pprint import pprint

import fiona
import numpy as np
import rasterio
import rasterio.mask
import tifffile as tf
from osgeo import gdal, osr
from sympy import Le, ln

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

"""def pb(kb_numerator, kb2, kb3, kb4, kb5, kb6, kb7, outputPath):
    
    # Albedo Weight 
    # via https://www.scielo.br/j/rbeaa/a/sX6cJjNXWMfHQ5p4h33B8Zz/?lang=en&format=pdf


    
    #kb_numerator = np.array(tf.imread(kb_numerator))
    #kb2 = np.array(tf.imread(kb2))
    #kb3 = np.array(tf.imread(kb3))
   # kb4 = np.array(tf.imread(kb4))
    #kb5 = np.array(tf.imread(kb5))
    #kb6 = np.array(tf.imread(kb6))
    #kb7 = np.array(tf.imread(kb7))

    pb = kb_numerator / (kb2 + kb3 + kb4 + kb5 + kb6 + kb7)

    savetif(pb, outputPath) 

    return pb"""

def toaplanetary(pb2,rb2, pb3, rb3, pb4, rb4,  pb5, rb5, pb6,rb6, pb7, rb7, outputPath):
    """band2 = np.array(tf.imread(band2))
    band3 = np.array(tf.imread(band3))
    band4 = np.array(tf.imread(band4))
    band5 = np.array(tf.imread(band5))
    band6 = np.array(tf.imread(band6))
    band7 = np.array(tf.imread(band7))"""


    #rb = np.array(tf.imread(rb))
    
    toaPlanet = pb2 * rb2 + pb3 * rb3 + pb4 * rb4 + pb5 * rb5 + pb6 * rb6 + pb7 * rb7
    savetif(toaPlanet, outputPath)
    
    return toaPlanet

def albedo(toaplanet, Toc, outputPath):
    """
    # Albedo 
    # Unitless or %
    # Range from 0 to 1 (ot 0 % to 100 %)

    toaplanet = Planetary Top Of Atmosphere Radiance  
    Toc = Atmospheric transmittance in the solar radiation domain
    """
    albedo = (toaplanet - 0.03) / (Toc ** 2)

    albedo[albedo < 0] = 0.1
    
    albedo[albedo > 1] = 0.99

    savetif(albedo, outputPath)
    return albedo


def longout(emisSurf, LST, outputPath):
    """
    # Outgoing Longwave Radiation [W/m2]
    # via Stephan-Boltzmann law https://doi.org/10.1016/j.jrmge.2016.10.004

    emisSurf = emissivity of surface [-]
    LST = Land Surface Temprature [˚C]
    """
    LST = LST + 273.15
    longOut = emisSurf * 5.6703 * 10.0 ** (-8.0) * LST ** 4
    savetif(longOut, outputPath)
    return longOut
######################################################################

def longin(emisAtm, LST, outputPath):
    """
    # Incoming Longwave Radiation [W/m2]
    # via Stephan-Boltzmann law https://doi.org/10.1016/j.jrmge.2016.10.004

    emis = emissivity of atm [-]
    LST = Land Surface Temprature [˚C]
    """
    LST = LST + 273.15
    longIn = emisAtm * 5.6703 * 10.0 ** (-8.0) * LST ** 4
    savetif(longIn, outputPath)
    return longIn


def shortout(albedo, shortin, outputPath):
    """
    # Outgoing Shortwave Radiation [W/m2]
    # via https://www.posmet.ufv.br/wp-content/uploads/2016/09/MET-479-Waters-et-al-SEBAL.pdf

    albedo = Albedo [-]
    shortin = Shortwave Incoming Radiation [W/m2]
    """
    shortOut =  albedo * shortin
    savetif(shortOut, outputPath)  
    return shortOut

def Rn(shortIn, shortOut, longIn, longOut, outputPath):
    """
    # Net Energy Bdget [W/m2]
    # via https://www.redalyc.org/journal/2736/273652409002/html/#redalyc_273652409002_ref4

    shortIn = Incoming Shortwave Radiation [W/m2]
    shortOut = Outgoing Shortwave Radiation [W/m2]
    longIn = Incoming Longwave Radiation [W/m2]
    longOut = Outgoing Longwave Radiation [W/m2]
    """
    Rn = shortIn - shortOut + longIn - longOut
    savetif(Rn, outputPath)
    return Rn

def ra(airDensity, LST, Ta, eo, es, psychro, Rn, OutputPath, cp = 0.001013):
    """
    # Aerodynamic Resistance
    # via https://www.posmet.ufv.br/wp-content/uploads/2016/09/MET-479-Waters-et-al-SEBAL.pdf

    airDensity = Density of Air [kg/m3]
    LST = Land Surface Temperature [˚C]
    Ta = Air Temperature [˚C]
    eo = Partial Water Vapour Pressure [kPa]
    es = Saturated vapour pressure [kPa]
    psychro = Psychrometric constant [kPa/˚C]
    Rn = Net Energy Budget [W/m2]
    cp = Specific heat at constant pressure [MJ/kg/°C]

    """
    ra = (airDensity * cp * ((LST - Ta) + ((eo - es) / psychro))) / Rn
    savetif(ra, OutputPath)
    return ra

def sensHFlux(airDens, LST, ra, Ta, outputPath, cp = 0.001013):
    """
    # Sensible Heat Flux
    # via

    LST = Land Surface Temperature [˚C]
    ra = Air Resistence [s/m]
    Ta = Air Temperature [˚C]
    cp = Specific heat at constant pressure [MJ/kg/°C]
    airDens = Density of Air [kg/m3]
    LST_K = Land Surface Temperature [K]
    Ta_K = Air Temperature [K]
    """
    
    LST_K = LST + 273.15
    Ta_K = Ta + 273.15

    H = (airDens * cp * (LST_K - Ta_K)) / ra
    savetif(H, outputPath)
    return H

def soilGFlux(LST, albedo, ndvi, Rn, outputPath):
    """
    # Soil/Ground Heat Flux [W/m2]
    # via Baasriansen, 2000 (BASTIAANSSEN, W. G. M. SEBAL - based sensible and latent heat fluxes in the irrigated Gediz Basin, Turkey. Journal of Hydrology, v.229, p.87-100, 2000.)

    LST = Land Surface Temperature [˚C]
    albedo = Albedo [-]
    ndvi = Normal Differential Vegetation Index [-]
    Rn - Net ENergy Budget [W/m-2]
    """
    G = LST / albedo * (0.0038 * albedo + 0.0074 * albedo ** 2) * (1 - 0.98 * ndvi ** 4) * Rn
    savetif(G, outputPath)
    return G


def evapoFraction(Ta_max, Ta, LST, outputPath):
    """
    # Fraction of Evapotranspiration
    # via https://www.sciencedirect.com/science/article/pii/S0309170811000145?via%3Dihub

    Ta_max = Maximum Air Temperature [˚C]
    Ta = Air Temperature [˚C]
    LST = Land Surface Temperature [˚C]
    outputPath = path to output directory
    """
    
    EF = (Ta_max/LST) / (Ta_max / Ta)
    savetif(EF, outputPath)
    return EF

def le(EF, Rn, G, outputPath):
    """
    # Latent HEat Flux [W/m2]
    # via Baasriansen, 2000 (BASTIAANSSEN, W. G. M. SEBAL - based sensible and latent heat fluxes in the irrigated Gediz Basin, Turkey. Journal of Hydrology, v.229, p.87-100, 2000.)

    EF = Frantion of Evaporation [-]
    Rn = Net Energy Budget [W/m2]
    G = Soil/Ground Heat Flux [W/m2]
    """
    LE = EF * (Rn - G)
    savetif(LE, outputPath)
    return LE

def bowenIndex(H, LE, outputPath):
    """
    # Bowen Index
    # via https://daac.ornl.gov/FIFE/Datasets/Surface_Flux/Bowen_Ratio_USGS.html
    H = Sensible Heat Flux [W/m2]
    LE = Latent Heat Flux [W/m2]
    """
    BI = H / LE
    BI[BI < 0] = 0.1
    BI[BI > 5] = 4.9
    savetif(BI, outputPath)
    return BI

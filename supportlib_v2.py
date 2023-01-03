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

############################################################################################################################################
#   GENERAL FUNCTIONS
############################################################################################################################################

def getfilepath(input_folder, suffix):
    
    """
    # Get a filtered list of paths to files containing a certain suffix 

    input_folder = folder contaning the input files
    suffix = file suffix to be filtered (eg. .TIF, .JSON)
    """

    pathListFolder = []
    for root, dirs, files in os.walk(input_folder, topdown=False):
        for name in files:
            if name.endswith(suffix):
                if name not in pathListFolder:
                    pathListFolder.append(os.path.join(root, name))
    return pathListFolder

######################################################################

def createmeteodict(csv_file):

    """
    # Create a dictionary of meteorology data
    # eg. {'date': [avg_temp': '22.60','max_temp': '28.4','min_temp': '13.3','relHum': '70.21','wind_sp': '0.83']}

    csv_file = file in CSV format containing meterological data with date in original format (YYYYMMDD)
    """
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

    """
    # Image clipping using a pre-prepared GeoPackage mask in the same coordinate system as the images

    maskPath = path to polygon mask
    inputBand = image to be clipped
    outImgPath = path to new cropped image

    """

    with fiona.open(maskPath, "r") as gpkg:
        shapes = [feature["geometry"] for feature in gpkg]

    with rasterio.open(inputBand) as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta
    
    out_meta.update({"driver": "GTiff", # output format GeoTiff
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform})

    with rasterio.open(outImgPath, "w", **out_meta) as dest:
        dest.write(out_image)
    return 

######################################################################

def savetif(img, outputPath, epsg = 32633):

    """
    # Save and trasnform cropped image

    img = input cropped image to be saved
    outputPath = place to save the image
    epsg = a code for coordinate system by https://epsg.io/ standard
    """

    new_dataset = rasterio.open(outputPath, "w", 
        driver = "GTiff",
        height = img.shape[0],
        width = img.shape[1],
        count = 1,
        nodata = -9999, # optinal value for nodata
        dtype = img.dtype,
        crs = epsg, # driver for coordinate system code
        
        # upper left, pixel size
        # rasterio.transform.from_origin(west, north, xsize, ysize) -> affine transformation 
        transform = rasterio.transform.from_origin(656268, 5503485, 30, 30))
        #transform = rasterio.transform.from_origin(606952.1, 5485756.2, 30, 30))
    new_dataset.write(img, 1)
    new_dataset.close()
    
    return

############################################################################################################################################
#   ALBEDO                                                                                                                                 #
# via https://www.scielo.br/j/rbeaa/a/sX6cJjNXWMfHQ5p4h33B8Zz/?lang=en&format=pdf                                                          #
############################################################################################################################################
                                                                                                                                    
def dr(sunDist):

    """
    # Correction of the eccentricity of the terrestrial orbit

    sunDist = distance Earth-Sun in AU extracted from image metadata (EARTH_SUN_DISTANCE)
    """
    dr = (1/sunDist) ** 2 
    return dr

############################################################################################################################################

def zenithAngle(sunElev):

    """
    # Calculation of sun zenith angle in radians

    sunElev = Sun Elevation angle extracted from image metadata (SUN_ELEVATION)
    """
    zenith_angle = ((90 - sunElev) * math.pi) / 180
    return zenith_angle

############################################################################################################################################

def lb_band(AddRad, MultRad, band, outputPath):
   
    """
    # Pixel Radiance

    AddRad = Additive Radiance Term from the Metadata (RADIANCE_ADD_BAND_X)
    MultRef =  Radiance Multiplicative Term from the Metadata (RADIANCE_MULT_BAND_X)
    band = Landsat band
    """
    band = np.array(tf.imread(band))
    lb = (MultRad * band) + AddRad
    
    savetif(lb, outputPath)
    return lb

############################################################################################################################################

def rb_band(AddRef, MultRef, band, dr, zenithAngle, outputPath):
   
    """
    # Band Reflectance [W/m2]

    AddRef = Reflectance Additive Term from the Metadata (REFLECTANCE_ADD_BAND_X)
    MultRef =  Reflectance Multiplicative Term from the Metadata (REFLECTANCE_MULT_BAND_X)
    band = Landsat band
    dr = Correction of the eccentricity of the terrestrial orbit
    zenithAngle = Sun zenith angle in radians
    """
    band = np.array(tf.imread(band))
    
    #refl = (MultRef * band) + AddRef
    #rb = refl / (math.cos(zenithAngle) * (dr))

    rb = ((MultRef * band) + AddRef) / (math.cos(zenithAngle)*(dr))

    savetif(rb, outputPath)
    return rb

######################################################################

def kb(rb, lb, dr, zenithAngle, outputPath):
    """
    # Solar Constant [W/m2]

    rb = Band reflectance
    lb = Radiance of each pixel [W/m2]   
    dr = correction of the eccentricity of the terrestrial orbit

    """
    
    #rb_band = np.array(tf.imread(rb))
    #lb_band = np.array(tf.imread(lb))
    
    kb = (math.pi * lb) / (rb * math.cos(zenithAngle) * dr)
    kb[kb < 0] = 0

    savetif(kb, outputPath) 
    return kb

######################################################################

def toaplanetary(pb2,rb2, pb3, rb3, pb4, rb4,  pb5, rb5, pb6,rb6, pb7, rb7, outputPath):
    
    """
    # pPanetary albedo (without atmospheric correction)
    # Unitless or %
    # Range from 0 to 1 (ot 0 % to 100 %)

    toaplanet = Planetary Top Of Atmosphere Radiance  
    pbx = weight of spectral band
    rbx = and reflectance
    """

    toaPlanet = (pb2 * rb2) + (pb3 * rb3) + (pb4 * rb4) + (pb5 * rb5) + (pb6 * rb6) + (pb7 * rb7)
    savetif(toaPlanet, outputPath)
    
    return toaPlanet

######################################################################
def toc(P, zenithAngle, tpw, Kt = 1.0): #Kt from 1.0 for clear air, 0.5 for extremely polluted air
    """
    # Atmospheric transmittance in the solar radiation domain
    # Unitless

    P = Atmospheric Pressure[kPa]
    tpw = Total Precipitable Water [kg/m2]
    Kt = Air Turbidity coefficient (Kt= 1.0 for clear air and Kt= 0.5 for extremely turbid or polluted air)
    """
    exp =(-(0.00146 * P)/(Kt * zenithAngle)) - 0.075 * ((tpw/zenithAngle) * 0.4)
   
   
    Toc = 0.35 + (0.627 ** exp)
    
    return Toc

######################################################################

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

############################################################################################################################################
#   ATMOSHEPRIC FUNCTIONS
############################################################################################################################################
def e0(T):
    """
    # Partial Water Vapour Pressure [kPa]
    # Tetens Formula via https://www.omnicalculator.com/chemistry/vapour-pressure-of-water

    T = Air Temperature [˚C]
    """
    e0 = 0.6108 * math.exp(((17.27 * T) / (T + 237.3)))
    return e0

######################################################################

def es(Tamax, Tamin, Ta):
    """
    # Saturated Vapour Pressure [kPa]
    # via https://www.fao.org/3/x0490e/x0490e07.htm#atmospheric%20pressure%20(p)
    
    Ta = Average Air Temperature (˚C)
    """
    
    es = (6.1078 * 10**(7.5 * Ta /(Ta + 237.3))) / 10
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

def psychroCons(P):
    """
    # Psychrometric constant [kPa/˚C]
    # via https://www.fao.org/3/x0490e/x0490e07.htm#psychrometric%20constant%20(g)

    P = Atmospheric Pressure [kPa]
    """
    y =  0.00065 * P
    return y

######################################################################

def densityair(P, Ta, RH, R = 278.05):
    """
    Density of Air [kg/m3]
    # via https://designbuilder.co.uk/helpv3.4/Content/Calculation_of_Air_Density.htm
    
    R  = Gas Constant (287.05 J/kg-K)
    Ta = Air Temperature in K (T [˚C] + 273.15)
    P = Standard Pressure [kPa]
    """
    P = P * 1000

    #saturated vapour pressure in Pa
    p1 = 6.1078 * (10**(7.5 * Ta /(Ta + 237.3)))
    
    # the water vapor pressure in Pa
    pv = p1 * RH

    #pressure of dry air
    pd = P - pv
    air_density = (pd / (R * (Ta + 273.15))) + (pv / (461.495  * (Ta + 273.15)))
  
    return air_density

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

def z0m(vegHeigth):

    """
    # Roughness length governing momentum transfer [m]

    vegHeigth = height of vegetation
    """
    z0m = 0.123 * vegHeigth
    return z0m

############################################################################################################################################
#   FLUXES
############################################################################################################################################

def soilGFlux(LST, albedo, ndvi, Rn, outputPath):
    """
    # Soil/Ground Heat Flux [W/m2]
    # via BASTIAANSSEN, 2000 (BASTIAANSSEN, W. G. M. SEBAL - based sensible and latent heat fluxes in the irrigated Gediz Basin, Turkey. Journal of Hydrology, v.229, p.87-100, 2000.)

    LST = Land Surface Temperature [˚C]
    albedo = Albedo [-]
    ndvi = Normal Differential Vegetation Index [-]
    Rn - Net ENergy Budget [W/m-2]
    """
    G = LST / albedo * (0.0038 * albedo + 0.0074 * albedo ** 2) * (1 - 0.98 * ndvi ** 4) * Rn
    
    #LST_K = LST - 273.15
    #G = (Rn * (-13.46 + 0.507 * (4 * np.exp(0.123*LST_K)))) + 0.0086

    savetif(G, outputPath)
    return G

######################################################################

def H(LE, Rn, G, outputPath):
    # gradient method

    #LE = np.array(tf.imread(LE))
    #Rn = np.array(tf.imread(Rn))
    #G = np.array(tf.imread(G))

    h = (Rn - G) - LE
    savetif(h, outputPath)

    return h

######################################################################

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

############################################################################################################################################
#   RADIATION
############################################################################################################################################

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

######################################################################

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

######################################################################

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

######################################################################

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

############################################################################################################################################
#   SURFACE
############################################################################################################################################

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

######################################################################

def LST(BT, emis, band, outputPath):
    """
    # Land Surface Temperature [˚C]
    # via Ugur Avdan and Gordana Jovanovska. “Automated
            Mapping of Land Surface Temperature Using
            LANDSAT 8 Satellite Data”, Journal of Sensors,
            Vol. 2016, Article ID 1480307, 2016.

    BT - Brightness temperature [K]
    emis - emissivity [-]
    band - Landsat band [-]
    """
    band = np.array(tf.imread(band))
    
    
    LST = (BT / (1 + ((0.0015 * BT)/1.4488) * np.log(emis)))
    
    savetif(LST, outputPath)
    return LST

######################################################################

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
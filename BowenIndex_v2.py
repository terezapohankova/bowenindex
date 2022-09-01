
from codecs import readbuffer_encode
import os, sys
from socketserver import DatagramRequestHandler
from sre_constants import IN
from pprint import pprint
import csv
#from BowenIndex import RadAddBandB2, RadAddBandB4_L1TP, RadAddBandB6_L1TP
import supportlib_v2
import numpy as np
import math
import json
from collections import defaultdict


# get data
INPUT_FOLDER = r'/home/tereza/ownCloud/skripty/BowenIndex/snimky_L9_testovaci'
OUTPUT_PATH = r'/home/tereza/Documents/testy_VSC'
MASK = r'/home/tereza/Documents/boundary/olomouc_32633.gpkg'
METEOROLOGY = r'/home/tereza/ownCloud/skripty/BowenIndex/snimky_L9_testovaci/weather.csv'


# constants
Z = 219 #elevetaion in m
CP = 0.001013 # cp specific heat at constant pressure, 1.013 10-3 [MJ kg-1 °C-1], https://www.fao.org/3/x0490e/x0490e07.htm

# folders for processed images
OUT_CLIP_FOLDER = os.path.join(OUTPUT_PATH,'clipped_bands')
OUT_VEG_INDEX_FOLDER = os.path.join(OUTPUT_PATH,'vegIndices')
OUT_VEG_HEIGHT_FOLDER = os.path.join(OUTPUT_PATH,'vegHeight')
OUT_FLUX_FOLDER = os.path.join(OUTPUT_PATH,'flux')
OUT_NET_FOLDER = os.path.join(OUTPUT_PATH,'netBudget')
OUT_ALBEDO_FOLDER = os.path.join(OUTPUT_PATH,'albedo')
OUT_LST_FOLDER = os.path.join(OUTPUT_PATH,'lst')
OUT_BOWEN_FOLDER = os.path.join(OUTPUT_PATH,'bowenindex')

os.makedirs(OUT_CLIP_FOLDER, exist_ok = True)
os.makedirs(OUT_VEG_INDEX_FOLDER, exist_ok = True)
os.makedirs(OUT_VEG_HEIGHT_FOLDER, exist_ok = True)
os.makedirs(OUT_FLUX_FOLDER, exist_ok = True)
os.makedirs(OUT_NET_FOLDER, exist_ok = True)
os.makedirs(OUT_ALBEDO_FOLDER, exist_ok = True)
os.makedirs(OUT_LST_FOLDER, exist_ok = True)
os.makedirs(OUT_BOWEN_FOLDER, exist_ok = True)

meteorologyDict = supportlib_v2.createmeteodict(METEOROLOGY) #{{'20220518': {'avg_temp': '14.25','max_temp': '19.8','min_temp': '8.7','relHum': '65.52','wind_sp': '1.25'}},
#pprint(meteorologyDict)
# e.g. meteorologyDict[date]['avg_temp']




# get paths to images and jsons
JSON_MTL_PATH = supportlib_v2.getfilepath(INPUT_FOLDER, 'MTL.json') #['root/snimky_L9_testovaci/LC09_L2SP_190025_20220518_20220520_02_T1/LC09_L2SP_190025_20220518_20220520_02_T1_MTL.json']
ORIGINAL_IMG = supportlib_v2.getfilepath(INPUT_FOLDER, '.TIF') #['root/snimky_L9_testovaci/18052022/LC09_L2SP_190025_20220518_20220518_02_T1_SZA.TIF']

# empty stuff
mtlJSONFile = {}
sensingDate = []
#nni knstanta
CLIPPED_IMG_PATHS = []

# load JSON MTL file with metadata into dictionary {sensingdate : {metadatafile}} for level 2 (level 2 MTL json includes level 1 MTL data)
for jsonFile in JSON_MTL_PATH:
    if 'L2SP' in jsonFile:
        loadJSON = supportlib_v2.load_json(jsonFile)
        sensDate = jsonFile.split('_')[5] # 20220518
        #imgLevel = jsonFile.split('_')[3] # 'L2SP'
        
        if sensDate not in sensingDate:
            sensingDate.append(sensDate)
            if not os.path.exists(os.path.join(OUTPUT_PATH, OUT_CLIP_FOLDER, sensDate)): 
                os.makedirs(os.path.join(OUTPUT_PATH, OUT_CLIP_FOLDER, sensDate))
        
        #for jsonFile in JSON_MTL_PATH:
            #if 'L2SP' in jsonFile:
                #loadJSON = supportlib_v2.load_json(jsonFile)
                #sensDate = jsonFile.split('_')[5] # 20220518
                #sensingDate.append(sensDate)
            
            
        mtlJSONFile[sensDate] = loadJSON
        
#pprint(((mtlJSONFile)))


# create output path for clipped images by pairing sensing date from JSON metadata file and sensing date on original images
#for inputBand in ORIGINAL_IMG:
#    for date in sensingDate:
#        if os.path.basename(inputBand).split('_')[3] == date: # if date on original input band equals date sensing date from json mtl, then append it to the list
#            CLIPPED_IMG_PATHS.append(os.path.join(OUTPUT_PATH, RESULT_FOLDERS[0], date, 'clipped_' + os.path.basename(inputBand)))


imgDict = {} # {sensingdate : {path : path, radiance : int ...} }
kbDict = {}
kbList = []

for inputBand in ORIGINAL_IMG:
    # if date on original input band equals date sensing date from json mtl, then append it to the list
    
    image_basename = os.path.basename(inputBand) # 'LC09_L1TP_190025_20220518_20220518_02_T1_B6_clipped.TIF'
    
    if 'B' in image_basename:
        image_name = image_basename.replace('.TIF','') #'LC09_L2SP_189026_20220612_20220614_02_T1_SR_B1'
    
        #pprint(imgBandList)
        
        date = image_basename.split('_')[3] # '20220612'
        
    
        #.split('_')[-1] - last splitted value which should be B1 - B10 
        band = image_basename.replace('.TIF','').split('_')[-1] # 'B1
    
        # from basename by splitting keep L1TP, by [:2] keep just L1
        image_level = image_basename.split('_')[1][:2] # 'L2'
        clippedImgPath = os.path.join(OUTPUT_PATH, OUT_CLIP_FOLDER, date, 'clipped_' + os.path.basename(inputBand)) # '/home/tereza/Documents/testy_VSC/clipped_bands/20220612/clipped_LC09_L2SP_189026_20220612_20220614_02_T1_SR_B1.TIF'
        
        # band_level_key must be unique
        band_level_key = f'{band}_{image_level}' # 'B4_L2'
        
        #supportlib_v2.clipimage(MASK, inputBand, clippedImgPath)

    if date not in imgDict:
        imgDict.setdefault(date, {})

    imgDict[date][band_level_key] = {
            'clipped_path' : clippedImgPath,
            'imageName' : image_name,
            'RADIANCE_ADD' : float(mtlJSONFile[date]['LANDSAT_METADATA_FILE']['LEVEL1_RADIOMETRIC_RESCALING'].get(f'RADIANCE_ADD_BAND_{band[1:]}') ), #[1:] - delete B
            'RADIANCE_MULT' : float(mtlJSONFile[date]['LANDSAT_METADATA_FILE']['LEVEL1_RADIOMETRIC_RESCALING'].get(f'RADIANCE_MULT_BAND_{band[1:]}') ),
            'KELVIN_CONS_1' : float(mtlJSONFile[date]['LANDSAT_METADATA_FILE']['LEVEL1_THERMAL_CONSTANTS'].get(f'K1_CONSTANT_BAND_{band[1:]}') or 0),
            'KELVIN_CONS_2' : float(mtlJSONFile[date]['LANDSAT_METADATA_FILE']['LEVEL1_THERMAL_CONSTANTS'].get(f'K2_CONSTANT_BAND_{band[1:]}') or 0),
            'REFLECTANCE_ADD': float(mtlJSONFile[date]['LANDSAT_METADATA_FILE']['LEVEL2_SURFACE_REFLECTANCE_PARAMETERS'].get(f'REFLECTANCE_ADD_BAND_{band[1:]}') or 0),
            'REFLECTANCE_MULT' : float(mtlJSONFile[date]['LANDSAT_METADATA_FILE']['LEVEL2_SURFACE_REFLECTANCE_PARAMETERS'].get(f'REFLECTANCE_MULT_BAND_{band[1:]}') or 0),
            'dES' : float(mtlJSONFile[date]['LANDSAT_METADATA_FILE']['IMAGE_ATTRIBUTES'].get('EARTH_SUN_DISTANCE')),
            'sunAzimuth' : float(mtlJSONFile[date]['LANDSAT_METADATA_FILE']['IMAGE_ATTRIBUTES'].get('SUN_AZIMUTH')),
            'sunElev' : float(mtlJSONFile[date]['LANDSAT_METADATA_FILE']['IMAGE_ATTRIBUTES'].get('SUN_ELEVATION')),
            }
    
    
# for each sensing date in, calculate meteorology values:
for date in sensingDate:

    
    #pprint(f' sensing date : {date}')
    
    e0 = supportlib_v2.e0(meteorologyDict[date]['avg_temp'])
    #pprint(f' e0 : {e0}')

    es = supportlib_v2.es(meteorologyDict[date]['max_temp'], meteorologyDict[date]['min_temp'])
    #pprint(f' es : {es}')

    p = supportlib_v2.atmPress(Z)
    #pprint(f' atmPress : {p}')

    psychro = supportlib_v2.psychroCons(p)
    #pprint(f' psychro : {psychro}')

    tpw = supportlib_v2.tpw(e0, p) 
    #pprint(f' tpw : {tpw}')

    toc = supportlib_v2.toc(p, Z, tpw)
    #pprint(f' toc : {toc}')

    atmEmis = supportlib_v2.atmemis(es, meteorologyDict[date]['avg_temp'])
    #pprint(f' atmEmis : {atmEmis}')

    airDensity = supportlib_v2.densityair(Z, meteorologyDict[date]['avg_temp'])
    #pprint(f' airDensity : {airDensity}')

    ndvi = supportlib_v2.ndvi(imgDict[date]['B5_L2']['clipped_path'], imgDict[date]['B4_L2']['clipped_path'], 
            os.path.join(OUTPUT_PATH, OUT_VEG_INDEX_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'ndvi' + '.TIF')))
    bt = supportlib_v2.bt(imgDict[date]['B10_L1']['KELVIN_CONS_1'], imgDict[date]['B10_L1']['KELVIN_CONS_2'], 
            imgDict[date]['B10_L1']['RADIANCE_ADD'], imgDict[date]['B10_L1']['RADIANCE_MULT'], imgDict[date]['B10_L1']['clipped_path'])

    pv = supportlib_v2.pv(ndvi)
    
    emisSurf = supportlib_v2.emis(imgDict[date]['B4_L2']['clipped_path'], ndvi, pv)

    lst = supportlib_v2.LST(bt, emisSurf, imgDict[date]['B10_L1']['clipped_path'], 
            os.path.join(OUTPUT_PATH, OUT_LST_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'lst' + '.TIF')))
    
    #calculate rb for each band from b2 to b7
    """for date, bandLevelDict in imgDict.items():
        for bandLevel, valuesDict in bandLevelDict.items():
            rb = supportlib_v2.rb_band(imgDict[date][bandLevel]['REFLECTANCE_ADD'], imgDict[date][bandLevel]['RADIANCE_MULT'], 
                    imgDict[date][bandLevel]['clipped_path'], imgDict[date][bandLevel]['sunElev'],
                    os.path.join(OUTPUT_PATH, OUT_LST_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'rb_' + bandLevel + '.TIF')))
            
            kb = supportlib_v2.kb(rb, Z, imgDict[date][bandLevel]['RADIANCE_ADD'], imgDict[date][bandLevel]['RADIANCE_MULT'], 
                    imgDict[date][bandLevel]['clipped_path'], imgDict[date]['B5_L2']['dES'],
                    os.path.join(OUTPUT_PATH, OUT_LST_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'kb_' + bandLevel + '.TIF')))
        
            #kbDict.setdefault(date, [])
            #pprint(kbDict[date].update('h'))

            # append all kb paths to list, get unique values by set and convert back to list
            kbList.append(os.path.join(OUTPUT_PATH, OUT_LST_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'kb_' + bandLevel + '.TIF')))
            kbList = (list(set(kbList)))"""

    rb2 = supportlib_v2.rb_band(imgDict[date]['B2_L2']['REFLECTANCE_ADD'], imgDict[date]['B2_L2']['REFLECTANCE_MULT'], 
                    imgDict[date]['B2_L2']['clipped_path'], imgDict[date]['B2_L2']['sunElev'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'rb_' + 'B2' + '.TIF')))
    rb3 = supportlib_v2.rb_band(imgDict[date]['B3_L2']['REFLECTANCE_ADD'], imgDict[date]['B3_L2']['REFLECTANCE_MULT'], 
                    imgDict[date]['B3_L2']['clipped_path'], imgDict[date]['B3_L2']['sunElev'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'rb_' + 'B3' + '.TIF')))
    rb4 = supportlib_v2.rb_band(imgDict[date]['B4_L2']['REFLECTANCE_ADD'], imgDict[date]['B4_L2']['REFLECTANCE_MULT'], 
                    imgDict[date]['B4_L2']['clipped_path'], imgDict[date]['B4_L2']['sunElev'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'rb_' + 'B4' + '.TIF')))
    rb5 = supportlib_v2.rb_band(imgDict[date]['B5_L2']['REFLECTANCE_ADD'], imgDict[date]['B5_L2']['REFLECTANCE_MULT'], 
                    imgDict[date]['B5_L2']['clipped_path'], imgDict[date]['B5_L2']['sunElev'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'rb_' + 'B5' + '.TIF')))
    rb6 = supportlib_v2.rb_band(imgDict[date]['B6_L2']['REFLECTANCE_ADD'], imgDict[date]['B6_L2']['REFLECTANCE_MULT'], 
                    imgDict[date]['B6_L2']['clipped_path'], imgDict[date]['B6_L2']['sunElev'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'rb_' + 'B6' + '.TIF')))
    rb7 = supportlib_v2.rb_band(imgDict[date]['B7_L2']['REFLECTANCE_ADD'], imgDict[date]['B7_L2']['REFLECTANCE_MULT'], 
                    imgDict[date]['B7_L2']['clipped_path'], imgDict[date]['B7_L2']['sunElev'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'rb_' + 'B7' + '.TIF')))
    
    
    kb2 = supportlib_v2.kb(rb2, Z, imgDict[date]['B2_L2']['RADIANCE_ADD'], imgDict[date]['B2_L2']['RADIANCE_MULT'], 
                    imgDict[date]['B2_L2']['clipped_path'], imgDict[date]['B2_L2']['dES'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'kb_' + 'B2' + '.TIF')))

    kb3 = supportlib_v2.kb(rb3, Z, imgDict[date]['B3_L2']['RADIANCE_ADD'], imgDict[date]['B3_L2']['RADIANCE_MULT'], 
                    imgDict[date]['B3_L2']['clipped_path'], imgDict[date]['B3_L2']['dES'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'kb_' + 'B3' + '.TIF')))

    kb4 = supportlib_v2.kb(rb4, Z, imgDict[date]['B4_L2']['RADIANCE_ADD'], imgDict[date]['B4_L2']['RADIANCE_MULT'], 
                    imgDict[date]['B4_L2']['clipped_path'], imgDict[date]['B4_L2']['dES'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'kb_' + 'B4' + '.TIF')))

    kb5 = supportlib_v2.kb(rb5, Z, imgDict[date]['B5_L2']['RADIANCE_ADD'], imgDict[date]['B5_L2']['RADIANCE_MULT'], 
                    imgDict[date]['B5_L2']['clipped_path'], imgDict[date]['B5_L2']['dES'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'kb_' + 'B5' + '.TIF')))

    kb6 = supportlib_v2.kb(rb6, Z, imgDict[date]['B6_L2']['RADIANCE_ADD'], imgDict[date]['B6_L2']['RADIANCE_MULT'], 
                    imgDict[date]['B6_L2']['clipped_path'], imgDict[date]['B6_L2']['dES'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'kb_' + 'B6' + '.TIF')))

    kb7 = supportlib_v2.kb(rb7, Z, imgDict[date]['B7_L2']['RADIANCE_ADD'], imgDict[date]['B7_L2']['RADIANCE_MULT'], 
                    imgDict[date]['B7_L2']['clipped_path'], imgDict[date]['B7_L2']['dES'],
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'kb_' + 'B7' + '.TIF')))

    
    """pb2 = supportlib_v2.pb(kb2, kb2, kb3, kb4, kb5, kb6, kb7, 
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'pb_' + 'B2' + '.TIF')))
    pb3 = supportlib_v2.pb(kb3, kb2, kb3, kb4, kb5, kb6, kb7, 
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'pb_' + 'B3' + '.TIF')))
    pb4 = supportlib_v2.pb(kb4, kb2, kb3, kb4, kb5, kb6, kb7, 
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'pb_' + 'B4' + '.TIF')))
    pb5 = supportlib_v2.pb(kb5, kb2, kb3, kb4, kb5, kb6, kb7, 
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'pb_' + 'B5' + '.TIF')))
    pb6 = supportlib_v2.pb(kb6, kb2, kb3, kb4, kb5, kb6, kb7, 
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'pb_' + 'B6' + '.TIF')))
    pb7 = supportlib_v2.pb(kb7, kb2, kb3, kb4, kb5, kb6, kb7, 
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'pb_' + 'B7' + '.TIF')))"""

    kbDterminator = (kb2 + kb3 + kb4 + kb5 + kb6 + kb7)
    pb2 = kb2 / kbDterminator   
    pb3 = kb3 / kbDterminator
    pb4 = kb4 / kbDterminator  
    pb5 = kb5 / kbDterminator 
    pb6 = kb6 / kbDterminator 
    pb7 = kb7 / kbDterminator

    toaPlanet = supportlib_v2.toaplanetary(pb2,rb2, pb3, rb3, pb4, rb4,  pb5, rb5, pb6,rb6, pb7, rb7, 
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'toaPlanet' + '.TIF')))

    albedo = supportlib_v2.albedo(toaPlanet, toc, 
                    os.path.join(OUTPUT_PATH, OUT_ALBEDO_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'albedo' + '.TIF')))

    longOut = supportlib_v2.longout(emisSurf, lst, os.path.join(OUTPUT_PATH, OUT_NET_FOLDER ,
                    os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'longOut' + '.TIF')))
    longIn = supportlib_v2.longin(atmEmis, lst, 
                    os.path.join(OUTPUT_PATH, OUT_NET_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'longIn' + '.TIF')))
    
    #shortnIn = supportlib_v2.shortIn(imgDict[date]['B2_L2']['sunElev'], imgDict[date]['B2_L2']['dES'], Z)
    shortIn = 650
    shortOut = supportlib_v2.shortout(albedo, shortIn, 
                    os.path.join(OUTPUT_PATH, OUT_NET_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'shortOut' + '.TIF')))

    rn = supportlib_v2.Rn(shortIn, shortOut, longIn, longOut, 
                    os.path.join(OUTPUT_PATH, OUT_NET_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'rn' + '.TIF')))

    ra = supportlib_v2.ra(airDensity, lst, meteorologyDict[date]['avg_temp'], e0, es, psychro, rn, 
                    os.path.join(OUTPUT_PATH, OUT_NET_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'ra' + '.TIF')))
    
    h = supportlib_v2.sensHFlux(airDensity, lst, ra, meteorologyDict[date]['avg_temp'], 
                    os.path.join(OUTPUT_PATH, OUT_FLUX_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'sensHeatFLux' + '.TIF')))
    
    g = supportlib_v2.soilGFlux(lst, albedo, ndvi, rn, 
                    os.path.join(OUTPUT_PATH, OUT_FLUX_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'SoilHeatFLux' + '.TIF')))

    ef = supportlib_v2.evapoFraction(meteorologyDict[date]['max_temp'], meteorologyDict[date]['avg_temp'], lst, 
                    os.path.join(OUTPUT_PATH, OUT_FLUX_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'ef' + '.TIF')))

    le = supportlib_v2.le(ef, rn, g, 
                    os.path.join(OUTPUT_PATH, OUT_FLUX_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'le' + '.TIF')))
    
    bowenIndex = supportlib_v2.bowenIndex(h, le, 
                    os.path.join(OUTPUT_PATH, OUT_BOWEN_FOLDER, os.path.basename(imgDict[date]['B5_L2']['clipped_path']).replace('B5.TIF', 'bowen' + '.TIF')))

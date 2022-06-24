# -*- coding: utf-8 -*-

#Sentinel-2 Imagery API
from sentinelsat import SentinelAPI

#Data Libraries
import geopandas as gpd
import pandas as pd
import fiona
from scipy import stats
import matplotlib.pyplot as plt

#Raster Libraries
import rasterio as rio
from rasterio.merge import merge
from rasterio.plot import show
from rasterstats import zonal_stats
import rasterio.plot as rplt
import numpy


#Py Libraries
import glob
import os
from zipfile import ZipFile
import pathlib
import warnings
import platform
import random
import progressbar
import sys

#Web Scraping Libraries
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException
from selenium.common.exceptions import SessionNotCreatedException


#Prompt user to provide working directory to files and return path to directory
def get_working_directory():
    working_directory_path = str(pathlib.Path(__file__).parent.absolute())
    os.chdir(working_directory_path)
    return working_directory_path

#Query sentinelsat API and downloaded imagery over LA and return a geodataframe of queried products(some hard coding required)
def get_sentinel_imagery(workspace_path):

    print('\nQuerying Sentinel API...')
    #create a bounding box for api query using the bounds of LA County (not including Southern Islands)
    #footprint = 'POLYGON((-118.9447 33.60,-118.9447 34.8233,-117.6463 34.8233,-117.6463 33.6,	-118.9447 33.6))'
    footprint = 'POLYGON((-118.30348139323311 34.11871855411627,-117.78085403190155 34.11871855411627,-117.78085403190155 34.50171775223244,-118.30348139323311 34.50171775223244,-118.30348139323311 34.11871855411627))'
    sent_user = 'final_project'
    sent_password = 'fppassword'
    api = SentinelAPI(sent_user, sent_password, 'https://scihub.copernicus.eu/dhus')
    print(api.api_url)

    #Query API looking for imagery that intersects with the footprint (LA), in June of 2020 (max veg growth),
    # cloud cover from 0-90 %, and with processing level 2A (bottom of atmosphere reflectance correction).
    products = api.query(footprint, area_relation='Intersects', date=('20210401', '20210429'),
                     platformname='Sentinel-2',
                     processinglevel='Level-2A',
                     cloudcoverpercentage=(0, 10))


    #save api query results as gdf and show products sorted by cloud cover percentange
    products_gdf = api.to_geodataframe(products)

    #drop overlapping images (with the same geometry)
    products_gdf = products_gdf.drop_duplicates(subset ='geometry',
                         keep='first')

    print('Query complete...')

    print('Downloading {} zipped Sentinel-2 images in Level-2C .SAFE format...'.format(len(products_gdf)))
    for product in products_gdf['uuid']:
        api.download(str(product))

    #Unzip .SAFE folders and move to Imagery folder
    print('\nDownload finished! Unzipping imagery and moving to Imagery folder...')
    zipped_files = glob.glob(os.path.join(workspace_path, '*.zip'))
    for file in zipped_files:
        file_path = os.path.join(workspace_path, str(file))
        if os.path.exists(file_path):
            with ZipFile(file_path, 'r') as zip:
                # Extract all the contents of zip file in different directory
                zip.extractall(os.path.join(workspace_path, 'Data/Imagery'))
        else:
             raise FileNotFoundError("Couldn't find Sentinel image. Download may have failed...")

    print("Imagery unzipped...")
    return products_gdf

#create Gtiffs of NDVI from imagery downloaded from API and write to individual images and a mosaic
def create_NDVI_mosaic(product_gdf, workspace_path):

    print('Beginning NDVI imagery mosaic...')
    imagery_path = workspace_path + '/Data/Imagery'
    count = 1
    #Go through each image in the unzipped imager folders and read into numpy arrays
    for product in product_gdf['title']:

        granule = '{}/{}.SAFE/GRANULE'.format(imagery_path, str(product))
        image = str(os.listdir(granule)[0])
        R10 = ('{}/{}/IMG_DATA/R10m'.format(granule, image))

        jp2s = os.listdir(R10)
        for jp2 in jp2s:
            if str(jp2).endswith('_B04_10m.jp2'):
                b4 = rio.open(R10 + '/' + str(jp2))
                red = b4.read()
            if str(jp2).endswith('_B08_10m.jp2'):
                b8 = rio.open(R10 + '/' + str(jp2))
                nir = b8.read()

        # Calculate NDVI into numpy array
        print('Calculating NDVI for image #{}...'.format(str(count)))
        numpy.seterr(divide='ignore', invalid='ignore')
        ndvi = ((nir.astype(float) - red.astype(float)) / (nir + red))

        # Define spatial characteristics of output object (basically same as input)
        kwargs = b8.meta

        # Update metadata (just in case)
        kwargs.update(
            driver='GTiff',
            dtype=rio.float32,
            nodata=0,
            count=1)

        #Write individual images to drive
        print('Writing NDVI image #{} to drive...'.format(str(count)))
        with rio.open((imagery_path + '/{}{}{}'.format('NDVI_', str(count),'.tif')), 'w', **kwargs) as dst:
            dst.write(ndvi.astype(rio.float32))

        count += 1

    # Make a search criteria to select the NDVI tifs
    search_criteria = "NDVI*.tif"
    q = os.path.join(imagery_path, search_criteria)

    # Get a list of rasterio image strings
    NDVI_imagery = glob.glob(q)
    NDVI_imagery_list = []
    for image in NDVI_imagery:
        tile = rio.open(str(image))
        NDVI_imagery_list.append(tile)

    # Merge images into one .tif
    print('Merging images into mosaic...')
    NDVI_mosaic, out_trans = merge(NDVI_imagery_list)

    # Write the mosaic to drive
    print('Writing mosaic to drive...')
    with rio.open(imagery_path + '/NDVI_mosaic.tif', 'w', **kwargs) as dst:
        dst.write(NDVI_mosaic.astype(rio.float32))
    print('Imagery prepared.')
    print('PART 1 COMPLETE...')

#limple function to return a list of columns from geodatafrmae
def list_columns(df):
    field_list = list(df)
    return field_list

# Open California School Campus Database (CSCD) and return public schools as Geodataframe, then project to (WGS84 UTM 11N), clip
# to imagery extent, drop duplicate CDS Codes, and filter unnecessary fields.
def CSCD_to_gdf(workspace_path, products_gdf):
    print('Preparing CSCD as GDF...')
    cscd_gdb = workspace_path + '/Data/CSCD_2021.gdb'
    schools_gdf = gpd.read_file(cscd_gdb, layer='Schools_Current_Stacked')
    schools_gdf = schools_gdf.drop(
        ['Ed_Type', 'Notes', 'Locations', 'GradesOffered', 'GradesServed', 'PartCount', 'CDS_uniq'],
        axis=1)
    products_gdf = products_gdf.to_crs(32611)  # EPSG code for WGS84 UTM Zone 11N
    schools_gdf = schools_gdf.to_crs(32611)
    schools_gdf = gpd.clip(schools_gdf, products_gdf)
    schools_gdf = schools_gdf.drop_duplicates(subset=['CDSCode', 'School'])
    print('CSCD prepared...')
    return schools_gdf

#Calculate average NDVI values for each school in the in the prepped CSCD gdf
def calculate_zonal_statistics(workspace_path, schools_gdf):
    ndvi = str(workspace_path) + '/Data/Imagery/NDVI_mosaic.tif'
    print('Calculating zonal statistics...')
    result = zonal_stats(schools_gdf, ndvi, stats=['mean'], geojson_out=True)
    geostats = gpd.GeoDataFrame.from_features(result)
    geostats = geostats.dropna(subset=['mean'])
    geostats = geostats.rename(columns={'mean': 'Mean_NDVI'})
    print('Zonal statistics calculated.')
    print('PART 2 COMPLETE...')
    return geostats

#Scrape California School Dashboard and dataframe containing (cds_code, percent of students on free or reduced lunch)
def scrape_CSD(workspace_path, cds_codes):

    #Set webdriver settings
    options = webdriver.ChromeOptions()
    options.add_argument(' — incognito')
    options.add_argument('-headless')
    options.add_argument('-no-sandbox')
    options.add_argument('-disable-dev-shm-usage')

    #Detect OS
    os = platform.system()

    if os == 'Darwin':
        try:
            driver_path = workspace_path + '/Selenium/Mac/0.90/chromedriver'
            wd = webdriver.Chrome(executable_path=driver_path, options=options)
        except SessionNotCreatedException:
            driver_path = workspace_path + '/Selenium/Mac/0.89/chromedriver'
            wd = webdriver.Chrome(executable_path=driver_path, options=options)
    elif os == 'Windows':
        try:
            driver_path = workspace_path + '/Selenium/Windows/0.90/chromedriver'
            wd = webdriver.Chrome(executable_path=driver_path, options=options)
        except SessionNotCreatedException:
            driver_path = workspace_path + '/Selenium/Windows/0.89/chromedriver'
            wd = webdriver.Chrome(executable_path=driver_path, options=options)
    elif os == 'Linux':
        try:
            driver_path = workspace_path + '/Selenium/Linux/0.90/chromedriver'
            wd = webdriver.Chrome(executable_path=driver_path, options=options)
        except SessionNotCreatedException:
            driver_path = workspace_path + '/Selenium/Linux/0.89/chromedriver'
            wd = webdriver.Chrome(executable_path=driver_path, options=options)
    else:
        raise RuntimeError('Unsupported operating system... Script only supports standard MacOS, Windows, and Linux systems')

    print('Scraping {} data points from CSD...'.format(len(cds_codes)))
    frl = [] #percent of students on free or reduced lunch
    bar = progressbar.ProgressBar(widgets=[
        progressbar.Percentage(),
        progressbar.Bar(),
    ], max_value=len(cds_codes)).start()
    for code in cds_codes:

      wd = webdriver.Chrome(executable_path=driver_path, options=options)
      wd.get('https://www.caschooldashboard.org/reports/' + str(code) + '/2020')

      delay = 2 # second
      try:
          source = WebDriverWait(wd, delay).until(EC.presence_of_element_located((By.ID, 'disadvantaged')))
      except TimeoutException:
        print("Loading took too much time!")

      tags = source.find_element_by_tag_name('p')
      frl.append(float(tags.text.rstrip('%')))
      bar += 1

    bar.finish()
    wd.quit()

    school_metric =  dict({'CDSCode': cds_codes, 'Percent_Disadvantaged': frl})
    school_metric_df = pd.DataFrame(data=school_metric)
    print('CSD scraped...')
    print('PART 3 COMPLETE...')
    return school_metric_df

def get_random_CDS_codes(school_stats_gdf, sample_n):
    cds_codes = school_stats_gdf['CDSCode'].tolist()
    cds_codes = random.sample(cds_codes, sample_n)
    print('CDS Code sampled...')
    return cds_codes

#Add scraped data to gdf containing school info and zonal statistics and provide option to save to GeoJSON
def merge_dataframes(workspace_path, school_stats_gdf, school_metric_df):
    final_gdf = pd.merge(school_stats_gdf, school_metric_df, on='CDSCode')
    final_gdf.to_file(str(workspace_path + '/Data/green_schools.geojson'), driver='GeoJSON')
    return final_gdf

#Run a linear regression analysis between % of students on free or reduced lunch (x) and mean campus NDVI values (y)
#Returns a dictonary containing the slope, intercept, r_value, p_value, and std_error
def linear_regression(final_gdf):
    x = final_gdf['Percent_Disadvantaged']
    y = final_gdf['Mean_NDVI']
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    stats_dict = {'slope' : slope, 'intercept' : intercept, 'r_value' : r_value, 'p_value' : p_value,'std_err' : std_err}
    return stats_dict

def show_NDVI_imagery(workspace_path, static):
    if static == True:
        src = rio.open(workspace_path + '/StaticData/NDVI_mosaic.tif')
    else:
        src = rio.open(workspace_path + '/Data/Imagery/NDVI_mosaic.tif')

    plt.imshow(src.read(1), cmap='RdYlGn')
    plt.pause(3)
    plt.close()

def show_regression_model(workspace_path, stats_dict, static):
    if static == True:
        filename = workspace_path + '/StaticData/green_schools.geojson'
    else:
        filename = workspace_path + '/Data/green_schools.geojson'

    final_gdf = gpd.read_file(filename)
    x = final_gdf['Percent_Disadvantaged']
    y = final_gdf['Mean_NDVI']
    plt.plot(x, y, 'o', label='schools')
    plt.plot(x, stats_dict['intercept'] + stats_dict['slope'] * x, color='red', label='Linear Regression')
    plt.xlabel("Percentge of Students on FRL")
    plt.ylabel('NDVI Value')
    plt.legend()
    plt.show(block=False)
    plt.pause(3)
    plt.close()

def main():

    args = sys.argv

    #user_input = str(input("Would you like to run the program in  'static' or 'default' mode?\n")).lower().strip()
    user_input = sys
    if len(args) == 1:
        static = False
        print('You selected to run the script in default mode...')
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=FutureWarning)
        workspace_path = get_working_directory()
        products_gdf = get_sentinel_imagery(workspace_path)
        create_NDVI_mosaic(products_gdf, workspace_path)
        school_stats_gdf = calculate_zonal_statistics(workspace_path, CSCD_to_gdf(workspace_path, products_gdf))
        sample_cds_codes = get_random_CDS_codes(school_stats_gdf, 100)
        school_metric_df = scrape_CSD(workspace_path, sample_cds_codes)
        final_gdf = merge_dataframes(workspace_path, school_stats_gdf, school_metric_df)
        print(final_gdf[['CDSCode', 'Mean_NDVI', 'Percent_Disadvantaged']])
        reg = linear_regression(final_gdf)
        print('\nr-value :' + str(reg['r_value']))
        print('Showing visualizations...')
        show_regression_model(workspace_path, reg, static)
        show_NDVI_imagery(workspace_path, static)
    elif str(args[1]) == '--static':
        static = True
        print('You selected to run the script in static mode...')
        workspace_path = get_working_directory()
        print(workspace_path)
        filename = workspace_path + '/StaticData/green_schools.geojson'
        green_schools_gdf = gpd.read_file(filename)
        print('Data:')
        green_schools_gdf = green_schools_gdf[green_schools_gdf['Percent_Disadvantaged'] != 0]
        print(green_schools_gdf[['CDSCode','Mean_NDVI','Percent_Disadvantaged']])
        analysis = linear_regression(green_schools_gdf)
        print('\nr-value: ' + str(analysis['r_value']))
        print('Showing visualizations...')
        show_regression_model(workspace_path, analysis, static)
        show_NDVI_imagery(workspace_path, static)

if __name__ == '__main__':
    main()

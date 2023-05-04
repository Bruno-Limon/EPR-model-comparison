# %% [markdown]
# ## **<font color="#84f745">0.1 IMPORTS & NOTEBOOK SETUP</font>**

# %%
# In case of running on Google Colab
%%capture
!apt-get install -qq curl g++ make
!curl -L http://download.osgeo.org/libspatialindex/spatialindex-src-1.8.5.tar.gz | tar xz
import os
os.chdir('spatialindex-src-1.8.5')
!./configure
!make
!make install
!pip install rtree
!ldconfig
!pip install scikit-mobility

!pip install geovoronoi
# Specifying pandas' version to overcome issue with 'TrajDataFrame' object has no attribute '_crs' during preprocessing phase
!pip install pandas==1.5.3

# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import geopandas as gpd
import random

import skmob
from skmob.tessellation import tilers
from skmob.models.epr import DensityEPR, SpatialEPR, Ditras
from skmob.models.markov_diary_generator import MarkovDiaryGenerator
from skmob.preprocessing import filtering, compression, detection, clustering
from skmob.measures.individual import jump_lengths, radius_of_gyration, uncorrelated_entropy, number_of_locations, number_of_visits, location_frequency
from skmob.measures.collective import visits_per_location
from skmob.utils.plot import *
from skmob.data.load import load_dataset, list_datasets
from shapely.geometry import MultiPolygon, Polygon, Point
from shapely.ops import unary_union

import sklearn
from sklearn.metrics import mean_squared_error

from scipy.stats import gaussian_kde, iqr
from scipy.spatial import Voronoi,voronoi_plot_2d
from geovoronoi import voronoi_regions_from_coords, points_to_coords

# %%
# Setting up plot style
sns.set_context(font_scale=2, rc={"font.size":10,"axes.titlesize":16,"axes.labelsize":14})
sns.set_style("whitegrid", {'grid.linestyle': '--', 'alpha': 0.25})
sns.set_style({'font.family':'serif', 'font.serif':'Computer Modern'})

# %%
np.seterr(divide = 'ignore') 

# %% [markdown]
# ## **<font color="#84f745">0.2 FUNCTIONS</font>**

# %%
# Function to obtain the distributions Mi, i in {1...6}
def compute_measures(Traj):
    m1 = jump_lengths(Traj, show_progress = False, merge=True)
    m1 = [m for m in m1 if m >= 1]

    m2 = list(radius_of_gyration(Traj, show_progress = False)['radius_of_gyration']) 
    m2 = [m for m in m2 if m >= 1]

    m3 = list(uncorrelated_entropy(Traj, show_progress = False)['uncorrelated_entropy']) 
    m3 = [m for m in m3 if m > 0]

    m4 = list(number_of_locations(Traj, show_progress = False)['number_of_locations']) 
    m4 = [m for m in m4 if m > 0]

    m5 = list(visits_per_location(Traj)['n_visits'])
    m5 = [m for m in m5 if m > 0]

    m6 = location_frequency(Traj, as_ranks = True, show_progress = False)
    m6 = [m for m in m6 if m > 0]

    list_measures = [m1, m2, m3, m4, m5, m6]
    list_min = []
    list_max = []
    list_avg = []
    list_std = []

    # Stats for previously computed measures
    for m in list_measures:
        list_min.append(round(min(m),4)) 
        list_max.append(round(max(m),4))
        list_avg.append(round(np.mean(m),4))
        list_std.append(round(np.std(m),4))

    measures_df = pd.DataFrame(list(zip(list_min, list_max, list_avg, list_std)), columns = ['Min', 'Max', 'Avg', 'Std'], index = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6'])
    print(measures_df.T)

    return list_measures

# %%
# Function to obtain distributions
def get_kde(dict_measures):
    dict_kde = {}
    dict_x = {}
    dict_y = {}

    # Using Kernel Density Estimation with np.geomspace to get bins in the log scale
    for model, measures in dict_measures.items():
        dict_kde[model] = [gaussian_kde(measures[i]) for i in range(6)]
        dict_x[model] = [np.geomspace(min(measures[i]), max(measures[i]), 10) for i in range(6)]
        dict_y[model] = [dict_kde[model][i](dict_x[model][i]) for i in range(6)]

    return dict_x, dict_y

# %%
# Function to build each model given an area and tessellation
def build_model(area_tessellation, n_agents, n_individuals):
    real_measures = dict_real_measures[area_tessellation]

    # SEPR MODEL
    print('Building S-EPR model for {}'.format(dict_area_tess_names[area_tessellation]))
    sepr_model = SpatialEPR()
    sepr_tdf = sepr_model.generate(start_date = start_time, 
                                   end_date = end_time, 
                                   spatial_tessellation = dict_tessellations[area_tessellation], 
                                   n_agents = n_agents,
                                   show_progress = True, 
                                   random_state = 42)
    print('\nS-EPR measures:')
    sepr_measures = compute_measures(sepr_tdf)

    # DEPR MODEL
    print('\nBuilding D-EPR model for {}'.format(dict_area_tess_names[area_tessellation]))
    depr_model = DensityEPR()
    depr_tdf = depr_model.generate(start_date = start_time, 
                                   end_date = end_time, 
                                   spatial_tessellation = dict_weighted_tess[area_tessellation], 
                                   relevance_column = 'relevance', 
                                   n_agents = n_agents,
                                   show_progress = True, 
                                   random_state = 42)   
    print('\nD-EPR measures:')
    depr_measures = compute_measures(depr_tdf)

    # DITRAS MODEL
    print('\nBuilding Ditras model for {}'.format(dict_area_tess_names[area_tessellation]))
    ctdf = compression.compress(dict_real_tdfs[area_tessellation])
    stdf = detection.stay_locations(ctdf)
    cstdf = clustering.cluster(stdf)

    mdg = MarkovDiaryGenerator()
    mdg.fit(cstdf, number_individuals, lid = 'cluster')

    ditras_model = Ditras(mdg)

    ditras_tdf = ditras_model.generate(start_date = start_time, 
                                       end_date = end_time, 
                                       spatial_tessellation = dict_weighted_tess[area_tessellation], 
                                       relevance_column = 'relevance', 
                                       n_agents = n_agents,
                                       od_matrix = None,
                                       show_progress = True, 
                                       random_state = 42)
    print('\nDitras measures:')
    ditras_measures = compute_measures(ditras_tdf)

    return real_measures, sepr_measures, depr_measures, ditras_measures

# %%
# Building models and plotting the results  
def plot_comparison(area_tessellation, n_agents, show_real):
    if area_tessellation in dict_area_tess_names.keys():
        real_measures, sepr_measures, depr_measures, ditras_measures = build_model(area_tessellation, n_agents, n_individuals = number_individuals)
    else:
        return print("Please input a valid area and tessellation pair, e.g. 'a1_t1' to work on NY State with squared tessellation")

    list_measures = ["Travel Distance (Jump Length)", 
                     "Radius of Gyration", 
                     "Uncorrelated Entropy", 
                     "Distinct visited Locations", 
                     "Visits per Location", 
                     "Location Frequency"]
    list_labels = ["Δr", "r_g", "S^{unc}", "N_u", "V_l", "L_i"]

    dict_models_measures = {"real": real_measures, "sepr": sepr_measures, "depr": depr_measures, "ditras": ditras_measures}
    
    # Getting x and y for each distribution
    dict_x, dict_y = get_kde(dict_models_measures)
    x_real, x_sepr, x_depr, x_ditras = [dict_x[m] for m in dict_models_measures.keys()]
    y_real, y_sepr, y_depr, y_ditras = [dict_y[m] for m in dict_models_measures.keys()]

    print('\n---------------------------------------- MEASURES DISTRIBUTIONS PLOT ----------------------------------------\n')
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(18, 9))

    k = 0
    for i in range(2):
        for j in range(3):
            # Calculating RMSE for each possible pair of the distributions
            rmse_sepr_depr = round(mean_squared_error(y_sepr[k], y_depr[k], squared = False),4)
            rmse_sepr_ditras = round(mean_squared_error(y_sepr[k], y_ditras[k], squared = False),4)
            rmse_depr_ditras = round(mean_squared_error(y_depr[k], y_ditras[k], squared = False),4)
            rmse_real_sepr = round(mean_squared_error(y_real[k], y_sepr[k], squared = False),4)
            rmse_real_depr = round(mean_squared_error(y_real[k], y_depr[k], squared = False),4)
            rmse_real_ditras = round(mean_squared_error(y_real[k], y_ditras[k], squared = False),4)

            ax[i,j].plot(x_sepr[k], y_sepr[k], marker = 's', linestyle = 'dotted', linewidth = 2, markersize = 5, label = 'S-EPR')
            ax[i,j].plot(x_depr[k], y_depr[k], marker = 'o', linestyle = 'dotted', linewidth = 2, markersize = 5, label = 'D-EPR')
            ax[i,j].plot(x_ditras[k], y_ditras[k], marker = '^', linestyle = 'dotted', linewidth = 2, markersize = 5, label = 'Ditras')
            ax[i,j].loglog()
            ax[i,j].set(title = list_measures[k])
            ax[i,j].set(ylabel = "$P({})$".format(list_labels[k]))
            ax[i,j].set(xlabel = "${}$".format(list_labels[k]))
            
            # Style for text box showing RMSE inside plot
            props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.3)

            # Showing the RMSE results, depending if the real measures are taken into account
            if show_real == True: 
                ax[i,j].plot(x_real[k], y_real[k], marker = 'D', linestyle = 'dotted', linewidth = 2, markersize = 5, color = 'black', label = 'Real')
                ax[i,j].text(-.6, .3, fontsize = 12, transform = ax[i,j].transAxes, bbox = props, s = ('RMSE M{} ${}$\n'
                                                                                                       '\n'
                                                                                                       'R/S:{}\n'
                                                                                                       'R/D:{}\n'
                                                                                                       'R/Di:{}\n'
                                                                                                       '\n'
                                                                                                       'S/D:{}\n'
                                                                                                       'S/Di:{}\n'
                                                                                                       'D/Di:{}').format(str(k+1),list_labels[k], 
                                                                                                                         rmse_real_sepr, 
                                                                                                                         rmse_real_depr, 
                                                                                                                         rmse_real_ditras,
                                                                                                                         rmse_sepr_depr,
                                                                                                                         rmse_sepr_ditras,
                                                                                                                         rmse_depr_ditras)) 
            else:
                ax[i,j].text(-.6, .4, fontsize = 12, transform = ax[i,j].transAxes, bbox = props, s = ('RMSE M{} ${}$\n''\nS/D:{}\n''S/Di:{}\n''D/Di:{}').format(str(k+1),
                                                                                                                                                                 list_labels[k],
                                                                                                                                                                 rmse_sepr_depr,
                                                                                                                                                                 rmse_sepr_ditras,
                                                                                                                                                                 rmse_depr_ditras)) 
            ax[i,j].legend()
            fig.tight_layout()
            k += 1

# %%
# Function to produce random points inside a polygon
def polygon_random_points(poly, num_points):
    min_x, min_y, max_x, max_y = poly.bounds
    points = []
    while len(points) < num_points:
        random_point = Point([random.uniform(min_x, max_x), random.uniform(min_y, max_y)])
        if (random_point.within(poly)):
            points.append([random_point.x, random_point.y])
    return np.array(points)

# Function to turn to geodataframe, used in get_voronoi_tessellation
def to_GeoDataFrame(region_polys):
    name=[]
    for i in range(1, len(region_polys) + 1):
        name.append('cell ' + str(i))
    gdf = gpd.GeoDataFrame(columns=['name','geometry'], crs={'init': 'epsg:4326'})
    gdf['name'] = name
    for index, row in gdf.iterrows():
        gdf.at[index, 'geometry'] = region_polys[index]
    return gdf

# Function to obtain voronoi tessellation from a set of points inside a polygon
def get_voronoi_tessellation(poly_ch, points):
    vor = Voronoi(points, qhull_options='Qbb Qc Qx')
    region_polys, region_pts = voronoi_regions_from_coords(points, poly_ch)
    tess_voronoi = to_GeoDataFrame(region_polys)
    return tess_voronoi

# %% [markdown]
# # **<font color="#ffb94f">1.0 PREPARING AREAS & TESSELLATIONS</font>**

# %% [markdown]
# The geographical areas to be considered are:
# 
# 
# 
# 1.   New York State
# 2.   San Francisco, California
# 3.   Austin, Texas
# 4.   Mexico City
# 
# 

# %% [markdown]
# The tessellations to be used are:
# 
# 1.   Squared
# 2.   Hexagonal
# 3.   Official
# 4.   Voronoi
# 
# Paired with the previously stated geographical areas, the possible combination of areas and tessellations can be identified by the form "a{}_t{}", with the blank spaces representing the number of area and tessellation respectively. For example, the hexagonal tessellation of Austin Texas shall be identified by "a3_t2". This form will be used throught the entirety of this notebook.
# 
# 

# %%
# Official administrative divisions to use as tessellations for each geographical area

# New York Counties 2011, provided within scikit mobility
url_area1 = 'https://raw.githubusercontent.com/scikit-mobility/tutorials/master/mda_masterbd2020/data/NY_counties_2011.geojson'
# San Francisco Neighborhoods Boundaries https://data.sfgov.org/Geographic-Locations-and-Boundaries/Analysis-Neighborhoods/p5b7-5n3h
url_area2 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/georef-us-SF_boundaries.geojson'
# Austin Open Data Texas - https://data.austintexas.gov/Locations-and-Maps/Counties/9pr5-nzce
url_area3 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/georef-us-austin-counties.geojson'
# División política municipal/Demarcación territorial - Vintage/Millésimé - Mexico https://data.opendatasoft.com/
url_area4 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/alcaldias_cdmx.zip'

# To obtain the basic shape and outline bounds for each area, some preprocessing steps were needed
# In the case of area 1 and 2 they were obtained using osmnx library and then imported to https://mapshaper.org/ to fix self-intersections
# For area 3 and 4 a unary union of the administrative division was necessary

area1_merged_url = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/georef-ny-flat.geojson'
area2_merged_url = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/georef-us-SF_flat.geojson'

area1_merged = gpd.read_file(area1_merged_url)
area2_merged = gpd.read_file(area2_merged_url)
area3_merged = gpd.read_file(url_area3) 
area3_merged = gpd.GeoSeries(unary_union(area3_merged['geometry']))
area4_merged = gpd.read_file(url_area4) 
area4_merged = gpd.GeoSeries(unary_union(area4_merged['geometry']))

# %%
# Creating squared, hexagonal and official tessellation geodataframes

%%capture
list_tessellations = ['squared', 'h3_tessellation', 'official']
list_official_tessellations = [url_area1, url_area2, url_area3, url_area4]
list_area_meter_squared = [40000, 1500, 20000, 6000]
list_area_meter_h3 = [20000, 500, 10000, 3000]

dict_tessellations = {}

for i, area in enumerate(list_official_tessellations, start = 1):
    for j, tessellation in enumerate(list_tessellations, start = 1):
        # Using the official tessellations as base shape
        base_gpd = gpd.read_file(area)
        base_gpd = base_gpd.explode()
        if j == 1:
            dict_tessellations['a{}_t{}'.format(i,j)] = tilers.tiler.get(tessellation, base_shape = base_gpd, meters = list_area_meter_squared[i-1])
            # Exploding geometries to end up with single polygons instead of multipolygons, which have issues later on with the skmob.tdf.mapping function
            # dict_tessellations['a{}_t{}'.format(i,j)] = dict_tessellations['a{}_t{}'.format(i,j)].explode()
        elif j == 2:
            dict_tessellations['a{}_t{}'.format(i,j)] = tilers.tiler.get(tessellation, base_shape = dict_tessellations['a{}_t{}'.format(i,j-1)], meters = list_area_meter_h3[i-1])
        elif j == 3:
            dict_tessellations['a{}_t{}'.format(i,j)] = gpd.read_file(list_official_tessellations[i-1])
            dict_tessellations['a{}_t{}'.format(i,j)] = dict_tessellations['a{}_t{}'.format(i,j)].explode()
            # Adding a tile_ID column, necesarry for the mapping function
            dict_tessellations['a{}_t{}'.format(i,j)].insert(0, 'tile_ID', range(0, len(dict_tessellations['a{}_t{}'.format(i,j)])))

        dict_tessellations['a{}_t{}'.format(i,j)] = dict_tessellations['a{}_t{}'.format(i,j)].set_crs({'init': 'epsg:4326'}, allow_override=True)

# %% [markdown]
# Voronoi Tessellations

# %%
# Creating polygon of city boundaries for each area using shapely's function unary_union
area1_boundary = unary_union(area1_merged.geometry)
area2_boundary = unary_union(area2_merged.geometry)
area3_boundary = unary_union(area3_merged.geometry)
area4_boundary = unary_union(area4_merged.geometry)

# %%
# Points to build voronoi tessellation for area 1
# Consisting of DMV office locaitons in NY state, obtained from https://data.ny.gov/Transportation/Department-of-Motor-Vehicle-DMV-Office-Locations/9upz-c7xg
url_points1 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/Points-NY.csv'
points1_df = pd.read_csv(url_points1)
to_drop = ['Office Name','Office Type','Public Phone Number','Public Phone Extension','Street Address Line 1','Street Address Line 2','Monday Beginning Hours','Monday Ending Hours',
            'Tuesday Beginning Hours','Tuesday Beginning Hours','Tuesday Ending Hours','Wednesday Beginning Hours','Wednesday Ending Hours','Thursday Beginning Hours','Thursday Ending Hours',
            'Friday Beginning Hours','Friday Ending Hours','Saturday Beginning Hours','Saturday Ending Hours']
points1_df.drop(to_drop, axis = 1, inplace = True)
points1_df = points1_df.dropna()

list_points = points1_df['Georeference'].tolist()
lat = []
lng = []
voronoi_points1 = []
for point in list_points:
    x = point.replace('POINT (','')
    x = x.replace(')','')
    latitude = x.split()[1]
    longitude = x.split()[0]
    lat.append(latitude)
    lng.append(longitude)

points1_df['lat'] = lat
points1_df['lng'] = lng
points1_df.drop('Georeference', axis = 1, inplace = True)
points1_df
points1_gdf = gpd.GeoDataFrame(points1_df, geometry=gpd.points_from_xy(points1_df.lng, points1_df.lat))
points1_gdf.crs = "EPSG:4326"

for lat, lng in zip(lat, lng):
    voronoi_points1.append([lng, lat])
voronoi_points1 = np.array(voronoi_points1)

# %%
# Points to build voronoi tessellation for area 2
# The dataset to use is obtained from https://data.sfgov.org/Geographic-Locations-and-Boundaries/ and represents healthcare facilities in san francisco
url_points2 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/Points-SF.csv'
points2_df = pd.read_csv(url_points2)
to_drop = ['OSHPD_ID','Facility Name','Facility Type','Services','Neighborhoods','SF Find Neighborhoods','Current Police Districts','Current Supervisor Districts','Analysis Neighborhoods']
points2_df.drop(to_drop, axis = 1, inplace = True)

list_points = points2_df['point'].tolist()
lat = []
lng = []
voronoi_points2 = []
for point in list_points:
    x = point.replace('POINT (','')
    x = x.replace(')','')
    latitude = x.split()[1]
    longitude = x.split()[0]
    lat.append(latitude)
    lng.append(longitude)

points2_df['lat'] = lat
points2_df['lng'] = lng
points2_df.drop('point', axis = 1, inplace = True)
points2_df
points2_gdf = gpd.GeoDataFrame(points2_df, geometry=gpd.points_from_xy(points2_df.lng, points2_df.lat))
points2_gdf.crs = "EPSG:4326"

for lat, lng in zip(lat, lng):
    voronoi_points2.append([lng, lat])
voronoi_points2 = np.array(voronoi_points2)

# %%
# Points to build voronoi tessellation for area 3
# In the case of area 3, a random approach to produce the points is used to see how it differs from actual points
voronoi_points3 = polygon_random_points(area3_boundary, 15)
points3_gdf = gpd.GeoDataFrame(geometry = gpd.points_from_xy(voronoi_points3[:,0], voronoi_points3[:,1]))

# %%
# Points to build voronoi tessellation for area 4
# The dataset to use is obtained from https://datos.cdmx.gob.mx/dataset/hospitales-y-centros-de-salud and represents healthcare facilities in mexico city
url_points4 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/Points-MX.csv'
points4_df = pd.read_csv(url_points4)
to_drop = ['direccion','colonia','alcaldia','titular','tel']
points4_df.drop(to_drop, axis = 1, inplace = True)

points4_gdf = gpd.GeoDataFrame(points4_df, geometry=gpd.points_from_xy(points4_df.longitud, points4_df.latitud))
points4_gdf.crs = "EPSG:4326"

voronoi_points4 = []
for lat, lng in zip(points4_df.latitud.to_list(), points4_df.longitud.to_list()):
    voronoi_points4.append([lng, lat])
voronoi_points4 = np.array(voronoi_points4)

# %%
# Producing Voronoi tessellations
%%capture
list_boundaries = [area1_boundary, area2_boundary, area3_boundary, area4_boundary]
list_voronoi_points = [voronoi_points1, voronoi_points2, voronoi_points3, voronoi_points4]

for i in range(1,len(list_official_tessellations)+1):
        dict_tessellations['a{}_t4'.format(i)] = get_voronoi_tessellation(list_boundaries[i-1], list_voronoi_points[i-1])
        dict_tessellations['a{}_t4'.format(i)] = dict_tessellations['a{}_t4'.format(i)].explode()
        dict_tessellations['a{}_t4'.format(i)]['tile_ID'] = [str(j) for j in range(1, len(dict_tessellations['a{}_t4'.format(i)])+1)]
        dict_tessellations['a{}_t4'.format(i)] = dict_tessellations['a{}_t4'.format(i)].set_crs({'init': 'epsg:4326'}, allow_override=True)

# %%
# Looking at number of tiles in each tessellation
for i in range(1,len(list_official_tessellations)+1):
    for j in range(1,len(list_official_tessellations)+1):
        print('Tiles in tessellation a{}_t{}: {}'.format(i, j, len(dict_tessellations['a{}_t{}'.format(i, j)])))   
        if j == 4 and i != 4: print('') 

# %% [markdown]
# Visualizing tessellations

# %%
background_color = 'gray'
tess_color = "tab20b"

fig, ax = plt.subplots(1,4, figsize=(32, 6)) 
fig.suptitle('Basic shapes of geographical areas', fontsize = 28)
np.vectorize(lambda ax:ax.axis('off'))(ax)

area1_merged.plot(ax = ax[0], color = background_color)
area2_merged.plot(ax = ax[1], color = background_color)
area3_merged.plot(ax = ax[2], color = background_color)
area4_merged.plot(ax = ax[3], color = background_color)

# %%
fig, ax = plt.subplots(1,4, figsize=(32, 6)) 
fig.suptitle('Squared tessellations', fontsize = 28)
np.vectorize(lambda ax:ax.axis('off'))(ax)

area1_merged.plot(ax = ax[0], color = background_color)
area2_merged.plot(ax = ax[1], color = background_color)
area3_merged.plot(ax = ax[2], color = background_color)
area4_merged.plot(ax = ax[3], color = background_color)

dict_tessellations['a1_t1'].plot(ax = ax[0], cmap = tess_color, alpha = .5)
dict_tessellations['a2_t1'].plot(ax = ax[1], cmap = tess_color, alpha = .5)
dict_tessellations['a3_t1'].plot(ax = ax[2], cmap = tess_color, alpha = .5)
dict_tessellations['a4_t1'].plot(ax = ax[3], cmap = tess_color, alpha = .5)

# %%
fig, ax = plt.subplots(1,4, figsize=(32, 6)) 
fig.suptitle('Hexagonal tessellations', fontsize = 28)
np.vectorize(lambda ax:ax.axis('off'))(ax)

area1_merged.plot(ax = ax[0], color = background_color)
area2_merged.plot(ax = ax[1], color = background_color)
area3_merged.plot(ax = ax[2], color = background_color)
area4_merged.plot(ax = ax[3], color = background_color)

dict_tessellations['a1_t2'].plot(ax = ax[0], cmap = tess_color, alpha = .5)
dict_tessellations['a2_t2'].plot(ax = ax[1], cmap = tess_color, alpha = .5)
dict_tessellations['a3_t2'].plot(ax = ax[2], cmap = tess_color, alpha = .5)
dict_tessellations['a4_t2'].plot(ax = ax[3], cmap = tess_color, alpha = .5)

# %%
fig, ax = plt.subplots(1,4, figsize=(32, 6)) 
fig.suptitle('Official tessellations', fontsize = 28)
np.vectorize(lambda ax:ax.axis('off'))(ax)

dict_tessellations['a1_t3'].plot(ax = ax[0], cmap = tess_color)
dict_tessellations['a2_t3'].plot(ax = ax[1], cmap = tess_color)
dict_tessellations['a3_t3'].plot(ax = ax[2], cmap = tess_color)
dict_tessellations['a4_t3'].plot(ax = ax[3], cmap = tess_color)

# %%
fig, ax = plt.subplots(1,4, figsize=(32, 6)) 
fig.suptitle('Voronoi tessellations', fontsize = 28)
np.vectorize(lambda ax:ax.axis('off'))(ax)

vor_point_color = 'white'
markersize = 10

dict_tessellations['a1_t4'].plot(ax = ax[0], cmap = tess_color)
dict_tessellations['a2_t4'].plot(ax = ax[1], cmap = tess_color)
dict_tessellations['a3_t4'].plot(ax = ax[2], cmap = tess_color)
dict_tessellations['a4_t4'].plot(ax = ax[3], cmap = tess_color)

points1_gdf.plot(ax = ax[0], color = vor_point_color, markersize = markersize)
points2_gdf.plot(ax = ax[1], color = vor_point_color, markersize = markersize)
points3_gdf.plot(ax = ax[2], color = vor_point_color, markersize = markersize)
points4_gdf.plot(ax = ax[3], color = vor_point_color, markersize = markersize)

# %% [markdown]
# Creating TDFs for each area
# 
# Using the following mobility traces datasets:
# 
# 1. Brightkite check-ins
# 2. San Francisco Taxicab traces
# 3. Gowalla check-ins
# 4. Mexico City Taxicab traces

# %%
# Creating trajectory dataframe for Area 1 - New York State
# Trajectories obtained from brightkite checkins dataset
url = "https://snap.stanford.edu/data/loc-brightkite_totalCheckins.txt.gz"
area1_df = pd.read_csv(url, sep='\t', header=0, nrows = 1000000, names=['user', 'check-in_time', "latitude", "longitude", "location id"])
area1_tdf = skmob.TrajDataFrame(area1_df, latitude='latitude', longitude='longitude', datetime='check-in_time', user_id='user')
print("Users: ", len(area1_tdf.uid.unique()))
print('Records: ', len(area1_tdf)) 
area1_tdf.head()

# %%
# Visualizing trajectories flowing through Area 1
map_f = folium.Map(location=[42.631610, -73.935242], zoom_start=7, tiles='cartodbdark_matter')

# Outlining official tessellation
for _, row in dict_tessellations['a1_t3'].iterrows():
    sim_geo = gpd.GeoSeries(row['geometry']).simplify(tolerance=0.001)
    geo_j = sim_geo.to_json()
    geo_j = folium.GeoJson(data=geo_j, style_function = lambda x: {'fillColor': '#98ff69', 'color': 'white'})
    geo_j.add_to(map_f)
    
area1_tdf.plot_trajectory(map_f = map_f, max_users = 20, weight = 1, max_points = None, opacity=0.9, start_end_markers = False)

# %%
# Creating trajectory dataframe for Area 2 -  San Francisco
# Trajectories obtained from taxicab dataset, retrieved from Crawdad repository (http://crawdad.org/epfl/mobility/20090224/index.html), using a sample of 6 taxis
url = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/new_abboip.txt'
url1 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/new_effneomi.txt'
url2 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/new_iawxben.txt'
url3 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/new_ofikco.txt'
url4 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/new_udwadla.txt'
url5 = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/new_upthin.txt'
list_url = [url, url1, url2, url3, url4, url5]

dataframes = []
for i, url in enumerate(list_url):
    df = pd.read_csv(url, sep=' ', names=['lat', 'lng', "occupancy", "datetime"])
    df['datetime'] = pd.to_datetime(df['datetime'], unit='s')
    df['uid'] = i+1
    dataframes.append(df)

area2_df = pd.concat(dataframes)
area2_tdf = skmob.TrajDataFrame(area2_df, latitude='lat', longitude='lng', datetime='datetime', user_id='uid')
print("Users: ", len(area2_tdf.uid.unique()))
print('Records: ', len(area2_tdf))
area2_tdf.head()

# %%
# Visualizing trajectories flowing through Area 2
map_f = folium.Map(location=[37.75134, -122.39488], zoom_start=12, tiles='cartodbdark_matter')

# Outlining official tessellation and adding name of neighborghoods as folium popup
for _, row in dict_tessellations['a2_t3'].iterrows():
    sim_geo = gpd.GeoSeries(row['geometry']).simplify(tolerance=0.001)
    geo_j = sim_geo.to_json()
    geo_j = folium.GeoJson(data=geo_j, style_function = lambda x: {'fillColor': '#98ff69', 'color': 'white'})
    folium.Popup(row['nhood']).add_to(geo_j)
    geo_j.add_to(map_f)
    
area2_tdf.plot_trajectory(map_f = map_f, hex_color = '#59fff7', max_users = 1, zoom = 12, weight = 1, max_points = 200, opacity=0.9, start_end_markers = False)

# %%
# Creating trajectory dataframe for Area 3 - Austin
# Trajectories obtained from Gowalla checkins dataset
url = "https://snap.stanford.edu/data/loc-gowalla_totalCheckins.txt.gz"
area3_df = pd.read_csv(url, sep='\t', header=0, nrows = 1000000, names=['user', 'check-in_time', "latitude", "longitude", "location id"])
area3_tdf = skmob.TrajDataFrame(area3_df, latitude='latitude', longitude='longitude', datetime='check-in_time', user_id='user')
print("Users: ", len(area3_tdf.uid.unique()))
print('Records: ', len(area3_tdf))
area3_tdf.head()

# %%
# Visualizing trajectories flowing through Area 3
map_f = folium.Map(location=[30.266666, -97.733330], zoom_start=8, tiles='cartodbdark_matter')

# Outlining official tessellation and adding name of counties as folium popup
for _, row in dict_tessellations['a3_t3'].iterrows():
    sim_geo = gpd.GeoSeries(row['geometry']).simplify(tolerance=0.001)
    geo_j = sim_geo.to_json()
    geo_j = folium.GeoJson(data=geo_j, style_function = lambda x: {'fillColor': '#98ff69', 'color': 'white'})
    folium.Popup(row['county_name']).add_to(geo_j)
    geo_j.add_to(map_f)
    
area3_tdf.plot_trajectory(map_f = map_f, max_users = 50, weight = 1, max_points = None, opacity=0.9, start_end_markers = False)

# %%
# Creating trajectory dataframe for Area 4 - Mexico City
# Trajectories obtained from Kaggle Taxi traces dataset. (https://www.kaggle.com/datasets/mnavas/taxi-routes-for-mexico-city-and-quito?select=mex_clean.csv)
url = 'https://raw.githubusercontent.com/Bruno-Limon/EPR-model-comparison/main/Data/mex_clean.csv'
area4_df = pd.read_csv(url)
mapping = {item:i for i, item in enumerate(area4_df["vendor_id"].unique())}
area4_df["uid"] = area4_df["vendor_id"].apply(lambda x: mapping[x])
area4_df.drop(['vendor_id', 'id', 'dropoff_datetime', 'dropoff_longitude', 'dropoff_latitude', 'store_and_fwd_flag', 'trip_duration', 'dist_meters', 'wait_sec'], axis=1, inplace=True)

area4_tdf = skmob.TrajDataFrame(area4_df, latitude='pickup_latitude', longitude='pickup_longitude', datetime='pickup_datetime', user_id='uid')
print("Users: ", len(area4_tdf.uid.unique()))
print('Records: ', len(area4_tdf))
area4_tdf.head()

# %%
# Visualizing trajectories flowing through Area 4
map_f = folium.Map(location=[19.2016874, -99.097369], zoom_start=10, tiles='cartodbdark_matter')

# Outlining official tessellation and adding name of neighborghoods as folium popup
for _, row in dict_tessellations['a4_t3'].iterrows():
    sim_geo = gpd.GeoSeries(row['geometry']).simplify(tolerance=0.001)
    geo_j = sim_geo.to_json()
    geo_j = folium.GeoJson(data=geo_j, style_function = lambda x: {'fillColor': '#98ff69', 'color': 'white'})
    folium.Popup(row['nomgeo']).add_to(geo_j)
    geo_j.add_to(map_f)
    
area4_tdf.plot_trajectory(map_f = map_f, max_users = 1, hex_color = '#59fff7', weight = 1, max_points = 200, opacity = 0.9, start_end_markers = False)

# %% [markdown]
# TDFS for the real traces dataframes

# %%
%%capture
list_areas_tdfs = [area1_tdf, area2_tdf, area3_tdf, area4_tdf]
list_tessellations = ['squared', 'hexagonal', 'official', 'voronoi']
dict_real_tdfs = {}
dict_weighted_tess = {}

for i in range(1, len(list_areas_tdfs)+1):
    for j in range(1, len(list_tessellations)+1):
        # Mapping area traces with tessellation
        dict_real_tdfs['a{}_t{}'.format(i,j)] = list_areas_tdfs[i-1].mapping(dict_tessellations['a{}_t{}'.format(i,j)], remove_na=True)
        
        # Counting visits per tile to use as parameter for weighted tessellation
        visits_per_tile = dict_real_tdfs['a{}_t{}'.format(i,j)].groupby("tile_ID", as_index=False).count()
        visits_per_tile = visits_per_tile[["tile_ID", "uid"]]
        visits_per_tile["relevance"] = visits_per_tile["uid"]
        visits_per_tile = visits_per_tile[["tile_ID", "relevance"]]

        dict_weighted_tess['a{}_t{}'.format(i,j)] = dict_tessellations['a{}_t{}'.format(i,j)].set_index("tile_ID").join(visits_per_tile.set_index("tile_ID"))
        dict_weighted_tess['a{}_t{}'.format(i,j)] = dict_weighted_tess['a{}_t{}'.format(i,j)].fillna(0)
        dict_weighted_tess['a{}_t{}'.format(i,j)] = dict_weighted_tess['a{}_t{}'.format(i,j)][dict_weighted_tess['a{}_t{}'.format(i,j)]["relevance"]>0]

# %%
# Creating a dictionary to map area_tessellation pairs with their names
dict_area_tess_names = {}
area_names = ['New York State', 'San Francisco', 'Houston', 'Mexico City']

for i, area in enumerate(area_names, start = 1):
    for j, tess in enumerate(list_tessellations, start = 1):
        dict_area_tess_names['a{}_t{}'.format(i,j)] = str(area + ' - ' + tess + ' tessellation')

# %%
# Calculating measures for the real TDFs
dict_real_measures = {}
show_df = True

for i in range(1, len(list_areas_tdfs)+1):
    for j in range(1, len(list_tessellations)+1):
        if show_df == True: print('\na{}_t{} measures'.format(i,j))
        dict_real_measures['a{}_t{}'.format(i,j)] = compute_measures(dict_real_tdfs['a{}_t{}'.format(i,j)])

# %% [markdown]
# # **<font color="#ffb94f">2.0 BUILDING MODELS & COMPARING DISTRIBUTIONS</font>**

# %% [markdown]
# The comparison among distributions is made on an area_tessellation basis, computing 6 measures of model realism, they are as follows:
# 
# 1. Travel Distance
# 2. Radius of Gyration
# 3. Uncorrelated Entropy
# 4. Distinct Visited Locations
# 5. Visits per Location
# 6. Location Frequency
# 
# Then, the distribution for each of these measures across all models is plotted, together with the RMSE of their pairs

# %%
# Denoting if the real TDF is taken into account into the comparison plot
show_real = False

# Parameters to use in EPR models
number_agents = 5000
number_individuals = 5000 # Used for Ditras Markov Diaries Generator
start_time = pd.to_datetime('2023/01/01 00:00:00')
end_time = pd.to_datetime('2023/01/15 00:00:00')

# %% [markdown]
# ## **<font color="#84f745">2.1 AREA 1 - NEW YORK STATE</font>**

# %% [markdown]
# ### **<font color="#4fffd9">2.1.1 AREA 1 - TESSELLATION 1</font>**

# %%
plot_comparison(area_tessellation = 'a1_t1', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.1.2 AREA 1 - TESSELLATION 2</font>**

# %%
plot_comparison(area_tessellation = 'a1_t2', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.1.3 AREA 1 - TESSELLATION 3</font>**

# %%
plot_comparison(area_tessellation = 'a1_t3', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.1.4 AREA 1 - TESSELLATION 4</font>**

# %%
plot_comparison(area_tessellation = 'a1_t4', n_agents = number_agents, show_real = show_real)

# %%


# %% [markdown]
# ## **<font color="#84f745">2.2 AREA 2 - SAN FRANCISCO</font>**

# %% [markdown]
# ### **<font color="#4fffd9">2.2.1 AREA 2 - TESSELLATION 1</font>**

# %%
plot_comparison(area_tessellation = 'a2_t1', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.2.2 AREA 2 - TESSELLATION 2</font>**

# %%
plot_comparison(area_tessellation = 'a2_t2', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.2.3 AREA 2 - TESSELLATION 3</font>**

# %%
plot_comparison(area_tessellation = 'a2_t3', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.2.4 AREA 2 - TESSELLATION 4</font>**

# %%
plot_comparison(area_tessellation = 'a2_t4', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ## **<font color="#84f745">2.3 AREA 3 - HOUSTON</font>**

# %% [markdown]
# ### **<font color="#4fffd9">2.3.1 AREA 3 - TESSELLATION 1</font>**

# %%
plot_comparison(area_tessellation = 'a3_t1', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.3.2 AREA 3 - TESSELLATION 2</font>**

# %%
plot_comparison(area_tessellation = 'a3_t2', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.3.3 AREA 3 - TESSELLATION 3</font>**

# %%
plot_comparison(area_tessellation = 'a3_t3', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.3.4 AREA 3 - TESSELLATION 4</font>**

# %%
plot_comparison(area_tessellation = 'a3_t4', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ## **<font color="#84f745">2.4 AREA 4 - MEXICO CITY</font>**

# %% [markdown]
# ### **<font color="#4fffd9">2.4.1 AREA 4 - TESSELLATION 1</font>**

# %%
plot_comparison(area_tessellation = 'a4_t1', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.4.2 AREA 4 - TESSELLATION 2</font>**

# %%
plot_comparison(area_tessellation = 'a4_t2', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.4.3 AREA 4 - TESSELLATION 3</font>**

# %%
plot_comparison(area_tessellation = 'a4_t3', n_agents = number_agents, show_real = show_real)

# %% [markdown]
# ### **<font color="#4fffd9">2.4.4 AREA 4 - TESSELLATION 4</font>**

# %%
plot_comparison(area_tessellation = 'a4_t4', n_agents = number_agents, show_real = show_real)



import os
import pandas as pd
import geopandas as gpd
import re
import xml.etree.ElementTree as et 

absolute_path = os.path.dirname(__file__)
INTEGRITY = False # TODO: Eigenlijk wil ik dit wel er in, maar er zijn drie overstorten die dubbel voorkomen. Dat moet beheerd worden om fouten te voorkomen.
scenarios = [ # Turn on or off the scenarios that must be ran
    'HIST',
    'GHG120'
    ]

# Read and process data
## Load data of lateral points and areas of Wageningen model, and assign Wagmod area code to laterals based on intersection.
gdf_laterals = (
    gpd.read_file(absolute_path + '/laterals_shapefile/' + 'laterals.shp')
    .set_index('CODE', verify_integrity = True)
    [['Opp_hactar','X','Y','geometry']]
)
gdf_laterals['type'] = 'laterals'

gdf_wagmod = (
    gpd.read_file(absolute_path + '/wageningenmodel_shapefile/' + 'NAvakken.shp')
    [['GAFIDEN','geometry']]
)

gdf_overstorten = gpd.read_file(
    absolute_path + '/overstortpunten_shapefile/'+'overstortlocaties.shp'
    ).set_index('RIODAT_ZRO', 
    verify_integrity=INTEGRITY 
    )[['geometry']]
gdf_overstorten['type'] = 'overstort'

# Create the shapefile containing all locations (laterals and overstorten) that can be read by the following script. This is a component of the final output.
pd.concat([gdf_laterals,gdf_overstorten], 
    verify_integrity=INTEGRITY 
    )[['geometry','type']].to_file(os.path.abspath(os.path.join(absolute_path, '../07_Rapportage/shp_afvoer.shp')))
    
gdf_laterals_joined = gpd.sjoin(gdf_laterals, gdf_wagmod, how='left')
gdf_missing = gdf_laterals_joined.loc[gdf_laterals_joined['index_right'].isna()].to_file(absolute_path + '/laterals_shapefile/unmatchedLaterals.shp')
#TODO: locaties die niet binnen een Wagmod gebied vallen moeten op e.o.a. manier gekoppeld worden aan Wagmod data. 
# Mijn voorkeur zou zijn een koppeltabel voor die 100 locaties.

for scenario in scenarios:
    
    ## Read Wageningen model output (discharges) per region, and multiply by the area belonging to each lateral
    wagmodFiles = [absolute_path + '/wageningenmodel_output/' + scenario+ '/' + region for region in os.listdir(absolute_path + '/wageningenmodel_output/' + scenario)]
    ls_wagmodoutput = []
    
    for region in wagmodFiles:
        metadata_region = pd.read_csv(region, nrows=1, skiprows=[0,2], header=None)[0][0]
        model = re.search('model:(.*) (?:van|met winterbui)', metadata_region).group(1)
                
        df_region = pd.read_csv(
                region, 
                header=4, 
                delim_whitespace=True, 
                usecols=['-I-','-QC-'],
                )

        df_region=df_region.set_index(
                pd.to_datetime(df_region['-I-'],
                origin=pd.Timestamp('2010-11-08 23:00:00'),
                unit='h')
            ).drop(labels=['-I-'], axis='columns').rename(columns={'-QC-':model})
        
        ls_wagmodoutput.append(df_region)

    # Combine data from all Wageningen models, transpose to later multiply with the areas, and finally multiply by 0.001 to adjust for units.
    df_wagmodoutput = pd.concat(ls_wagmodoutput, axis=1, verify_integrity=True).T * 0.001 

    ### Multiplication
    output_laterals = gdf_laterals_joined.join(df_wagmodoutput, on='GAFIDEN', how='left')
    output_laterals = output_laterals.iloc[:,7:].multiply(output_laterals['Opp_hactar'],axis='index').T

    ## Read overstort schatter output (xml files) 
    NS = '{http://www.wldelft.nl/fews/PI}' # xml default namespace

    #overstortFiles = os.listdir(absolute_path + '/overstortschatter_output/')
    overstortFiles = [absolute_path + '/overstortschatter_output/' + scenario+ '/' + region for region in os.listdir(absolute_path + '/overstortschatter_output/' + scenario)]
    ls_multipleyears = []

    for overstort in overstortFiles:
        xtree = et.parse(overstort)
        xroot = xtree.getroot()
        ls_overstort = []

        for i in xroot.iter(NS+'series'):
            locId = i.find('./'+NS+'header/'+NS+'locationId').text
            events = i.findall('./'+NS+'event')
            ls=[item.attrib for item in events]
            df=pd.DataFrame(ls).head()
            df_resetIndex=(
                df.set_index(
                pd.DatetimeIndex(pd.to_datetime(df['date']+' '+df['time'], format="%Y-%m-%d %H:%M:%S")),
                verify_integrity=True
                ).drop(labels=['date','time'], axis='columns')
                .astype('float')
                .rename({'value':locId}, axis='columns')
            )
            ls_overstort.append(df_resetIndex)
        pd_oneyear=pd.concat(ls_overstort, axis='columns',verify_integrity=True)
        ls_multipleyears.append(pd_oneyear)
    df_overstorten = pd.concat(ls_multipleyears, axis='index', verify_integrity=True) / 3600

    df_combined = pd.concat([output_laterals,df_overstorten], axis='columns').fillna(value=0)
    df_combined.to_csv(os.path.abspath(os.path.join(absolute_path, '../07_Rapportage/'+'afvoer_' + scenario +'.csv')), sep = ';')
   
    ## RWZI output # Left out for now, because there are no relevant RWZIs in the pilot area


    # Format output


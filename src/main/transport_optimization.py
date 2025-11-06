"""
@authors:
 - Claire Halloran
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk
Contains code originally written by Leander MÃ¼ller, RWTH Aachen University

Calculates the cost-optimal commodity transportation strategy to the nearest 
demand center.
Calculates cost of pipeline transport and demand profile based on optimal size.
"""
import geopandas as gpd
import geopy.distance
import numpy as np
import pandas as pd
from shapely.geometry import Point
import shapely.wkt

from functions import CRF, cheapest_trucking_strategy, h2_conversion_stand, \
                            cheapest_pipeline_strategy, calculate_trucking_costs, \
                            calculate_nh3_pipeline_costs
from utils import check_folder_exists

def main():
    plant_type = snakemake.wildcards.plant_type
    tech_params_filepath = snakemake.input.technology_parameters
    demand_params_filepath = snakemake.input.demand_parameters
    country_params_filepath = snakemake.input.country_parameters
    transport_params_filepath = snakemake.input.transport_parameters
    pipeline_params_filepath = snakemake.input.pipeline_parameters
    if plant_type == "hydrogen":
        conversion_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/conversion_parameters.xlsx'

    infra_data = pd.read_excel(tech_params_filepath,
                           sheet_name='Infra',
                           index_col='Infrastructure')
    
    demand_center_list = pd.read_excel(demand_params_filepath,
                                    sheet_name='Demand centers',
                                    index_col='Demand center',)
    demand_centers = demand_center_list.index

    country_params = pd.read_excel(country_params_filepath,
                                        index_col='Country')

    hexagons = gpd.read_file(snakemake.input.hexagons)
    if hexagons.crs is None:
        hexagons = hexagons.set_crs(epsg=4326)
    elif hexagons.crs.to_epsg() != 4326:
        hexagons = hexagons.to_crs(epsg=4326)
    len_hexagons = len(hexagons)

    needs_pipeline_construction = snakemake.config["transport"]["pipeline_construction"]
    needs_road_construction = snakemake.config["transport"]["road_construction"]

    long_road_capex = infra_data.at['Long road','CAPEX']
    short_road_capex = infra_data.at['Short road','CAPEX']
    road_opex = infra_data.at['Short road','OPEX']

    # Prices from the country excel file
    currency = snakemake.config["currency"]
    elec_price = country_params[f'Electricity price ({currency}/kWh)'].iloc[0]
    heat_price = country_params[f'Heat price ({currency}/kWh)'].iloc[0]
    plant_interest_rate = country_params['Plant interest rate'].iloc[0]
    infrastructure_interest_rate = country_params['Infrastructure interest rate'].iloc[0]
    infrastructure_lifetime = country_params['Infrastructure lifetime (years)'].iloc[0]

    check_folder_exists("resources")
    
    # Calculate cost of hydrogen state conversion and transportation for demand
    # Loop through all demand centers
    for demand_center in demand_centers:
        print(f"\nOptimisation for {demand_center} begins...\n")
        # Demand location based variables
        demand_center_lat = demand_center_list.loc[demand_center,'Lat [deg]']
        demand_center_lon = demand_center_list.loc[demand_center,'Lon [deg]']
        demand_location = Point(demand_center_lon, demand_center_lat)
        annual_demand_quantity = demand_center_list.loc[demand_center,'Annual demand [kg/a]']
        demand_state = demand_center_list.loc[demand_center,'Demand state']

        # Storage arrays containing zeros to be filled with costs and states or
        # left as zeros, as necessary
        road_construction_costs = np.zeros(len_hexagons)
        trucking_states = np.zeros(len_hexagons,dtype='<U10')
        trucking_costs = np.zeros(len_hexagons)
        pipeline_costs = np.zeros(len_hexagons)

        if demand_state not in ['500 bar', 'LH2', 'NH3', 'LOHC']:
            raise NotImplementedError(f'{demand_state} demand not supported.')
        
        # Loop through all hexagons
        for i in range(len_hexagons):
            print(f"Currently optimising {i+1} of {len_hexagons} hexagons...")
            dist_to_road = hexagons['road_dist'][i]
            hex_geometry = hexagons['geometry'][i]
            
            # Calculating distance to demand
            poly = shapely.wkt.loads(str(hex_geometry))
            center = poly.centroid
            demand_coords = (demand_center_lat, demand_center_lon)
            hexagon_coords = (center.y, center.x)
            dist_to_demand = geopy.distance.geodesic(demand_coords, hexagon_coords).km

            #!!! maybe this is the place to set a restriction based on distance to demand center-- for all hexagons with a distance below some cutoff point
            # label demand location under consideration
            # Different calculations dependent if in demand location or not

            # If the hexagon contains the demand location:
            if hex_geometry.contains(demand_location) == True:
                # Calculate cost of converting hydrogen to a demand state for local demand (i.e. no transport)
                # State changed to "None" to label the hexagon that contains the demand center
                # Leave road construction cost as 0, as no road is needed
                if plant_type == "hydrogen":
                    if demand_state == 'NH3':
                        trucking_costs[i]=pipeline_costs[i]=h2_conversion_stand(demand_state+'_load',
                                                annual_demand_quantity,
                                                elec_price,
                                                heat_price,
                                                plant_interest_rate,
                                                conversion_params_filepath,
                                                currency
                                                )[2]/annual_demand_quantity
                        trucking_states[i] = "None"
                    else:
                        trucking_costs[i]=pipeline_costs[i]=h2_conversion_stand(demand_state,
                                                annual_demand_quantity,
                                                elec_price,
                                                heat_price,
                                                plant_interest_rate,
                                                conversion_params_filepath,
                                                currency
                                                )[2]/annual_demand_quantity
                        trucking_states[i] = "None"
                # Leave the costs as 0 for ammonia, just change state
                elif plant_type == "ammonia":
                    trucking_states[i] = "None"
            # Else, if the hexagon does not contain the demand center:
            else:
                # Calculate costs of constructing a pipeline to the hexagon if allowed
                if needs_pipeline_construction== True:
                    if plant_type == "hydrogen":
                        pipeline_costs[i] =\
                            cheapest_pipeline_strategy(demand_state,
                                                    annual_demand_quantity,
                                                    dist_to_demand,
                                                    elec_price,
                                                    heat_price,
                                                    infrastructure_interest_rate,
                                                    conversion_params_filepath,
                                                    pipeline_params_filepath, 
                                                    currency)[0]
                    elif plant_type == "ammonia":
                        pipeline_costs[i] =\
                            calculate_nh3_pipeline_costs(dist_to_demand,
                                                        annual_demand_quantity,
                                                        elec_price,
                                                        pipeline_params_filepath,
                                                        infrastructure_interest_rate,
                                                        currency)[0]
                else:
                    pipeline_costs[i] = np.nan
                
                # Calculate the cost of constructing a road to the hexagon if needed
                # and leave cost as 0 if road construction is not allowed or if distance
                # to road is 0
                if needs_road_construction == True:
                    # If the distance is more than 0 and less than 10, use short road costs
                    if dist_to_road!=0 and dist_to_road<10:
                        road_construction_costs[i] = (dist_to_road * 
                                                      short_road_capex * 
                                                      CRF(infrastructure_interest_rate, 
                                                          infrastructure_lifetime) + 
                                                        dist_to_road *
                                                        road_opex)/annual_demand_quantity
                    # Else, when the distance is more than 10, use long road costs
                    else:
                        road_construction_costs[i] = (dist_to_road * 
                                                      long_road_capex * 
                                                      CRF(infrastructure_interest_rate, 
                                                          infrastructure_lifetime) + 
                                                        dist_to_road *
                                                        road_opex)/annual_demand_quantity
                # Else if hexagon cannot reach a road (road construction is not allowed 
                # and distance to road is > 0), trucking is not viable:
                elif needs_road_construction == False and dist_to_road>0:
                    trucking_costs[i]=trucking_states[i] = np.nan
                    continue
                
                # Store cheapest trucking strategy values
                if plant_type == "hydrogen":
                    trucking_costs[i], trucking_states[i] =\
                        cheapest_trucking_strategy(demand_state,
                                                    annual_demand_quantity,
                                                    dist_to_demand,
                                                    elec_price,
                                                    heat_price,
                                                    infrastructure_interest_rate,
                                                    conversion_params_filepath,
                                                    transport_params_filepath,
                                                    currency)
                elif plant_type == "ammonia":
                    trucking_costs[i] = calculate_trucking_costs(demand_state,
                                                        dist_to_demand, 
                                                        annual_demand_quantity,
                                                        infrastructure_interest_rate, 
                                                        transport_params_filepath,
                                                        currency)/annual_demand_quantity
                    trucking_states[i] = "NH3"

        print("\nOptimisation complete.\n")
        # Hexagon file updated with each demand center's costs and states
        hexagons[f'{demand_center} road construction costs'] = road_construction_costs
        if plant_type == "hydrogen":
            hexagons[f'{demand_center} trucking transport and conversion costs'] = trucking_costs
            hexagons[f'{demand_center} pipeline transport and conversion costs'] = pipeline_costs
        elif plant_type == "ammonia":
            hexagons[f'{demand_center} trucking transport costs'] = trucking_costs
            hexagons[f'{demand_center} pipeline transport costs'] = pipeline_costs
        hexagons[f'{demand_center} trucking state'] = trucking_states

    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')

if __name__ == "__main__":
    main()
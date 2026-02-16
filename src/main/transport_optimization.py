"""
@authors:
 - Claire Halloran
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk
Contains code originally written by Leander MÃ¼ller, RWTH Aachen University

Calculates the cost-optimal commodity transportation strategy to the nearest 
demand center.
Calculates the cost-optimal transportation strategy of feedstock to the 
hexagon being optimised (Copper only).
Calculates cost of pipeline transport and demand profile based on optimal size
(Hydrogen and Ammonia only).
"""
import geopandas as gpd
import geopy.distance
import numpy as np
import pandas as pd
from shapely.geometry import Point
import shapely.wkt

from functions import cheapest_trucking_strategy, h2_conversion_stand, \
                    cheapest_pipeline_strategy, calculate_trucking_costs, \
                    calculate_nh3_pipeline_costs, \
                    calculate_conversion_facility_costs, \
                    calculate_construction_costs, \
                    compute_geodesic_distance_matrix, \
                    nearest_hexagon_to_location, allocate_feedstock_sources
from utils import check_folder_exists

def main():
    plant_type = snakemake.wildcards.plant_type
    tech_params_filepath = snakemake.input.technology_parameters
    demand_params_filepath = snakemake.input.demand_parameters
    country_params_filepath = snakemake.input.country_parameters
    transport_params_filepath = snakemake.input.transport_parameters

    if plant_type == 'hydrogen' or plant_type == 'ammonia':
        pipeline_params_filepath = snakemake.input.pipeline_parameters
        if plant_type == "hydrogen":
            conversion_params_filepath = snakemake.input.conversion_parameters
        needs_pipeline_construction = snakemake.config["transport"]["pipeline_construction"]
        
    infra_data = pd.read_excel(tech_params_filepath,
                               sheet_name='Infra',
                               index_col='Infrastructure')
    
    demand_center_list = pd.read_excel(demand_params_filepath,
                                       sheet_name='Demand centers',
                                       index_col='Demand center',)
    demand_centers = demand_center_list.index

    country_params = pd.read_excel(country_params_filepath, index_col='Country')
    country_series = country_params.iloc[0]

    hexagons = gpd.read_file(snakemake.input.hexagons)
    if hexagons.crs is None:
        hexagons = hexagons.set_crs(epsg=4326)
    elif hexagons.crs.to_epsg() != 4326:
        hexagons = hexagons.to_crs(epsg=4326)
    len_hexagons = len(hexagons)

    needs_road_construction = snakemake.config["transport"]["road_construction"]

    long_road_capex = infra_data.at['Long road','CAPEX']
    short_road_capex = infra_data.at['Short road','CAPEX']
    long_road_opex = infra_data.at['Long road','OPEX']
    short_road_opex = infra_data.at['Short road','OPEX']

    # Prices from the country excel file
    currency = snakemake.config["currency"]
    elec_price = country_series[f'Electricity price ({currency}/kWh)']
    heat_price = country_series[f'Heat price ({currency}/kWh)']
    plant_interest_rate = float(country_series['Plant interest rate'])
    infrastructure_interest_rate = float(country_series['Infrastructure interest rate'])
    infrastructure_lifetime = float(country_series['Infrastructure lifetime (years)'])

    check_folder_exists("resources")

    # Convert demand locations into a geodataframe
    demand_points_gdf = gpd.GeoDataFrame(demand_center_list, geometry=[Point(xy) for xy in zip(demand_center_list['Lon [deg]'], demand_center_list['Lat [deg]'])]).set_crs(epsg=4326)
    # Calculate distances to demand
    hexagon_to_demand_distance_matrix = compute_geodesic_distance_matrix(hexagons, demand_points_gdf)
    # Find closest hexagon to demand
    demand_points_gdf["nearest hexidx"] = [nearest_hexagon_to_location(idx, hexagon_to_demand_distance_matrix) for idx in demand_points_gdf.index]

    if plant_type == 'copper':
        size_offgrid_feedstocks = snakemake.config["size_offgrid_feedstocks"]
        conversion_params_filepath = snakemake.input.conversion_parameters
        feedstock_center_list = pd.read_excel(demand_params_filepath,
                                              sheet_name='Feedstock (kg.year-1)',
                                              index_col='Feedstock center')
        
        # Convert feedstock locations into a geodataframe
        feedstock_points_gdf = gpd.GeoDataFrame(feedstock_center_list, geometry=[Point(xy) for xy in zip(feedstock_center_list['Lon [deg]'], feedstock_center_list['Lat [deg]'])]).set_crs(epsg=4326)
        # Calculate distances to feedstock
        hexagon_to_feedstock_distance_matrix = compute_geodesic_distance_matrix(hexagons, feedstock_points_gdf)
        # Find closest hexagon to feedstock
        feedstock_points_gdf["nearest hexidx"] = [nearest_hexagon_to_location(idx, hexagon_to_feedstock_distance_matrix) for idx in feedstock_points_gdf.index]

        if size_offgrid_feedstocks:
            fs_output_path = f"resources/feedstocks_transport_{snakemake.wildcards.country}_{snakemake.wildcards.plant_type}.geojson"
            feedstock_points_gdf.to_file(fs_output_path, driver='GeoJSON', encoding='utf-8')
            
    # Calculate road construction costs
    if needs_road_construction:
        hexagons_road_construction = hexagons.apply(calculate_construction_costs,
                                                    args=(infrastructure_interest_rate,
                                                            infrastructure_lifetime,
                                                            long_road_capex,
                                                            short_road_capex,
                                                            long_road_opex,
                                                            short_road_opex,
                                                            'road_dist'),
                                                    axis=1)
        
    # Calculate cost of hydrogen state conversion and transportation for demand
    # Loop through all demand centers
    for demand_center in demand_centers:
        print(f"\nOptimisation for {demand_center} begins\n")
        # Demand location based variables
        demand_center_lat = demand_center_list.loc[demand_center,'Lat [deg]']
        demand_center_lon = demand_center_list.loc[demand_center,'Lon [deg]']
        demand_location = Point(demand_center_lon, demand_center_lat)
        annual_demand_quantity = demand_center_list.loc[demand_center,'Annual demand [kg/a]']
        demand_state = demand_center_list.loc[demand_center,'Demand state']
        demand_hix = demand_points_gdf["nearest hexidx"][demand_center]

        transport_params = pd.read_excel(transport_params_filepath,
                                            sheet_name = demand_state,
                                            index_col = 'Parameter'
                                            ).squeeze('columns')
        
        # Storage arrays containing zeros to be filled with costs and states or
        # left as zeros, as necessary
        road_construction_costs = np.zeros(len_hexagons)
        trucking_states = np.zeros(len_hexagons,dtype='<U100')
        trucking_costs = np.zeros(len_hexagons)
        if plant_type == 'copper':
            facility_annual_cost_per_kg = np.zeros(len_hexagons)
            facility_annual_capex = np.zeros(len_hexagons)
            facility_annual_opex = np.zeros(len_hexagons)
            demand_trucking_costs = np.zeros(len_hexagons)
            feedstocks_trucking_costs = np.zeros(len_hexagons)
            feedstock_costs = np.zeros(len_hexagons)
            feedstock_quantities = np.zeros(len_hexagons)
            diesel_costs = np.zeros(len_hexagons)
            diesel_demands = np.zeros(len_hexagons)
        elif plant_type == 'hydrogen' or plant_type == 'ammonia':
            pipeline_costs = np.zeros(len_hexagons)

        if plant_type == 'hydrogen':
            if demand_state not in ['500 bar', 'LH2', 'NH3', 'LOHC', "CuAnode"]:
                raise NotImplementedError(f'{demand_state} demand not supported.')
        elif plant_type == 'copper':
            if demand_state not in ["CuAnode", "CuCathode"]:
                raise NotImplementedError(f'{demand_state} demand not supported.')
            
            conversion_params = pd.read_excel(conversion_params_filepath,
                                    sheet_name = demand_state,
                                    index_col = 'Parameter').squeeze('columns')
            
            cucon_transport_params = pd.read_excel(transport_params_filepath,
                                            sheet_name = "CuConcentrate",
                                            index_col = 'Parameter'
                                            ).squeeze('columns')

            feedstock_quantity = annual_demand_quantity / float(conversion_params.loc["Efficiency (kg product / kg feedstock)"])
            cost_of_feedstock = (feedstock_quantity * float(conversion_params.loc[f"Feedstock price ({currency} per kg feedstock)"]))

            if feedstock_quantity > feedstock_points_gdf["Annual capacity [kg/a]"].sum():
                raise NotImplementedError(f'Not enough feedstock demand to meet {demand_center} {demand_state}')

            # Calculating diesel costs and demand
            diesel_unit_demand = conversion_params.loc["Diesel demand (kWh per kg product)"]
            diesel_demand = diesel_unit_demand * annual_demand_quantity
            diesel_cost = diesel_demand * country_series[f'Diesel price ({currency}/kWh)']
            
        # Loop through all hexagons
        for i in range(len_hexagons):
            print(f"Currently optimising {i+1} of {len_hexagons} hexagons")
            hex_geometry = hexagons['geometry'][i]
            
            # Calculating distance to demand
            poly = shapely.wkt.loads(str(hex_geometry))
            center = poly.centroid
            demand_coords = (demand_center_lat, demand_center_lon)
            hexagon_coords = (center.y, center.x)
            dist_to_demand = geopy.distance.geodesic(demand_coords, hexagon_coords).km
            
            if plant_type == "copper":
                feedstock_costs[i] = cost_of_feedstock / annual_demand_quantity
                feedstock_quantities[i] = feedstock_quantity
                # Allocate feedstock sources
                feedstock_sources = allocate_feedstock_sources(feedstock_points_gdf,
                                                                hexagon_to_feedstock_distance_matrix,
                                                                i,
                                                                feedstock_quantity)

                # Facility costs need same calculation for all hexagons
                if hex_geometry.contains(demand_location) == True or snakemake.config['restrict_hexagons'] == False:
                    facility_annual_capex[i], facility_annual_opex[i] = calculate_conversion_facility_costs(annual_demand_quantity,
                                                                                                    plant_interest_rate,
                                                                                                    conversion_params,
                                                                                                    currency)
                    facility_annual_cost_per_kg[i] = facility_annual_capex[i] + facility_annual_opex[i]

                # Store diesel information
                diesel_costs[i] = diesel_cost / annual_demand_quantity
                diesel_demands[i] = diesel_demand / annual_demand_quantity
            
            # Different calculations dependent if in demand location or not
            if hex_geometry.contains(demand_location) == True:
                if plant_type == 'copper':
                    feedstock_road_construction_costs = 0
                    feedstocks_trucking_cost = 0
                    # Road construction needed if feedstock is not in the hexagon
                    for f in feedstock_sources.index:
                        feedstock_hix = feedstock_points_gdf["nearest hexidx"][f]
                        if i != feedstock_hix:
                            # Build road for feedstock hexagon
                            feedstock_road_construction_costs += hexagons_road_construction[feedstock_hix]
                            
                            # Calculate road transport from feedstock to hexagon
                            source_quantity = feedstock_sources.loc[f, "Feedstock used [kg/year]"]
                            feedstocks_trucking_cost += calculate_trucking_costs("CuConcentrate",
                                                                                    hexagon_to_feedstock_distance_matrix.loc[i, f],
                                                                                    source_quantity,
                                                                                    infrastructure_interest_rate,
                                                                                    cucon_transport_params,
                                                                                    currency)
                            
                    feedstocks_trucking_costs[i] = feedstocks_trucking_cost/annual_demand_quantity
                    road_construction_costs[i] = feedstock_road_construction_costs/ annual_demand_quantity
                    demand_trucking_costs[i] = np.nan
                    trucking_states[i] = f"{demand_center} hexagon"
                elif plant_type == "hydrogen":
                    # Calculate cost of converting hydrogen to a demand state for local demand (i.e. no transport)
                    # State changed to label the hexagon with which demand center it contains
                    # Leave road construction cost as 0, as no road is needed
                    if demand_state == 'NH3':
                        trucking_costs[i]=pipeline_costs[i]=h2_conversion_stand(demand_state+'_load',
                                                                                annual_demand_quantity,
                                                                                elec_price,
                                                                                heat_price,
                                                                                plant_interest_rate,
                                                                                conversion_params_filepath,
                                                                                currency
                                                                                )[2]/annual_demand_quantity
                        trucking_states[i] = f"{demand_center} hexagon"
                    else:
                        trucking_costs[i]=pipeline_costs[i]=h2_conversion_stand(demand_state,
                                                                                annual_demand_quantity,
                                                                                elec_price,
                                                                                heat_price,
                                                                                plant_interest_rate,
                                                                                conversion_params_filepath,
                                                                                currency
                                                                                )[2]/annual_demand_quantity
                        trucking_states[i] = f"{demand_center} hexagon"
                elif plant_type == "ammonia":
                    # Leave the costs as 0 for ammonia, just change state
                    trucking_states[i] = f"{demand_center} hexagon"
            # Else, if the hexagon does not contain the demand center:
            else:
                if plant_type == 'copper':
                    # Only optimises the demand center hexagon
                    if snakemake.config['restrict_hexagons']:
                        facility_annual_cost_per_kg[i] = np.nan
                        facility_annual_capex[i] = np.nan
                        facility_annual_opex[i] = np.nan

                        road_construction_costs[i] = np.nan
                        trucking_costs[i] = np.nan
                        trucking_states[i] = np.nan
                        feedstocks_trucking_costs[i] = np.nan
                        demand_trucking_costs[i] = np.nan

                        feedstock_quantities[i] = np.nan
                        feedstock_costs[i] = np.nan
                        diesel_costs[i] = np.nan
                        diesel_demands[i] = np.nan
                        continue

                    # Calculate the cost of constructing a road to the hexagon if needed
                    if needs_road_construction:
                        # Build road for hexagon and demand hexagon
                        demand_road_construction = (hexagons_road_construction[i]
                                                    + hexagons_road_construction[demand_hix])
                        
                        # Calculating feedstock road construction costs
                        feedstock_road_construction = 0
                        for f in feedstock_sources.index:
                            feedstock_hix = feedstock_points_gdf["nearest hexidx"][f]
                            
                            if feedstock_hix != i:
                                # Build road for feedstock hexagon
                                feedstock_road_construction += hexagons_road_construction[feedstock_hix]
                        
                        # Total road construction costs
                        road_construction_costs[i] = (demand_road_construction + 
                                                    feedstock_road_construction)/annual_demand_quantity
                    elif hexagon_to_demand_distance_matrix.loc[i, demand_center]>0:
                        # Else if hexagon cannot reach a road (road construction is not allowed 
                        # and distance to road is > 0), trucking is not viable:
                        trucking_costs[i]=np.nan
                        continue
                    
                    # Calculate trucking costs to demand
                    demand_trucking_costs[i] = calculate_trucking_costs(demand_state,
                                                                        hexagon_to_demand_distance_matrix.loc[i, demand_center],
                                                                        annual_demand_quantity,
                                                                        infrastructure_interest_rate,
                                                                        transport_params,
                                                                        currency)/annual_demand_quantity

                    # Calculate trucking costs from feedstock
                    feedstocks_trucking_cost = 0
                    for f in feedstock_sources.index:
                        source_quantity = feedstock_sources.loc[f, "Feedstock used [kg/year]"]
                        feedstock_hix = feedstock_points_gdf["nearest hexidx"][f]

                        if feedstock_hix != i:
                            feedstocks_trucking_cost += calculate_trucking_costs("CuConcentrate",
                                                                                hexagon_to_feedstock_distance_matrix.loc[i, f],
                                                                                source_quantity,
                                                                                infrastructure_interest_rate,
                                                                                transport_params,
                                                                                currency)
                    feedstocks_trucking_costs[i] = feedstocks_trucking_cost/annual_demand_quantity
                    
                    # Storing trucking totals
                    trucking_costs[i] = road_construction_costs[i] + feedstocks_trucking_costs[i] + demand_trucking_costs[i]
                    trucking_states[i] = demand_state
                elif plant_type == 'hydrogen':
                    # Only optimises the demand center hexagon
                    if snakemake.config['restrict_hexagons']:
                        road_construction_costs[i] = np.nan
                        trucking_costs[i] = np.nan
                        pipeline_costs[i] = np.nan
                        trucking_states[i] = np.nan
                        continue

                    # Calculate costs of constructing a pipeline to the hexagon if allowed
                    if needs_pipeline_construction== True:
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
                    else:
                        pipeline_costs[i] = np.nan

                    # Calculate the cost of constructing a road to the hexagon if needed
                    if needs_road_construction:
                        # Road costs for hexagon and demand hexagon
                        road_construction_costs[i] = (hexagons_road_construction[i]
                                                    + hexagons_road_construction[demand_hix])/annual_demand_quantity
                    elif hexagon_to_demand_distance_matrix.loc[i, demand_center]>0:
                        # Else if hexagon cannot reach a road (road construction is not allowed 
                        # and distance to road is > 0), trucking is not viable:
                        trucking_costs[i]=trucking_states[i] = np.nan

                    # Store cheapest trucking strategy values
                    trucking_costs[i], trucking_states[i] =\
                        cheapest_trucking_strategy(demand_state,
                                                    annual_demand_quantity,
                                                    dist_to_demand,
                                                    elec_price,
                                                    heat_price,
                                                    infrastructure_interest_rate,
                                                    conversion_params_filepath,
                                                    transport_params,
                                                    currency)
                elif plant_type == 'ammonia':
                    # Only optimises the demand center hexagon
                    if snakemake.config['restrict_hexagons']:
                        road_construction_costs[i] = np.nan
                        trucking_costs[i] = np.nan
                        pipeline_costs[i] = np.nan
                        trucking_states[i] = np.nan
                        continue

                    # Calculate costs of constructing a pipeline to the hexagon if allowed
                    if needs_pipeline_construction== True:
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
                    if needs_road_construction:
                        # Road costs for hexagon and demand hexagon
                        road_construction_costs[i] = (hexagons_road_construction[i]
                                                    + hexagons_road_construction[demand_hix])/annual_demand_quantity
                    elif hexagon_to_demand_distance_matrix.loc[i, demand_center]>0:
                        # Else if hexagon cannot reach a road (road construction is not allowed 
                        # and distance to road is > 0), trucking is not viable:
                        trucking_costs[i]=trucking_states[i] = np.nan

                    # Store cheapest trucking strategy values
                    trucking_costs[i] = calculate_trucking_costs(demand_state,
                                                                dist_to_demand, 
                                                                annual_demand_quantity,
                                                                infrastructure_interest_rate, 
                                                                transport_params,
                                                                currency)/annual_demand_quantity
                    trucking_states[i] = "NH3"

        print("\nOptimisation complete.\n")
        # Hexagon file updated with each demand center's costs and states
        if plant_type == 'copper':
            # Facility costs
            hexagons[f'{demand_center} annual facility costs ({currency}/kg/year)'] = facility_annual_cost_per_kg
            hexagons[f'{demand_center} annual facility capex ({currency}/kg/year)'] = facility_annual_capex
            hexagons[f'{demand_center} annual facility opex ({currency}/kg/year)'] = facility_annual_opex

            # Transport related
            hexagons[f'{demand_center} road construction costs ({currency}/kg/year)'] = road_construction_costs
            hexagons[f'{demand_center} trucking transport costs ({currency}/kg/year)'] = trucking_costs
            hexagons[f'{demand_center} feedstock trucking transport costs ({currency}/kg/year)'] = feedstocks_trucking_costs
            hexagons[f'{demand_center} product trucking transport costs ({currency}/kg/year)'] = demand_trucking_costs
            hexagons[f'{demand_center} trucking state'] = trucking_states

            # Feedstock related
            hexagons[f'{demand_center} feedstock quantity (kg/year)'] = feedstock_quantities 
            hexagons[f'{demand_center} feedstock cost ({currency}/kg/year)'] = feedstock_costs
            
            # Diesel related
            hexagons[f'{demand_center} diesel cost ({currency}/kg/year)'] = diesel_costs
            hexagons[f'{demand_center} diesel usage (kWh/kg/year)'] = diesel_demands
        else:
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
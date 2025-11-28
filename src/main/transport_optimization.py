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

from functions import cheapest_trucking_strategy, h2_conversion_stand, \
                            cheapest_pipeline_strategy, calculate_trucking_costs, \
                            calculate_nh3_pipeline_costs, mineral_conversion_facility, \
                            calculate_construction_costs, geodesic_matrix, find_nearest_hex, \
                            determine_feedstock_sources, calculate_train_costs
from utils import check_folder_exists

def main():
    file = open("file.txt", "a")
    plant_type = snakemake.wildcards.plant_type
    tech_params_filepath = snakemake.input.technology_parameters
    demand_params_filepath = snakemake.input.demand_parameters
    country_params_filepath = snakemake.input.country_parameters
    transport_params_filepath = snakemake.input.transport_parameters

    if plant_type == 'minx':
        rail_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/rail_parameters.xlsx'
        conversion_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/conversion_parameters.xlsx'
        feedstock_center_list = pd.read_excel(demand_params_filepath,
                                    sheet_name='Feedstock (kg.year-1)',
                                      index_col='Feedstock center')
        needs_rail_construction = snakemake.config["transport"]["rail_construction"]
    else:
        pipeline_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/pipeline_parameters.xlsx'
        if plant_type == "hydrogen":
            conversion_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/conversion_parameters.xlsx'
        needs_pipeline_construction = snakemake.config["transport"]["pipeline_construction"]
        
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

    needs_road_construction = snakemake.config["transport"]["road_construction"]

    long_road_capex = infra_data.at['Long road','CAPEX']
    short_road_capex = infra_data.at['Short road','CAPEX']
    long_road_opex = infra_data.at['Long road','OPEX']
    short_road_opex = infra_data.at['Short road','OPEX']
    if plant_type == 'minx':
        long_rail_capex = infra_data.at['Long rail','CAPEX']
        short_rail_capex = infra_data.at['Short rail','CAPEX']
        long_rail_opex = infra_data.at['Long rail','OPEX']
        short_rail_opex = infra_data.at['Short rail','OPEX']

    # Prices from the country excel file
    currency = snakemake.config["currency"]
    elec_price = country_params[f'Electricity price ({currency}/kWh)'].iloc[0]
    heat_price = country_params[f'Heat price ({currency}/kWh)'].iloc[0]
    plant_interest_rate = float(country_params['Plant interest rate'].iloc[0])
    infrastructure_interest_rate = float(country_params['Infrastructure interest rate'].iloc[0])
    infrastructure_lifetime = float(country_params['Infrastructure lifetime (years)'].iloc[0])

    check_folder_exists("resources")

    # Convert demand locations into a geodataframe
    demand_points_gdf = gpd.GeoDataFrame(demand_center_list, geometry=[Point(xy) for xy in zip(demand_center_list['Lon [deg]'], demand_center_list['Lat [deg]'])]).set_crs(epsg=4326)
    # Calculate distances to demand
    hexagon_to_demand_distance_matrix = geodesic_matrix(hexagons, demand_points_gdf)
    # Find closest hexagon to demand
    demand_points_gdf["nearest hexidx"] = [find_nearest_hex(idx, hexagon_to_demand_distance_matrix) for idx in demand_points_gdf.index]

    if plant_type == 'minx':
        # Convert feedstock locations into a geodataframe
        feedstock_points_gdf = gpd.GeoDataFrame(feedstock_center_list, geometry=[Point(xy) for xy in zip(feedstock_center_list['Lon [deg]'], feedstock_center_list['Lat [deg]'])]).set_crs(epsg=4326)
        # Calculate distances to feedstock
        hexagon_to_feedstock_distance_matrix = geodesic_matrix(hexagons, feedstock_points_gdf)
        # Find closest hexagon to feedstock
        feedstock_points_gdf["nearest hexidx"] = [find_nearest_hex(idx, hexagon_to_feedstock_distance_matrix) for idx in feedstock_points_gdf.index]

        # Calculate rail construction costs
        if needs_rail_construction:
            hexagons_rail_construction = hexagons.apply(calculate_construction_costs,
                                                                args=(infrastructure_interest_rate,
                                                                    infrastructure_lifetime,
                                                                    long_rail_capex,
                                                                    short_rail_capex,
                                                                    long_rail_opex,
                                                                    short_rail_opex,
                                                                    'rail'),
                                                                axis=1)
            
    # Calculate road construction costs
    if needs_road_construction:
        hexagons_road_construction = hexagons.apply(calculate_construction_costs,
                                                    args=(infrastructure_interest_rate,
                                                            infrastructure_lifetime,
                                                            long_road_capex,
                                                            short_road_capex,
                                                            long_road_opex,
                                                            short_road_opex,
                                                            'road'),
                                                    axis=1)
        
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
        demand_hix = demand_points_gdf["nearest hexidx"][demand_center]

        # Storage arrays containing zeros to be filled with costs and states or
        # left as zeros, as necessary
        road_construction_costs = np.zeros(len_hexagons)
        trucking_states = np.zeros(len_hexagons,dtype='<U10')
        trucking_costs = np.zeros(len_hexagons)
        if plant_type == 'minx':
            facility_annual_cost_per_kg = np.zeros(len_hexagons)
            facility_annual_capex = np.zeros(len_hexagons)
            facility_annual_opex = np.zeros(len_hexagons)
            demand_trucking_costs = np.zeros(len_hexagons)
            feedstocks_trucking_costs = np.zeros(len_hexagons)
            rail_construction_costs = np.zeros(len_hexagons)
            demand_train_costs = np.zeros(len_hexagons)
            feedstocks_train_costs = np.zeros(len_hexagons)
            train_costs = np.zeros(len_hexagons)
        else:
            pipeline_costs = np.zeros(len_hexagons)

        if plant_type == 'hydrogen':
            if demand_state not in ['500 bar', 'LH2', 'NH3', 'LOHC', "CuAnode"]:
                raise NotImplementedError(f'{demand_state} demand not supported.')
        elif plant_type == 'minx':
            if demand_state not in ["CuAnode", "CuCathode"]:
                raise NotImplementedError(f'{demand_state} demand not supported.')
            else:
                conversion_params = pd.read_excel(conversion_params_filepath,
                                    sheet_name=demand_state,
                                    index_col='Parameter',)

            # -- BELOW WOULD CAUSE ISSUES FOR FUTURE FILES, DOES THIS NEED TO BE DONE?
            # skip if equal to zero or NA
            if ((annual_demand_quantity==0.0) or (np.isnan(annual_demand_quantity))):
                continue
            feedstock_quantity = annual_demand_quantity / float(conversion_params.loc["Efficiency (kg product / kg feedstock)"])
            cost_of_feedstock = (feedstock_quantity * float(conversion_params.loc["Feedstock price (euros per kg feedstock)"]))

            if feedstock_quantity > feedstock_points_gdf["Annual capacity [kg/a]"].sum():
                raise NotImplementedError(f'Not enough feedstock demand to meet {demand_center} {demand_state}')

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
            
            if plant_type == "minx":
                # determine feedstock sources
                feedstock_sources, feedstock_ranked_idxs = determine_feedstock_sources(feedstock_points_gdf,
                                                                                        hexagon_to_feedstock_distance_matrix,
                                                                                        i,
                                                                                        feedstock_quantity,file)
                # --Below is misc and is just added to the hexagon and not used anywhere else
                # feedstock_locs[i] = json.dumps(list(feedstock_ranked_idxs.values))

                # Facility costs need same calculation for all hexagons
                if hex_geometry.contains(demand_location) == True or snakemake.config['restrict_hexagons'] == False:
                    facility_annual_capex[i], facility_annual_opex[i] = mineral_conversion_facility(demand_state,
                                                                                                    annual_demand_quantity,
                                                                                                    plant_interest_rate,
                                                                                                    conversion_params_filepath)
                    facility_annual_cost_per_kg[i] = facility_annual_capex[i] + facility_annual_opex[i]

            # If the hexagon contains the demand location:
            if hex_geometry.contains(demand_location) == True:
                if plant_type == 'minx':
                    feedstock_road_construction_costs = 0
                    feedstocks_trucking_cost = 0
                    feedstocks_train_cost = 0
                    # Only road construction needed is if feedstock is not in the hexagon
                    for f in feedstock_sources.index:
                        feedstock_hix = feedstock_points_gdf["nearest hexidx"][f]
                        if i != feedstock_hix:
                            # build road for feedstock hexagon
                            feedstock_road_construction_costs += hexagons_road_construction[feedstock_hix]
                            
                            # calculate road transport from feedstock
                            source_quantity = feedstock_sources.loc[f, "Feedstock used [kg/year]"]
                            feedstocks_trucking_cost += calculate_trucking_costs("CuConcentrate",
                                                                                    hexagon_to_feedstock_distance_matrix.loc[i, f],
                                                                                    source_quantity,
                                                                                    infrastructure_interest_rate,
                                                                                    transport_params_filepath,
                                                                                    currency)
                            
                            for f in feedstock_sources.index:
                                source_quantity = feedstock_sources.loc[f, "Feedstock used [kg/year]"]
                                feedstocks_train_cost += calculate_train_costs("CuConcentrate",
                                                                                hexagon_to_feedstock_distance_matrix.loc[i, f],
                                                                                source_quantity,
                                                                                infrastructure_interest_rate,
                                                                                rail_params_filepath,
                                                                                currency)
                    
                    feedstocks_train_costs[i] = feedstocks_train_cost/annual_demand_quantity        
                    feedstocks_trucking_costs[i] = feedstocks_trucking_cost/annual_demand_quantity
                    road_construction_costs[i] = feedstock_road_construction_costs/ annual_demand_quantity
                    demand_trucking_costs[i] = np.nan

                    # Calculate rail transport from feedstock
                    feedstocks_train_cost = 0
                    for f in feedstock_sources.index:
                        source_quantity = feedstock_sources.loc[f, "Feedstock used [kg/year]"]
                        feedstocks_train_cost += calculate_train_costs("CuConcentrate",
                                                                        hexagon_to_feedstock_distance_matrix.loc[i, f],
                                                                        source_quantity,
                                                                        infrastructure_interest_rate,
                                                                        rail_params_filepath,
                                                                        currency)
                    feedstocks_train_costs[i] = feedstocks_train_cost/annual_demand_quantity
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
                if plant_type == 'minx':
                    # Forces only consideration of the demand center hexagon
                    if snakemake.config['restrict_hexagons']:
                        facility_annual_cost_per_kg[i] = np.nan
                        facility_annual_capex[i] = np.nan
                        facility_annual_opex[i] = np.nan

                        # Road costs
                        road_construction_costs[i] = np.nan
                        trucking_costs[i] = np.nan
                        feedstocks_trucking_costs[i] = np.nan
                        demand_trucking_costs[i] = np.nan
                        
                        # Rail costs
                        rail_construction_costs[i] = np.nan
                        demand_train_costs[i] = np.nan
                        feedstocks_train_costs[i] = np.nan
                        train_costs[i] = np.nan
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
                                # build road for feedstock hexagon
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
                                                                        transport_params_filepath,
                                                                        currency)/annual_demand_quantity

                    # Calculate trucking costs from feedstock
                    feedstocks_trucking_cost = 0
                    for f in feedstock_sources.index:
                        source_quantity = feedstock_sources.loc[f, "Feedstock used [kg/year]"]
                        # file.write(f"{source_quantity}\n\n")
                        feedstock_hix = feedstock_points_gdf["nearest hexidx"][f]

                        # if feedstock_hix != i:
                        feedstocks_trucking_cost += calculate_trucking_costs("CuConcentrate",
                                                                            hexagon_to_feedstock_distance_matrix.loc[i, f],
                                                                            source_quantity,
                                                                            infrastructure_interest_rate,
                                                                            transport_params_filepath,
                                                                            currency)
                    feedstocks_trucking_costs[i] = feedstocks_trucking_cost/annual_demand_quantity
                    
                    # Storing trucking totals
                    trucking_costs[i] = road_construction_costs[i] + feedstocks_trucking_costs[i] + demand_trucking_costs[i]

                    # Calculate the cost of constructing a rail to the hexagon if needed
                    if needs_rail_construction:
                        # Build rail for hexagon and demand hexagon
                        demand_rail_construction = (hexagons_rail_construction[i]
                                                    + hexagons_rail_construction[demand_hix])
                        
                        # Calculating feedstock rail construction costs
                        feedstock_rail_construction = 0
                        for f in feedstock_sources.index:
                            feedstock_hix = feedstock_points_gdf["nearest hexidx"][f]

                            if feedstock_hix != i:
                                # build rail for feedstock hexagon
                                feedstock_rail_construction += hexagons_rail_construction[feedstock_hix]
                        
                        # Total rail construction costs
                        rail_construction_costs[i] = (demand_rail_construction + 
                                                    feedstock_rail_construction)/annual_demand_quantity
                    elif hexagon_to_demand_distance_matrix.loc[i, demand_center]>0:
                        # Else if hexagon cannot reach a rail (rail construction is not allowed 
                        # and distance to rail is > 0), train is not viable:
                        train_costs[i] = np.nan
                        continue
                    
                    # Calculate rail transport to demand
                    demand_train_costs[i] = calculate_train_costs(demand_state,
                                                                hexagon_to_demand_distance_matrix.loc[i, demand_center],
                                                                annual_demand_quantity,
                                                                infrastructure_interest_rate,
                                                                rail_params_filepath,
                                                                currency)/annual_demand_quantity

                    # Calculate rail transport from feedstock
                    feedstocks_train_cost = 0
                    for f in feedstock_sources.index:
                        source_quantity = feedstock_sources.loc[f, "Feedstock used [kg/year]"]

                        feedstocks_train_cost += calculate_train_costs("CuConcentrate",
                                                                        hexagon_to_feedstock_distance_matrix.loc[i, f],
                                                                        source_quantity,
                                                                        infrastructure_interest_rate,
                                                                        rail_params_filepath,
                                                                        currency)
                    feedstocks_train_costs[i] = feedstocks_train_cost/annual_demand_quantity

                    # Storing train totals
                    train_costs[i] = rail_construction_costs[i] + feedstocks_train_costs[i] + demand_train_costs[i]
                elif plant_type == 'hydrogen':
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
                                                    transport_params_filepath,
                                                    currency)
                elif plant_type == 'ammonia':
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
                                                                transport_params_filepath,
                                                                currency)/annual_demand_quantity
                    trucking_states[i] = "NH3"

        print("\nOptimisation complete.\n")
        # Hexagon file updated with each demand center's costs and states
        if plant_type == 'minx':
            # Facility costs
            # Below this with up to the space have been calculated
            hexagons[f'{demand_center} annual facility costs ({currency}/kg/year)'] = facility_annual_cost_per_kg
            hexagons[f'{demand_center} annual facility capex ({currency}/kg/year)'] = facility_annual_capex
            hexagons[f'{demand_center} annual facility opex ({currency}/kg/year)'] = facility_annual_opex

            # Road costs
            hexagons[f'{demand_center} road construction costs ({currency}/kg/year)'] = road_construction_costs
            hexagons[f'{demand_center} trucking transport costs ({currency}/kg/year)'] = trucking_costs
            hexagons[f'{demand_center} feedstock trucking transport costs ({currency}/kg/year)'] = feedstocks_trucking_costs
            hexagons[f'{demand_center} product trucking transport costs ({currency}/kg/year)'] = demand_trucking_costs

            # Rail costs
            hexagons[f'{demand_center} rail construction costs ({currency}/kg/year)'] = rail_construction_costs
            hexagons[f'{demand_center} train transport costs ({currency}/kg/year)'] = train_costs  # cost of rail construction, supply conversion, train transport, and demand conversion
            hexagons[f'{demand_center} feedstock train transport costs ({currency}/kg/year)'] = feedstocks_train_costs  # cost of road construction, supply conversion, trucking transport, and demand conversion
            hexagons[f'{demand_center} product train transport costs ({currency}/kg/year)'] = demand_train_costs  # cost of road construction, supply conversion, trucking transport, and demand conversion

            # --- This below has already been calculated above, you just have to create the arrays to store in the info in.
            # --- Also, do we need this information?
            # hexagons[f'{demand_center} feedstock quantity (kg/year)'] = feedstock_quantities 
            # hexagons[f'{demand_center} feedstock cost ({currency}/kg/year)'] = feedstock_costs / product_quantity
            # hexagons[f'{demand_center} feedstock locs'] = feedstock_locs
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
    file.close()

if __name__ == "__main__":
    main()
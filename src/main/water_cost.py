"""
@authors: 
 - Claire Halloran
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk

Water costs for commodity production in each hexagon.
"""
import geopandas as gpd
import numpy as np
import pandas as pd

def main():
    print("Calculations begin...\n")
    hexagons = gpd.read_file(str(snakemake.input.hexagons))
    len_hexagons = len(hexagons)
    water_data = pd.read_excel(str(snakemake.input.technology_parameters), sheet_name='Water', index_col='Parameter').squeeze("columns")
    plant_type = str(snakemake.wildcards.plant_type)
    currency = snakemake.config["currency"]

    water_transport_cost = water_data[f'Water transport cost ({currency}/100 km/m3)']
    water_spec_cost = water_data[f'Water specific cost ({currency}/m3)']
    if plant_type == "minx":
        conversion_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/conversion_parameters.xlsx'
        demand_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/demand_parameters.xlsx'
        demand_center_list = pd.read_excel(demand_params_filepath,
                                    sheet_name='Demand centers',
                                    index_col='Demand center',)
        demand_centers = demand_center_list.index

        # Water cost related variables
        freshwater_costs = np.zeros(len_hexagons)
        freshwater_cost_per_kg_product = np.zeros(len_hexagons)
        ocean_water_costs = np.zeros(len_hexagons)
        ocean_water_cost_per_kg_product = np.zeros(len_hexagons)
        water_usages_L = np.zeros(len_hexagons)
        lowest_cost = np.zeros(len_hexagons)
    else:
        country_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/country_parameters.xlsx'
        country_params = pd.read_excel(country_params_filepath, index_col='Country')

        electricity_demand_h2o_treatment = water_data['Freshwater treatment electricity demand (kWh/m3)']
        electricity_demand_ocean_h2o_treatment = water_data['Ocean water treatment electricity demand (kWh/m3)']
        elec_price = country_params[f'Electricity price ({currency}/kWh)'].iloc[0]
        water_demand = water_data['Water demand  (L/kg of commodity)']

        # Water cost related variables
        h2o_costs_dom_water_bodies = np.zeros(len_hexagons)
        h2o_costs_ocean = np.zeros(len_hexagons)
        min_h2o_costs = np.zeros(len_hexagons)
    
    if plant_type == "minx":
        print(f"\Calculating water costs for {demand_center}...")
        for demand_center in demand_centers:
            annual_demand_quantity = demand_center_list.loc[demand_center,'Annual demand [kg/a]']
            demand_state = demand_center_list.loc[demand_center,'Demand state']
            conversion_params = pd.read_excel(conversion_params_filepath,
                                             sheet_name = demand_state,
                                             index_col = 'Parameter'
                                             ).squeeze('columns')
            water_intensity = conversion_params["Water intensity (L/kg product)"]

            # Calculating water costs for each hexagon
            for i in range(len_hexagons):
                print(f"Calculating water costs for {i+1} of {len_hexagons}...")
                waterbody_dist = hexagons['waterbody_dist'][i]
                waterway_dist = hexagons['waterway_dist'][i]
                ocean_dist = hexagons['ocean_dist'][i]

                water_usages_L[i] = annual_demand_quantity * water_intensity
                freshwater_costs[i] = (water_spec_cost * water_usages_L[i]
                                            + water_transport_cost * water_usages_L[i] * np.min([waterbody_dist, waterway_dist])/100)
                freshwater_cost_per_kg_product[i] = freshwater_costs[i] / annual_demand_quantity
                
                ocean_water_costs[i] = (water_spec_cost * water_usages_L[i]
                                        + water_transport_cost * water_usages_L[i] * ocean_dist/100)
                ocean_water_cost_per_kg_product[i] = ocean_water_costs[i] / annual_demand_quantity

                lowest_cost[i] = np.min(ocean_water_cost_per_kg_product[i], freshwater_cost_per_kg_product[i])

            hexagons[f'{demand_center} water usage (L)': water_usages_L]
            hexagons[f'{demand_center} freshwater cost ({currency}/kg/year)': freshwater_cost_per_kg_product]
            hexagons[f'{demand_center} ocean water cost ({currency}/kg/year)': ocean_water_cost_per_kg_product]
            hexagons[f'{demand_center} lowest water cost ({currency}/kg/year)': lowest_cost]
    else:
        # Calculating water costs for each hexagon
        for i in range(len_hexagons):
            print(f"Calculating water costs for {i+1} of {len_hexagons}...")
            waterbody_dist = hexagons['waterbody_dist'][i]
            waterway_dist = hexagons['waterway_dist'][i]
            ocean_dist = hexagons['ocean_dist'][i]
            
            h2o_costs_dom_water_bodies[i] =(water_spec_cost 
                                                + (water_transport_cost/100)
                                                * min(waterbody_dist, waterway_dist) 
                                                + electricity_demand_h2o_treatment
                                                * elec_price
                                            ) * water_demand/1000
            
            h2o_costs_ocean[i] =(water_spec_cost 
                                    + (water_transport_cost/100)
                                    * ocean_dist
                                    + electricity_demand_ocean_h2o_treatment
                                    * elec_price
                                ) * water_demand/1000
            
            min_h2o_costs[i] = min(h2o_costs_dom_water_bodies[i], h2o_costs_ocean[i])
        
        hexagons['Ocean water cost'] = h2o_costs_ocean
        hexagons['Freshwater cost'] = h2o_costs_dom_water_bodies
        hexagons['Lowest water cost'] = min_h2o_costs

    print("\nCalculations complete.\n")

    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')

if __name__ == "__main__":
    main()
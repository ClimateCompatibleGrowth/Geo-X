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
    tech_params_filepath = str(snakemake.input.technology_parameters)
    country_params_filepath = str(snakemake.input.country_parameters)

    water_data = pd.read_excel(tech_params_filepath, sheet_name='Water', index_col='Parameter').squeeze("columns")
    country_params = pd.read_excel(country_params_filepath, index_col='Country')

    # Water cost for each hexagon for each kg commodity produced

    # Water cost related variables
    h2o_costs_dom_water_bodies = np.empty(len(hexagons))
    h2o_costs_ocean = np.empty(len(hexagons))
    min_h2o_costs = np.empty(len(hexagons))
    
    currency = snakemake.config["currency"]
    electricity_demand_h2o_treatment = water_data['Freshwater treatment electricity demand (kWh/m3)']
    electricity_demand_ocean_h2o_treatment = water_data['Ocean water treatment electricity demand (kWh/m3)']
    water_transport_costs = water_data[f'Water transport cost ({currency}/100 km/m3)']
    water_spec_cost = water_data[f'Water specific cost ({currency}/m3)']
    water_demand = water_data['Water demand  (L/kg of commodity)']
    elec_price = country_params[f'Electricity price ({currency}/kWh)'].iloc[0]
    
    # Loop through all hexagons
    # Calculating water costs for each hexagon
    for i in range(len(hexagons)):
        print(f"Calculating water costs for {i+1} of {len(hexagons)}...")
        waterbody_dist = hexagons['waterbody_dist'][i]
        waterway_dist = hexagons['waterway_dist'][i]
        ocean_dist = hexagons['ocean_dist'][i]
        
        h2o_costs_dom_water_bodies[i] =(water_spec_cost 
                                            + (water_transport_costs/100)
                                            * min(waterbody_dist, waterway_dist) 
                                            + electricity_demand_h2o_treatment
                                            * elec_price
                                        ) * water_demand/1000
        
        h2o_costs_ocean[i] =(water_spec_cost 
                                + (water_transport_costs/100)
                                * ocean_dist
                                + electricity_demand_ocean_h2o_treatment
                                * elec_price
                            ) * water_demand/1000
        
        min_h2o_costs[i] = min(h2o_costs_dom_water_bodies[i], h2o_costs_ocean[i])

    print("\nCalculations complete.\n")
    hexagons['Ocean water costs'] = h2o_costs_ocean
    hexagons['Freshwater costs'] = h2o_costs_dom_water_bodies
    hexagons['Lowest water cost'] = min_h2o_costs

    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')

if __name__ == "__main__":
    main()
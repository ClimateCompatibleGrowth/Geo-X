"""
@authors: 
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk

Add attributes to hex file for cost of each component.
"""
import geopandas as gpd
import pandas as pd

from functions import CRF

if __name__ == "__main__":
    # Load hexagons
    hexagons = gpd.read_file(snakemake.input.hexagons)

    generators = snakemake.config['generators_dict']
    plant_type = snakemake.wildcards.plant_type

    # Load necessary parameters
    demand_excel_path = snakemake.input.demand_parameters
    demand_parameters = pd.read_excel(demand_excel_path, index_col='Demand center')
    country_excel_path = snakemake.input.country_parameters
    country_parameters = pd.read_excel(country_excel_path, index_col='Country')
    country_series = country_parameters.iloc[0]
    if plant_type == "hydrogen":
        # Battery
        storage_csv_path = 'parameters/basic_h2_plant/storage_units.csv'
        storage_parameters = pd.read_csv(storage_csv_path, index_col='name')
        # H2 storage
        stores_csv_path = 'parameters/basic_h2_plant/stores.csv'
        stores_parameters = pd.read_csv(stores_csv_path, index_col='name')
        # Electrolyzer
        links_csv_path = 'parameters/basic_h2_plant/links.csv'
        links_parameters = pd.read_csv(links_csv_path, index_col='name')
        # Generators
        generators_csv_path = 'parameters/basic_h2_plant/generators.csv'
        generators_parameters = pd.read_csv(generators_csv_path, index_col='name')
    elif plant_type == "ammonia":
        # H2 storage and battery
        stores_csv_path = 'parameters/basic_nh3_plant/stores.csv'
        stores_parameters = pd.read_csv(stores_csv_path, index_col='name')
        # Electrolyzer
        links_csv_path = 'parameters/basic_nh3_plant/links.csv' 
        links_parameters = pd.read_csv(links_csv_path, index_col='name')
        # Generators 
        generators_csv_path = 'parameters/basic_nh3_plant/generators.csv'
        generators_parameters = pd.read_csv(generators_csv_path, index_col='name')

    demand_centers = demand_parameters.index

    # For each demand center, get costs for each component
    for demand_center in demand_centers:
        print(f"\nCalculating for {demand_center} begins...")
  
        # Store CRF from plant data
        interest_plant = country_series['Plant interest rate']
        lifetime_plant = country_series['Plant lifetime (years)']
        crf_plant = CRF(interest_plant, lifetime_plant)
        
        # Work out the cost for each component using the data for the country 
        # you are looking at
        # Battery
        if plant_type == "hydrogen":
            capital_cost_battery = storage_parameters.loc['Battery', 'capital_cost']
        elif plant_type == "ammonia":
            capital_cost_battery = stores_parameters.loc['Battery', 'capital_cost']
        hexagons[f'{demand_center} battery costs'] = \
            hexagons[f'{demand_center} battery capacity'] *\
                capital_cost_battery * crf_plant
        hexagons[f'{demand_center} LC - battery costs portion'] = \
            hexagons[f'{demand_center} battery costs']/ \
                demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
        
        # Electrolyzer
        capital_cost_electrolyzer = links_parameters.loc['Electrolysis', 'capital_cost']
        hexagons[f'{demand_center} electrolyzer costs'] = \
            hexagons[f'{demand_center} electrolyzer capacity'] *\
                capital_cost_electrolyzer * crf_plant
        hexagons[f'{demand_center} LC - electrolyzer portion'] = \
            hexagons[f'{demand_center} electrolyzer costs']/ \
                demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
        
        # H2 Storage
        if plant_type == "hydrogen":
            capital_cost_h2_storage = stores_parameters.loc['Compressed H2 Store', 'capital_cost']
        else:
            capital_cost_h2_storage = stores_parameters.loc['CompressedH2Store', 'capital_cost']
        hexagons[f'{demand_center} H2 storage costs'] = \
            hexagons[f'{demand_center} H2 storage capacity'] *\
                capital_cost_h2_storage * crf_plant
        hexagons[f'{demand_center} LC - H2 storage portion'] = \
            hexagons[f'{demand_center} H2 storage costs']/ \
                demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
        
        # HB - ammonia only
        if plant_type == "ammonia":
            capital_cost_hb = links_parameters.loc['HB', 'capital_cost']
            hexagons[f'{demand_center} HB costs'] = \
                hexagons[f'{demand_center} HB capacity'] *\
                    capital_cost_hb * crf_plant
            hexagons[f'{demand_center} LC - HB portion'] = \
                hexagons[f'{demand_center} HB costs']/ \
                    demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

        # Work out CRF, then work out the cost for each generator using the data 
        # for the country you are looking at
        for generator in generators:
            generator_capitalized = generator.capitalize()
            interest_generator = country_series[f'{generator_capitalized} interest rate']
            lifetime_generator = country_series[f'{generator_capitalized} lifetime (years)']
            crf_generator = CRF(interest_generator, lifetime_generator)
            capital_cost_generator = generators_parameters.loc[f'{generator_capitalized}', 'capital_cost']
            hexagons[f'{demand_center} {generator} costs'] = \
                hexagons[f'{demand_center} {generator} capacity'] * capital_cost_generator * crf_generator
            hexagons[f'{demand_center} LC - {generator} portion'] = \
                hexagons[f'{demand_center} {generator} costs']/ \
                    demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
            
    print("\nCalculations complete.\n")
    hexagons.to_file(snakemake.output[0], driver='GeoJSON', encoding='utf-8')
    hexagons.to_csv(snakemake.output[1], encoding='latin-1')
"""
@authors: 
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Mulako Mukelabai, University of Oxford, mulako.mukelabai@eng.ox.ac.uk


Add attributes to hex file for cost of each component.
"""
import geopandas as gpd
import pandas as pd

from functions import CRF, annualise_capex_with_replacements

GRID_CONSTRUCTION_CLASSIFICATION = "INFRA"


def _canonicalize_component_columns(hexagons: pd.DataFrame, currency: str) -> pd.DataFrame:
    """
    Convert active component output columns to the canonical decision-CSV suffixes.

    Only `LC - ...` component columns are renamed. Aggregate LCOP, energy, water,
    and other decision metrics keep their established human-readable names.
    Legacy component suffixes are removed from the returned frame.
    """
    atomic_suffix = f" cost ({currency}/kg/year)"
    aggregate_suffix = f" total cost ({currency}/kg/year)"
    canonical_map = {}
    for col in hexagons.columns:
        if not (isinstance(col, str) and " LC - " in col):
            continue
        if col.endswith(aggregate_suffix):
            canonical_map[col] = col[: -len(aggregate_suffix)] + " total_cost_eur_per_kg_year"
        elif col.endswith(atomic_suffix):
            canonical_map[col] = col[: -len(atomic_suffix)] + " cost_eur_per_kg_year"

    if canonical_map:
        hexagons = hexagons.rename(columns=canonical_map)

    legacy_portion_suffix = " " + "portion"
    legacy_component_suffix = " component" + "_cost_eur_per_kg_year"
    drop_cols = [
        c for c in hexagons.columns
        if isinstance(c, str) and (c.endswith(legacy_portion_suffix) or c.endswith(legacy_component_suffix))
    ]
    if drop_cols:
        hexagons = hexagons.drop(columns=drop_cols)

    return hexagons

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
        generators_csv_path = 'parameters/basic_cu_plant/generators.csv'
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
    elif plant_type == "copper":
        # Rectifier and Inverter
        links_csv_path = 'parameters/basic_cu_plant/links.csv'
        links_parameters = pd.read_csv(links_csv_path, index_col='name')
        # Battery
        storage_csv_path = 'parameters/basic_cu_plant/storage_units.csv'
        storage_parameters = pd.read_csv(storage_csv_path, index_col='name')
        # Generators
        generators_csv_path = 'parameters/basic_h2_plant/generators.csv'
        generators_parameters = pd.read_csv(generators_csv_path, index_col='name')
        # Currency
        currency = snakemake.config["currency"]

    demand_centers = demand_parameters.index

    # For each demand center, get costs for each component
    for demand_center in demand_centers:
        print(f"\nCalculating for {demand_center} begins")
  
        # Store CRF from plant data
        interest_plant = country_series['Plant interest rate']
        lifetime_plant = country_series['Plant lifetime (years)']
        crf_plant = CRF(interest_plant, lifetime_plant)
        
        # Work out the cost for each component using the data for the country 
        # you are looking at
        if plant_type == 'copper':
            if GRID_CONSTRUCTION_CLASSIFICATION not in {"ENERGY", "INFRA"}:
                raise ValueError(
                    f"issue=6 file=src/main/costs_by_component.py context=invalid_grid_construction_classification missing=[{GRID_CONSTRUCTION_CLASSIFICATION}]"
                )
            annual_demand = demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
            battery_interest = country_series['Battery interest rate']
            battery_lifetime = country_series['Battery lifetime (years)']
            plant_lifetime = country_series['Plant lifetime (years)']
            system_types = ['offgrid','hybrid']
            for system_type in system_types:
                # Rectifier
                rectifier_capacity_col = f'{demand_center} {system_type} rectifier capacity (MW)'
                if rectifier_capacity_col not in hexagons.columns:
                    raise ValueError(f"issue=6 file=src/main/costs_by_component.py context=rectifier_capacity missing={rectifier_capacity_col}")
                capital_cost_rectifier = links_parameters.loc['Rectifier', 'capital_cost']
                hexagons[f'{demand_center} {system_type} rectifier costs'] = \
                    hexagons[rectifier_capacity_col] *\
                        capital_cost_rectifier * crf_plant
                hexagons[f'{demand_center} LC - {system_type} rectifier costs cost ({currency}/kg/year)'] = \
                    hexagons[f'{demand_center} {system_type} rectifier costs']/ \
                        annual_demand
                
                # Inverter
                inverter_capacity_col = f'{demand_center} {system_type} inverter capacity (MW)'
                if inverter_capacity_col not in hexagons.columns:
                    raise ValueError(f"issue=6 file=src/main/costs_by_component.py context=inverter_capacity missing={inverter_capacity_col}")
                capital_cost_inverter = links_parameters.loc['Inverter', 'capital_cost']
                hexagons[f'{demand_center} {system_type} inverter costs'] = \
                    hexagons[inverter_capacity_col] *\
                        capital_cost_inverter * crf_plant
                hexagons[f'{demand_center} LC - {system_type} inverter costs cost ({currency}/kg/year)'] = \
                    hexagons[f'{demand_center} {system_type} inverter costs']/ \
                        annual_demand
                
                # Battery
                battery_capacity_col = f'{demand_center} {system_type} battery capacity (MW)'
                if battery_capacity_col not in hexagons.columns:
                    raise ValueError(f"issue=9 file=src/main/costs_by_component.py context=battery_capacity missing={battery_capacity_col}")
                capital_cost_battery = storage_parameters.loc['Battery', 'capital_cost']
                battery_annualised_per_mw = annualise_capex_with_replacements(
                    capital_cost_battery,
                    battery_interest,
                    plant_lifetime,
                    battery_lifetime,
                )
                hexagons[f'{demand_center} {system_type} battery costs'] = \
                    hexagons[battery_capacity_col] * battery_annualised_per_mw
                hexagons[f'{demand_center} LC - {system_type} battery costs cost ({currency}/kg/year)'] = \
                    hexagons[f'{demand_center} {system_type} battery costs']/ \
                        annual_demand

                # Work out CRF, then work out the cost for each generator using the data 
                # for the country you are looking at
                for generator in generators:
                    objective_component_col = (
                        f'{demand_center} {system_type} {generator} objective component cost ({currency}/kg/year)'
                    )
                    if objective_component_col not in hexagons.columns:
                        raise ValueError(
                            f"issue=6 file=src/main/costs_by_component.py context=missing_objective_based_renewables missing={[objective_component_col]}"
                        )
                    hexagons[f'{demand_center} {system_type} {generator} costs'] = (
                        hexagons[objective_component_col] * annual_demand
                    )
                    hexagons[f'{demand_center} LC - {system_type} {generator} cost ({currency}/kg/year)'] = \
                        hexagons[objective_component_col]
                
                if system_type == 'hybrid':
                    required_hybrid_cols = [
                        f'{demand_center} hybrid grid purchase costs ({currency}/kg/year)',
                        f'{demand_center} hybrid grid capacity charge ({currency}/kg/year)',
                        f'{demand_center} hybrid grid fixed charge ({currency}/kg/year)',
                        f'{demand_center} hybrid grid construction charge ({currency}/kg/year)',
                    ]
                    missing_hybrid_cols = [c for c in required_hybrid_cols if c not in hexagons.columns]
                    if missing_hybrid_cols:
                        raise ValueError(f"issue=6 file=src/main/costs_by_component.py context=hybrid_grid_components missing={missing_hybrid_cols}")
                    hexagons[f'{demand_center} LC - hybrid grid purchase costs cost ({currency}/kg/year)'] = \
                        hexagons[f'{demand_center} hybrid grid purchase costs ({currency}/kg/year)']
                    hexagons[f'{demand_center} LC - hybrid grid capacity charge cost ({currency}/kg/year)'] = \
                        hexagons[f'{demand_center} hybrid grid capacity charge ({currency}/kg/year)']
                    hexagons[f'{demand_center} LC - hybrid grid fixed charge cost ({currency}/kg/year)'] = \
                        hexagons[f'{demand_center} hybrid grid fixed charge ({currency}/kg/year)']
                    hexagons[f'{demand_center} LC - hybrid grid construction charge cost ({currency}/kg/year)'] = \
                        hexagons[f'{demand_center} hybrid grid construction charge ({currency}/kg/year)']
                    # Legacy aggregate column retained for backwards compatibility.
                    hexagons[f'{demand_center} LC - {system_type} grid costs total cost ({currency}/kg/year)'] = (
                        hexagons[f'{demand_center} LC - hybrid grid purchase costs cost ({currency}/kg/year)']
                        + hexagons[f'{demand_center} LC - hybrid grid capacity charge cost ({currency}/kg/year)']
                        + hexagons[f'{demand_center} LC - hybrid grid fixed charge cost ({currency}/kg/year)']
                        + hexagons[f'{demand_center} LC - hybrid grid construction charge cost ({currency}/kg/year)']
                    )
                    hexagons[f'{demand_center} {system_type} grid costs'] = \
                        hexagons[f'{demand_center} LC - {system_type} grid costs total cost ({currency}/kg/year)'] * annual_demand

        else:
            # Battery
            if plant_type == "hydrogen":
                capital_cost_battery = storage_parameters.loc['Battery', 'capital_cost']
            elif plant_type == "ammonia":
                capital_cost_battery = stores_parameters.loc['Battery', 'capital_cost']
            hexagons[f'{demand_center} battery costs'] = \
                hexagons[f'{demand_center} battery capacity'] *\
                    capital_cost_battery * crf_plant
            hexagons[f'{demand_center} LC - battery costs cost ({currency}/kg/year)'] = \
                hexagons[f'{demand_center} battery costs']/ \
                    demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
            
            # Electrolyzer
            if plant_type == 'hydrogen' or plant_type == 'ammonia':
                electrolyzer_annual_cost_col = f'{demand_center} electrolyzer annual costs'
                if electrolyzer_annual_cost_col not in hexagons.columns:
                    raise ValueError(
                        f"issue=10 file=src/main/costs_by_component.py context=electrolyzer_annual_costs missing=['{electrolyzer_annual_cost_col}']"
                    )
                hexagons[f'{demand_center} electrolyzer costs'] = \
                    hexagons[electrolyzer_annual_cost_col]
                hexagons[f'{demand_center} LC - electrolyzer cost ({currency}/kg/year)'] = \
                    hexagons[f'{demand_center} electrolyzer costs']/ \
                        demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

            # H2 Storage
            if plant_type == 'hydrogen' or plant_type == 'ammonia':
                if plant_type == "hydrogen":
                    capital_cost_h2_storage = stores_parameters.loc['Compressed H2 Store', 'capital_cost']
                else:
                    capital_cost_h2_storage = stores_parameters.loc['CompressedH2Store', 'capital_cost']
                hexagons[f'{demand_center} H2 storage costs'] = \
                    hexagons[f'{demand_center} H2 storage capacity'] *\
                        capital_cost_h2_storage * crf_plant
                hexagons[f'{demand_center} LC - H2 storage cost ({currency}/kg/year)'] = \
                    hexagons[f'{demand_center} H2 storage costs']/ \
                        demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
            
            # HB and NH3 storage - ammonia only
            if plant_type == "ammonia":
                capital_cost_hb = links_parameters.loc['HB', 'capital_cost']
                hexagons[f'{demand_center} HB costs'] = \
                    hexagons[f'{demand_center} HB capacity'] *\
                        capital_cost_hb * crf_plant
                hexagons[f'{demand_center} LC - HB cost ({currency}/kg/year)'] = \
                    hexagons[f'{demand_center} HB costs']/ \
                        demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
                
                capital_cost_nh3_storage = links_parameters.loc['Ammonia', 'capital_cost']
                hexagons[f'{demand_center} NH3 storage costs'] = \
                    hexagons[f'{demand_center} NH3 storage capacity'] *\
                        capital_cost_nh3_storage * crf_plant
                hexagons[f'{demand_center} LC - NH3 storage cost ({currency}/kg/year)'] = \
                    hexagons[f'{demand_center} NH3 storage costs']/ \
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
                hexagons[f'{demand_center} LC - {generator} cost ({currency}/kg/year)'] = \
                    hexagons[f'{demand_center} {generator} costs']/ \
                        demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
            
    hexagons = _canonicalize_component_columns(hexagons, currency)

    print("\nCalculations complete.\n")
    hexagons.to_file(snakemake.output[0], driver='GeoJSON', encoding='utf-8')
    hexagons.to_csv(snakemake.output[1], encoding='latin-1')

"""
@authors:
 - Claire Halloran
 - Samiyha Naqvi
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk
 - Lukas Schirren, Imperial College London, lukas.schirren@imperial.ac.uk
 - Mulako Mukelabai, University of Oxford, mulako.mukelabai@eng.ox.ac.uk
Includes code from Nicholas Salmon, University of Oxford, for optimizing
the commodity's plant capacity.

Optimisation of the commodity production plant, accounting for e.g., demand 
schedule, water constraints, and renewable energy availability.
"""
import atlite
from copy import deepcopy
import geopandas as gpd
import logging
import numpy as np
import pandas as pd
import pyomo.environ as pm
import pypsa
from scipy.constants import physical_constants
import xarray as xr
import glob

from functions import (
    CRF,
    annualise_capex_with_hourly_replacements,
)
from network import Network

def hydropower_potential_with_capacity(flowrate, head, capacity, eta):
    '''
    Calculate the hydropower potential considering the capacity limit.

    ...
    Parameters
    ----------
    flowrate : float
        flowrate calculated with runoff multiplied by the hydro-basin area, 
        in cubic meters per hour.
    head : float
        height difference at the hydropower site, in meters.
    capacity : float
        maximum hydropower capacity in Megawatts (MW).
    eta : float
        efficiency of the hydropower plant.

    Returns
    -------
    capacity_factor : xarray DataArray
        limited potential divided by the capacity.
    '''
    # kg/m3; Density of water
    rho = 997
    
    # m/s2; Based on the CODATA constants 2018
    g = physical_constants['standard acceleration of gravity'][0]
    
    # Transform flowrate per hour into flowrate per second
    # (if Atlite is updated to 0.4.0 the term (1000/24) can be removed)
    Q = (flowrate/(1000/24)) / 3600

    potential = (eta * rho * g * Q * head) / (1000 * 1000) # MW
    limited_potential = xr.where(potential > capacity, capacity, potential)
    capacity_factor = limited_potential / capacity

    return capacity_factor

def calculate_grid_construction(hexagon, infra_data, country_series):
    '''
    Calculate the costs associated with constructing a grid in a hexagon.

    ...
    Parameters
    ----------
    hexagon : pandas Series
        contains geospatial information regarding a hexagonally shaped area.
    infra_data : pandas DataFrame
        contains data from the 'Infra' sheet in the technology parameters 
        excel document.
    country_series : pandas Series
        contains data from country parameters excel sheet.

    Returns
    -------
    grid_construction_cost : float
        cost of constructing a grid in a hexagon.
    '''
    grid_capex = infra_data.at['Grid','CAPEX']
    grid_opex = infra_data.at['Grid','OPEX']
    infrastructure_interest_rate = float(country_series['Infrastructure interest rate'])
    infrastructure_lifetime = float(country_series['Infrastructure lifetime (years)']) 

    if hexagon['grid_dist']==0:
        grid_construction_cost = 0.
    else:
        grid_construction_cost = (hexagon['grid_dist'] 
                                   * grid_capex 
                                   * CRF(infrastructure_interest_rate,
                                         infrastructure_lifetime)
                                   + hexagon['grid_dist'] * grid_opex)
    
    return grid_construction_cost

def calculate_grid_costs(demand, elec_kWh_per_kg, country_series, currency):
    '''
    Calculate the costs and capacities associated with using a grid to supply
    the energy required.

    ...
    Parameters
    ----------
    demand : pandas DataFrame
        dataframe of commodity demand in kg in frequency configured.
    elec_kWh_per_kg : float
        electricity demand (kWh per kg product) based on demand state.
    country_series : pandas Series
        contains data from country parameters excel sheet.
    currency : string
        unit of currency that is used in the parameter files.
    
    Returns
    -------
    grid_energy_cost : float
        cost for energy from grid to fulfil demand required.
    lcoe : float
        levelized cost of energy in currency unit/ kWh.
    grid_capacity : float
        optimal grid capacity in MW.
    '''
    grid_capacity = (elec_kWh_per_kg
                        * demand).max().iloc[0]/1000 # in MW if hourly data
    
    grid_energy_cost = ((elec_kWh_per_kg * demand * country_series[f"Electricity price ({currency}/kWh)"]).sum().iloc[0]
                        + (grid_capacity * country_series[f"Grid connection cost ({currency}/kW)"] * 1000)
                        + country_series[f"Electricity fixed charge ({currency}/year)"])
    
    lcoe = grid_energy_cost / (elec_kWh_per_kg * demand).sum().iloc[0]  # currency/kWh
    
    return grid_energy_cost, lcoe, grid_capacity

def _get_objective_snapshot_weights(network):
    """
    Return the snapshot weights that contribute to the PyPSA objective value.

    Parameters
    ----------
    network : pypsa.Network
        Solved network whose objective is being decomposed.

    Returns
    -------
    pandas.Series or pandas.Index
        Snapshot weights aligned with the dispatch time series used in the
        network objective.
    """
    snapshot_weightings = network.snapshot_weightings
    if isinstance(snapshot_weightings, pd.DataFrame):
        if 'objective' in snapshot_weightings.columns:
            return snapshot_weightings['objective']
        if 'generators' in snapshot_weightings.columns:
            return snapshot_weightings['generators']
        return snapshot_weightings.iloc[:, 0]
    return snapshot_weightings


def _calculate_renewable_objective_costs(network, generators, annual_demand_quantity):
    """
    Decompose renewable objective costs into per-kilogram component charges.

    Parameters
    ----------
    network : pypsa.Network
        Solved network whose objective is being decomposed.
    generators : list[str]
        Generator names from the scenario config, stored in lowercase.
    annual_demand_quantity : float
        Annual production quantity used to normalize annual costs.
    Returns
    -------
    dict[str, float]
        Per-kilogram annualized objective contribution for each configured
        renewable generator.
    """
    snapshot_weights = _get_objective_snapshot_weights(network)
    component_costs = {}
    for generator in generators:
        generator_name = generator.capitalize()
        annual_component_cost = (
            network.generators.at[generator_name, 'capital_cost'] * network.generators.p_nom_opt[generator_name]
            + network.generators.at[generator_name, 'marginal_cost']
            * (network.generators_t.p[generator_name] * snapshot_weights).sum()
        )
        component_costs[generator] = annual_component_cost / annual_demand_quantity

    return component_costs


def _calculate_link_operating_hours(network, link_name):
    """
    Calculate weighted annual operating hours for a solved link dispatch series.

    Parameters
    ----------
    network : pypsa.Network
        Solved network containing the link dispatch results.
    link_name : str
        Name of the link whose operating hours should be measured.

    Returns
    -------
    float
        Weighted annual operating hours, counting each active weighted snapshot
        once when the link dispatch is non-zero.
    """
    snapshot_weights = _get_objective_snapshot_weights(network)
    active = network.links_t.p0[link_name].abs() > 1e-9
    return float(snapshot_weights[active].sum())


def _calculate_product_output_kg(network_class):
    """
    Calculate solved annual product output in kilograms from the plant load.

    Parameters
    ----------
    network_class : Network
        Solved plant network wrapper.

    Returns
    -------
    float
        Weighted annual product output in kilograms.
    """
    if network_class.type == "hydrogen":
        load_name = "Hydrogen demand"
        specific_energy_mwh_per_tonne = 39.4
    elif network_class.type == "ammonia":
        load_name = "Ammonia demand"
        specific_energy_mwh_per_tonne = 6.25
    else:
        return 0.0

    snapshot_weights = _get_objective_snapshot_weights(network_class.n)
    if hasattr(network_class.n.loads_t, "p") and load_name in network_class.n.loads_t.p.columns:
        load_series = network_class.n.loads_t.p[load_name].abs()
    else:
        load_series = network_class.n.loads_t.p_set[load_name].abs()

    annual_energy_mwh = float((load_series * snapshot_weights).sum())
    return annual_energy_mwh / specific_energy_mwh_per_tonne * 1000.0


def _calculate_electrolyser_annual_cost(network_class, country_series):
    """
    Calculate replacement-adjusted annual electrolyser cost for a solved plant.

    Parameters
    ----------
    network_class : Network
        Solved plant network wrapper.
    country_series : pandas.Series
        Country-level plant financing inputs.

    Returns
    -------
    tuple[float, float]
        Replacement-adjusted annual electrolyser cost and weighted annual
        operating hours.
    """
    link_name = 'Electrolysis'
    p_nom_opt = float(network_class.n.links.p_nom_opt[link_name])
    annual_operating_hours = _calculate_link_operating_hours(network_class.n, link_name)
    solved_product_output_kg = _calculate_product_output_kg(network_class)
    base_annual_cost = p_nom_opt * float(network_class.n.links.at[link_name, 'capital_cost'])

    if (
        network_class.type not in {"hydrogen", "ammonia"}
        or link_name not in network_class.n.links.index
        or network_class.electrolyser_raw_capital_cost is None
        or p_nom_opt <= 0
        or annual_operating_hours <= 0
        or solved_product_output_kg <= 0
    ):
        return base_annual_cost, annual_operating_hours

    replacement_inputs = network_class.electrolyser_stack_replacement_inputs or {}
    replacement_fraction = float(replacement_inputs.get("replacement_fraction", 0.0))
    stack_lifetime_hours = float(replacement_inputs.get("stack_lifetime_hours", 0.0))
    if replacement_fraction <= 0 or stack_lifetime_hours <= 0:
        return base_annual_cost, annual_operating_hours

    raw_capex = p_nom_opt * float(network_class.electrolyser_raw_capital_cost)
    replacement_annual_cost = annualise_capex_with_hourly_replacements(
        raw_capex,
        float(country_series['Plant interest rate']),
        float(country_series['Plant lifetime (years)']),
        annual_operating_hours,
        stack_lifetime_hours,
        replacement_fraction,
    )
    return replacement_annual_cost, annual_operating_hours


def estimate_offgrid_pop(hexagons, country_series):
    """
    Estimation of offgrid population.

    Parameters
    ----------
    hexagons : pandas DataFrame
        contains geospatial information on country of interest split into 
        hexagons.
    country_series : pandas Series
        contains data from country parameters excel sheet.

    Returns
    -------
    hexagons : pandas DataFrame
        updated DataFrame containing geospatial information on country of 
        interest split into hexagons.
    """
    hexagons = hexagons.copy()
    # Calculate the off-grid population for distance_to_grid == 0 (40% of the population)
    hexagons['offgrid_population'] = hexagons['pop'].where(hexagons['grid_dist'] == 0, 0) * (100-country_series["Near grid uptake (%)"])/100

    # Calculate the total offgrid population based on the national average
    total_population = hexagons['pop'].sum()
    total_offgrid_population = total_population * (100-country_series["Electricity Access (%)"])/100

    # Sum of offgrid population for areas with distance_to_grid == 0
    offgrid_population_assigned = hexagons['offgrid_population'].sum()

    # Calculate the remaining offgrid population to be distributed
    remaining_offgrid_population = total_offgrid_population - offgrid_population_assigned

    # Get the total population for areas with distance_to_grid > 0
    remaining_population = hexagons.loc[hexagons['grid_dist'] > 0, 'pop'].sum()

    # Calculate the offgrid population for areas with distance_to_grid > 0
    hexagons.loc[hexagons['grid_dist'] > 0, 'offgrid_population'] = \
        hexagons['pop'] * (remaining_offgrid_population / remaining_population)
    
    hexagons["offgrid_households"] = hexagons["offgrid_population"] / country_series["Household size"]

    return hexagons

def _nh3_pyomo_constraints(n, snapshots):
    """Includes a series of additional constraints which make the ammonia plant work as needed:
    i) Battery sizing
    ii) Ramp hard constraints down (Cannot be violated)
    iii) Ramp hard constraints up (Cannot be violated)
    iv) Ramp soft constraints down
    v) Ramp soft constraints up
    (iv) and (v) just softly suppress ramping so that the model doesn't 'zig-zag', which looks a bit odd on operation.
    Makes very little difference on LCOA. """
    timestep = int(snakemake.config['freq'].split('H')[0])
    # The battery constraint is built here - it doesn't need a special function because it doesn't depend on time
    n.model.battery_interface = pm.Constraint(
        rule=lambda model: n.model.link_p_nom['BatteryInterfaceIn'] ==
                        n.model.link_p_nom['BatteryInterfaceOut'] /
                        n.links.efficiency["BatteryInterfaceOut"])

    # Constrain the maximum discharge of the H2 storage relative to its size
    time_step_cycle = 4/8760*timestep*0.5  # Factor 0.5 for 3 hour time step, 0.5 for oversized storage
    n.model.cycling_limit = pm.Constraint(
        rule=lambda model: n.model.link_p_nom['BatteryInterfaceOut'] ==
                        n.model.store_e_nom['CompressedH2Store'] * time_step_cycle)

    # The HB Ramp constraints are functions of time, so we need to create some pyomo sets/parameters to represent them.
    n.model.t = pm.Set(initialize=n.snapshots)
    n.model.HB_max_ramp_down = pm.Param(initialize=n.links.loc['HB'].ramp_limit_down)
    n.model.HB_max_ramp_up = pm.Param(initialize=n.links.loc['HB'].ramp_limit_up)

    # Using those sets/parameters, we can now implement the constraints...
    logging.warning('Pypsa has been overridden - Ramp rates on NH3 plant are included')
    n.model.NH3_pyomo_overwrite_ramp_down = pm.Constraint(n.model.t, rule=_nh3_ramp_down)
    n.model.NH3_pyomo_overwrite_ramp_up = pm.Constraint(n.model.t, rule=_nh3_ramp_up)
    # n.model.NH3_pyomo_penalise_ramp_down = pm.Constraint(n.model.t, rule=_penalise_ramp_down)
    # n.model.NH3_pyomo_penalise_ramp_up = pm.Constraint(n.model.t, rule=_penalise_ramp_up)

def _nh3_ramp_down(model, t):
    """Places a cap on how quickly the ammonia plant can ramp down"""
    timestep = int(snakemake.config['freq'].split('H')[0])
    if t == model.t.at(1):

        old_rate = model.link_p['HB', model.t.at(-1)]
    else:
        # old_rate = model.link_p['HB', t - 1]
        old_rate = model.link_p['HB', t - pd.Timedelta(timestep, unit = 'H')]

    return old_rate - model.link_p['HB', t] <= \
        model.link_p_nom['HB'] * model.HB_max_ramp_down
    # Note 20 is the UB of the size of the ammonia plant; essentially if x = 0 then the constraint is not active


def _nh3_ramp_up(model, t):
    """Places a cap on how quickly the ammonia plant can ramp up."""
    timestep = int(snakemake.config['freq'].split('H')[0])
    if t == model.t.at(1):
        old_rate = model.link_p['HB', model.t.at(-1)]
    else:
        # old_rate = model.link_p['HB', t - 1]
        old_rate = model.link_p['HB', t - pd.Timedelta(timestep, unit = 'H')]


    return model.link_p['HB', t] - old_rate <= \
        model.link_p_nom['HB'] * model.HB_max_ramp_up()

if __name__ == "__main__":
    country_params_filepath = snakemake.input.country_parameters
    demand_params_filepath = snakemake.input.demand_parameters
    country_params = pd.read_excel(country_params_filepath, index_col='Country')
    country_series = country_params.iloc[0]
    demand_params = pd.read_excel(demand_params_filepath, sheet_name='Demand centers', index_col='Demand center')
    demand_centers = demand_params.index

    years = snakemake.config["years_to_check"]
    weather_year = int(snakemake.wildcards.weather_year)
    end_weather_year = weather_year + years
    start_date = f'{weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'
    solver = snakemake.config['solver']
    generators = snakemake.config['generators_dict']
    hexagons = gpd.read_file(snakemake.input.hexagons)
    # Get a uniform capacity layout for all grid cells. https://atlite.readthedocs.io/en/master/ref_api.html
    cutout_filepath = f'cutouts/{snakemake.wildcards.country}_{snakemake.wildcards.weather_year}.nc'

    # Checking that the weather file is in the 'cutouts' folder
    try:
        with open(cutout_filepath) as f:
            print("Weather file found")
    except FileNotFoundError as e:
        e.strerror = "Missing file: Weather file not found. Please run weather file rule or add weather file to cutouts folder."
        raise e
    cutout = atlite.Cutout(cutout_filepath)
    layout = cutout.uniform_layout()
    profiles = {}
    len_hexagons = len(hexagons)
    has_water_limit = snakemake.config['water']['has_limit']
    water_limit_amount = snakemake.config['water']['annual_limit']
    freq = snakemake.config['freq']
    plant_type = snakemake.wildcards.plant_type
    
    if plant_type == 'copper':
        needs_grid_construction = snakemake.config["grid_construction"]
        tech_params_filepath = snakemake.input.technology_parameters
        infra_data = pd.read_excel(tech_params_filepath,
                                    sheet_name='Infra',
                                    index_col='Infrastructure')
        conversion_params_filepath = snakemake.input.conversion_parameters
        currency = snakemake.config["currency"]
        # Whether to include sizing of off-grid system for feedstock
        size_offgrid_feedstocks = snakemake.config['size_offgrid_feedstocks']
        feedstock_output = False

        if needs_grid_construction:
            # Calculate grid contruction costs
            hexagons_grid_construction = hexagons.apply(calculate_grid_construction,
                                                        args=[infra_data, country_series],
                                                        axis=1)
                
    if "hydro" in generators:
        # Creating hydropower generator profile
        print("Creating hydropower profile for network")
        location_hydro = gpd.read_file(f'data/{snakemake.wildcards.country}/{snakemake.wildcards.country}_hydropower_dams.gpkg')
        location_hydro['lat'] = location_hydro.geometry.y
        location_hydro['lon'] = location_hydro.geometry.x
        hydrobasins = gpd.read_file(glob.glob(f'data/{snakemake.wildcards.country}/*.shp')[0])
        
        runoff = cutout.hydro(
            plants=location_hydro,
            hydrobasins= hydrobasins,
            # Normalize output per unit area
            per_unit=True
            ).resample(time=freq).mean()
        
        # Efficiency of hydropower plant
        eta = snakemake.config['efficiency']['hydro']

        capacity_factor = xr.apply_ufunc(
            hydropower_potential_with_capacity,
            runoff,
            xr.DataArray(location_hydro['head'].values, dims=['plant']),
            xr.DataArray(location_hydro['capacity'].values, dims=['plant']),
            eta,
            vectorize=True,
            # Dask for parallel computation
            dask='parallelized',
            output_dtypes=[float]
            )

        hydro_hex_mapping = gpd.sjoin(location_hydro, hexagons, how='left', 
                                      predicate='within')
        num_time_steps = len(capacity_factor.time)

        hydro_profile = xr.DataArray(
            data=np.zeros((len_hexagons, num_time_steps)),
            dims=['hexagon', 'time'],
            coords={'hexagon': np.arange(len_hexagons), 
                    'time': capacity_factor.time}
            )

        for hex_index in range(len_hexagons):
            plants_in_hex = hydro_hex_mapping[hydro_hex_mapping['index_right'] == hex_index].index.tolist()
            if len(plants_in_hex) > 0:
                hex_capacity_factor = capacity_factor.sel(plant=plants_in_hex)
                plant_capacities = xr.DataArray(location_hydro.loc[plants_in_hex]['capacity'].values, dims=['plant'])

                weights = plant_capacities / plant_capacities.sum()
                weighted_avg_capacity_factor = (hex_capacity_factor * weights).sum(dim='plant')
                hydro_profile.loc[hex_index] = weighted_avg_capacity_factor
        
        profiles["hydro"] = hydro_profile

    if "geothermal" in generators:
        # Creating geothermal generator profile
        print("Creating geothermal profile for network")
        geo_locations = gpd.read_file(f'data/{snakemake.wildcards.country}/{snakemake.wildcards.country}_geothermal_plants.gpkg')
        geo_hex_mapping = gpd.sjoin(geo_locations, hexagons, how='left', 
                                      predicate='within')
        
        geothermal_profile = xr.DataArray(
            data=np.zeros((len_hexagons)),
            dims=['hexagon'],
            coords={'hexagon': np.arange(len_hexagons)}
            )
        
        for hex_index in range(len_hexagons):
            plants_in_hex = geo_hex_mapping[geo_hex_mapping['index_right'] == hex_index].index.tolist()
            if len(plants_in_hex) > 0:
                hex_capacity_factor = snakemake.config['efficiency']['geothermal']
                plant_capacities = xr.DataArray(geo_locations.loc[plants_in_hex]['capacity'].values, dims=['plant'])

                weights = plant_capacities / plant_capacities.sum()
                weighted_avg_capacity_factor = (hex_capacity_factor * weights).sum(dim='plant')
                geothermal_profile.loc[hex_index] = weighted_avg_capacity_factor
        
        profiles["geothermal"] = geothermal_profile
    
    if "nuclear" in generators:
        # Creating nuclear generator profile
        print("Creating nuclear profile for network")
        nuclear_locations = gpd.read_file(f'data/{snakemake.wildcards.country}/{snakemake.wildcards.country}_nuclear_plants.gpkg')
        nuclear_hex_mapping = gpd.sjoin(nuclear_locations, hexagons, how='left', 
                                      predicate='within')
        
        nuclear_profile = xr.DataArray(
            data=np.zeros((len_hexagons)),
            dims=['hexagon'],
            coords={'hexagon': np.arange(len_hexagons)}
            )
        
        for hex_index in range(len_hexagons):
            plants_in_hex = nuclear_hex_mapping[nuclear_hex_mapping['index_right'] == hex_index].index.tolist()
            if len(plants_in_hex) > 0:
                hex_capacity_factor = snakemake.config['efficiency']['nuclear']
                plant_capacities = xr.DataArray(nuclear_locations.loc[plants_in_hex]['capacity'].values, dims=['plant'])

                weights = plant_capacities / plant_capacities.sum()
                weighted_avg_capacity_factor = (hex_capacity_factor * weights).sum(dim='plant')
                nuclear_profile.loc[hex_index] = weighted_avg_capacity_factor
        
        profiles["nuclear"] = nuclear_profile
        
    if "solar" in generators:
        # Creating solar generator profile
        solar_profile = cutout.pv(panel= snakemake.config["panel"],
                                orientation='latitude_optimal',
                                layout = layout,
                                shapes = hexagons,
                                per_unit = True).resample(time=freq).mean()
        solar_profile = solar_profile.rename(dict(dim_0='hexagon'))
        profiles["solar"] = solar_profile
    if "wind" in generators:
        # Creating wind generator profile
        wind_profile = cutout.wind(turbine = snakemake.config["turbine"],
                                    layout = layout,
                                    shapes = hexagons,
                                    per_unit = True).resample(time=freq).mean()
        wind_profile = wind_profile.rename(dict(dim_0='hexagon'))
        profiles["wind"] = wind_profile

    # Loop through all demand centers
    for demand_center in demand_centers:
        print(f"\nOptimisation for {demand_center} begins")
        if plant_type == "copper":
            generators_list = [gen for gen in snakemake.config['generators_dict']]
            # Offgrid storage variables
            offgrid_E_costs = np.zeros(len_hexagons)
            offgrid_lcoes = np.zeros(len_hexagons)
            offgrid_battery_capacities= np.zeros(len_hexagons)
            offgrid_rectifier_capacities= np.zeros(len_hexagons)
            offgrid_inverter_capacities= np.zeros(len_hexagons)
            offgrid_generators_capacities = deepcopy(snakemake.config['generators_dict'])
            offgrid_generators_component_costs = {
                gen: [np.nan] * len_hexagons for gen in snakemake.config['generators_dict']
            }

            # Hybrid storage variables
            hybrid_E_costs = np.zeros(len_hexagons)
            hybrid_lcoes = np.zeros(len_hexagons)
            hybrid_battery_capacities= np.zeros(len_hexagons)
            hybrid_rectifier_capacities= np.zeros(len_hexagons)
            hybrid_inverter_capacities= np.zeros(len_hexagons)
            hybrid_grid_capacities= np.zeros(len_hexagons)
            hybrid_grid_purchase_costs = np.zeros(len_hexagons)
            hybrid_grid_capacity_charge_costs = np.zeros(len_hexagons)
            hybrid_grid_fixed_charge_costs = np.zeros(len_hexagons)
            hybrid_grid_construction_costs = np.zeros(len_hexagons)
            hybrid_generators_capacities = deepcopy(snakemake.config['generators_dict'])
            hybrid_generators_component_costs = {
                gen: [np.nan] * len_hexagons for gen in snakemake.config['generators_dict']
            }

            # Grid storage variables
            grid_construction_costs = np.zeros(len_hexagons)
            ongrid_totE_costs = np.zeros(len_hexagons)
            grid_lcoes = np.zeros(len_hexagons)
            grid_capacities = np.zeros(len_hexagons)

            demand_state = demand_params.loc[demand_center,'Demand state']
            conversion_params = pd.read_excel(conversion_params_filepath,
                                    sheet_name=demand_state,
                                    index_col='Parameter').squeeze('columns')
            electricity_demand = conversion_params['Electricity demand (kWh per kg product)']
        else:
            lcs = np.zeros(len_hexagons)
            generators_capacities = deepcopy(snakemake.config['generators_dict'])
            electrolyzer_capacities= np.zeros(len_hexagons)
            electrolyzer_annual_costs = np.zeros(len_hexagons)
            battery_capacities = np.zeros(len_hexagons)
            h2_storages= np.zeros(len_hexagons)

            if plant_type == "ammonia":
                nh3_storages = np.zeros(len_hexagons)
                hb_capacities = np.zeros(len_hexagons)
        
        annual_demand_quantity = float(demand_params.loc[demand_center,'Annual demand [kg/a]'])
        total_demand = annual_demand_quantity*years

        # Loop through all hexagons
        for i in range(len_hexagons):
            print(f"\nCurrently optimising {i+1} of {len_hexagons} hexagons")
            generators = deepcopy(snakemake.config['generators_dict'])
            
            # Creating demand schedule
            demand_schedule_index = pd.date_range(start_date, end_date, freq = "1H", inclusive='left')
            freq_demand_quantity = total_demand/demand_schedule_index.size
            demand_schedule = pd.DataFrame(freq_demand_quantity, index=demand_schedule_index, columns=['Demand'])
            demand_resampled_schedule = demand_schedule.resample(freq).mean()
            transport_state = hexagons[f'{demand_center} trucking state'][i]

            # Collecting information that is needed for Network set up.
            for gen in generators:
                potential = profiles[gen].sel(hexagon=i)

                gen_capacity = snakemake.config['gen_capacity'][gen]
                max_capacity = hexagons.loc[i,gen]*gen_capacity

                generators[gen].append(potential)
                generators[gen].append(max_capacity)

            if plant_type == 'copper':
                # Only optimises the demand center hexagon
                if snakemake.config['restrict_hexagons'] and transport_state != f"{demand_center} hexagon":
                    for gen in generators:
                        offgrid_generators_capacities[gen].append(np.nan)
                    offgrid_battery_capacities[i] = np.nan
                    offgrid_rectifier_capacities[i] = np.nan
                    offgrid_inverter_capacities[i] = np.nan
                    offgrid_E_costs[i] = np.nan
                    offgrid_lcoes[i] = np.nan

                    for gen in generators:
                        hybrid_generators_capacities[gen].append(np.nan)
                    hybrid_battery_capacities[i] = np.nan
                    hybrid_grid_capacities[i] = np.nan
                    hybrid_rectifier_capacities[i] = np.nan
                    hybrid_inverter_capacities[i] = np.nan
                    hybrid_E_costs[i] = np.nan
                    hybrid_lcoes[i] = np.nan
                    hybrid_grid_purchase_costs[i] = np.nan
                    hybrid_grid_capacity_charge_costs[i] = np.nan
                    hybrid_grid_fixed_charge_costs[i] = np.nan
                    hybrid_grid_construction_costs[i] = np.nan

                    grid_construction_costs[i] = np.nan
                    ongrid_totE_costs[i] = np.nan
                    grid_capacities[i] = np.nan
                    grid_lcoes[i] = np.nan
                    continue

                # Grid-only energy optimisation
                if needs_grid_construction:
                    grid_construction_costs[i] = hexagons_grid_construction[i]/total_demand

                (ongrid_E_cost,
                grid_lcoes[i],
                grid_capacities[i]) = calculate_grid_costs(demand_resampled_schedule,
                                                electricity_demand,
                                                country_series,
                                                currency)
                
                ongrid_totE_costs[i] = (ongrid_E_cost + hexagons_grid_construction[i])/total_demand
                

                # Set up the network
                network_class = Network(plant_type, generators, pypsa.Network())
                network_class.set_network(demand_resampled_schedule, country_series,
                                          demand_state, electricity_demand)
                network_class.update_generators(country_series)

                # Offgrid-only optimisation of facility
                # Solve offgrid 
                network_class.n.lopf(solver_name=solver,
                                    solver_options = {'OutputFlag': 0},
                                    pyomo=False,
                                    )

                offgrid_E_costs[i] = network_class.n.objective
                offgrid_lcoes[i] = offgrid_E_costs[i] / (network_class.n.generators_t.p.sum().sum()*1000)  # total system cost / energy produced (currency/kWh)
                offgrid_renewable_components = _calculate_renewable_objective_costs(
                    network_class.n,
                    generators_list,
                    annual_demand_quantity,
                )
                for gen, per_kg_cost in offgrid_renewable_components.items():
                    offgrid_generators_component_costs[gen][i] = per_kg_cost
                for gen in generators:
                    offgrid_generators_capacities[gen].append(network_class.n.generators.p_nom_opt[gen.capitalize()])
                offgrid_battery_capacities[i] = network_class.n.storage_units.p_nom_opt['Battery']
                offgrid_rectifier_capacities[i] = network_class.n.links.p_nom_opt['Rectifier']
                offgrid_inverter_capacities[i] = network_class.n.links.p_nom_opt['Inverter']

                # Hybrid energy optimisation of facility
                network_class.add_grid(country_series, currency)

                # Solve hybrid
                network_class.n.lopf(solver_name=solver,
                                    solver_options = {'OutputFlag': 0},
                                    pyomo=False,
                                    )
                
                hybrid_E_costs[i] = network_class.n.objective + country_series[f"Electricity fixed charge ({currency}/year)"]
                # Hybrid LCOE = total system cost / energy produced (currency/kWh)
                hybrid_lcoes[i] = hybrid_E_costs[i] / (network_class.n.generators_t.p.sum().sum()*1000)
                hybrid_renewable_components = _calculate_renewable_objective_costs(
                    network_class.n,
                    generators_list,
                    annual_demand_quantity,
                )
                for gen, per_kg_cost in hybrid_renewable_components.items():
                    hybrid_generators_component_costs[gen][i] = per_kg_cost
                for gen in generators:
                        hybrid_generators_capacities[gen].append(network_class.n.generators.p_nom_opt[gen.capitalize()])
                hybrid_battery_capacities[i] = network_class.n.storage_units.p_nom_opt['Battery']
                hybrid_rectifier_capacities[i] = network_class.n.links.p_nom_opt['Rectifier']
                hybrid_inverter_capacities[i] = network_class.n.links.p_nom_opt['Inverter']
                hybrid_grid_capacities[i] = network_class.n.generators.p_nom_opt['Grid']
                if hasattr(network_class.n.snapshot_weightings, 'columns'):
                    if 'objective' in network_class.n.snapshot_weightings.columns:
                        snapshot_weights = network_class.n.snapshot_weightings['objective']
                    elif 'generators' in network_class.n.snapshot_weightings.columns:
                        snapshot_weights = network_class.n.snapshot_weightings['generators']
                    else:
                        snapshot_weights = network_class.n.snapshot_weightings.iloc[:, 0]
                else:
                    snapshot_weights = network_class.n.snapshot_weightings
                grid_dispatch = network_class.n.generators_t.p['Grid']
                hybrid_grid_purchase_annual_cost = (
                    (grid_dispatch * snapshot_weights).sum()
                    * network_class.n.generators.at['Grid', 'marginal_cost']
                )
                hybrid_grid_capacity_annual_cost = (
                    network_class.n.generators.p_nom_opt['Grid']
                    * network_class.n.generators.at['Grid', 'capital_cost']
                )
                hybrid_grid_fixed_annual_cost = country_series[f"Electricity fixed charge ({currency}/year)"]
                hybrid_grid_purchase_costs[i] = hybrid_grid_purchase_annual_cost / annual_demand_quantity
                hybrid_grid_capacity_charge_costs[i] = hybrid_grid_capacity_annual_cost / annual_demand_quantity
                hybrid_grid_fixed_charge_costs[i] = hybrid_grid_fixed_annual_cost / annual_demand_quantity
                hybrid_grid_construction_costs[i] = grid_construction_costs[i]

                # Calculate offgrid system for feedstock centers if neccessary
                if (size_offgrid_feedstocks) and feedstock_output != True:
                    feedstock_points_gdf = gpd.read_file(f"resources/feedstocks_transport_{snakemake.wildcards.country}_{snakemake.wildcards.plant_type}.geojson")
                    feedstocks_amount= len(feedstock_points_gdf)
                    
                    fs_lcoes = np.empty(feedstocks_amount)
                    fs_generators_capacities = deepcopy(snakemake.config['generators_dict'])
                    fs_battery_capacities = np.empty(feedstocks_amount)
                    
                    feedstock_index = pd.date_range(start_date, end_date, freq=freq, inclusive='left')

                    for c in feedstock_points_gdf.index:
                        hex_id = feedstock_points_gdf.loc[c, "nearest hexidx"]
                        hourly_demand = (feedstock_points_gdf.loc[c, "Electrical Energy (MWh/year)"]*1000)/(365*24)  # kwh
                        feedstock_demand_schedule = pd.DataFrame(hourly_demand,
                                                        index=feedstock_index, columns = ['Demand'])
                        fs_demand_resampled = feedstock_demand_schedule.resample(freq).mean()

                        # Set up the network
                        fs_network_class = Network(plant_type, generators, pypsa.Network())
                        fs_network_class.set_network(fs_demand_resampled, 
                                                     country_series, 
                                                     "Feedstock",
                                                     1.0,)
                        fs_network_class.update_generators(country_series)

                        # Solve feedstocks offgrid system
                        fs_network_class.n.lopf(solver_name=solver,
                                            solver_options = {'OutputFlag': 0},
                                            pyomo=False,
                                            )
                        
                        # LCOE = total system cost / energy produced (currency/kWh)
                        fs_lcoes[c] = fs_network_class.n.objective / (fs_network_class.n.generators_t.p.sum().sum()*1000) 
                        fs_battery_capacities[c] = fs_network_class.n.storage_units.p_nom_opt['Battery']
                        for gen in generators:
                            fs_generators_capacities[gen].append(fs_network_class.n.generators.p_nom_opt[gen.capitalize()])
                        
                    # Store off-grid feedstock capacities and cost
                    for gen, capacities in fs_generators_capacities.items():
                        feedstock_points_gdf[f'Offgrid {gen} capacity (MW)'] = capacities
                    feedstock_points_gdf['Offgrid battery capacity (MW)'] = fs_battery_capacities
                    feedstock_points_gdf[f"Levlised cost of energy ({currency}/kg/year)"] = fs_lcoes
                    
                    # Save feedstocks information to a separate file
                    fs_output_path = f"results/feedstocks_{snakemake.wildcards.country}_{snakemake.wildcards.plant_type}.geojson"
                    feedstock_points_gdf.to_file(fs_output_path, driver='GeoJSON', encoding='utf-8')
                    feedstock_output = True
            else:
                # Only optimises the demand center hexagon
                if snakemake.config['restrict_hexagons'] and transport_state != f"{demand_center} hexagon":
                    for gen in generators:
                            generators_capacities[gen].append(np.nan)
                    electrolyzer_capacities[i] = np.nan
                    electrolyzer_annual_costs[i] = np.nan
                    battery_capacities[i] = np.nan
                    h2_storages[i] = np.nan
                    lcs[i] = np.nan
                    
                    if plant_type == 'ammonia':
                        nh3_storages[i] = np.nan
                        hb_capacities[i] = np.nan
                    continue

                # Check for water constraint before any solving occurs
                if has_water_limit == True:
                    if plant_type == "hydrogen":
                        # Check if hydrogen demand can be met based on water 
                        # availability kg H2 per cubic meter of water
                        water_constraint =  annual_demand_quantity > water_limit_amount * 111.57
                    elif plant_type == "ammonia":
                        # Total hydrogen demand in kg
                        # Converts kg ammonia to kg H2
                        total_hydrogen_demand = annual_demand_quantity * 17 / 3  
                        # Check if hydrogen demand can be met based on water 
                        # availability kg H2 per cubic meter of water
                        water_constraint = total_hydrogen_demand > water_limit_amount * 111.57
                        
                    # Note: this constraint is purely stoichiometric
                    # - more water may be needed for cooling or other processes
                    if water_constraint == True:
                        print('Not enough water to meet demand!')
                        lcs[i] = np.nan
                        electrolyzer_capacities[i] = np.nan
                        electrolyzer_annual_costs[i] = np.nan
                        battery_capacities[i] = np.nan
                        h2_storages[i] = np.nan
                        for gen in generators:
                            generators_capacities[gen].append(np.nan)
                        if plant_type == "ammonia": 
                            nh3_storages[i] = np.nan
                            hb_capacities[i] = np.nan
                        continue
                
                # Set up the network
                network_class = Network(plant_type, generators)
                network_class.set_network(demand_resampled_schedule, country_series)
                network_class.update_generators(country_series)

                if plant_type == "hydrogen":
                    # Solve hydrogen offgrid
                    network_class.n.lopf(solver_name=solver,
                        solver_options = {'LogToConsole':0, 'OutputFlag':0},
                        pyomo=False,
                        )
                    
                    h2_storages[i] = network_class.n.stores.e_nom_opt['Compressed H2 Store']
                    battery_capacities[i] = network_class.n.storage_units.p_nom_opt['Battery']
                elif plant_type == "ammonia":
                    # Solve ammonia offgrid
                    network_class.n.lopf(solver_name=solver,
                        solver_options={'LogToConsole': 0, 'OutputFlag': 0},
                        pyomo=True,
                        extra_functionality=_nh3_pyomo_constraints)
                    
                    h2_storages[i] = network_class.n.stores.e_nom_opt['CompressedH2Store']
                    nh3_storages[i] = network_class.n.stores.e_nom_opt['Ammonia']
                    hb_capacities[i] = network_class.n.links.p_nom_opt['HB']
                    battery_capacities[i] = network_class.n.stores.e_nom_opt['Battery']

                electrolyzer_capacities[i] = network_class.n.links.p_nom_opt['Electrolysis']
                electrolyzer_annual_costs[i], _ = _calculate_electrolyser_annual_cost(
                    network_class,
                    country_series,
                )
                base_electrolyzer_annual_cost = (
                    electrolyzer_capacities[i]
                    * network_class.n.links.at['Electrolysis', 'capital_cost']
                )
                lcs[i] = (
                    network_class.n.objective
                    + (electrolyzer_annual_costs[i] - base_electrolyzer_annual_cost)
                ) / total_demand
                
                for gen in generators:
                    generators_capacities[gen].append(network_class.n.generators.p_nom_opt[gen.capitalize()])
    
        print("\nOptimisation complete.\n")    
        # Updating results in hexagon file
        if plant_type == 'copper':
            # Off-grid costs
            for gen, capacities in offgrid_generators_capacities.items():
                hexagons[f'{demand_center} offgrid {gen} capacity (MW)'] = capacities
            hexagons[f'{demand_center} offgrid rectifier capacity (MW)'] = offgrid_rectifier_capacities
            hexagons[f'{demand_center} offgrid inverter capacity (MW)'] = offgrid_inverter_capacities
            hexagons[f'{demand_center} offgrid battery capacity (MW)'] = offgrid_battery_capacities
            hexagons[f'{demand_center} offgrid total energy cost ({currency}/kg/year)'] = offgrid_E_costs / annual_demand_quantity
            hexagons[f'{demand_center} offgrid lcoe ({currency}/kWh)'] = offgrid_lcoes
            for gen, component_costs in offgrid_generators_component_costs.items():
                hexagons[f'{demand_center} offgrid {gen} objective component cost ({currency}/kg/year)'] = component_costs
            
            # Hybrid-grid costs
            for gen, capacities in hybrid_generators_capacities.items():
                hexagons[f'{demand_center} hybrid {gen} capacity (MW)'] = capacities
            hexagons[f'{demand_center} hybrid rectifier capacity (MW)'] = hybrid_rectifier_capacities
            hexagons[f'{demand_center} hybrid inverter capacity (MW)'] = hybrid_inverter_capacities
            hexagons[f'{demand_center} hybrid battery capacity (MW)'] = hybrid_battery_capacities
            hexagons[f'{demand_center} hybrid grid capacity (MW)'] = hybrid_grid_capacities
            hexagons[f'{demand_center} hybrid total energy cost ({currency}/kg/year)'] = (hybrid_E_costs + grid_construction_costs) / annual_demand_quantity
            hexagons[f'{demand_center} hybrid grid purchase costs ({currency}/kg/year)'] = hybrid_grid_purchase_costs
            hexagons[f'{demand_center} hybrid grid capacity charge ({currency}/kg/year)'] = hybrid_grid_capacity_charge_costs
            hexagons[f'{demand_center} hybrid grid fixed charge ({currency}/kg/year)'] = hybrid_grid_fixed_charge_costs
            hexagons[f'{demand_center} hybrid grid construction charge ({currency}/kg/year)'] = hybrid_grid_construction_costs
            hexagons[f'{demand_center} hybrid lcoe ({currency}/kWh)'] = hybrid_lcoes
            for gen, component_costs in hybrid_generators_component_costs.items():
                hexagons[f'{demand_center} hybrid {gen} objective component cost ({currency}/kg/year)'] = component_costs
            
            # Grid costs
            hexagons[f'{demand_center} grid construction ({currency}/kg/year)'] = grid_construction_costs
            hexagons[f'{demand_center} grid total energy cost ({currency}/kg/year)'] = ongrid_totE_costs 
            hexagons[f'{demand_center} grid capacity (MW)'] = grid_capacities
            hexagons[f'{demand_center} grid lcoe ({currency}/kWh)'] = grid_lcoes
        else:
            for gen, capacities in generators_capacities.items():
                hexagons[f'{demand_center} {gen} capacity'] = capacities
            hexagons[f'{demand_center} electrolyzer capacity'] = electrolyzer_capacities
            hexagons[f'{demand_center} electrolyzer annual costs'] = electrolyzer_annual_costs
            hexagons[f'{demand_center} battery capacity'] = battery_capacities
            hexagons[f'{demand_center} H2 storage capacity'] = h2_storages
            hexagons[f'{demand_center} production cost'] = lcs
            if plant_type == "ammonia":        
                hexagons[f'{demand_center} NH3 storage capacity'] = nh3_storages
                hexagons[f'{demand_center} HB capacity'] = hb_capacities

    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')

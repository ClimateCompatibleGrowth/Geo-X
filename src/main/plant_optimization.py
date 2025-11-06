"""
@authors:
 - Claire Halloran
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk
 - Lukas Schirren, Imperial College London, lukas.schirren@imperial.ac.uk
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
from scipy.constants import physical_constants
import xarray as xr
import glob

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

def get_h2_results(n, generators):
    '''
    Get final results from network optimisation

    ...
    Parameters
    ----------
    n : network Object
        network.
    generators : dictionary
        contains types of generators that this plant uses.

    Returns
    -------
    lc : float
        levelized cost per kg commodity.
    generator_capactities : dictionary
        contains each generator with their optimal capacity in MW.
    electrolyzer_capacity : float
        optimal electrolyzer capacity in MW.
    battery_capacity : float
        optimal battery storage capacity in MW/MWh (1 hour batteries).
    h2_storage : float
        optimal hydrogen storage capacity in MWh.
    '''
    generator_capacities = {}
    # Convert back to kg H2
    lc = n.objective/(n.loads_t.p_set.sum()[0]/39.4*1000)
    for generator in generators:
            generator_capacities[generator] = n.generators.p_nom_opt[generator.capitalize()]
    electrolyzer_capacity = n.links.p_nom_opt['Electrolysis']
    battery_capacity = n.storage_units.p_nom_opt['Battery']
    h2_storage = n.stores.e_nom_opt['Compressed H2 Store']
    
    return lc, generator_capacities, electrolyzer_capacity, battery_capacity, h2_storage

def get_nh3_results(n, generators):
    '''
    Get final results from network optimisation

    ...
    Parameters
    ----------
    n : network Object
        network.
    generators : list
        contains types of generators that this plant uses.

    Returns
    -------
    lc : float
        levelized cost per kg commodity.
    generator_capactities : dictionary
        contains each generator with their optimal capacity in MW.
    electrolyzer_capacity : float
        optimal electrolyzer capacity in MW.
    battery_capacity : float
        optimal battery storage capacity in MW/MWh (1 hour batteries).
    h2_storage : float
        optimal hydrogen storage capacity in MWh.
    nh3_storage : float
        optimal ammonia storage capacity in MWh.
    hb_capacity : float
        optimal Haber-Bosch capacity in MWh.
    '''
    generator_capacities = {}
    lc = n.objective / ((n.loads_t.p_set['Ammonia demand'] * n.snapshot_weightings[
        'objective']).sum() / 6.25 * 1000)  # convert back to kg NH3
    for generator in generators:
            generator_capacities[generator] = n.generators.p_nom_opt[generator.capitalize()]
    electrolyzer_capacity = n.links.p_nom_opt['Electrolysis']
    battery_capacity = n.stores.e_nom_opt['Battery']
    h2_storage = n.stores.e_nom_opt['CompressedH2Store']
    # !!! need to save ammonia storage capacity as well
    nh3_storage = n.stores.e_nom_opt['Ammonia']
    hb_capacity = n.links.p_nom_opt['HB']

    return lc, generator_capacities, electrolyzer_capacity, battery_capacity, h2_storage, nh3_storage, hb_capacity

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
    """Places a cap on how quickly the ammonia plant can ramp down"""
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
    demand_params = pd.read_excel(demand_params_filepath, index_col='Demand center')
    demand_centers = demand_params.index

    years = snakemake.config["years_to_check"]
    weather_year = int(snakemake.wildcards.weather_year)
    end_weather_year = weather_year + years
    start_date = f'{weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'
    solver = snakemake.config['solver']
    generators = snakemake.config['generators_dict']
    hexagons = gpd.read_file(snakemake.input.hexagons)
    needs_pipeline_construction = snakemake.config['transport']['pipeline_construction']

    # Get a uniform capacity layout for all grid cells. https://atlite.readthedocs.io/en/master/ref_api.html
    cutout_filepath = f'cutouts/{snakemake.wildcards.country}_{snakemake.wildcards.weather_year}.nc'

    # Checking that the weather file is in the 'cutouts' folder
    try:
        with open(cutout_filepath) as f:
            print("Weather file found")
    except FileNotFoundError as e:
        e.strerror = "Weather file not found. Please run weather file rule or add weather file to cutouts folder. Missing file"
        raise e
    
    cutout = atlite.Cutout(cutout_filepath)
    layout = cutout.uniform_layout()
    profiles = []
    len_hexagons = len(hexagons)
    has_water_limit = snakemake.config['water']['has_limit']
    water_limit_amount = snakemake.config['water']['annual_limit']
    freq = snakemake.config['freq']
    plant_type = snakemake.wildcards.plant_type
    
    # Creating hydropower generator profile
    if "hydro" in generators:
        print("Creating hydropower profile for network")
        location_hydro = gpd.read_file(f'data/{snakemake.wildcards.country}/hydro/{snakemake.wildcards.country}_hydropower_dams.gpkg')
        location_hydro['lat'] = location_hydro.geometry.y
        location_hydro['lon'] = location_hydro.geometry.x
        hydrobasins = gpd.read_file(glob.glob(f'data/{snakemake.wildcards.country}/hydro/*.shp')[0])
        
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

    # Creating geothermal generator profile
    if "geothermal" in generators:
        print("Creating geothermal profile for network")
        geo_locations = gpd.read_file(f'data/{snakemake.wildcards.country}/geothermal/{snakemake.wildcards.country}_geothermal_plants.gpkg')
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
    
    # Creating solar generator profile
    if "solar" in generators:
        solar_profile = cutout.pv(panel= snakemake.config["panel"],
                                orientation='latitude_optimal',
                                layout = layout,
                                shapes = hexagons,
                                per_unit = True).resample(time=freq).mean()
        solar_profile = solar_profile.rename(dict(dim_0='hexagon'))
    
    # Creating wind generator profile
    if "wind" in generators:
        wind_profile = cutout.wind(turbine = snakemake.config["turbine"],
                                    layout = layout,
                                    shapes = hexagons,
                                    per_unit = True).resample(time=freq).mean()
        wind_profile = wind_profile.rename(dict(dim_0='hexagon'))

    # Adding generators into profiles list
    for gen in generators:
        # profile = f'{gen}_profile'
        if gen == 'solar':
            profiles.append(solar_profile)
        elif gen == 'wind':
            profiles.append(wind_profile)
        elif gen == 'hydro':
            profiles.append(hydro_profile)
        elif gen == 'geothermal':
            profiles.append(geothermal_profile)

    # Loop through all demand centers
    for demand_center in demand_centers:
        print(f"\nOptimisation for {demand_center} begins...")
        # Store results
        lcs = np.zeros(len_hexagons)
        generators_capacities = deepcopy(snakemake.config['generators_dict'])
        electrolyzer_capacities= np.zeros(len_hexagons)
        battery_capacities = np.zeros(len_hexagons)
        h2_storages= np.zeros(len_hexagons)

        # Ammonia extra storage variables
        if plant_type == "ammonia":
            nh3_storages = np.zeros(len_hexagons)
            hb_capacities = np.zeros(len_hexagons)

        annual_demand_quantity = float(demand_params.loc[demand_center,'Annual demand [kg/a]'])

        # Loop through all hexagons
        for i in range(len_hexagons):
            print(f"\nCurrently optimising {i+1} of {len_hexagons} hexagons...")
            trucking_state = hexagons.loc[i, f'{demand_center} trucking state']
            pipeline_transport_costs = hexagons.loc[i, f'{demand_center} pipeline transport costs']
            gen_index = 0
            generators = deepcopy(snakemake.config['generators_dict'])
            
            # Creating demand schedule
            demand_schedule_index = pd.date_range(start_date, end_date, freq = "1H", inclusive='left')
            freq_demand_quantity = (annual_demand_quantity*years)/demand_schedule_index.size
            demand_schedule = pd.DataFrame(freq_demand_quantity, index=demand_schedule_index, columns=['Demand'])
            demand_resampled_schedule = demand_schedule.resample(freq).mean()
            
            # Collecting information that is needed for Network set up.
            for gen in generators:
                potential = profiles[gen_index].sel(hexagon=i)

                gen_capacity = snakemake.config['gen_capacity'][gen]
                if gen == "wind":
                    max_capacity = hexagons.loc[i,'theo_turbines']*gen_capacity
                elif gen == "solar":
                    max_capacity = hexagons.loc[i,'theo_pv']*gen_capacity
                # -- Need to change solar and wind names in data-prep to have the below only
                else:
                    max_capacity = hexagons.loc[i,gen]*gen_capacity

                generators[gen].append(potential)
                generators[gen].append(max_capacity)
                gen_index += 1

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
            network_class.set_generators_in_network(country_series)

            if plant_type == "hydrogen":
                # Solve
                network_class.n.lopf(solver_name=solver,
                    solver_options = {'LogToConsole':0, 'OutputFlag':0},
                    pyomo=False,
                    )
                
                h2_storages[i] = network_class.n.stores.e_nom_opt['CompressedH2Store']
            elif plant_type == "ammonia":
                # Solve 
                network_class.n.lopf(solver_name=solver,
                    solver_options={'LogToConsole': 0, 'OutputFlag': 0},
                    pyomo=True,
                    extra_functionality=_nh3_pyomo_constraints)
                
                h2_storages[i] = network_class.n.stores.e_nom_opt['CompressedH2Store']
                # !!! need to save ammonia storage capacity as well
                nh3_storages[i] = network_class.n.stores.e_nom_opt['Ammonia']
                hb_capacities[i] = network_class.n.links.p_nom_opt['HB']

            lcs[i] = network_class.n.objective/annual_demand_quantity*years
            electrolyzer_capacities[i] = network_class.n.links.p_nom_opt['Electrolysis']
            battery_capacities[i] = network_class.n.stores.e_nom_opt['Battery']

            for gen in generators:
                generators_capacities[gen].append(network_class.n.generators.p_nom_opt[gen.capitalize()])
        
        print("\nOptimisation complete.\n")        
        # Updating results in hexagon file
        for gen, capacities in generators_capacities.items():
            hexagons[f'{demand_center} {gen} capacity'] = capacities
        hexagons[f'{demand_center} electrolyzer capacity'] = electrolyzer_capacities
        hexagons[f'{demand_center} battery capacity'] = battery_capacities
        hexagons[f'{demand_center} H2 storage capacity'] = h2_storages
        hexagons[f'{demand_center} production cost'] = lcs
        if plant_type == "ammonia":        
            hexagons[f'{demand_center} NH3 storage capacity'] = nh3_storages
            hexagons[f'{demand_center} HB capacity'] = hb_capacities

    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')
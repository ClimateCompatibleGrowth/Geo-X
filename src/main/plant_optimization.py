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

def get_demand_schedule(quantity, start_date, end_date, transport_state, transport_params_filepath, freq):
    '''
    Calculates hourly commodity demand for truck shipment and pipeline transport.

    ...
    Parameters
    ----------
    quantity : float
        annual amount of commodity to transport in kilograms.
    start_date : string
        start date for demand schedule in the format YYYY-MM-DD.
    end_date : string
        end date for demand schedule in the format YYYY-MM-DD.
    transport_state : string
        state commodity is transported in, one of '500 bar', 'LH2', 'LOHC', 
        or 'NH3'.
    transport_params_filepath : string
        path to transport_parameters.xlsx file.
    freq : string
        frequency that the schedules should be in.

    Returns
    -------
    trucking_hourly_demand_schedule : pandas DataFrame
        hourly demand profile for trucking transport.
    pipeline_hourly_demand_schedule : pandas DataFrame
        hourly demand profile for pipeline transport.
    '''
    # Schedule for pipeline
    index = pd.date_range(start_date, end_date, freq = freq)
    pipeline_hourly_quantity = quantity/index.size
    pipeline_hourly_demand_schedule = pd.DataFrame(pipeline_hourly_quantity, index=index,  columns=['Demand'])
    # Resample pipeline schedule
    pipeline_demand_resampled_schedule = pipeline_hourly_demand_schedule.resample(freq).mean()

    # NaN is where there is no road access and no construction so hexagon 
    # is infeasible for trucking
    if pd.isnull(transport_state):
        trucking_hourly_demand_schedule = np.nan
    # If demand center is in hexagon
    elif transport_state=="None":
        # A demand schedule is required
        # Even though no trucking will happen, we keep the name for consistency
        # as only one hexagon will have no trucking happening in it
        annual_deliveries = 365*24
        trucking_hourly_demand = quantity/annual_deliveries
        index = pd.date_range(start_date, end_date, periods=annual_deliveries)
        trucking_demand_schedule = pd.DataFrame(trucking_hourly_demand, index=index, columns=['Demand'])
        # First create hourly schedule
        trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)
        # Then resample to desired frequency using mean
        trucking_demand_resampled_schedule = trucking_hourly_demand_schedule.resample(freq).mean()
    else:
        transport_params = pd.read_excel(transport_params_filepath,
                                            sheet_name = transport_state,
                                            index_col = 'Parameter'
                                            ).squeeze('columns')

        truck_capacity = transport_params['Net capacity (kg of commodity)']

        # Schedule for trucking
        annual_deliveries = quantity/truck_capacity
        quantity_per_delivery = quantity/annual_deliveries
        index = pd.date_range(start_date, end_date, periods=annual_deliveries)
        trucking_demand_schedule = pd.DataFrame(quantity_per_delivery, index=index, columns=['Demand'])
        # First create hourly schedule
        trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)
        # Then resample to desired frequency using mean
        trucking_demand_resampled_schedule = trucking_hourly_demand_schedule.resample(freq).mean()

    return trucking_demand_resampled_schedule, pipeline_demand_resampled_schedule

def get_water_constraint(n, demand_profile, water_limit): 
    '''
    Calculates the water constraint.

    ...
    Parameters
    ----------
    n : network Object
        network.
    demand_profile : pandas DataFrame
        hourly dataframe of commodity demand in kg.
    water_limit : float
        annual limit on water available for electrolysis in hexagon, in cubic meters.

    Returns
    -------
    water_constraint : boolean
        whether there is a water constraint or not.
    '''
    if n.type == "hydrogen":
        # Total hydrogen demand in kg
        total_hydrogen_demand = demand_profile['Demand'].sum()
        # Check if hydrogen demand can be met based on hexagon water availability
        # kg H2 per cubic meter of water
        water_constraint =  total_hydrogen_demand <= water_limit * 111.57
    elif n.type == "ammonia":
        # Total ammonia demand in kg
        total_ammonia_demand = (
                    (n.loads_t.p_set['Ammonia demand'] * n.snapshot_weightings['objective']).sum() / 6.25 * 1000)
        # Total hydrogen demand in kg
        # Converts kg ammonia to kg H2
        total_hydrogen_demand = total_ammonia_demand * 17 / 3  
        # Check if hydrogen demand can be met based on hexagon water availability
        # kg H2 per cubic meter of water
        water_constraint = total_hydrogen_demand <= water_limit * 111.57
        
        # Note: this constraint is purely stoichiometric
        # - more water may be needed for cooling or other processes
    
    return water_constraint

def get_generator_profile(generator, cutout, layout, hexagons, freq):
    '''
    Determines the generation profile of the specified generator in the cutout based on weather data.

    ...
    Parameters
    ----------
    generator : string
        the name of the generator type to be used.
        (i.e., "solar" for pv, "wind" for wind turbines)
    cutout : atlite cutout object
        a spatial and temporal subset of ERA-5 weather data consisting of grid 
        cells.
        https://atlite.readthedocs.io/en/latest/introduction.html
    layout : xarray DataArray
        the capacity to be built in each of the grid_cells.
        https://atlite.readthedocs.io/en/master/ref_api.html
    hexagons : geodataframe
        hexagon file outputted after the water_cost file is ran.
    freq : string
        frequency that the profiles should be in.
    
    Returns
    -------
    profile : xarray DataArray
        A profile where weather data has been converted into a generation
        time-series.
        https://atlite.readthedocs.io/en/master/ref_api.html
    '''
    if generator == "solar":
        # The panel string should be in the config file as well in case people want to change that in the prep and main.
        # Alycia to double-check that the CSi is 1MW
        profile = cutout.pv(panel= str(snakemake.config["panel"]),
                            orientation='latitude_optimal',
                            layout = layout,
                            shapes = hexagons,
                            per_unit = True).resample(time=freq).mean()
        profile = profile.rename(dict(dim_0='hexagon'))
    elif generator == "wind":
        # This string is what we should need to put in the config file (turbine) for both data prep, replacing the constant 4, replacing here.
        profile = cutout.wind(turbine = str(snakemake.config["turbine"]),
                            layout = layout,
                            shapes = hexagons,
                            per_unit = True).resample(time=freq).mean()
        profile = profile.rename(dict(dim_0='hexagon'))
    
    return profile

def solve_model(network_class, solver):
    '''
    Solves model using the provided solver.

    ...
    Parameters
    ----------
    n : network Object
        network.
    solver : string
        name of solver to be used.
    '''
    if network_class.type == "hydrogen":
        network_class.n.lopf(solver_name=solver,
            solver_options = {'LogToConsole':0, 'OutputFlag':0},
            pyomo=False,
            )
    elif network_class.type == "ammonia":
        network_class.n.lopf(solver_name=solver,
            solver_options={'LogToConsole': 0, 'OutputFlag': 0},
            pyomo=True,
            extra_functionality=_nh3_pyomo_constraints,
            )

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
    timestep = int(snakemake.config['freq'][0])
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
    timestep = int(snakemake.config['freq'][0])
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
    timestep = int(snakemake.config['freq'][0])
    if t == model.t.at(1):
        old_rate = model.link_p['HB', model.t.at(-1)]
    else:
        # old_rate = model.link_p['HB', t - 1]
        old_rate = model.link_p['HB', t - pd.Timedelta(timestep, unit = 'H')]


    return model.link_p['HB', t] - old_rate <= \
        model.link_p_nom['HB'] * model.HB_max_ramp_up()

if __name__ == "__main__":
    transport_params_filepath = str(snakemake.input.transport_parameters)
    country_params_filepath = str(snakemake.input.country_parameters)
    demand_params_filepath = str(snakemake.input.demand_parameters )
    country_params = pd.read_excel(country_params_filepath, index_col='Country')
    country_series = country_params.iloc[0]
    demand_params = pd.read_excel(demand_params_filepath, index_col='Demand center')
    demand_centers = demand_params.index

    weather_year = int(snakemake.wildcards.weather_year)
    end_weather_year = int(snakemake.wildcards.weather_year)+int(snakemake.config["years_to_check"])
    start_date = f'{weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'
    solver = str(snakemake.config['solver'])
    generators = dict(snakemake.config['generators_dict'])
    hexagons = gpd.read_file(str(snakemake.input.hexagons))
    pipeline_construction = True # snakemake config

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
    water_limit = bool(snakemake.config['water_limit'])
    freq = str(snakemake.config['freq'])
    
    # Creating hydropower generator profile
    if "hydro" in generators:
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
        eta = 0.75

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
        num_hexagons = len(hexagons)
        num_time_steps = len(capacity_factor.time)

        hydro_profile = xr.DataArray(
            data=np.zeros((num_hexagons, num_time_steps)),
            dims=['hexagon', 'time'],
            coords={'hexagon': np.arange(num_hexagons), 
                    'time': capacity_factor.time}
            )

        for hex_index in range(num_hexagons):
            plants_in_hex = hydro_hex_mapping[hydro_hex_mapping['index_right'] == hex_index].index.tolist()
            if len(plants_in_hex) > 0:
                hex_capacity_factor = capacity_factor.sel(plant=plants_in_hex)
                plant_capacities = xr.DataArray(location_hydro.loc[plants_in_hex]['capacity'].values, dims=['plant'])

                weights = plant_capacities / plant_capacities.sum()
                weighted_avg_capacity_factor = (hex_capacity_factor * weights).sum(dim='plant')
                hydro_profile.loc[hex_index] = weighted_avg_capacity_factor

    # Adding generators into profiles list
    for gen in generators:
         if gen == "hydro":
            profiles.append(hydro_profile)         
         else:
            profiles.append(get_generator_profile(gen, cutout, layout, hexagons, freq))
    
    times = profiles[0].time
    plant_type = str(snakemake.wildcards.plant_type)

    # Loop through all demand centers
    for demand_center in demand_centers:
        print(f"\nOptimisation for {demand_center} begins...")
        # Store trucking results
        trucking_lcs = np.zeros(len_hexagons)
        t_generators_capacities = deepcopy(snakemake.config['generators_dict'])
        t_electrolyzer_capacities= np.zeros(len_hexagons)
        t_battery_capacities = np.zeros(len_hexagons)
        t_h2_storages= np.zeros(len_hexagons)
        
        # Store pipeline variables
        pipeline_lcs = np.zeros(len_hexagons)
        p_generators_capacities = deepcopy(snakemake.config['generators_dict'])
        p_electrolyzer_capacities= np.zeros(len_hexagons)
        p_battery_capacities = np.zeros(len_hexagons)
        p_h2_storages= np.zeros(len_hexagons)

        # Ammonia extra storage variables
        if plant_type == "ammonia":
            t_nh3_storages = np.zeros(len_hexagons)
            p_nh3_storages = np.zeros(len_hexagons)
            t_hb_capacities = np.zeros(len_hexagons)
            p_hb_capacities = np.zeros(len_hexagons)

        annual_demand_quantity = demand_params.loc[demand_center,'Annual demand [kg/a]']

        # Loop through all hexagons
        for i in range(len_hexagons):
            print(f"\nCurrently optimising {i+1} of {len(hexagons)} hexagons...")
            trucking_state = hexagons.loc[i, f'{demand_center} trucking state']
            gen_index = 0
            generators = deepcopy(snakemake.config['generators_dict'])
            
            # Get the demand schedule for both pipeline and trucking transport
            trucking_demand_schedule, pipeline_demand_schedule =\
                get_demand_schedule(annual_demand_quantity,
                                start_date,
                                end_date,
                                trucking_state,
                                transport_params_filepath,
                                freq)
            
            # Get the max capacity for each generation type
            for gen in generators:
                if plant_type == "hydrogen":
                    potential = profiles[gen_index].sel(hexagon = i)
                elif plant_type == "ammonia":
                    potential = profiles[gen_index].sel(hexagon=i, time=trucking_demand_schedule.index)
                # -- Eventually make a for loop - we can change the theo_turbines name to be Wind
                gen_capacity = int(snakemake.config['gen_capacity'][f'{gen}'])
                if gen == "wind":
                    # -- We'll need to remove this hard-coded 4 eventually CONFIG FILE - 4 MW turbine in spatial data prep
                    max_capacity = hexagons.loc[i,'theo_turbines']*gen_capacity
                elif gen == "solar":
                    max_capacity = hexagons.loc[i,'theo_pv']*gen_capacity
                ##### elif added newly for hydro, need to double check the column name
                elif gen == "hydro":
                    max_capacity = hexagons.loc[i,'hydro']*gen_capacity
                # -- Eventually move loops to something like this so we don't have ifs - max_capacity = hexagons.loc[i, gen] * SNAKEMAKE_CONFIG_GEN_SIZE
                
                generators[gen].append(potential)
                generators[gen].append(max_capacity)
                gen_index += 1

            # For each transport type, set up the network and solve
            transport_types = ["trucking", "pipeline"]
            for transport in transport_types:
                network = Network(plant_type, generators)

                # If trucking is viable
                if transport == "trucking" and pd.isnull(trucking_state) == False:
                    network.set_network(trucking_demand_schedule, times, country_series)

                    # Check for water constraint before any solving occurs
                    if water_limit != False:
                        water_constraint = get_water_constraint(network, trucking_demand_schedule, water_limit)
                        if water_constraint == False:
                            print('Not enough water to meet demand!')
                            trucking_lcs[i], generators_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], \
                            t_h2_storages[i] = np.nan
                            for gen, capacity in generators_capacities.items():
                                t_generators_capacities[gen].append(np.nan)
                            if plant_type == "ammonia": 
                                t_nh3_storages[i], t_hb_capacities[i] = np.nan
                            continue
                    
                    network.set_generators_in_network(country_series)
                    solve_model(network, solver)

                    if plant_type == "hydrogen":
                        trucking_lcs[i], generators_capacities, \
                        t_electrolyzer_capacities[i], t_battery_capacities[i], \
                        t_h2_storages[i] = get_h2_results(network.n, generators)
                    elif plant_type == "ammonia":
                        trucking_lcs[i], generators_capacities, \
                        t_electrolyzer_capacities[i], t_battery_capacities[i], \
                        t_h2_storages[i], t_nh3_storages[i], t_hb_capacities[i]  = get_nh3_results(network.n, generators)
                    
                    for gen, capacity in generators_capacities.items():
                        t_generators_capacities[gen].append(capacity)

                # If the hexagon has no viable trucking state (i.e., no roads reach it), set everything to nan.
                elif transport == "trucking" and pd.isnull(trucking_state) == True:
                    trucking_lcs[i], generators_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], \
                    t_h2_storages[i] = np.nan
                    for gen, capacity in generators_capacities.items():
                        t_generators_capacities[gen].append(np.nan)
                    if plant_type == "ammonia":
                        t_nh3_storages[i], t_hb_capacities[i] = np.nan

                # For pipeline, set it up with pipeline demand schedule if construction is true
                else:
                    if pipeline_construction == True:
                        network.set_network(pipeline_demand_schedule, times, country_series)

                        # Check for water constraint before any solving occurs
                        if water_limit != False:
                            water_constraint = get_water_constraint(network, pipeline_demand_schedule, water_limit)
                            if water_constraint == False:
                                print('Not enough water to meet demand!')
                                pipeline_lcs[i], generators_capacities, p_electrolyzer_capacities[i], p_battery_capacities[i], \
                                p_h2_storages[i] = np.nan
                                for gen, capacity in generators_capacities.items():
                                    p_generators_capacities[gen].append(np.nan)
                                if plant_type == "ammonia": 
                                    p_nh3_storages[i], p_hb_capacities[i] = np.nan
                                continue

                        network.set_generators_in_network(country_series)
                        solve_model(network, solver)

                        if plant_type == "hydrogen":
                            pipeline_lcs[i], generators_capacities, \
                            p_electrolyzer_capacities[i], p_battery_capacities[i], \
                            p_h2_storages[i] = get_h2_results(network.n, generators)
                        elif plant_type == "ammonia":
                            pipeline_lcs[i], generators_capacities, \
                            p_electrolyzer_capacities[i], p_battery_capacities[i], \
                            p_h2_storages[i], p_nh3_storages[i], p_hb_capacities[i] = get_nh3_results(network.n, generators)

                        for gen, capacity in generators_capacities.items():
                            p_generators_capacities[gen].append(capacity)

                    # If construction is false, you can't transport it
                    # Everything gets nan UNLESS in the demand centre hexagon
                    else:
                        # -- Demand location has trucking state as None, this is an easy check
                        if trucking_state == "None":
                            network.set_network(pipeline_demand_schedule, times, country_series)

                            # Check for water constraint before any solving occurs
                            if water_limit != False:
                                water_constraint = get_water_constraint(network, pipeline_demand_schedule, water_limit)
                                if water_constraint == False:
                                    print('Not enough water to meet demand!')
                                    pipeline_lcs[i], generators_capacities, p_electrolyzer_capacities[i], p_battery_capacities[i], \
                                    p_h2_storages[i] = np.nan
                                    for gen, capacity in generators_capacities.items():
                                        p_generators_capacities[gen].append(np.nan)
                                    if plant_type == "ammonia": 
                                        p_nh3_storages[i], p_hb_capacities[i] = np.nan
                                    continue

                            network.set_generators_in_network(country_series)
                            solve_model(network, solver)

                            if plant_type == "hydrogen":
                                pipeline_lcs[i], generators_capacities, \
                                p_electrolyzer_capacities[i], p_battery_capacities[i], \
                                p_h2_storages[i] = get_h2_results(network.n, generators)
                            elif plant_type == "ammonia":
                                pipeline_lcs[i], generators_capacities, \
                                p_electrolyzer_capacities[i], p_battery_capacities[i], \
                                p_h2_storages[i], p_nh3_storages[i], p_hb_capacities[i] = get_nh3_results(network.n, generators)

                            for gen, capacity in generators_capacities.items():
                                p_generators_capacities[gen].append(capacity)
                        else:
                            pipeline_lcs[i], p_electrolyzer_capacities[i], p_battery_capacities[i], p_h2_storages[i] = np.nan
                            for gen in p_generators_capacities:
                                p_generators_capacities[gen].append(np.nan)
                            if plant_type == "ammonia": 
                                p_nh3_storages[i], p_hb_capacities[i] = np.nan
        
        print("\nOptimisation complete.\n")        
        # Updating trucking-based results in hexagon file
        for gen, capacities in t_generators_capacities.items():
            hexagons[f'{demand_center} trucking {gen} capacity'] = capacities
        hexagons[f'{demand_center} trucking electrolyzer capacity'] = t_electrolyzer_capacities
        hexagons[f'{demand_center} trucking battery capacity'] = t_battery_capacities
        hexagons[f'{demand_center} trucking H2 storage capacity'] = t_h2_storages
        hexagons[f'{demand_center} trucking production cost'] = trucking_lcs
        if plant_type == "ammonia":        
            hexagons[f'{demand_center} trucking NH3 storage capacity'] = t_nh3_storages
            hexagons[f'{demand_center} trucking HB capacity'] = t_hb_capacities

        # Updating pipeline-based results in hexagon file
        for gen, capacities in p_generators_capacities.items():
            hexagons[f'{demand_center} pipeline {gen} capacity'] = capacities
        hexagons[f'{demand_center} pipeline electrolyzer capacity'] = p_electrolyzer_capacities
        hexagons[f'{demand_center} pipeline battery capacity'] = p_battery_capacities
        hexagons[f'{demand_center} pipeline H2 storage capacity'] = p_h2_storages
        hexagons[f'{demand_center} pipeline production cost'] = pipeline_lcs
        if plant_type == "ammonia":  
            hexagons[f'{demand_center} pipeline NH3 storage capacity'] = p_nh3_storages
            hexagons[f'{demand_center} pipeline HB capacity'] = p_hb_capacities

    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')

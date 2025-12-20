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
import pypsa
from scipy.constants import physical_constants
import xarray as xr
import glob

from functions import CRF
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

def calculate_grid_construction(hexagon, infra_data, country_params):
    grid_capex = infra_data.at['Grid','CAPEX']
    grid_opex = infra_data.at['Grid','OPEX']
    infrastructure_interest_rate = float(country_params['Infrastructure interest rate'].iloc[0])
    infrastructure_lifetime = float(country_params['Infrastructure lifetime (years)'].iloc[0]) 

    if hexagon['grid_dist']==0:
        grid_construction_costs = 0.
    else:
        grid_construction_costs = (hexagon['grid_dist'] 
                                   * grid_capex 
                                   * CRF(infrastructure_interest_rate,
                                         infrastructure_lifetime)
                                   + hexagon['grid_dist'] * grid_opex)
    return grid_construction_costs

def grid_connection(demand, elec_kWh_per_kg, product_quantity, country_series):
    
    grid_capacity = (elec_kWh_per_kg
                        * demand).max().iloc[0]/1000 # in MW if hourly data
    
    grid_energy_cost = ((elec_kWh_per_kg * demand * country_series[f"Electricity price ({currency}/kWh)"]).sum().iloc[0]
                        + (grid_capacity * country_series[f"Grid connection cost ({currency}/kW)"] * 1000)
                        + country_series[f"Electricity fixed charge ({currency}/year)"])
    
    grid_energy_cost_per_kg = grid_energy_cost / product_quantity
    lcoe = grid_energy_cost / (elec_kWh_per_kg * demand).sum().iloc[0]  # currency/kWh
    return grid_energy_cost, lcoe, grid_energy_cost_per_kg, grid_capacity

def estimate_offgrid_pop(df, elec_access, near_grid_uptake, avg_household):
    """
    The most basic estimate. Multiply the population per hexagon by the
    national electricity access rate.

    Parameters
    ----------
    population : TYPE
        DESCRIPTION.
    grid_dist : TYPE
        DESCRIPTION.
    elec_access : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    df = df.copy()
    # Calculate the off-grid population for distance_to_grid == 0 (40% of the population)
    df['offgrid_population'] = df['pop'].where(df['grid_dist'] == 0, 0) * (100-near_grid_uptake)/100

    # Calculate the total offgrid population based on the national average
    total_population = df['pop'].sum()
    total_offgrid_population = total_population * (100-elec_access)/100

    # Sum of offgrid population for areas with distance_to_grid == 0
    offgrid_population_assigned = df['offgrid_population'].sum()

    # Calculate the remaining offgrid population to be distributed
    remaining_offgrid_population = total_offgrid_population - offgrid_population_assigned

    # Get the total population for areas with distance_to_grid > 0
    remaining_population = df.loc[df['grid_dist'] > 0, 'pop'].sum()

    # Calculate the offgrid population for areas with distance_to_grid > 0
    df.loc[df['grid_dist'] > 0, 'offgrid_population'] = \
        df['pop'] * (remaining_offgrid_population / remaining_population)
    
    df["offgrid_households"] = df["offgrid_population"] / avg_household
    return df

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
    
    if plant_type == 'copper':
        needs_grid_construction = snakemake.config["grid_construction"]
        tech_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/technology_parameters.xlsx'
        infra_data = pd.read_excel(tech_params_filepath,
                                    sheet_name='Infra',
                                    index_col='Infrastructure')
        conversion_params_filepath = f'parameters/{snakemake.wildcards.country}/{snakemake.wildcards.plant_type}/conversion_parameters.xlsx'
        currency = snakemake.config["currency"]
        # Whether to include sizing of off-grid system for feedstock
        size_offgrid_feedstocks = snakemake.config['size offgrid feedstocks']
        # Either False or Float (indicating % of hexagon offgrid population to consider, or absolute households)
        community_energy_access = snakemake.config['community energy access']
        feedstock_output = False

        # estimate offgrid population
        if community_energy_access != False:
            hexagons = estimate_offgrid_pop(hexagons,
                                            country_params["Electricity Access (%)"],
                                            country_params["Near grid uptake (%)"],
                                            country_params["Household size"])
        
        # Calculate grid contruction costs
        if needs_grid_construction:
            hexagons_grid_construction = hexagons.apply(calculate_grid_construction,
                                                                    args=[infra_data, country_params],
                                                                    axis=1)
                
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
        if plant_type == "copper":
            # Offgrid storage variables
            offgrid_E_costs = np.zeros(len_hexagons)
            offgrid_lcoes = np.zeros(len_hexagons)
            offgrid_battery_capacities= np.zeros(len_hexagons)
            offgrid_generators_capacities = deepcopy(snakemake.config['generators_dict'])

            # Hybrid storage variables
            hybrid_E_costs = np.zeros(len_hexagons)
            hybrid_lcoes = np.zeros(len_hexagons)
            hybrid_battery_capacities= np.zeros(len_hexagons)
            hybrid_grid_capacities= np.zeros(len_hexagons)
            hybrid_generators_capacities = deepcopy(snakemake.config['generators_dict'])


            # Grid storage variables
            grid_construction_costs = np.zeros(len_hexagons)
            ongrid_totE_costs = np.zeros(len_hexagons)
            grid_lcoes = np.zeros(len_hexagons)
            grid_capacities = np.zeros(len_hexagons)
        else:
            lcs = np.zeros(len_hexagons)
            generators_capacities = deepcopy(snakemake.config['generators_dict'])
            electrolyzer_capacities= np.zeros(len_hexagons)
            battery_capacities = np.zeros(len_hexagons)
            h2_storages= np.zeros(len_hexagons)

        # Ammonia extra storage variables
        if plant_type == "ammonia":
            nh3_storages = np.zeros(len_hexagons)
            hb_capacities = np.zeros(len_hexagons)
        
        if plant_type =="copper":
            demand_state = demand_params.loc[demand_center,'Demand state']
            conversion_params = pd.read_excel(conversion_params_filepath,
                                    sheet_name=demand_state,
                                    index_col='Parameter').squeeze('columns')
            
            electricity_demand = conversion_params['Electricity demand (kWh per kg product)']

        annual_demand_quantity = float(demand_params.loc[demand_center,'Annual demand [kg/a]'])
        total_demand = annual_demand_quantity*years

        # Loop through all hexagons
        for i in range(len_hexagons):
            print(f"\nCurrently optimising {i+1} of {len_hexagons} hexagons...")
            gen_index = 0
            generators = deepcopy(snakemake.config['generators_dict'])
            
            # Creating demand schedule
            demand_schedule_index = pd.date_range(start_date, end_date, freq = "1H", inclusive='left')
            freq_demand_quantity = total_demand/demand_schedule_index.size
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

            if plant_type == 'copper':
                # Setting the energy access connections variable
                if (community_energy_access != False) and (community_energy_access <= 1.0):
                    # Use a percentage of households
                    energy_access_connections = hexagons.loc[i, "offgrid_households"] * community_energy_access
                elif (community_energy_access != False) and (community_energy_access > 1.0):
                    # Use community_energy_access as absolute number of households
                    energy_access_connections = community_energy_access
                else:
                    # Set as false - do not include
                    energy_access_connections = False
                    
                # Grid only energy optimisation
                if needs_grid_construction:
                    grid_construction_costs[i] = hexagons_grid_construction[i]/total_demand

                # Grid-only power
                (ongrid_E_cost,
                grid_lcoes[i],
                grid_energy_cost_per_kg,
                grid_capacities[i]) = grid_connection(demand_resampled_schedule,
                                                electricity_demand,
                                                total_demand,
                                                country_series)
                
                ongrid_totE_costs[i] = (ongrid_E_cost + hexagons_grid_construction[i])/total_demand
                

                # Set up the network
                network_class = Network(plant_type, generators, pypsa.Network())
                network_class.set_network(demand_resampled_schedule, country_series,
                                          demand_state,electricity_demand)
                network_class.update_generators(country_series)

                # Offgrid only optimisation of facility
                # Facility energy optimisation
                if demand_state == 'CuConcentrate':
                    offgrid_E_costs[i] = 0
                    offgrid_lcoes[i] = 0
                    offgrid_battery_capacities[i]= 0
                    for gen in generators:
                        offgrid_generators_capacities[gen].append(0)
                else:
                    if energy_access_connections != False:
                        network_class.add_community_energy_demand(energy_access_connections,
                                                                  f'data/community_elec_access_profile_{snakemake.wildcards.country}.csv')

                    network_class.n.lopf(solver_name=solver,
                                        solver_options = {'OutputFlag': 0},
                                        pyomo=False,
                                        )
                    
                    offgrid_E_costs[i] = network_class.n.objective
                    offgrid_lcoes[i] = offgrid_E_costs[i] / (network_class.n.generators_t.p.sum().sum()*1000)  # total system cost / energy produced (currency/kWh)
                    for gen in generators:
                        offgrid_generators_capacities[gen].append(network_class.n.generators.p_nom_opt[gen.capitalize()])
                    offgrid_battery_capacities[i] = network_class.n.storage_units.p_nom_opt['Battery']

                # Hybrid energy optimisation of facility
                network_class.add_grid(country_series, currency)

                network_class.n.lopf(solver_name=solver,
                                    solver_options = {'OutputFlag': 0},
                                    pyomo=False,
                                    )
                
                hybrid_E_costs[i] = network_class.n.objective + country_series[f"Electricity fixed charge ({currency}/year)"]
                # Hybrid LCOE = total system cost / energy produced (currency/kWh)
                hybrid_lcoes[i] = hybrid_E_costs[i] / (network_class.n.generators_t.p.sum().sum()*1000)
                for gen in generators:
                        hybrid_generators_capacities[gen].append(network_class.n.generators.p_nom_opt[gen.capitalize()])
                hybrid_battery_capacities[i] = network_class.n.storage_units.p_nom_opt['Battery']
                hybrid_grid_capacities[i] = network_class.n.generators.p_nom_opt['Grid']

                # Calculate offgrid system for feedstock centers if neccessary -- the results don't match geominx yet
                if (size_offgrid_feedstocks) and feedstock_output != True:
                    feedstock_points_gdf = gpd.read_file(f"resources/feedstocks_transport_{snakemake.wildcards.country}_{snakemake.wildcards.plant_type}.geojson")
                    
                    feedstocks_amount= len(feedstock_points_gdf)
                    # Off-grid costs
                    fs_lcoms = np.empty(feedstocks_amount)
                    fs_generators_capacities = deepcopy(snakemake.config['generators_dict'])
                    fs_battery_capacities = np.empty(feedstocks_amount)
                    fs_generators = deepcopy(snakemake.config['generators_dict'])
                    
                    feedstock_index = pd.date_range(start_date, end_date, freq=freq, inclusive='left')
                    
                    for c in feedstock_points_gdf.index:
                        hex_id = feedstock_points_gdf.loc[c, "nearest hexidx"]
                        hourly_demand = (feedstock_points_gdf.loc[c, "Electrical Energy (MWh/year)"]*1000)/(365*24)  # kwh
                        feedstock_demand_schedule = pd.DataFrame(hourly_demand,
                                                        index=feedstock_index, columns = ['Demand'])
                        fs_demand_resampled = feedstock_demand_schedule.resample(freq).mean()
                        
                        # Collecting information that is needed for Network set up.
                        gen_index=0
                        for gen in fs_generators:
                            potential = profiles[gen_index].sel(hexagon=c)

                            gen_capacity = snakemake.config['gen_capacity'][gen]
                            if gen == "wind":
                                max_capacity = hexagons.loc[c,'theo_turbines']*gen_capacity
                            elif gen == "solar":
                                max_capacity = hexagons.loc[c,'theo_pv']*gen_capacity
                            # -- Need to change solar and wind names in data-prep to have the below only
                            else:
                                max_capacity = hexagons.loc[c,gen]*gen_capacity

                            fs_generators[gen].append(potential)
                            fs_generators[gen].append(max_capacity)
                            gen_index += 1

                        # Set up the network
                        fs_network_class = Network(plant_type, fs_generators, pypsa.Network())
                        # -- The below elec kWh per kg is hardcoded...
                        fs_network_class.set_network(fs_demand_resampled, country_series,
                                                    elec_kWh_per_kg=1.0, feedstocks_opt=True )
                        fs_network_class.update_generators(country_series)

                        fs_network_class.n.lopf(solver_name=solver,
                                            solver_options = {'OutputFlag': 0},
                                            pyomo=False,
                                            )
                        
                        fs_lcoms[c] = fs_network_class.n.objective / (fs_network_class.n.generators_t.p.sum().sum()*1000)  # total system cost / energy produced (currency/kWh)
                        fs_battery_capacities[c] = fs_network_class.n.storage_units.p_nom_opt['Battery']
                        for gen in fs_generators:
                            fs_generators_capacities[gen].append(fs_network_class.n.generators.p_nom_opt[gen.capitalize()])
                        
                    # Store off-grid feedstock capacities and cost
                    for gen, capacities in fs_generators_capacities.items():
                        feedstock_points_gdf[f'Offgrid {gen} capacity (MW)'] = capacities
                    feedstock_points_gdf['Offgrid battery capacity (MW)'] = fs_battery_capacities
                    feedstock_points_gdf["Levlised cost of energy (euros/kg/year)"] = fs_lcoms
                    
                    # Output feedstocks information
                    fs_output_path = f"results/feedstocks_lc_{snakemake.wildcards.country}_{snakemake.wildcards.plant_type}.geojson"
                    feedstock_points_gdf.to_file(fs_output_path, driver='GeoJSON', encoding='utf-8')
                    feedstock_output = True
            else:
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
                network_class.update_generators(country_series)

                if plant_type == "hydrogen":
                    # Solve
                    network_class.n.lopf(solver_name=solver,
                        solver_options = {'LogToConsole':0, 'OutputFlag':0},
                        pyomo=False,
                        )
                    
                    h2_storages[i] = network_class.n.stores.e_nom_opt['Compressed H2 Store']
                    battery_capacities[i] = network_class.n.storage_units.p_nom_opt['Battery']
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
                    battery_capacities[i] = network_class.n.stores.e_nom_opt['Battery']

                lcs[i] = network_class.n.objective/total_demand
                electrolyzer_capacities[i] = network_class.n.links.p_nom_opt['Electrolysis']
                
                for gen in generators:
                    generators_capacities[gen].append(network_class.n.generators.p_nom_opt[gen.capitalize()])
    
        print("\nOptimisation complete.\n")        
        # Updating results in hexagon file
        if plant_type == 'copper':
            # off-grid costs
            for gen, capacities in offgrid_generators_capacities.items():
                hexagons[f'{demand_center} offgrid {gen} capacity (MW)'] = capacities
            hexagons[f'{demand_center} offgrid battery capacity (MW)'] = offgrid_battery_capacities
            hexagons[f'{demand_center} offgrid total energy cost ({currency}/kg/year)'] = offgrid_E_costs / annual_demand_quantity
            hexagons[f'{demand_center} offgrid lcoe ({currency}/kWh)'] = offgrid_lcoes
            
            # hybrid-grid costs
            for gen, capacities in hybrid_generators_capacities.items():
                hexagons[f'{demand_center} hybrid {gen} capacity (MW)'] = capacities
            hexagons[f'{demand_center} hybrid battery capacity (MW)'] = hybrid_battery_capacities
            hexagons[f'{demand_center} hybrid grid capacity (MW)'] = hybrid_grid_capacities
            hexagons[f'{demand_center} hybrid total energy cost ({currency}/kg/year)'] = (hybrid_E_costs + grid_construction_costs) / annual_demand_quantity
            hexagons[f'{demand_center} hybrid lcoe ({currency}/kWh)'] = hybrid_lcoes
            
            # grid costs
            hexagons[f'{demand_center} grid construction ({currency}/kg/year)'] = grid_construction_costs
            hexagons[f'{demand_center} grid total energy cost ({currency}/kg/year)'] = ongrid_totE_costs 
            hexagons[f'{demand_center} grid capacity (MW)'] = grid_capacities
            hexagons[f'{demand_center} grid lcoe ({currency}/kWh)'] = grid_lcoes

        else:
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
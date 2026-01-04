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
import matplotlib.pyplot as plt
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

# Helper: Add hydrogen demand flexibility buffer for load shifting (DSM)
def add_h2_demand_flexibility_storage(n, flex_hours, *,
                                     power_fraction=1.0,
                                     energy_fraction=1.0,
                                     cyclic=True,
                                     name_prefix="H2DemandFlex"):
    """Add a simple demand-flexibility buffer for hydrogen demand.

    This implements *load shifting* (not shedding): a Store on a dedicated DSM bus
    plus two Links (charge/discharge) connected to the hydrogen demand bus.

    Intuition:
      - Charging increases net demand in low-cost hours (pre-consumption).
      - Discharging reduces net demand in high-cost hours (postponed demand).
      - With `e_cyclic=True`, the buffer must end the optimisation horizon at the
        same state as it started (net-zero over the horizon).

    Parameters
    ----------
    n : pypsa.Network
        The network to modify (in-place).
    flex_hours : float
        Approximate hours of demand that can be shifted (energy buffer size is
        derived from the mean hydrogen demand).
    power_fraction : float
        Charging/discharging power limit as a fraction of peak hydrogen demand.
        (1.0 ~ allow shifting at up to peak demand rate)
    energy_fraction : float
        Energy buffer scaling as a fraction of mean_demand * flex_hours.
    cyclic : bool
        Whether the DSM store is cyclic (recommended for pure shifting).
    name_prefix : str
        Prefix for component names.
    """
    try:
        flex_hours = float(flex_hours)
    except Exception:
        return

    if flex_hours <= 0:
        return

    # Need at least one load time series to infer sizing.
    if getattr(n, "loads", None) is None or n.loads.empty:
        return
    if not hasattr(n, "loads_t") or not hasattr(n.loads_t, "p_set") or n.loads_t.p_set.empty:
        return

    # Identify the main hydrogen demand load and its bus.
    demand_load_name = "Hydrogen demand" if "Hydrogen demand" in n.loads.index else n.loads.index[0]
    demand_bus = n.loads.loc[demand_load_name, "bus"]

    # Size relative to the actual load time series (units are whatever the H2 bus uses, typically MW_H2).
    total_load = n.loads_t.p_set.sum(axis=1)
    if total_load.empty:
        return

    mean_load = float(total_load.mean())
    peak_load = float(total_load.max())

    # Derive buffer size (energy) and power limits (rate) from demand.
    e_nom = max(1e-6, mean_load * flex_hours * float(energy_fraction))
    p_nom = max(1e-6, peak_load * float(power_fraction))

    dsm_bus = f"{demand_bus}__dsm"
    if dsm_bus not in n.buses.index:
        # Keep carrier consistent where possible
        carrier = ""
        try:
            if "carrier" in n.buses.columns and demand_bus in n.buses.index:
                carrier = n.buses.loc[demand_bus, "carrier"]
        except Exception:
            carrier = ""
        n.add("Bus", dsm_bus, carrier=carrier)

    store_name = f"{name_prefix} Store"
    link_charge = f"{name_prefix} Charge"
    link_discharge = f"{name_prefix} Discharge"

    # Avoid accidental duplication if called twice
    if store_name in getattr(n, "stores", pd.DataFrame()).index:
        return

    n.add(
        "Store",
        store_name,
        bus=dsm_bus,
        e_nom=e_nom,
        e_initial=0.0,
        e_cyclic=bool(cyclic),
        marginal_cost=0.0,
        capital_cost=0.0,
    )

    # Power limits for charging/discharging (Store itself has no power cap in PyPSA).
    n.add(
        "Link",
        link_charge,
        bus0=demand_bus,
        bus1=dsm_bus,
        p_nom=p_nom,
        efficiency=1.0,
    )
    n.add(
        "Link",
        link_discharge,
        bus0=dsm_bus,
        bus1=demand_bus,
        p_nom=p_nom,
        efficiency=1.0,
    )

    logging.info(
        f"Added H2 demand flexibility buffer: e_nom={e_nom:.3f}, p_nom={p_nom:.3f}, flex_hours={flex_hours}"
    )


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

def plot_dispatch(n, time, demand, hex, tran):
    print("Plotting dispatch...")
    fig, ax = plt.subplots(figsize=(6, 3))
    
    p_by_carrier = n.generators_t.p.loc[time]
    # print("Generators shape:", p_by_carrier.shape)
    # # Debug: inspect index info
    # idx = p_by_carrier.index
    # print("Index dtype:", idx.dtype)
    # print("Index length:", len(idx))
    # if len(idx) > 0:
    #     print("Index range:", idx.min(), "->", idx.max())
    #     print("First 5 index values:", idx[:5])
    
    # print("Time type:", type(time))

    if not n.storage_units.empty:
    
        sto = (
            n.storage_units_t.p.loc[time]
        )
        
        p_by_carrier = pd.concat([p_by_carrier, sto], axis=1)
    if not n.stores.empty:
        sto = (
            n.stores_t.p.loc[time]
            )
        p_by_carrier = pd.concat([p_by_carrier, sto], axis=1)
    
    p_by_carrier.where(p_by_carrier > 0).loc[time].plot.area(
        ax=ax,
        linewidth=0,
    )
    print(n.generators_t.p.shape) 
    print(n.generators_t.p.columns)
    print(type(time), repr(time))
    print(p_by_carrier.shape, p_by_carrier.head())
    
    charge = p_by_carrier.where(p_by_carrier < 0).dropna(how="all", axis=1).loc[time]
    
    if not charge.empty:
        charge.plot.area(
            ax=ax,
            linewidth=0,
        )
    
    n.loads_t.p_set.sum(axis=1).loc[time].plot(ax=ax, c="k")
    
    plt.legend(loc=(1.05, 0.2))
    ax.set_ylabel("MW")
    ax.set_ylim(-1500, 2000)
    fig.savefig(f'resources\{demand}_{hex}_temporal_{tran}_June.png', bbox_inches='tight')
    plt.close()

def plot_dispatch_2(n, time, demand, hex, tran):
    """
    Produces one plot per carrier. Each plot shows stacked positive (generation/discharge)
    and stacked negative (consumption/charging) time series for the original components
    (generators, storage_units, stores, and split links). Legend entries are the original
    series names (component names) rather than only the carrier name.

    Use: replace calls to plot_dispatch_2(...) with plot_dispatch_by_series(...)
    """
    p = n.generators_t.p.copy() if hasattr(n, "generators_t") else pd.DataFrame()
    if not getattr(n, "storage_units", pd.DataFrame()).empty and hasattr(n.storage_units_t, "p"):
        p = pd.concat([p, n.storage_units_t.p], axis=1)
    if not getattr(n, "stores", pd.DataFrame()).empty and hasattr(n.stores_t, "p"):
        p = pd.concat([p, n.stores_t.p], axis=1)

    # Include link flows split into positive on receiving carrier and negative on sending carrier
    print("hasattr links_t", hasattr(n, "links_t"))
    #print("attribute links_t:", getattr(n, "links_t", None)). 
    print("hasattr links_t.p0", hasattr(n.links_t, "p0"))
    print("attribute links_t.p0:", getattr(n.links_t, "p0", None))
    if hasattr(n, "links_t") and hasattr(n.links_t, "p0") and not n.links_t.p0.empty:
        links_p = n.links_t.p0.copy()
        print("Links header", links_p.head())
        link_contrib = pd.DataFrame(index=links_p.index)
        for ln in links_p.columns:
            s = links_p[ln].fillna(0)
            bus_from = None
            bus_to = None
            try:
                if ln in n.links.index and "bus0" in n.links.columns:
                    bus_from = n.links.loc[ln, "bus0"]
            except Exception:
                bus_from = None
            try:
                if ln in n.links.index and "bus1" in n.links.columns:
                    bus_to = n.links.loc[ln, "bus1"]
            except Exception:
                bus_to = None

            carrier_from = "unknown"
            carrier_to = "unknown"
            try:
                if bus_from in getattr(n, "buses", pd.DataFrame()).index:
                    carrier_from = n.buses.loc[bus_from, "carrier"]
            except Exception:
                carrier_from = "unknown"
            try:
                if bus_to in getattr(n, "buses", pd.DataFrame()).index:
                    carrier_to = n.buses.loc[bus_to, "carrier"]
            except Exception:
                carrier_to = carrier_from if carrier_from != "unknown" else "unknown"

            link_contrib[f"{carrier_from}__{ln}"] = s      # positive on receiving carrier
            # link_contrib[f"{carrier_from}__{ln}"] = -s   # negative on sending carrier

        p = pd.concat([p, link_contrib], axis=1)

    # ensure datetime index
    if not isinstance(p.index, pd.DatetimeIndex):
        try:
            p.index = pd.to_datetime(p.index)
        except Exception:
            print("Could not convert generators index to DatetimeIndex; aborting carrier-by-series plots.")
            return

    # prepare carrier map: column -> carrier
    carrier_map = {}
    for col in p.columns:
        if col in getattr(n, "generators", pd.DataFrame()).index:
            carrier_map[col] = n.generators.loc[col, 'carrier']
        elif hasattr(n, 'storage_units') and col in getattr(n, "storage_units", pd.DataFrame()).index:
            carrier_map[col] = n.storage_units.loc[col, 'carrier']
        elif hasattr(n, 'stores') and col in getattr(n, "stores", pd.DataFrame()).index:
            carrier_map[col] = n.stores.loc[col, 'carrier']
        else:
            # columns created earlier for links use a "{carrier}__{ln}" prefix
            if "__" in col:
                carrier_map[col] = col.split("__", 1)[0]
            else:
                carrier_map[col] = 'other'

    pd.DataFrame.from_dict(carrier_map, orient='index', columns=['carrier']).to_csv(f"resources/{demand}_{hex}_{tran}.csv") 

    # total load series for the black line
    total_load = n.loads_t.p_set.sum(axis=1)

    # time slicing
    if isinstance(time, slice):
        start = pd.to_datetime(time.start)
        stop = pd.to_datetime(time.stop)
        mask = (p.index >= start) & (p.index < stop)
        if not mask.any():
            print(f"No timestamps in requested slice {start} -> {stop}; skipping carrier plots.")
            return

    # iterate carriers
    for carrier in sorted(set(carrier_map.values())):
        cols = [c for c, ca in carrier_map.items() if ca == carrier]
        if not cols:
            continue
        carrier_df = p[cols]

        # slice if requested
        try:
            carrier_slice = carrier_df.loc[time]
            load_slice = total_load.loc[time] if not total_load.empty else pd.Series(dtype=float)
        except Exception:
            # fallback boolean mask for slice compatibility
            if isinstance(time, slice):
                carrier_slice = carrier_df[(carrier_df.index >= start) & (carrier_df.index < stop)]
                load_slice = total_load[(total_load.index >= start) & (total_load.index < stop)]
            else:
                carrier_slice = carrier_df
                load_slice = total_load

        if carrier_slice.empty:
            print(f"Carrier '{carrier}' has no data in slice; skipping.")
            continue

        pos = carrier_slice.where(carrier_slice > 0).fillna(0)
        neg = carrier_slice.where(carrier_slice < 0).fillna(0)

        fig, ax = plt.subplots(figsize=(8, 3.5))
        # stacked positive components (each series keeps its own name in the legend)
        if not pos.empty and pos.sum().sum() != 0:
            pos.plot.area(ax=ax, linewidth=0)

        # stacked negative components (these will appear below zero)
        if not neg.empty and neg.sum().sum() != 0:
            neg.plot.area(ax=ax, linewidth=0)

        # plot total load line
        if not load_slice.empty:
            load_slice.plot(ax=ax, color="k", linewidth=1, label="Total load")

        ax.set_title(f"Carrier: {carrier} — {demand} hex {hex} ({tran})")
        ax.set_ylabel("MW")
        plt.legend(loc=(1.01, 0.2))
        ax.set_ylim(-1500, 2000)
        plt.tight_layout()
        fname = f"resources/{demand}_{hex}_carrier_{carrier}_{tran}.png"
        fig.savefig(fname, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved carrier plot: {fname}")

def save_operational_csv(n, time, demand, hex, tran):
    """
    Save operational time-series from network `n` for the given `time` window to a single CSV.

    Signature matches plot_dispatch_2: (n, time, demand, hex, tran).
    Produces: resources/{demand}_{hex}_ops_{tran}.csv
    """
    #import os
    #os.makedirs("resources", exist_ok=True)

    # collect time-series containers (prefer common attribute names)
    ts_frames = []

    if hasattr(n, "generators_t") and hasattr(n.generators_t, "p"):
        ts_frames.append(n.generators_t.p.copy())

    if hasattr(n, "storage_units_t") and hasattr(n.storage_units_t, "p"):
        ts_frames.append(n.storage_units_t.p.copy())

    if hasattr(n, "stores_t") and hasattr(n.stores_t, "p"):
        ts_frames.append(n.stores_t.p.copy())

    # links time-series: prefer .p then .p0 (PyPSA versions differ)
    links_p = None
    links_p1 = None
    if hasattr(n, "links_t"):
        if hasattr(n.links_t, "p"):
            links_p = getattr(n.links_t, "p")
        elif hasattr(n.links_t, "p0"):
            links_p = getattr(n.links_t, "p0")
            links_p1 = getattr(n.links_t, "p1")
    if links_p is not None:
        try:
            ts_frames.append(links_p.copy())
            ts_frames.append(links_p1.copy())
        except Exception:
            pass

    # merge all series (align on index)
    if ts_frames:
        ops = pd.concat(ts_frames, axis=1)
    else:
        ops = pd.DataFrame()

    # add loads / total load if available
    total_load = None
    if hasattr(n, "loads_t") and hasattr(n.loads_t, "p_set"):
        try:
            total_load = n.loads_t.p_set.sum(axis=1).rename("total_load")
            # join will align indices
            ops = ops.join(total_load, how="outer")
        except Exception:
            pass

    # ensure datetime index where possible
    if not isinstance(ops.index, pd.DatetimeIndex):
        try:
            ops.index = pd.to_datetime(ops.index)
        except Exception:
            # leave as-is if conversion fails
            pass

    # slice selection: support slice(...) or label-based selection
    if isinstance(time, slice):
        start = pd.to_datetime(time.start)
        stop = pd.to_datetime(time.stop)
        ops_slice = ops[(ops.index >= start) & (ops.index < stop)]
    else:
        try:
            ops_slice = ops.loc[time]
        except Exception:
            # fallback -> return full ops
            ops_slice = ops

    if ops_slice.empty:
        print(f"No operational data to save for slice {time}.")
        return

    fname = f"resources/{demand}_{hex}_ops_{tran}.csv"
    # write CSV (preserve index as first column)
    ops_slice.to_csv(fname, index=True)
    print(f"Saved operational CSV: {fname}")

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
    profiles = {}
    len_hexagons = len(hexagons)
    has_water_limit = snakemake.config['water']['has_limit']
    water_limit_amount = snakemake.config['water']['annual_limit']
    freq = snakemake.config['freq']
    plant_type = snakemake.wildcards.plant_type

    generic_generator_cfg = snakemake.config.get("generic_generator", {})
    if "existing" in generators:
        if "marginal_cost" not in generic_generator_cfg:
            raise KeyError(
                "config.yaml requires generic_generator.marginal_cost when 'existing' is enabled in generators_dict."
            )
        existing_marginal_cost = float(generic_generator_cfg["marginal_cost"])

        missing_cols = [c for c in ["hydro", "nuclear"] if c not in hexagons.columns]
        if missing_cols:
            raise KeyError(
                f"Hexagon file must contain {missing_cols} for existing generator p_nom = hydro + nuclear."
            )

        existing_capacity_by_hex = (
            pd.to_numeric(hexagons["hydro"], errors="coerce").fillna(0.0)
            + pd.to_numeric(hexagons["nuclear"], errors="coerce").fillna(0.0)
        )
    else:
        existing_marginal_cost = None
        existing_capacity_by_hex = None

    
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
        profiles["solar"] = solar_profile.rename(dict(dim_0='hexagon'))
        pd.DataFrame(profiles["solar"]).to_csv(f"resources/{snakemake.wildcards.country}_solar_profile_{snakemake.wildcards.weather_year}.csv")
    # Creating wind generator profile
    if "wind" in generators:
        wind_profile = cutout.wind(turbine = snakemake.config["turbine"],
                                    layout = layout,
                                    shapes = hexagons,
                                    per_unit = True).resample(time=freq).mean()
        profiles["wind"] = wind_profile.rename(dict(dim_0='hexagon'))
        pd.DataFrame(profiles["wind"]).to_csv(f"resources/{snakemake.wildcards.country}_wind_profile_{snakemake.wildcards.weather_year}.csv")
    
    if "existing" in generators:
        existing_time_index = pd.date_range(start_date, end_date, freq=freq, inclusive='left')
        profiles["existing"] = xr.DataArray(
            data=np.ones((len_hexagons, len(existing_time_index))),
            dims=["hexagon", "time"],
            coords={"hexagon": np.arange(len_hexagons), "time": existing_time_index},
        )

    for gen in generators:
        if gen not in profiles:
            raise ValueError(f"No profile was prepared for generator '{gen}'.")

    # Adding generators into profiles list
    # for gen in generators:
    #     # profile = f'{gen}_profile'
    #     if gen == 'solar':
    #         profiles.append(solar_profile)
    #     elif gen == 'wind':
    #         profiles.append(wind_profile)
    #     elif gen == 'hydro':
    #         profiles.append(hydro_profile)
    #     elif gen == 'geothermal':
    #         profiles.append(geothermal_profile)


    # Loop through all demand centers
    for demand_center in demand_centers:
        print(f"\nOptimisation for {demand_center} begins...")
        # Store results
        lcs = np.zeros(len_hexagons)
        generators_capacities = deepcopy(snakemake.config['generators_dict'])
        electrolyzer_capacities= np.zeros(len_hexagons)
        battery_capacities = np.zeros(len_hexagons)
        h2_storages= np.zeros(len_hexagons)
        existing_generation = np.zeros(len_hexagons) if "existing" in snakemake.config['generators_dict'] else None


        # Ammonia extra storage variables
        if plant_type == "ammonia":
            nh3_storages = np.zeros(len_hexagons)
            hb_capacities = np.zeros(len_hexagons)

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
                potential = profiles[gen].sel(hexagon=i)

                if gen == "wind":
                    gen_capacity = snakemake.config['gen_capacity'][gen]
                    max_capacity = hexagons.loc[i,'theo_turbines']*gen_capacity
                elif gen == "solar":
                    gen_capacity = snakemake.config['gen_capacity'][gen]
                    max_capacity = hexagons.loc[i,'theo_pv']*gen_capacity
                elif gen == "existing":
                    max_capacity = float(existing_capacity_by_hex.iloc[i])
                else:
                    gen_capacity = snakemake.config['gen_capacity'][gen]
                    max_capacity = hexagons.loc[i,gen]*gen_capacity

                generators[gen].append(potential)
                generators[gen].append(max_capacity)
                if gen == "existing":
                    generators[gen].append(existing_marginal_cost)

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
                    if existing_generation is not None:
                        existing_generation[i] = np.nan
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
                dsm_enabled = bool(snakemake.config.get("h2_demand_flex", False))
                dsm_hours = float(snakemake.config.get("h2_demand_flex_hours", 0) or 0)
                if dsm_enabled and dsm_hours > 0:
                    add_h2_demand_flexibility_storage(
                        network_class.n,
                        dsm_hours,
                        power_fraction=float(snakemake.config.get("h2_demand_flex_power_fraction", 1.0) or 1.0),
                        energy_fraction=float(snakemake.config.get("h2_demand_flex_energy_fraction", 1.0) or 1.0),
                        cyclic=True,
                        name_prefix="H2DemandFlex",
                        )
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
            if existing_generation is not None:
                if "Existing" in network_class.n.generators_t.p.columns:
                    snapshot_weight = 8760 / len(network_class.n.snapshots)
                    existing_generation[i] = (
                        network_class.n.generators_t.p["Existing"].sum() * snapshot_weight
                    )
                else:
                    existing_generation[i] = 0.0
        
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
        if existing_generation is not None:
            hexagons[f'{demand_center} existing generation'] = existing_generation


    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')
"""
@authors:
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk
Includes code from Nicholas Salmon, University of Oxford, for setting up the 
network.

Class representing a PyPSA Network. Contains functions used for set up.
"""
import numpy as np
import pandas as pd
import pypsa

from functions import CRF

class Network:
    """
    A class representing a PyPSA Network.

    ...
    Attributes
    ----------
    type : string
        type of commodity plant.
    generators : dictionary
        keys - generator types, values - list of respective potential and 
        maximum capacity.
    n : pypsa Network Object
        network. Default is None.

    Methods
    -------
    set_network(demand_profile, times, country_series):
        sets up the network.
    set_generators_in_network(country_series):
        sets provided generator in the network.
    _create_override_components():
        set up new component attributes as required.
    """
    def __init__(self, type, generators, n=None):
        """
        Provide a network object.
        """
        self.type = type
        self.generators = generators
        self.n = n
    
    def set_network(self, demand_profile, country_series, demand_state=None, elec_kWh_per_kg=None):
        '''
        Sets up the network.
        
        ...
        Parameters
        ----------
        demand_profile : pandas DataFrame
            dataframe of commodity demand in kg in frequency configured.
        country_series : pandas Series
            interest rate and lifetime information.
        '''
        # Set standard Network, if none provided
        if self.n == None:
            self.n = pypsa.Network(override_component_attrs=self._create_override_components())
        
        # Set the time values for the network
        self.n.set_snapshots(demand_profile.index)
        demand_profile['weights'] = 8760 / len(self.n.snapshots)
        self.n.snapshot_weightings = demand_profile['weights']

        if self.type == "hydrogen":
            # Import the design of the H2 plant into the network
            self.n.import_from_csv_folder("parameters/basic_h2_plant")

            # Import demand profile
            # Note: All flows are in MW or MWh, conversions for hydrogen done 
            # using HHVs. Hydrogen HHV = 39.4 MWh/t
            self.n.add('Load',
                'Hydrogen demand',
                bus='Hydrogen',
                p_set=demand_profile['Demand']/1000*39.4,
                )

            for item in [self.n.links, self.n.stores, self.n.storage_units]:
                item.capital_cost = item.capital_cost * CRF(country_series['Plant interest rate'], 
                                                            country_series['Plant lifetime (years)'])
        elif self.type == "ammonia":
            # Import the design of the NH3 plant into the network
            self.n.import_from_csv_folder("parameters/basic_nh3_plant")

            # Import demand profile
            # Note: All flows are in MW or MWh, conversions for ammonia done 
            # using HHVs. Ammonia HHV = 6.25 MWh/t
            self.n.add('Load',
                'Ammonia demand',
                bus='Ammonia',
                p_set=demand_profile['Demand']/1000*6.25,
                )
            
            for item in [self.n.links, self.n.stores]:
                item.capital_cost = item.capital_cost \
                                    * CRF(country_series['Plant interest rate'],
                                        country_series['Plant lifetime (years)'])
            # Stops pointless cycling through storage
            self.n.links.loc['HydrogenCompression', 'marginal_cost'] = 0.0001
        elif self.type == "copper":
            # Import the design of the Cu plant into the network
            self.n.import_from_csv_folder("parameters/basic_cu_plant")

            # Import demand profile
            # Note: All flows are in MW or MWh, conversions for Concentrate done using 0.717 kWh per kg
            # Add the load
            self.n.add('Load',
                f'{demand_state} demand',
                bus = 'AC_bus',
                p_set = (demand_profile["Demand"] * elec_kWh_per_kg) / 1000,
                )
    
            # Update and set capital costs
            self.n.storage_units.loc["Battery", "capital_cost"] *= CRF(country_series['Battery interest rate'],
                                                                  country_series['Battery lifetime (years)'])
            self.n.links.loc["Inverter", "capital_cost"] *= CRF(country_series['Plant interest rate'],
                                                                country_series['Plant lifetime (years)'])
            self.n.links.loc["Rectifier", "capital_cost"] *= CRF(country_series['Plant interest rate'],
                                                                country_series['Plant lifetime (years)'])
                                                                   
    def add_community_energy_demand(self, energy_access_connections, filepath):
        '''
        Adds community demand as a load.
        
        ...
        Parameters
        ----------
        energy_access_connections : float
                number of energy access connections required.
        filepath : string
            pathway to the community electricity access profile.
        '''
        # Community energy demand
        # Basic model assumes connected to same single AC bus as facility
        energy_access_df = pd.read_csv(filepath, index_col=0, parse_dates=True)
        energy_access_demand = energy_access_df * energy_access_connections
        self.n.add('Load',
                    'Community Demand',
                    bus = 'AC_bus',
                    p_set = energy_access_demand['MW'])
    
    def add_grid(self, country_series, currency):
        '''
        Adds grid as a generator.

        ...
        Parameters
        ----------
        country_series : pandas Series
            electricity price and grid connection cost information.
        currency : string
        unit of currency that is used in the parameter files.

        '''
        self.n.add("Generator",
                   "Grid",
                   bus="AC_bus",
                   p_nom_extendable=True,
                   marginal_cost=country_series[f"Electricity price ({currency}/kWh)"]*1000,
                   capital_cost=country_series[f"Grid connection cost ({currency}/kW)"]*1000)

    def update_generators(self, country_series):
        '''
        Updates provided generators in the network.

        ...
        Parameters
        ----------
        country_series : pandas Series
            interest rate and lifetime information.
        '''
        # Send the generator data to the network
        for gen, gen_list in self.generators.items():
            if gen == "geothermal":
                self.n.generators.loc[gen.capitalize(), "p_max_pu"] = gen_list[0]
            else:
                self.n.generators_t.p_max_pu[gen.capitalize()] = gen_list[0]

            # Specify maximum capacity based on land use
            self.n.generators.loc[gen.capitalize(),'p_nom_max'] = gen_list[1]

            # Specify technology-specific and country-specific WACC and lifetime here
            self.n.generators.loc[gen.capitalize(),'capital_cost'] *= \
                CRF(country_series[f'{gen.capitalize()} interest rate'], 
                    country_series[f'{gen.capitalize()} lifetime (years)'])

    def _create_override_components(self):
        """Set up new component attributes as required"""
        # Modify the capacity of a link so that it can attach to 2 buses.
        override_component_attrs = pypsa.descriptors.Dict(
            {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
        )

        override_component_attrs["Link"].loc["bus2"] = [
            "string",
            np.nan,
            np.nan,
            "2nd bus",
            "Input (optional)",
        ]
        override_component_attrs["Link"].loc["efficiency2"] = [
            "static or series",
            "per unit",
            1.0,
            "2nd bus efficiency",
            "Input (optional)",
        ]
        override_component_attrs["Link"].loc["p2"] = [
            "series",
            "MW",
            0.0,
            "2nd bus output",
            "Output",
        ]
        
        return override_component_attrs
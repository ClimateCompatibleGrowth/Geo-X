"""
@authors:
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk
Includes code from Nicholas Salmon, University of Oxford, for setting up the 
network.

Class representing a PyPSA Network. Contains functions used for set up.
"""
import numpy as np
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
        Provide an network object.
        """
        self.type = type
        self.generators = generators
        self.n = n

    def set_network(self, demand_profile, country_series):
        '''
        Sets up the network.
        
        ...
        Parameters
        ----------
        demand_profile : pandas DataFrame
            dataframe of commodity demand in kg in frequency configured.
        country_series: pandas Series
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

    def set_generators_in_network(self, country_series):
        '''
        Sets provided generator in the network.

        ...
        Parameters
        ----------
        country_series: pandas Series
            interest rate and lifetime information.
        '''
        for gen, gen_list in self.generators.items():
            gen_name = gen.capitalize()

            # Allow generators created programmatically (not only from CSV templates)
            if gen_name not in self.n.generators.index:
                self.n.add(
                    "Generator",
                    gen_name,
                    bus="Power",
                    carrier=gen,
                    p_nom=0.0,
                    p_nom_extendable=True,
                    capital_cost=0.0,
                    marginal_cost=0.0,
                )

            if gen == "geothermal":
                self.n.generators.loc[gen_name, "p_max_pu"] = gen_list[0]
            else:
                self.n.generators_t.p_max_pu[gen_name] = gen_list[0]

            if gen == "existing":
                if len(gen_list) < 3:
                    raise ValueError(
                        "Generator 'existing' expects [profile, p_nom, marginal_cost]."
                    )
                fixed_capacity = float(gen_list[1])
                marginal_cost = float(gen_list[2])

                self.n.generators.loc[gen_name, "p_nom_extendable"] = False
                self.n.generators.loc[gen_name, "p_nom"] = fixed_capacity
                self.n.generators.loc[gen_name, "p_nom_max"] = fixed_capacity
                self.n.generators.loc[gen_name, "capital_cost"] = 0.0
                self.n.generators.loc[gen_name, "marginal_cost"] = marginal_cost
                continue

            self.n.generators.loc[gen_name, 'p_nom_max'] = gen_list[1]
            self.n.generators.loc[gen_name, 'capital_cost'] *= \
                CRF(country_series[f'{gen_name} interest rate'],
                    country_series[f'{gen_name} lifetime (years)'])


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
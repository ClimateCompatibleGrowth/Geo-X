Introduction
============

Geo-X is a geospatial modelling framework for analysing the cost of commodity production and supply.

It combines geospatial data, `ERA5 <https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview>`_ weather data (via the `Climate Data Store API <https://cds.climate.copernicus.eu/how-to-api>`_), and energy system optimisation using `PyPSA <https://pypsa.org/>`_ to estimate the levelised cost of supplying commodities at specified demand locations.

What Geo-X does
---------------

Geo-X calculates the cost of commodity production across many potential locations and identifies least-cost supply options.

It includes:

- Optimisation of production plant design
- Modelling of renewable energy generation
- Assessment of water availability and supply costs
- Estimation of storage and transport costs
- Conversion between commodity states (where relevant)
- Spatially resolved cost estimates

Results can be compared with current or projected market prices to assess economic competitiveness.

Supported commodities
---------------------

Geo-X currently supports:

- **Green hydrogen**

  Produced using off-grid renewable energy systems

- **Green ammonia**

  Produced from hydrogen using off-grid renewable energy systems

- **Green copper**

  Modelled under multiple electricity supply options:

  - Off-grid renewable systems
  - Hybrid systems (grid + on-site generation)
  - Fully grid-connected systems

Hydrogen and ammonia are restricted to off-grid renewable production, while copper can incorporate grid electricity depending on the scenario.

Electricity generation
----------------------

Electricity supply is configurable within the model.

Supported generation technologies include:

- Solar
- Wind
- Hydropower
- Geothermal
- Nuclear

Users can select which technologies to include depending on the scenario.

How it works
------------

Geo-X evaluates a set of spatial locations (represented as `H3 hexagons <https://h3geo.org/>`_) across a region.

For each location, it:

- Uses ERA5 weather data (via the Climate Data Store) to estimate renewable generation
- Assesses water availability and supply costs for production
- Designs an optimal production system using PyPSA
- Calculates transport costs to demand centres
- Combines all components to estimate total supply cost

The workflow is managed using Snakemake, which resolves dependencies and runs the required steps automatically.

Input data
----------

Geo-X requires three main types of input:

- Spatial data (hexagons and infrastructure distances)
- Techno-economic parameters (costs, efficiencies, demand)
- A :doc:`configuration file </configuration>` defining the scenario

Additional data are required when including hydropower, geothermal, or nuclear.

Current limitations
-------------------

- Assumes greenfield development (no existing infrastructure) for solar and wind; other technologies rely on user-provided data
- Limited representation of grid electricity (primarily in copper scenarios)
- Transport is limited to land-based options
- Terrain and routing complexity are not explicitly modelled
- Results depend on input data resolution and assumptions
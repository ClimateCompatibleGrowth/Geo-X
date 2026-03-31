Configuration
=============

Geo-X uses a ``config.yaml`` file to define the scenario and control how the workflow is run.

This file specifies the country, weather data, technologies, and key model settings used in the analysis.

Scenario settings
-----------------

The ``scenario`` section defines the core setup for a run:

.. code-block:: yaml

   scenario:
     country: ['NA']
     weather_year: [2023]
     plant_type: 'ammonia'

- ``country``: `ISO 2-letter country code(s) <https://www.iban.com/country-codes>`_
- ``weather_year``: Start year for weather data (ERA5)
- ``plant_type``: Commodity to model (``hydrogen``, ``ammonia``, or ``copper``)

Multiple countries or years can be provided as lists. Required input data must exist for each combination.

Weather data
------------

.. code-block:: yaml

   years_to_check: 1
   freq: '1H'

- ``years_to_check``: Number of years of weather data to include
- ``freq``: Temporal resolution of the data (e.g. ``1H`` for hourly, ``3H`` for three-hourly)

Generators
----------

.. code-block:: yaml

   generators_dict: {'solar': [], 'wind': [], 'hydro': [], 'geothermal': [], 'nuclear': []}
   panel: 'CSi'
   turbine: 'NREL_ReferenceTurbine_2020ATB_4MW'

The ``generators_dict`` defines which generation technologies are included in the optimisation.

Supported options are:

- ``solar``
- ``wind``
- ``hydro``
- ``geothermal``
- ``nuclear``

Solar panel and wind turbine types can be specified using:

- ``panel``: solar technology
- ``turbine``: wind turbine model

Installed capacity assumptions are defined as:

.. code-block:: yaml

   gen_capacity:
       solar: 1
       wind: 4
       hydro: 1
       geothermal: 1
       nuclear: 1

- Values for solar and wind can be adjusted
- Values for hydro, geothermal, and nuclear should not be changed

Generator efficiency
--------------------

.. code-block:: yaml

   efficiency:
       hydro: 0.75
       geothermal: 0.9
       nuclear: 0.33

Defines average efficiencies for generator types where relevant.

Solver
------

.. code-block:: yaml

   solver: 'gurobi'

Specifies the optimisation solver to use. This must match an installed solver compatible with PyPSA.

Water constraints
-----------------

.. code-block:: yaml

   water:
       has_ocean_dists: True
       has_limit: False
       annual_limit:

- ``has_ocean_dists``: Set to ``True`` for coastal countries to include ocean distance calculations
- ``has_limit``: Enable water availability constraints (for hydrogen and ammonia)
- ``annual_limit``: Annual water availability (m³), required if ``has_limit`` is ``True``

Transport options
-----------------

.. code-block:: yaml

   transport:
       road_construction: True
       pipeline_construction: True

- ``road_construction``: Enable or disable road infrastructure
- ``pipeline_construction``: Enable pipeline infrastructure (hydrogen and ammonia only)

General settings
----------------

.. code-block:: yaml

   currency: 'euros'
   restrict_hexagons: False

- ``currency``: Must match the currency used in parameter files
- ``restrict_hexagons``: If ``True``, results are only produced for hexagons containing demand centres (all others are set to ``null``)

Copper-specific settings
------------------------

These options are only used when ``plant_type`` is ``copper``:

.. code-block:: yaml

   grid_construction: True
   size_offgrid_feedstocks: False
   community_energy_access: False

- ``grid_construction``: Allow grid infrastructure to be built
- ``size_offgrid_feedstocks``: Output a file with off-grid feedstock cost information
- ``community_energy_access``:
  
  - ``False``: Not considered  
  - ``float``: Percentage of population or number of households to include  

Notes
-----

- All required input files must match the selected country, year, and plant type
- Changes to scenario settings (e.g. country, weather year, or plant type) will trigger Snakemake to rerun relevant steps; other changes may not be picked up automatically. If needed, use ``--forceall`` to force a full rerun.
- Ensure consistency between this file and the parameter files
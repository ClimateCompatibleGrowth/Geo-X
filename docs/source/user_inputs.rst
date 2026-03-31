User inputs
===========

Geo-X and Geo-X-data-prep both require user-provided input data.

Geo-X uses prepared spatial hexagon files together with techno-economic parameters and scenario settings. Geo-X-data-prep uses raw geospatial datasets to create these hexagon files.

Geo-X input data
-----------------

Geo-X requires three main types of input:

- Spatial input data
- Techno-economic parameters
- A configuration file

Spatial input data
~~~~~~~~~~~~~~~~~~

Geo-X uses prepared hexagon files for the area of interest.

These files must contain the spatial parameters required by the model, including distances to infrastructure and other location-specific attributes.

The final hexagon file must be placed in the ``data`` folder of the Geo-X repository using the required country code structure and can be created using the :doc:`helper tool <helper_tool>`.

Techno-economic parameter files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Parameter files are organised by country and commodity. Each country folder contains subfolders for each commodity (e.g. ``hydrogen``, ``ammonia``, ``copper``), which include Excel and CSV files defining costs, technical assumptions, and demand.

The main files include:

- ``basic_[commodity]_plant``: Defines plant design parameters and available generators
- ``country_parameters.xlsx``: Country-specific costs (e.g. electricity, interest rates, lifetimes)
- ``demand_parameters.xlsx``: Locations and demand at each demand centre
- ``technology_parameters.xlsx``: Technology-specific assumptions (e.g. water costs, infrastructure)
- ``transport_parameters.xlsx``: Road transport costs and characteristics
- ``pipeline_parameters.xlsx``: Pipeline costs and capacities (hydrogen and ammonia)
- ``conversion_parameters.xlsx``: Conversion between commodity states (hydrogen and copper)

Not all files are required for every commodity. For example, pipeline parameters are not used for copper.

All parameter values should be checked and adjusted to match the scenario. Example templates can be found in the `parameters folder <https://github.com/ClimateCompatibleGrowth/Geo-X/tree/main/parameters>`_.

Configuration file
~~~~~~~~~~~~~~~~~~

Geo-X uses ``config.yaml`` to define the scenario.

This file specifies the country, weather year, commodity type, generator settings, transport options, and other modelling assumptions. Detailed information is available in the :doc:`/configuration` page.

Geo-X-data-prep input data
--------------------------

Geo-X-data-prep is used to create Geo-X-ready spatial input files from raw geospatial datasets.

It requires a set of raw geospatial datasets, including:

- Country boundary data
- Ocean and sea boundary data
- OpenStreetMap layers
- Land cover data

Depending on the analysis, additional datasets may also be needed for:

- Hydropower
- Geothermal
- Nuclear
- Grid infrastructure
- Slope exclusion

The exact files and naming conventions are described in the Geo-X-data-prep repository `README <https://github.com/ClimateCompatibleGrowth/Geo-X-data-prep/blob/main/README.md>`_.

Notes
-----

- Required files depend on the country being analysed
- Additional files are needed for some technologies and scenarios
- File names and folder structure must match the conventions used by the repositories

.. code-block:: console

   # conda environments:
   geox
   prep
   spider
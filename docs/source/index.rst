===================
Geo-X documentation
===================

Geo-X is an `open-source <https://opensource.com/resources/what-open-source>`_ tool for analysing where commodities can be produced most cost-effectively.

It combines geospatial data, weather data, and energy system optimisation to estimate production costs across different locations.

Geo-X calculates the levelised cost of production at different locations, taking into account energy supply, production, storage, transport, and conversion where needed. Results are reported at chosen demand centres.


Geo-X currently supports:

- Green hydrogen
- Green ammonia
- Green copper

The model evaluates many possible production locations and identifies the lowest-cost supply options. These results can be compared with current or future market prices to assess competitiveness.

Electricity supply is configurable. The model currently supports:

- Solar
- Wind
- Hydropower
- Geothermal
- Nuclear

Geo-X uses a Snakemake workflow, making it reproducible and easy to run across different scenarios and regions.

Origin
------
Geo-X was developed from the `GeoH2 model <https://www.sciencedirect.com/science/article/pii/S2215016124001146>`_.

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Getting Started

   introduction
   installation
   helper_tool
   
.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: User guide

   snakemake
   configuration
   user_inputs
   running_model
   support

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Reference

   license
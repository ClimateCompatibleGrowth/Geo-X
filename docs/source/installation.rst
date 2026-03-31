Installation
============

This page describes how to install Geo-X and its dependencies.

Geo-X is run using `Snakemake <https://snakemake.readthedocs.io/>`_, with a Python environment managed using `mamba <https://mamba.readthedocs.io/>`_. It also requires access to climate data and an optimisation solver.

Prerequisites
-------------

Before installing Geo-X, ensure you have the following:

- `Git <https://git-scm.com/>`_ for cloning the repository
- `mamba <https://mamba.readthedocs.io/>`_ (or conda) for environment management
- A `Climate Data Store (CDS) API key <https://cds.climate.copernicus.eu/how-to-api>`_
- An optimisation solver compatible with PyPSA (e.g. CBC or Gurobi)

Clone the repository
--------------------

Clone the Geo-X repository from GitHub:

.. code-block:: console

   git clone https://github.com/ClimateCompatibleGrowth/Geo-X.git

To move into the repository directory:

.. code-block:: console

   cd Geo-X

Create the environment
----------------------

Geo-X provides an ``environment.yaml`` file with all required dependencies.

Create a new environment with:

.. code-block:: console

   mamba env create -f environment.yaml

To activate the environment:

.. code-block:: console

   mamba activate geox

Set up the Climate Data Store API
---------------------------------

Geo-X uses ERA5 weather data for modelling renewable energy generation.

To enable data download:

1. Register for an `ECMWF account <https://accounts.ecmwf.int/auth/realms/ecmwf/login-actions/registration?client_id=cds&tab_id=-DMamVkY_24>`_
2. Generate an `API key <https://cds.climate.copernicus.eu/how-to-api>`_
3. Add the key to your local configuration file (``.cdsapirc``)

.. note::

   If weather data downloads fail, check that your CDS API key is correctly configured and that you have accepted the necessary licences.

Install an optimisation solver
-------------------------------

Geo-X requires a solver to run the plant optimisation step.

Supported options include:

- CBC (free and open source)
- Gurobi (commercial, free academic licences available)

Install your chosen solver following its documentation, and ensure it is accessible from your environment.

Computational requirements
--------------------------

Geo-X is designed to run on standard machines, but requirements depend on the size of the study area and the amount of weather data used.

As a general guide:

- CPU: Multi-core support is useful, but not necessary (up to 4 cores can be used)
- Memory: A few GB is sufficient for small to medium cases
- Disk space: Weather data can require several GB depending on spatial and temporal resolution

For large regions or multi-year analyses, higher memory and storage may be required.

Version-specific note
---------------------

If you are using a version of the repository from commit ``ba7abe3`` (2 October 2025) or later, weather data must be downloaded again on first use.

This is due to changes in:

``src/prep/get_weather_data.py``

Next steps
----------

Once installation is complete, you can:

- Prepare input data (see :doc:`/user_inputs`)
- Configure your scenario (see :doc:`/configuration`)
- Run the workflow using Snakemake
Helper tool
===========

Geo-X-data-prep provides spatial data preparation tools for Geo-X users.

Geo-X requires hexagon-based spatial input files with associated parameters for each country. 

This helper tool converts raw geospatial inputs into Geo-X-ready files by working with the `Slope-Exclusion <https://github.com/ClimateCompatibleGrowth/Slope-Exclusion/tree/main>`_, `GLAES <https://github.com/FZJ-IEK3-VSA/glaes/tree/master/>`_, and `SPIDER <https://github.com/carderne/ccg-spider/tree/main>`_ repositories.

It produces final hexagon files for use in Geo-X, along with optional supporting files for technologies such as hydropower, geothermal, and nuclear.

Installation
------------

Clone the repository with submodules:

.. code-block:: console

   git clone --recurse-submodules https://github.com/ClimateCompatibleGrowth/Geo-X-data-prep.git

Move into the repository directory and create the main environment:

.. code-block:: console

   cd Geo-X-data-prep
   mamba env create -f environment.yaml
   mamba activate prep

After installation, deactivate the environment before creating the SPIDER environment:

.. code-block:: console

   mamba deactivate

Install the SPIDER environment from the ``ccg-spider/prep`` folder:

.. code-block:: console

   cd ccg-spider/prep
   mamba create -n spider
   mamba activate spider
   mamba install pip gdal
   pip install -e .

Return to the top-level repository directory.

Preparing input data
--------------------

Geo-X-data-prep requires a set of raw geospatial datasets (e.g. country boundaries, OpenStreetMap data, and land cover data).

Additional datasets may be needed depending on the technologies being included (e.g. hydropower, geothermal, nuclear, or grid data).

See the :doc:`user inputs </user_inputs>` page and the repository `README <https://github.com/ClimateCompatibleGrowth/Geo-X-data-prep/blob/main/README.md>`_ for full details on required files and naming conventions.

Workflow overview
-----------------

Geo-X-data-prep is organised into three main stages:

1. Initial data preparation before SPIDER
2. Running SPIDER to generate hexagon files
3. Final data preparation after SPIDER

Optional preprocessing steps are available for slope exclusion and grid data.

Notes
-----

- Country names must match those used in the Natural Earth dataset
- Some steps require switching between environments (``prep`` and ``spider``)
- Some steps must be run separately for each country
- Intermediate files can be large; cleaning up unused files after runs is recommended

Support
-------

Additional support information is available in the :doc:`support section </support>`.

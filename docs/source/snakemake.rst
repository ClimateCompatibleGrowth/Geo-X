Snakemake
=========

Geo-X uses `Snakemake <https://snakemake.readthedocs.io/>`_ to manage the workflow.

Snakemake handles dependencies between workflow steps, so you can run a single rule and it will automatically execute any required upstream steps first. This makes the workflow easier to reproduce and reduces the need to run each step manually.

How the workflow is organised
-----------------------------

Geo-X is built around a set of Snakemake rules, each of which performs a specific task in the modelling workflow.

The main stages are:

- Downloading weather data and creating cutouts
- Running transport, water, and plant optimisation steps
- Combining results into total cost outputs
- Producing maps and component breakdowns

You do not usually need to run these steps one by one. Snakemake will run the correct sequence of rules based on the output you request.

Running Geo-X
-------------

Geo-X is run by calling Snakemake with a target rule (the output you want to generate).

Common examples include:

.. code-block:: console

   snakemake -j [NUMBER OF CORES TO BE USED] run_weather

.. code-block:: console

   snakemake -j [NUMBER OF CORES TO BE USED] optimise_all

.. code-block:: console

   snakemake -j [NUMBER OF CORES TO BE USED] map_all

The number of cores must be specified when running any rule. Geo-X currently supports up to 4 cores.

If input files are changed after a run has completed, the same command can be run again. Snakemake will only rerun the steps that are needed to update the outputs.

See :doc:`/running_model` for a full example.

Main rules
----------

The most commonly used rules are:

- ``run_weather``: downloads weather data for the selected country and year
- ``optimise_all``: runs the full cost optimisation workflow
- ``map_all``: generates cost maps and component breakdowns

Lower-level rules are also available for specific tasks, such as transport optimisation, water cost calculations, or generating intermediate outputs.

Input files and configuration
-----------------------------

Snakemake uses the scenario settings in ``config.yaml`` together with the input files stored in the ``data`` and ``parameters`` folders.

Before running the workflow, make sure that:

- The required hexagon file is present
- All parameter files are in the correct country and commodity folders
- Required weather files are available in ``cutouts``
- Any optional technology files are present if required

Common outputs
--------------

Depending on the rule that is run, Geo-X will produce files in folders such as:

- ``cutouts``
- ``resources``
- ``results``
- ``plots``

These include intermediate geospatial files, final cost outputs, and visualisations.

Notes
-----

- Some rules depend on weather data being available before they can run.
- Weather data files may need to be downloaded again after updates to the data preparation workflow.
- Large runs may take several hours depending on the country, the number of demand centres, and the temporal resolution of the weather data (e.g. hourly or 3-hourly).
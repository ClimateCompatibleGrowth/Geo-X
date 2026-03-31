Support
=======

This page provides support for Geo-X and the helper tool, :doc:`Geo-X-data-prep </helper_tool>`.

Geo-X
-----

Common issues when running Geo-X:

Input data incorrect
~~~~~~~~~~~~~~~~~~~~

- Check that all input files are present, correctly named, and follow the required structure
- Check that parameter values are valid and consistent with the selected country and commodity

Terminal seems stuck and not showing updates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Check your system's activity monitor or task manager to confirm that Python is still running
- Some steps may take time without producing terminal output

Weather data download fails
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the ``run_weather`` rule fails:

- Check that your CDS API key is correctly configured (``.cdsapirc``)
- Check you have accepted the ERA5 dataset licence
- Confirm you have an active internet connection

Model does not run or fails during optimisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Check that a compatible solver (e.g. CBC or Gurobi) is installed and any required licences are configured
- Check the solver name in ``config.yaml`` matches your installation
- Check that all required input files are present and correctly formatted

Missing or incorrect results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Check that the hexagon input file is correctly placed in ``data/[ISO CODE]/``, where ``ISO CODE`` is the country's 2-letter ISO code
- Check that parameter files match the selected country and commodity
- Confirm that configuration settings (e.g. ``plant_type``) are consistent with input data
- If demand is set too high relative to available resources, the model may return zero or null results

Snakemake automatic rerun issue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- If Snakemake does not rerun after changes to inputs or configuration, use ``--forceall`` to trigger a full rerun

Geo-X-data-prep
---------------

Common issues when preparing spatial data:

Input data incorrect
~~~~~~~~~~~~~~~~~~~~

- Check that all required datasets are present and correctly placed in the ``data`` folder
- Check that file names match the expected naming conventions

Terminal seems stuck and not showing updates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Check your system's activity monitor or task manager to confirm that Python is still running
- Some steps may take time without producing terminal output

Building wheel error
~~~~~~~~~~~~~~~~~~~~

- If SPIDER installation fails when building the wheel, install the missing dependency using ``conda install fiona``
- Rerun ``pip install -e .`` after installing dependencies

Environment issues
~~~~~~~~~~~~~~~~~~

- Check the correct environment is activated (``prep`` or ``spider``)
- Check that required packages (e.g. GDAL) are installed
- If issues persist, recreate the environment using the appropriate YAML file

Getting help
------------

If you are unable to resolve an issue:

- Check the repository documentation and README files
- Search existing issues on GitHub
- Open a new issue with:

  - A clear description of the problem  
  - The command used  
  - Any error messages  
  - Your system setup (OS, environment, solver)
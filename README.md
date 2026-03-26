# Geo-X
**Geospatial analysis of commodity production costs**

> [!IMPORTANT] 
> Using a version of this repository from and after **commit ba7abe3** (committed on October 2nd 2025) requires weather data to be re-collected again when used for the first time. This is due to changes made in `src/prep/get_weather_data.py`.

Geo-X calculates the **levelised cost (LC)** of commodity production at specified demand locations, accounting for **production, storage, transport, and where relevant conversion costs** across all potential production locations.

These results can be compared with current or projected regional prices to assess economic competitiveness.

### Implemented commodities

Geo-X currently supports the following commodities:
- **Green hydrogen**
- **Green ammonia**
- **Green copper**

Copper production is modelled under **off-grid, hybrid, and grid-connected** electricity supply scenarios, which are reported separately. Hydrogen and ammonia production, by contrast, are restricted to **off-grid renewable electricity**.

The codebase is written in a **modular and extensible** manner, allowing additional commodities to be implemented in the future.

### Example use cases

The repository includes illustrative case studies:
- **Namibia**: hydrogen and ammonia production
- **Zambia**: copper production

These should be changed for your run(s).
Parameter references for the hydrogen case are provided at the end of this file.  
Because the framework is generalised, Geo-X can be applied to a wide range of geographic regions beyond these examples.

Geo-X builds upon:
- A preliminary model by **Leander Müller** (CC-BY-4.0):  
  https://github.com/leandermue/GEOH2
- Code developed by **Nick Salmon** (MIT licence):  
  https://github.com/nsalmon11/LCOH_Optimisation
___
# Contents
- [Setup instructions](#setup-instructions)
- [Preparing input data](#preparing-input-data)
- [Running Geo-X with Snakemake](#running-geo-x-with-snakemake)
- [Definitions of all rules](#definitions-of-all-rules)
- [Limitations](#limitations)
- [Citation](#citation)
- [Case study parameters](#case-study-parameters)

# Setup instructions

Follow the steps below before running Geo-X.

Geo-X is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/), which automates workflow execution using **rules**. Each rule corresponds to a specific part of the codebase and performs a distinct function.

This installation requires **Git** for version control. If Git is not installed, follow these instructions:
https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

## 1) Clone the repository

Clone the Geo-X repository using Git. 

`... % git clone https://github.com/ClimateCompatibleGrowth/Geo-X.git`

## 2) Set up the environment

The Python package requirements are specified in the `environment.yaml` file.

Install these requirements in a new environment using the `mamba` package and environment manager (installation instructions are available [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)):

` .../geo-x % mamba env create -f environment.yaml`

Activate the newly created environment:

`.../geo-x % mamba activate geox`

## 3) Get a Climate Data Store (CDS) API key

To simulate commodities produced using renewable energy, historical climate data are required as an input.

The `get_weather_data` rule can be used to download the relevant data from the ERA5 reanalysis dataset using [Atlite](https://atlite.readthedocs.io/en/latest/) to create a weather cutout.

To enable this step, register for and configure your CDS API key as described on the [Climate Data Store website](https://cds.climate.copernicus.eu/how-to-api).

## 4) Install a solver

To run the `plant_optimization` rule, an optimisation solver must be installed. 

You can use any solver compatible with [PyPSA](https://docs.pypsa.org/v0.26.0/installation.html), such as:
- **Cbc**: free and open source
- **Gurobi**: commercial, with free academic licences available

Install your chosen solver following the instructions provided in the solver's documentation. 
___

# Preparing input data

Before running Geo-X, you must prepare the required input data.

Geo-X uses **three types of input data**:
1. Spatial input data  
2. Techno-economic and plant parameters  
3. A configuration file  

If you intend to include **hydropower** or **geothermal** or **nuclear** as generation technologies, additional input files are required, as described below.

## 1) Prepare input hexagons 

> [!NOTE]
> This step (and any optional steps, if required) must be completed for **each country** you wish to analyse.

Spatial input data are represented as a set of hexagons using the [H3 standard](https://h3geo.org/).

To analyse a different area of interest, the hexagon file must be replaced with one that follows the same structure and conventions as the provided example. Both example hexagon files might require updating at the time of use.

A full walkthrough for creating these hexagons, including the required tools, is available in the [Geo-X-data-prep](https://github.com/ClimateCompatibleGrowth/Geo-X-data-prep) repository.

### Required hexagon attributes

Each hexagon file must contain the following attributes:

  - `waterbody_dist`: Distance to selected waterbodies in area of interest
  - `waterway_dist`: Distance to selected waterways
  - `road_dist`: Distance to road network
  - `theo_pv`: Theoretical potential of standardised PV plants (can be calculated with [GLAES](https://github.com/FZJ-IEK3-VSA/glaes))
  - `theo_wind`: Theoretical potential of standardised wind turbines (can be calculated with [GLAES](https://github.com/FZJ-IEK3-VSA/glaes))

### Additional attributes

For **coastal** countries, the following attribute is required:
  - `ocean_dist`: Distance to ocean coastline

For **copper** analyses, the following additional attributes are required:

  - `grid_dist`: Distance to the electricity grid network
  - `pop`: Population within the hexagon (required when `community_energy_access` is due to be set to a float value in the config file)

The [ccg-spider](https://github.com/carderne/ccg-spider) repository can be used to combine multiple spatial data into standard H3 hexagons.

### File placement

Once the hexagon file has been created, place it in the following location:

`data/[COUNTRY ISO CODE]/hex_final_[COUNTRY ISO CODE].geojson`

> [!IMPORTANT]  
> `COUNTRY ISO CODE` must be the country’s [ISO](https://www.iban.com/country-codes) standard 2-letter abbreviation.

### (Optional) Hydropower files 

Hydropower input files can be obtained by:
- Running the [Geo-X-data-prep](https://github.com/ClimateCompatibleGrowth/Geo-X-data-prep) repository 
- Downloading data from the [HydroBASINS section of HydroSHEDS](https://www.hydrosheds.org/products/hydrobasins), choose a level that suits the country being analysed

Place the following files inside `data/[COUNTRY ISO CODE]/`, one file must be renamed to match the required format:
- `[COUNTRY ISO CODE]_hydropower_dams.gpkg`
- All the files from the HydroBASINS download (or the Shapefile cannot be used)

### (Optional) Geothermal GeoPackage file

The geothermal GeoPackage file can be obtained by running the [Geo-X-data-prep](https://github.com/ClimateCompatibleGrowth/Geo-X-data-prep) repository.

Place the following file inside `data/[COUNTRY ISO CODE]/`, it must be renamed to match the required format:
- `[COUNTRY ISO CODE]_geothermal_plants.gpkg`

### (Optional) Nuclear GeoPackage file

The nuclear GeoPackage file can be obtained by running the [Geo-X-data-prep](https://github.com/ClimateCompatibleGrowth/Geo-X-data-prep) repository.

Place the following file inside `data/[COUNTRY ISO CODE]/`, it must be renamed to match the required format:
- `[COUNTRY ISO CODE]_nuclear_plants.gpkg`

## 2) Prepare input parameter Excel files

Next, prepare the **techno-economic input parameter spreadsheets**.

Required input parameters include the spatial area of interest, total annual demand for the commodity, and prices and cost of capital for infrastructure investments.

These values can be either **current conditions** or **projected values** for a single snapshot in time.

All parameter values used by the model are specified in a set of Excel and CSV files located in the `parameters` folder. All files should be checked to ensure they are correctly set up before running the model.

### Parameter file overview

- **Basic plant parameters**
  The `basic_[COMMODITY ABBREVIATION]_plant` folder contains several CSV files defining global parameters for optimizing the plant design.
    - All power units must be in MW
    - All energy units must be in MWh
    - The `generators.csv` file **must only include** the generators that you wish to use in the optimisation

For further details on these parameters, refer to the [PyPSA documentation](https://docs.pypsa.org/v0.26.0/components.html).

> [!IMPORTANT]
> `COMMODITY ABBREVIATION` can be `h2`, `nh3` or `cu`, corresponding to the currently implemented commodities.

- **Conversion parameters**
  The `conversion_parameters.xlsx` file defines parameters related to converting between different states of the commodity.
  This file is required for the `hydrogen` and `copper` folders.

- **Country parameters**
  The`country_parameters.xlsx` file includes country- and technology-specific parameters such as:
  - Interest rates
  - Heat and electricity costs
  - Asset lifetimes.

  Additional notes:
  - Interest rates should be expressed as decimals (e.g. 5% as `0.05`).
  - Asset lifetimes should be in years.
  - The file may contain columns for generators that are not being considered.

- **Demand parameters**
  The `demand_parameters.xlsx` file defines the set of demand centers.

  For each demand center, the following must be specified:
  - Latitude and longitude
  - Annual demand
  - Commodity state at the demand location

- **Pipeline parameters**
  The `pipeline_parameters.xlsx` file contains price, capacity, and lifetime data for different pipeline sizes.
  This file is required for the `hydrogen` and `ammonia` folders.

- **Technology parameters**
  The `technology_parameters.xlsx` file includes parameters related to:
  - Water costs and/or usage
  - Road infrastructure

- **Transport parameters**
  The `transport_parameters.xlsx` file defines parameters for road transport of the commodity, including:
  - Truck speed
  - Transport cost
  - Asset lifetime
  - Transport capacity

> [!IMPORTANT]
> All parameter files must be stored in a folder structure organised by **country** and **commodity**. 
> Each commodity must be placed in a sub-folder named after the commodity, within a folder named using the corresponding **country’s ISO standard 2-letter abbreviation**.
> As currently implemented, the commodity must be one of `hydrogen`, `ammonia`, or `copper`.
> For the illustrative case studies:
> - Namibia uses the folder `NA`, with sub-folders `hydrogen` and `ammonia`
> - Zambia uses the folder `ZM`, with the sub-folder `copper`

## 3) Modify the config file

High-level workflow settings are controlled in the configuration file: `config.yaml`.

### Wildcards

Wildcards are defined in the `scenario` section of the config file and determine which data are used in the workflow.

They can be changed to match the case you're analysing:

- `country`: ISO standard 2-letter country abbreviation
- `weather_year`: A 4-digit year from 1940 onwards included in the ERA5 dataset
- `plant_type`: The commodity to be produced (`hydrogen`, `ammonia`, or `copper`)

### Weather data

The amount of weather data used in the analysis is controlled using the following parameters:

- `years_to_check`: The number of years of weather data to use (must be an integer value)

For example, if you want weather data to start in 2015 and span five years, set  
`weather_year` to `2015` and `years_to_check` to `5`.

The temporal resolution of the data used in the optimisation can be set using:

- `freq`: Data frequency (e.g., `1H` for hourly, `3H` for three-hourly)

### Generators

The renewable generation technologies considered for plant construction are defined in the `generators_dict` section.

Currently supported generator types are:
- `solar`
- `wind`
- `hydro`
- `geothermal`
- `nuclear`

Ensure that all generators you wish to consider are included in the dictionary, and remove any that are not required.

Generator-specific technology assumptions can be modified as follows:
- `panel`: Solar panel technology (used if `solar` is included)
- `turbine`: Wind turbine technology (used if `wind` is included)

Installed capacity assumptions are defined in the `gen_capacity` section:
- `solar` and `wind`: Values can be adjusted as needed for the analysis
- `hydro`, `geothermal` and `nuclear`: Values should not be changed

### Efficiency of generator plants

Average generation efficiencies can be specified for:
- `hydro`
- `geothermal`
- `nuclear`

Efficiencies should be provided as float values and adjusted as necessary.

### Other configuration options

- `solver`: Specify the name of the optimisation solver to be used

#### Water constraints
The `water` section controls whether water scarcity is considered:

- `has_ocean_dists`: Set to `True` if country is coastal, otherwise `False` for a landlocked country
- `has_limit`: Set to `True` to include a water availability constraint. Only required for `hydrogen` and `ammonia`
- `annual_limit`: Annual water availability limit in cubic meters (required if `has_limit` is set to `True`)

#### Transport infrastructure
The `transport` section controls whether infrastructure construction is permitted:

- `pipeline_construction`: Enable or disable pipeline construction. Only required for `hydrogen` and `ammonia`
- `road_construction`: Enable or disable road construction

#### Currency
Specify the `currency` used across all parameter files.

For example, in any case study’s `technology_parameters.xlsx` file, the **Water** sheet contains a row labelled `Water specific cost (euros/m3)`. If you change the currency in the parameter files, you must also update this value to ensure consistency across the model.

#### Spatial filtering
- `restrict_hexagons`: If set to `True`, results will only be produced for hexagon(s) containing demand centres. All other hexagons will be set to `null`.

### Copper-only configuration options

The following settings are only relevant when `plant_type` is set to `copper`:

- `grid_construction`: If `True`, grid infrastructure may be constructed to supply electricity

- `size_offgrid_feedstocks`: If `True`, a file containing costs related to an off-grid feedstock system will be outputted to the `results/` folder

- `community_energy_access`: Either `False` or a float value:
  - If `False`, community energy access is not considered
  - If a float is provided, it represents either:
    - The percentage of the off-grid population within a hexagon to consider
    - An absolute number of households

> [!NOTE]  
> `country` and `weather_year` may be provided as lists if multiple countries or years are being analysed.  
> You must ensure that all required input files for each country and year are present in the correct locations, as described earlier.
___

# Running Geo-X with Snakemake

Once the repository setup is complete and all input data have been prepared, you are ready to run Geo-X!

This repository uses [Snakemake](https://snakemake.readthedocs.io/en/stable/) to automate the workflow.
For a gentle introduction to Snakemake, see [Getting Started with Snakemake](https://carpentries-incubator.github.io/workflows-snakemake/) on The Carpentries Incubator.

Geo-X is executed by specifying a **rule name** when running Snakemake. Snakemake will automatically resolve dependencies and execute all required rules defined in the `Snakefile` to produce the requested outputs.

### Computational resources

When running any Snakemake rule, you must specify the number of CPU cores to use; `NUMBER OF CORES TO BE USED`. The maximum is 4 cores.

### Expected run times

Run times vary depending on the rule being executed:

- `get_weather_data`: Depending on country size and internet connection, this rule may take from a few minutes to several hours. Ensure that sufficient disk space is available, as weather data files can be several GB in size.

- `optimize_plant`: Depending on country size and the number of demand centres, this rule may take from several minutes to several hours.

- All other rules: Usually complete within a few seconds.


## 1) Get weather data (if needed)

If you already have weather data for your area of interest, you can proceed directly to the next step.
Otherwise, run the weather data rule to download the required data from the CDS API:

```
snakemake -j [NUMBER OF CORES TO BE USED] run_weather
```

## 2) Run optimization or mapping rules

Next, run the Geo-X optimisation workflow.

The rules in this section execute the **entire workflow**, so you do not need to run each individual rule step-by-step. If any input files are modified after a completed run, the same command can be executed again and Snakemake will automatically re-run only the steps required to update the results.

Before running these rules, ensure that the required weather file(s) are present in the `cutouts` directory. Each file must be named `[COUNTRY ISO CODE]_[WEATHER YEAR].nc` for each country being analysed.

### Run total commodity cost optimisation

To calculate the total commodity cost for all scenarios, run:
```
snakemake -j [NUMBER OF CORES TO BE USED] optimise_all
```

### Generate component breakdown and cost maps

To generate a further breakdown of the energy component and spatial maps of commodity costs for all scenarios, run:
```
snakemake -j [NUMBER OF CORES TO BE USED] map_all
```
___
# Definitions of all rules

All rules are documented below for completeness; however, **you do not need to run each rule individually**.

In most cases, you can simply run one of the optimisation or mapping rules described in the previous section. Snakemake will automatically ensure that all required rules are executed in the correct order.

> [!NOTE]
> `PLANT TYPE` refers to the commodity being produced (currently either `hydrogen`, `ammonia`, or `copper`).

### `get_weather_data` rule
Downloads weather data for the area of interest to support generation optimisation.
```
snakemake -j [NUMBER OF CORES TO BE USED] cutouts/[COUNTRY ISO CODE]_[WEATHER YEAR].nc
```

### `optimize_transport` rule

Calculates the cost of the optimal commodity transportation strategy from each hexagon to each demand centre, using available transport options.

This rule uses parameters from:
- `technology_parameters.xlsx`
- `demand_parameters.xlsx`
- `country_parameters.xlsx`

For hydrogen and copper, the parameters from `conversion_parameters.xlsx` are also included.
```
snakemake -j [NUMBER OF CORES TO BE USED] resources/hex_transport_[COUNTRY ISO CODE]_[PLANT TYPE].geojson
```

### `calculate_water_costs` rule

Calculates water supply costs from water sources for commodity production in each hexagon.

This rule uses parameters from:
- `technology_parameters.xlsx`
- `country_parameters.xlsx`
```
snakemake -j [NUMBER OF CORES TO BE USED] resources/hex_water_[COUNTRY ISO CODE]_[PLANT TYPE].geojson
```

### `optimize_plant` rule

Designs a production plant to meet the commodity demand profile for each demand centre, accounting for all available transportation options.

Ensure that:
- Plant parameters are specified in `basic_[COMMODITY ABBREVIATION]_plant`
- Demand centres are defined in `demand_parameters.xlsx`
```
snakemake -j [NUMBER OF CORES TO BE USED] resources/hex_lc_[COUNTRY ISO CODE]_[WEATHER YEAR]_[PLANT TYPE].geojson
```

### `calculate_total_cost` rule

Combines calculated componenets to identify the lowest-cost supply option for each demand centre.
```
snakemake -j [NUMBER OF CORES TO BE USED] results/hex_total_cost_[COUNTRY ISO CODE]_[WEATHER YEAR]_[PLANT TYPE].geojson
```

### `calculate_cost_components` rule

Calculates the cost contribution of components within the energy component for each hexagon.
```
snakemake -j [NUMBER OF CORES TO BE USED] results/hex_cost_components_[COUNTRY ISO CODE]_[WEATHER YEAR]_[PLANT TYPE].geojson
```

### `calculate_map_costs` rule
Generates spatial visualisations of different cost components per kilogram of commodity.
```
snakemake -j [NUMBER OF CORES TO BE USED] plots/[COUNTRY ISO CODE]_[WEATHER YEAR]_[PLANT TYPE]
```
___

# Limitations

Geo-X has several important limitations that should be considered when interpreting results.

The model primarily considers greenfield generation, meaning that it does not account for existing generation assets.  
Electricity supply is assumed to be based on newly constructed infrastructure, and any excess electricity generation is assumed to be fully curtailed.

Grid electricity is not considered for most commodities, except where explicitly enabled (e.g. grid-connected and hybrid scenarios for copper) and still takes a broad approach with calculating the costs associated.

While the design of the commodity plant is convex and therefore guaranteed to find the global optimum solution if it exists, the selection of the trucking strategy is greedy to avoid the long computation times and potential computational intractability associated with a mixed-integer optimization problem.

Currently, only land transport is considered in the model. 
To calculate the cost of commodity production for export, any additional costs for conversion and transport via ship or undersea pipeline must be added in post-processing.

Transport costs are calculated from the center of the hexagon to the demand center. 
When using large hexagon sizes, this assumption may over- or underestimate transportation costs significantly. 
Additionally, only path length is considered when calculating the cost of road and pipeline construction. 
Additional costs due to terrain are not considered.
___

# Citation

If you decide to use Geo-X, please kindly cite us using the following: 

*Halloran, C., Leonard, A., Salmon, N., Müller, L., & Hirmer, S. (2024). 
GeoH2 model: Geospatial cost optimization of green hydrogen production including storage and transportation. 
MethodsX, 12, 102660. https://doi.org/10.1016/j.mex.2024.102660.*

```commandline
@article{Halloran_GeoH2_model_Geospatial_2024,
author = {Halloran, Claire and Leonard, Alycia and Salmon, Nicholas and Müller, Leander and Hirmer, Stephanie},
doi = {10.1016/j.mex.2024.102660},
journal = {MethodsX},
month = jun,
pages = {102660},
title = {{GeoH2 model: Geospatial cost optimization of green hydrogen production including storage and transportation}},
volume = {12},
year = {2024}
}
```
___

# Case study parameters

This repository includes sample parameters for a hydrogen production case in Namibia.
References for these parameters are included in the tables below for reference.
For the results of this case, please refer to the GeoH2 model MethodsX article: https://doi.org/10.1016/j.mex.2024.102660. 

**Green hydrogen plant parameters:**

| Hardware                   | Parameter             | Value     | Units                | Ref.                                                                                                                            |
|----------------------------|-----------------------|-----------|----------------------|---------------------------------------------------------------------------------------------------------------------------------|
| Solar photovoltaic         | Capex                 | 1,470,000 | €/MW                 | [Allington et al., 2021](https://doi.org/10.1016/j.dib.2022.108021)                                                             |
| Wind turbines              | Capex                 | 1,580,000 | €/MW                 | [Allington et al., 2021](https://doi.org/10.1016/j.dib.2022.108021)                                                             |
| Hydrogen electrolysis      | Capex                 | 1,250,000 | €/MW                 | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                           |
| Hydrogen electrolysis      | Efficiency            | 0.59      | MWh H2/MWh el        | [Taibi et al., 2020](https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2020/Dec/IRENA_Green_hydrogen_cost_2020.pdf)  |
| Hydrogen compression       | Isentropic efficiency | 0.051     | MWh el/MWh H2        | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                           |
| Hydrogen storage unloading | Efficiency            | 1         | MWh H2/MWh H2-stored | Assumption                                                                                                                      |
| Battery                    | Capex                 | 95,000    | €/MW                 | [BloombergNEF, 2022](https://about.bnef.com/blog/lithium-ion-battery-pack-prices-rise-for-first-time-to-an-average-of-151-kwh/) |
| Hydrogen storage           | Capex                 | 21,700    | €/MWh                | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                           |

**Conversion parameters:**

| Process                | Parameter                    | Value       | Units             | Ref.                                                                                                                                             |
|------------------------|------------------------------|-------------|-------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| 500 bar compression    | Heat capacity                | 0.0039444   | kWh/kg/K          | Kurzweil and Dietlmeier, 2016                                                                                                                    |
| 500 bar compression    | Input temperature            | 298.15      | K                 | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| 500 bar compression    | Input pressure               | 25          | bar               | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| 500 bar compression    | Isentropic exponent          | 1.402       |                   | Kurzweil and Dietlmeier, 2016                                                                                                                    |
| 500 bar compression    | Isentropic efficiency        | 0.8         |                   | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| 500 bar compression    | Compressor lifetime          | 15          | years             | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar compression    | Compressor capex coefficient | 40,035      | €/kg H2/day       | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar compression    | Compressor opex              | 4           | % capex/year      | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| Hydrogen liquification | Electricity demand           | 9.93        | kWh/kg H2         | [Ausfelder and Dura](https://dechema.de/dechema_media/Downloads/Positionspapiere/2021_DEC_P2X_III_Technischer_Anhang.pdf)                        |
| Hydrogen liquification | Capex quadratic coefficient  | -0.0002     | €/(kg H2)^2       | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| Hydrogen liquification | Capex linear coefficient     | 1,781.9     | €/kg  H2          | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| Hydrogen liquification | Capex constant               | 300,000,000 | €                 | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| Hydrogen liquification | Opex                         | 8           | % capex/year      | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| Hydrogen liquification | Plant lifetime               | 20          | years             | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| LOHC hydrogenation     | Electricity demand           | 0.35        | kWh/kg H2         | [Andersson and Grönkvist, 2019](https://doi.org/10.1016/j.ijhydene.2019.03.063)                                                                  |
| LOHC hydrogenation     | Heat demand                  | -9          | kWh/kg H2         | [Hydrogenious, 2022](https://hydrogenious.net/how/)                                                                                              |
| LOHC hydrogenation     | Capex coefficient            | 0.84        | kWh/kg H2/year    | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC hydrogenation     | Opex                         | 4           | % capex/year      | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC hydrogenation     | Plant lifetime               | 25          | years             | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC hydrogenation     | Carrier costs                | 2           | €/kg carrier      | [Clark, 2020](https://hydrogenious.net/lohc-global-hydrogen-opportunity/)                                                                        |
| LOHC hydrogenation     | Carrier ratio                | 16.1        | kg carrier/kg  H2 | [Arlt and Obermeier, 2017](https://www.encn.de/fileadmin/user_upload/Studie_Wasserstoff_und_Schwerlastverkehr_WEB.pdf)                           |
| LOHC dehydrogenation   | Electricity demand           | 0.35        | kWh/kg H2         | [Andersson and Grönkvist, 2019](https://doi.org/10.1016/j.ijhydene.2019.03.063)                                                                  |
| LOHC dehydrogenation   | Heat demand                  | 12          | kWh/kg H2         | [Hydrogenious, 2022](https://hydrogenious.net/how/)                                                                                              |
| LOHC dehydrogenation   | Capex coefficient            | 2.46        | kWh/kg H2         | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC dehydrogenation   | Opex                         | 4           | % capex/year      | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC dehydrogenation   | Plant lifetime               | 25          | years             | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia synthesis      | Electricity demand           | 2.809       | kWh/kg H2         | [IEA, 2021](https://www.iea.org/reports/ammonia-technology-roadmap)                                                                              |
| Ammonia synthesis      | Capex coefficient            | 0.75717     | kWh/g H2/year     | [IEA, 2021](https://www.iea.org/reports/ammonia-technology-roadmap)                                                                              |
| Ammonia synthesis      | Opex                         | 1.5         | % capex/year      | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia synthesis      | Plant lifetime               | 25          | years             | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia cracking       | Heat demand                  | 4.2         | kWh/kg H2         | [Andersson and Grönkvist, 2019](https://doi.org/10.1016/j.ijhydene.2019.03.063)                                                                  |
| Ammonia cracking       | Capex coefficient            | 17,262,450  | kWh/g H2/hour     | [Cesaro et al., 2021](https://doi.org/10.1016/j.apenergy.2020.116009)                                                                            |
| Ammonia cracking       | Opex                         | 2           | % capex/year      | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| Ammonia cracking       | Plant lifetime               | 25          | years             | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |

**Trucking parameters:**

| Hardware                 | Parameter                  | Value   | Units        | Ref.                                                                                                                                             |
|--------------------------|----------------------------|---------|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| All trucks               | Average truck speed        | 70      | km/h         | Assumption                                                                                                                                       |
| All trucks               | Working hours              | 24      | h/day        | Assumption                                                                                                                                       |
| All trucks               | Diesel price               | 1.5     | €/L          | Assumption                                                                                                                                       |
| All trucks               | Driver wage                | 2.85    | €/h          | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| All trucks               | Working days               | 365     | days/year    | Assumption                                                                                                                                       |
| All trucks               | Max driving distance       | 160,000 | km/year      | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                                                                            |
| All trucks               | Truck capex                | 160,000 | €            | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| All trucks               | Truck Opex                 | 12      | % capex/year | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| All trucks               | Diesel consumption         | 35      | L/100 km     | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| All trucks               | Truck lifetime             | 8       | years        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| All trucks               | Trailer lifetime           | 12      | years        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| 500 bar hydrogen trailer | Trailer capex              | 660,000 | €            | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar hydrogen trailer | Trailer opex               | 2       | % capex/year | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar hydrogen trailer | Trailer capacity           | 1,100   | kg	H2        | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| 500 bar hydrogen trailer | Loading and unloading time | 1.5     | hours        | [Cerniauskas, 2021](https://juser.fz-juelich.de/record/906356)                                                                                   |
| Liquid hydrogen trailer  | Trailer capex              | 860,000 | €            | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| Liquid hydrogen trailer  | Trailer opex               | 2       | % capex/year | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| Liquid hydrogen trailer  | Trailer capacity           | 4,300   | kg H2        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| Liquid hydrogen trailer  | Loading and unloading time | 3       | hours        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| LOHC trailer             | Trailer capex              | 660,000 | €            | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| LOHC trailer             | Trailer opex               | 2       | % capex/year | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| LOHC trailer             | Trailer capacity           | 1,800   | kg H2        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| LOHC trailer             | Loading and unloading time | 1.5     | hours        | [Reuss et al., 2017](https://doi.org/10.1016/j.apenergy.2017.05.050)                                                                             |
| Ammonia trailer          | Trailer capex              | 210,000 | €            | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia trailer          | Trailer opex               | 2       | % capex/year | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia trailer          | Trailer capacity           | 2,600   | kg H2        | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |
| Ammonia trailer          | Loading and unloading time | 1.5     | hours        | [IEA, 2020](https://iea.blob.core.windows.net/assets/29b027e5-fefc-47df-aed0-456b1bb38844/IEA-The-Future-of-Hydrogen-Assumptions-Annex_CORR.pdf) |

**Road parameters:**

| Road length         | Parameter | Value      | Units     | Ref.                                                                  |
|---------------------|-----------|------------|-----------|-----------------------------------------------------------------------|
| Short road (<10 km) | Capex     | 626,478.45 | €/km      | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219) |
| Long road (>10 km)  | Capex     | 481,866.6  | €/km      | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219) |
| All roads           | Opex      | 7,149.7    | €/km/year | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219) |

**Pipeline parameters:**

| Pipeline size   | Parameter           | Value     | Units        | Ref.                                                                                             |
|-----------------|---------------------|-----------|--------------|--------------------------------------------------------------------------------------------------|
| All pipelines   | Opex                | 1.25      | % capex/year | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| All pipelines   | Availability        | 95        | %            | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)                            |
| All pipelines   | Pipeline lifetime   | 42.5      | years        | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| All pipelines   | Compressor lifetime | 24        | years        | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| All pipelines   | Electricity demand  | 0.000614  | kWh/kg H2/km | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Large pipeline  | Maximum capacity    | 13        | GW           | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Large pipeline  | Pipeline capex      | 2,800,000 | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Large pipeline  | Compressor capex    | 620,000   | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Medium pipeline | Maximum capacity    | 4.7       | GW           | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Medium pipeline | Pipeline capex      | 2,200,000 | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Medium pipeline | Compressor capex    | 310,000   | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Small pipeline  | Maximum capacity    | 1.2       | GW           | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Small pipeline  | Pipeline capex      | 90,000    | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |
| Small pipeline  | Compressor capex    | 90,000    | €/km         | [Jens et al., 2021](https://ehb.eu/files/downloads/European-Hydrogen-Backbone-April-2021-V3.pdf) |

**Water parameters:**

| Type        | Parameter                    | Value | Units              | Ref.                                                                                                                                                                              |
|-------------|------------------------------|-------|--------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Freshwater  | Treatment electricity demand | 0.4   | kWh/m^3 water      | [US Dept. of Energy, 2016](https://betterbuildingssolutioncenter.energy.gov/sites/default/files/Primer%20on%20energy%20efficiency%20in%20water%20and%20wastewater%20plants_0.pdf) |
| Ocean water | Treatment electricity demand | 3.7   | kWh/m^3 water      | [Patterson et al., 2019](https://doi.org/10.1073/pnas.1902335116)                                                                                                                 |
| All water   | Transport cost               | 0.1   | €/100 km/m^3 water | [Zhou and Tol, 2005](https://doi.org/10.1029/2004WR003749)                                                                                                                        |
| All water   | Water specific cost          | 1.25  | €/m^3 water        | [Wasreb, 2019](https://wasreb.go.ke/wasrebsystems/tariffs/about-us.html)                                                                                                          |
| All water   | Water demand                 | 21    | L water/kg H2      | [Taibi et al., 2020](https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2020/Dec/IRENA_Green_hydrogen_cost_2020.pdf)                                                    |

**Country-specific parameters:**

| Country | Parameter                    | Value   | Units | Ref.                                                                               |
|---------|------------------------------|---------|-------|------------------------------------------------------------------------------------|
| Namibia | Electricity price            | 0.10465 | €/kWh | [GlobalPetrolPrices.com]({https://www.globalpetrolprices.com/electricity_prices/}) |
| Namibia | Heat price                   | 0.02    | €/kWh | Assumption                                                                         |
| Namibia | Solar interest rate          | 6       | %     | Assumption                                                                         |
| Namibia | Solar lifetime               | 20      | years | Assumption                                                                         |
| Namibia | Wind interest rate           | 6       | %     | Assumption                                                                         |
| Namibia | Wind lifetime                | 20      | years | Assumption                                                                         |
| Namibia | Plant interest rate          | 6       | %     | Assumption                                                                         |
| Namibia | Plant lifetime               | 20      | years | Assumption                                                                         |
| Namibia | Infrastructure interest rate | 6       | %     | Assumption                                                                         |
| Namibia | Infrastructure lifetime      | 50      | years | [Müller et al., 2022](https://doi.org/10.1016/j.apenergy.2023.121219)              |

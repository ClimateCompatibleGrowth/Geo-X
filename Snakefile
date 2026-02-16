configfile: 'config.yaml'

wildcard_constraints:
    # ISO alpha-2 country code
    country = '[A-Z]{2}',
    # ERA5 weather year (1940-2023)
    weather_year = '(19[4-9]\d|20[0-1]\d|202[0-3])',

def get_transport_inputs(wildcards):
    input_dict = {'hexagons': f'data/{wildcards.country}/hex_final_{wildcards.country}.geojson', 
                'technology_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/technology_parameters.xlsx',
                'demand_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/demand_parameters.xlsx',
                'country_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/country_parameters.xlsx',
                'transport_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/transport_parameters.xlsx',
                'conversion_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/conversion_parameters.xlsx',
                'pipeline_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/pipeline_parameters.xlsx'
                }
        
    if wildcards.plant_type == 'copper':
        input_dict.pop("pipeline_parameters")
    elif wildcards.plant_type == 'ammonia':
        input_dict.pop("conversion_parameters")

    return input_dict

def get_water_inputs(wildcards):
    input_dict = {'hexagons': f'resources/hex_transport_{wildcards.country}_{wildcards.plant_type}.geojson', 
                'technology_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/technology_parameters.xlsx',
                'demand_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/demand_parameters.xlsx',
                'country_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/country_parameters.xlsx',
                'conversion_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/conversion_parameters.xlsx',
                }
        
    if wildcards.plant_type == 'copper':
        input_dict.pop("country_parameters")
    elif wildcards.plant_type == 'ammonia' or wildcards.plant_type == 'hydrogen':
        input_dict.pop("conversion_parameters")
        input_dict.pop("demand_parameters")

    return input_dict

def get_lc_inputs(wildcards):
    input_dict = {'hexagons': f'resources/hex_water_{wildcards.country}_{wildcards.plant_type}.geojson', 
                'demand_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/demand_parameters.xlsx',
                'country_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/country_parameters.xlsx',
                'transport_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/transport_parameters.xlsx',
                'technology_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/technology_parameters.xlsx',
                'conversion_parameters': f'parameters/{wildcards.country}/{wildcards.plant_type}/conversion_parameters.xlsx',
                }
    
    if wildcards.plant_type == 'copper':
        input_dict['plant_folder'] = 'parameters/basic_cu_plant'
    elif wildcards.plant_type == 'ammonia':
        input_dict['plant_folder'] = 'parameters/basic_nh3_plant'
        input_dict.pop("conversion_parameters")
        input_dict.pop("technology_parameters")
    elif wildcards.plant_type == 'hydrogen':
        input_dict['plant_folder'] = 'parameters/basic_h2_plant'
        input_dict.pop("conversion_parameters")
        input_dict.pop("technology_parameters")

    return input_dict

# Bulk run rule to run all countries and years listed in config file
rule optimise_all:
   input:
        expand('results/hex_total_cost_{country}_{weather_year}_{plant_type}.geojson',
        **config['scenario']
        ),
        
# Bulk run rule to map all countries and years listed in config file
rule map_all:
   input:
       expand('plots/{country}_{weather_year}_{plant_type}',
        **config['scenario']
        ),
        
# Clearer rule to run weather
rule run_weather:
    input:
        expand('cutouts/{country}_{weather_year}.nc',
        **config['scenario']
        ),


rule get_weather_data:
    input:
        hexagons = 'data/{country}/hex_final_{country}.geojson'
    output:
        'cutouts/{country}_{weather_year}.nc',
    script:
        'src/prep/get_weather_data.py'

rule optimize_transport:
    input:
        unpack(get_transport_inputs)
    output:
        'resources/hex_transport_{country}_{plant_type}.geojson'
    script:
        'src/main/transport_optimization.py'

rule calculate_water_costs:
    input:
        unpack(get_water_inputs)
    output:
        'resources/hex_water_{country}_{plant_type}.geojson'
    script:
        'src/main/water_cost.py'

rule optimize_plant:
    input:
        unpack(get_lc_inputs)
    output:
        'resources/hex_lc_{country}_{weather_year}_{plant_type}.geojson'
    script:
        'src/main/plant_optimization.py'

rule calculate_total_costs:
    input:
        hexagons = 'resources/hex_lc_{country}_{weather_year}_{plant_type}.geojson',
        demand_parameters = 'parameters/{country}/{plant_type}/demand_parameters.xlsx'
    output:
        'results/hex_total_cost_{country}_{weather_year}_{plant_type}.geojson'
    script:
        'src/main/total_costs.py'

rule calculate_cost_components:
    input:
        hexagons = 'results/hex_total_cost_{country}_{weather_year}_{plant_type}.geojson',
        demand_parameters = 'parameters/{country}/{plant_type}/demand_parameters.xlsx',
        country_parameters = 'parameters/{country}/{plant_type}/country_parameters.xlsx'
    output:
        'results/hex_cost_components_{country}_{weather_year}_{plant_type}.geojson',
        'results/hex_cost_components_{country}_{weather_year}_{plant_type}.csv'
    script:
        'src/main/costs_by_component.py'

rule calculate_map_costs:
    input:
        hexagons = 'results/hex_cost_components_{country}_{weather_year}_{plant_type}.geojson',
        demand_parameters = 'parameters/{country}/{plant_type}/demand_parameters.xlsx'
    output:
        directory('plots/{country}_{weather_year}_{plant_type}')
    script:
        'src/main/map_costs.py'
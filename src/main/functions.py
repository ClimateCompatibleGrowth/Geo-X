"""
@authors: 
 - Claire Halloran
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk

This script contains functions needed by the following files:
 - transport_optimization.py
 - network.py
"""
import math
import numpy as np
import pandas as pd

def CRF(interest, lifetime):
    '''
    Calculates the capital recovery factor of a capital investment.

    Parameters
    ----------
    interest : float
        interest rate.
    lifetime : float or integer
        lifetime of asset.

    Returns
    -------
    CRF : float
        present value factor.
    '''
    interest = float(interest)
    lifetime = float(lifetime)

    CRF = (((1 + interest)**lifetime) * interest)/(((1 + interest)**lifetime) - 1)
    
    return CRF


def calculate_trucking_costs(transport_state, distance, quantity, interest, 
                             transport_params_filepath, currency):
    '''
    Calculates the annual cost of transporting the commodity by truck.

    Parameters
    ----------
    transport_state : string
        state commodity is transported in, one of '500 bar', 'LH2', 'LOHC', or
        'NH3'.
    distance : float
        distance between commodity production site and demand site.
    quantity : float
        annual amount of commodity to transport.
    interest : float
        interest rate on capital investments.
    transport_params_filepath : string
        path to transport_parameters.xlsx file.
    currency : string
        type of currency that is used in the parameter files.

    Returns
    -------
    annual_costs : float
        annual cost of commodity transport with specified method.
    '''
    daily_quantity = quantity/365

    transport_params = pd.read_excel(transport_params_filepath,
                                         sheet_name = transport_state,
                                         index_col = 'Parameter'
                                         ).squeeze('columns')

    average_truck_speed = transport_params['Average truck speed (km/h)']
    working_hours = transport_params['Working hours (h/day)']
    diesel_price = transport_params[f'Diesel price ({currency}/L)']
    costs_for_driver = transport_params[f'Costs for driver ({currency}/h)']
    working_days = transport_params['Working days (per year)']
    max_driving_dist = transport_params['Max driving distance (km/a)']

    spec_capex_truck = transport_params[f'Spec capex truck ({currency})']
    spec_opex_truck = transport_params['Spec opex truck (% of capex/a)']
    diesel_consumption = transport_params['Diesel consumption (L/100 km)']
    truck_lifetime = transport_params['Truck lifetime (a)']

    spec_capex_trailor = transport_params[f'Spec capex trailer ({currency})']
    spec_opex_trailor =transport_params['Spec opex trailer (% of capex/a)']
    net_capacity = transport_params['Net capacity (kg of commodity)']
    trailor_lifetime = transport_params['Trailer lifetime (a)']
    loading_unloading_time = transport_params['Loading unloading time (h)']

    # Calculate deliveries needed per day
    amount_deliveries_needed = daily_quantity/net_capacity
    # Calculate how many deliveries each truck can do per day
    deliveries_per_truck = working_hours/(loading_unloading_time +
                                          (2 * distance/average_truck_speed))
    
    # Deliveries per day / Deliveries per truck = Trucks per day
    # 0.5 is used to round up so full demand is met
    trailors_needed = round((amount_deliveries_needed/
                             deliveries_per_truck) + 0.5)
    total_drives_day = round(amount_deliveries_needed + 0.5)
    if transport_state == 'NH3':
        trucks_needed = trailors_needed
    else:
        trucks_needed = max(round((total_drives_day * 2 * distance *
                                        working_days/max_driving_dist) + 0.5),
                                            trailors_needed)
    # Get the capex of all the trucks and trailors needed
    capex_trucks = trucks_needed * spec_capex_truck
    capex_trailor = trailors_needed * spec_capex_trailor
    # Get fuel costs and wages
    if amount_deliveries_needed < 1:
        # In the lines below, 365 refers to  days to spread it over the year
        # and 100 is there because diesel_consumption is in liters/100km
        fuel_costs = (amount_deliveries_needed * 2 *
                        distance * 365/100) * diesel_consumption * diesel_price
        wages = amount_deliveries_needed * ((distance/average_truck_speed) *
                    2 + loading_unloading_time) * working_days * costs_for_driver
    else:
        fuel_costs = (round(amount_deliveries_needed + 0.5) *
                        2 * distance * 365/100) * diesel_consumption * diesel_price
        wages = round(amount_deliveries_needed + 0.5) * (
                    (distance/average_truck_speed) * 2 + loading_unloading_time
                    ) * working_days * costs_for_driver
    # Get total annual costs including capex, fuel costs, wages
    annual_costs = (capex_trucks * CRF(interest, truck_lifetime) + capex_trailor *
                        CRF(interest, trailor_lifetime)) +\
                            capex_trucks * spec_opex_truck + capex_trailor *\
                                spec_opex_trailor + fuel_costs + wages
    
    return annual_costs


def h2_conversion_stand(final_state, quantity, electricity_costs, heat_costs, 
                        interest, conversion_params_filepath, currency):
    '''
    Calculates the annual cost and electricity and heating demand for converting 
    hydrogen to a given state.

    Parameters
    ----------
    final_state : string
        final state to convert hydrogen to, one of 'standard condition', '500 bar',
        'LH2', 'LOHC_load', 'LOHC_unload', 'NH3_load', or 'NH3_unload'.
    quantity : float
        annual quantity of hydrogen to convert in kg.
    electricity_costs : float
        unit price for electricity.
    heat_costs : float
        unit costs for heat.
    interest : float
        interest rate applicable to hydrogen converter investments.
    conversion_params_filepath : string
        path to conversion parameters excel sheet.
    currency : string
        type of currency that is used in the parameter files.

    Returns
    -------
    elec_demand : float
        annual electricity demand.
    heat_demand : float
        annual heat demand.
    annual_costs : float
        annual hydrogen conversion costs.
    '''
    daily_throughput = quantity/365
    
    if final_state == 'standard condition':
        elec_demand = 0.0
        heat_demand = 0.0
        annual_costs = 0.0
    else:
        conversion_params = pd.read_excel(conversion_params_filepath,
                                             sheet_name = final_state,
                                             index_col = 'Parameter'
                                             ).squeeze('columns')
        
        if final_state == '500 bar':
            cp = conversion_params['Heat capacity']
            tein = conversion_params['Input temperature (K)']
            pein = conversion_params['Input pressure (bar)']
            k = conversion_params['Isentropic exponent']
            n_isentrop = conversion_params['Isentropic efficiency']
                        
            compressor_lifetime = conversion_params['Compressor lifetime (a)']
            capex_coef = conversion_params[f'Compressor capex coefficient ({currency} per kilograms H2 per day)']
            opex_compressor = conversion_params['Compressor opex (% capex)']
            elec_demand_per_kg_h2 = (cp * tein * (((500/pein)**((k - 1)/k)) - 1))/n_isentrop
            elec_demand = elec_demand_per_kg_h2 * quantity
            heat_demand = 0 

            capex_compressor = capex_coef * ((daily_throughput)**0.6038)

            annual_costs = (capex_compressor * CRF(interest, compressor_lifetime)) +\
                                (capex_compressor * opex_compressor) +\
                                    elec_demand * electricity_costs +\
                                        heat_demand * heat_costs
            
        elif final_state == 'LH2':
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            capex_quadratic_coef = conversion_params[f'Capex quadratic coefficient ({currency} (kg H2)-2)']
            capex_linear_coef = conversion_params[f'Capex linear coefficient ({currency} per kg H2)']
            capex_constant = conversion_params[f'Capex constant ({currency})']
            opex_liquid_plant = conversion_params['Opex (% of capex)']
            liquid_plant_lifetime = conversion_params['Plant lifetime (a)']
            
            heat_demand = 0
            elec_demand = electricity_unit_demand * quantity
            capex_liquid_plant = capex_quadratic_coef * (daily_throughput**2) +\
                                    capex_linear_coef * daily_throughput +\
                                        capex_constant

            annual_costs = (capex_liquid_plant * CRF(interest, liquid_plant_lifetime)) +\
                                (capex_liquid_plant * opex_liquid_plant) +\
                                    elec_demand * electricity_costs +\
                                        heat_demand * heat_costs
            
        elif final_state == 'LOHC_load':
            # In this conversion you "load" the hydrogen molecule to a carrier liquid.
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            heat_unit_demand = conversion_params['Heat demand (kWh per kg H2)']
            capex_coef = conversion_params[f'Capex coefficient ({currency} per kilograms H2 per year)']
            opex_hydrogenation = conversion_params['Opex (% of capex)']
            hydrogenation_lifetime = conversion_params['Hydrogenation lifetime (a)']
            costs_carrier = conversion_params[f'Carrier costs ({currency} per kg carrier)']
            ratio_carrier = conversion_params['Carrier ratio (kg carrier: kg hydrogen)']
            
            elec_demand = electricity_unit_demand * quantity 
            heat_demand = heat_unit_demand * quantity              
            capex_hydrogenation = capex_coef * quantity

            annual_costs = (capex_hydrogenation + costs_carrier *
                                ratio_carrier * daily_throughput) *\
                                    CRF(interest, hydrogenation_lifetime) +\
                                        capex_hydrogenation * opex_hydrogenation +\
                                            elec_demand * electricity_costs +\
                                                heat_demand * heat_costs
            
        elif final_state == 'LOHC_unload':
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            heat_unit_demand = conversion_params['Heat demand (kWh per kg H2)']
            capex_coef = conversion_params[f'Capex coefficient ({currency} per kilograms H2 per year)']
            opex_dehydrogenation = conversion_params['Opex (% of capex)']
            dehydrogenation_lifetime = conversion_params['Hydrogenation lifetime (a)']
            
            elec_demand = electricity_unit_demand * quantity 
            heat_demand = heat_unit_demand * quantity
            capex_dehydrogenation = capex_coef * quantity
            
            annual_costs = (capex_dehydrogenation *
                                CRF(interest, dehydrogenation_lifetime)) +\
                                    (capex_dehydrogenation * opex_dehydrogenation) +\
                                        elec_demand * electricity_costs +\
                                            heat_demand * heat_costs
            
        elif final_state == 'NH3_load':
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            heat_unit_demand = conversion_params['Heat demand (kWh per kg H2)']
            capex_coefficient = conversion_params[f'Capex coefficient ({currency} per annual g H2)']
            opex_NH3_plant = conversion_params['Opex (% of capex)']
            NH3_plant_lifetime = conversion_params['Plant lifetime (a)']
            
            
            elec_demand = electricity_unit_demand * quantity
            heat_demand = heat_unit_demand * quantity
            capex_NH3_plant = capex_coefficient * quantity

            annual_costs = capex_NH3_plant * CRF(interest, NH3_plant_lifetime) +\
                                capex_NH3_plant * opex_NH3_plant +\
                                    elec_demand * electricity_costs +\
                                        heat_demand * heat_costs
            
        elif final_state == 'NH3_unload':
            electricity_unit_demand = conversion_params['Electricity demand (kWh per kg H2)']
            heat_unit_demand = conversion_params['Heat demand (kWh per kg H2)']
            capex_coefficient = conversion_params[f'Capex coefficient ({currency} per hourly g H2)']
            opex_NH3_plant = conversion_params['Opex (% of capex)']
            NH3_plant_lifetime = conversion_params['Plant lifetime (a)']
            
            elec_demand = electricity_unit_demand * quantity
            heat_demand = heat_unit_demand * quantity

            capex_NH3_plant = capex_coefficient * ((quantity/1000/365/24) ** 0.7451)    

            annual_costs = capex_NH3_plant *\
                                CRF(interest, NH3_plant_lifetime) +\
                                    capex_NH3_plant * opex_NH3_plant +\
                                        elec_demand * electricity_costs +\
                                            heat_demand * heat_costs
        else:
            raise NotImplementedError(f'Conversion costs for {final_state} not currently supported.')
        
    return elec_demand, heat_demand, annual_costs


def cheapest_trucking_strategy(final_state, quantity, distance, 
                                elec_costs, heat_costs, interest,
                                conversion_params_filepath,
                                transport_params_filepath, currency):
    '''
    Calculates the lowest-cost state to transport hydrogen by truck.

    Parameters
    ----------
    final_state : string
        final state for hydrogen demand.
    quantity : float
        annual demand for hydrogen in kg.
    distance : float
        distance to transport hydrogen.
    elec_costs : float
        cost per kWh of electricity for that country.
    heat_costs : float
        cost per kWh of heat.
    interest : float
        interest on conversion and trucking capital investments 
        (not including roads).
    conversion_params_filepath : string
        path to conversion parameters excel sheet.
    transport_params_filepath : string
        path to transport parameters excel sheet. 
    currency : string
        type of currency that is used in the parameter files.
    
    Returns
    -------
    costs_per_unit : float
        storage, conversion, and transport costs for the cheapest trucking option.
    cheapest_option : string
        the lowest-cost state in which to transport hydrogen by truck.
    '''
    if final_state == '500 bar':
        dist_costs_500bar = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('500 bar', distance, quantity, interest, transport_params_filepath, currency)
    elif final_state == 'NH3':
        dist_costs_500bar = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('500 bar',distance,quantity,interest,transport_params_filepath, currency) +\
                    h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]
    else:  
        dist_costs_500bar = h2_conversion_stand('500 bar', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('500 bar',distance,quantity,interest,transport_params_filepath, currency) +\
                    h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]
    
    if final_state == 'LH2':
        dist_costs_lh2 = h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('LH2',distance, quantity,interest,transport_params_filepath, currency)
    elif final_state == 'NH3':
        dist_costs_lh2 = h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('LH2',distance,quantity,interest,transport_params_filepath, currency) +\
                    h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]
    else:
        dist_costs_lh2 = h2_conversion_stand('LH2', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('LH2',distance, quantity,interest,transport_params_filepath, currency) +\
                    h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]
    
    if final_state == 'NH3':
        dist_costs_nh3 = h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('NH3',distance, quantity, interest,transport_params_filepath, currency)
        dist_costs_lohc = h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('LOHC',distance, quantity, interest,transport_params_filepath, currency) +\
                    h2_conversion_stand('LOHC_unload', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                        h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]
    else:
        dist_costs_nh3 = h2_conversion_stand('NH3_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('NH3',distance, quantity,interest,transport_params_filepath, currency) +\
                    h2_conversion_stand('NH3_unload', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                        h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]
        dist_costs_lohc = h2_conversion_stand('LOHC_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                calculate_trucking_costs('LOHC',distance, quantity,interest,transport_params_filepath, currency) +\
                    h2_conversion_stand('LOHC_unload', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2] +\
                        h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]

    lowest_cost = np.nanmin([dist_costs_500bar, dist_costs_lh2, dist_costs_lohc, dist_costs_nh3])
    
    # Allocating cheapest option
    if dist_costs_500bar == lowest_cost:
        cheapest_option = '500 bar'
    elif dist_costs_lh2 == lowest_cost:
        cheapest_option = 'LH2'
    elif dist_costs_lohc == lowest_cost:
        cheapest_option = 'LOHC'
    elif dist_costs_nh3 == lowest_cost: 
         cheapest_option = 'NH3'
    
    costs_per_unit = lowest_cost/quantity
    
    return costs_per_unit, cheapest_option
   
    
def cheapest_pipeline_strategy(final_state, quantity, distance, 
                                elec_costs, heat_costs, interest, 
                                conversion_params_filepath,
                                pipeline_params_filepath,
                                currency,
                                elec_cost_grid = 0.0):
    '''
    Calculates the lowest-cost way to transport hydrogen via pipeline

    Parameters
    ----------
    final_state : string
        final state for hydrogen demand.
    quantity : float
        annual demand for hydrogen in kg.
    distance : float
        distance to transport hydrogen.
    elec_costs : float
        cost per kWh of electricity for that country.
    heat_costs : float
        cost per kWh of heat.
    interest : float
        interest on pipeline capital investments.
    conversion_params_filepath: string
        path to conversion parameters excel sheet.
    pipeline_params_filepath : string
        path to pipeline parameters excel sheet.
    currency : string
        type of currency that is used in the parameter files.
    elec_cost_grid : float
        grid electricity costs that pipeline compressors pay. Default 0.0

    Returns
    -------
    costs_per_unit : float
        storage, conversion, and transport costs for the cheapest option.
    cheapest_option : string
        the lowest-cost state in which to transport hydrogen by truck.
    '''
    if final_state == 'NH3':
        dist_costs_pipeline = pipeline_costs(distance,quantity, elec_cost_grid, pipeline_params_filepath, interest, currency)[0] +\
                h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]  
    else:
        dist_costs_pipeline = pipeline_costs(distance, quantity, elec_cost_grid, pipeline_params_filepath, interest, currency)[0] +\
                h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]

    costs_per_unit = dist_costs_pipeline/quantity
    cheapest_option = pipeline_costs(distance, quantity, elec_cost_grid, pipeline_params_filepath, interest, currency)[1] 

    return costs_per_unit, cheapest_option


def pipeline_costs(distance, quantity, elec_cost, pipeline_params_filepath, 
                   interest, currency):
    '''
    Calculates the annualized cost of building a pipeline.

    Parameters
    ----------
    distance : float
        distance from production site to demand site in km.
    quantity : float
        annual quantity of hydrogen demanded in kg.
    elec_cost : float
        price of electricity along pipeline.
    pipeline_params_filepath: string
        path to conversion parameters excel sheet.
    interest : float
        interest rate on capital investments.
    currency : string
        type of currency that is used in the parameter files.

    Returns
    -------
    annual_costs : float
        annual costs for pipeline.
     : string
        size of pipeline to build.
    '''
    all_parameters = pd.read_excel(pipeline_params_filepath,
                                   sheet_name='All',
                                    index_col = 'Parameter'
                                    ).squeeze('columns')
    opex = all_parameters['Opex (% of capex)']
    availability = all_parameters['Availability']
    lifetime_pipeline = all_parameters['Pipeline lifetime (a)']
    lifetime_compressors = all_parameters['Compressor lifetime (a)']
    electricity_demand = all_parameters['Electricity demand (kWh/kg*km)']
    large_max_capacity = all_parameters['Large pipeline max capacity (GW)']
    med_max_capacity = all_parameters['Medium pipeline max capacity (GW)']
    sml_max_capacity = all_parameters['Small pipeline max capcity (GW)']

    # 33.333 (kWh/kg) is the energy density of hydrogen
    # 8760 are hours in a year
    large_max_throughput = (((large_max_capacity * (10**6))/33.333)) * 8760 * availability
    med_max_throughput = (((med_max_capacity * (10**6))/33.333)) * 8760 * availability
    sml_max_throughput = (((sml_max_capacity * (10**6))/33.333)) * 8760 * availability

    if quantity <= sml_max_throughput:
        pipeline_type = 'Small'
        
    elif quantity > sml_max_throughput and quantity <= med_max_throughput:
        pipeline_type = 'Medium'
    
    elif quantity > med_max_throughput and quantity <= large_max_throughput:
        pipeline_type = 'Large'

    else:
        return np.nan,'No Pipeline big enough'
    
    pipeline_parameters = pd.read_excel(pipeline_params_filepath,
                                   sheet_name=pipeline_type,
                                    index_col = 'Parameter'
                                    ).squeeze('columns')
    capex_pipeline = pipeline_parameters[f'Pipeline capex ({currency})']
    capex_compressor = pipeline_parameters[f'Compressor capex ({currency})']
    
    capex_annual = ((capex_pipeline * distance) *
                        CRF(interest, lifetime_pipeline)) +\
                            ((capex_compressor * distance) *\
                                CRF(interest, lifetime_compressors))
    opex_annual = opex * (capex_pipeline + capex_compressor) * distance
    electricity_costs = electricity_demand * distance * quantity * elec_cost

    annual_costs = capex_annual + opex_annual + electricity_costs

    return annual_costs, f"{pipeline_type} Pipeline"


def calculate_nh3_pipeline_costs(distance, quantity, elec_cost, 
                                 pipeline_params_filepath, interest, currency):
    '''
    Calculates the annualized cost of building a pipeline.

    Parameters
    ----------
    distance : float
        distance from production site to demand site in km.
    quantity : float
        annual quantity of ammonia demanded in kg.
    elec_cost : float
        price of electricity along pipeline.
    pipeline_params_filepath: string
        path to conversion parameters excel sheet.
    interest : float
        interest rate on capital investments.
    currency : string
        type of currency that is used in the parameter files.

    Returns
    -------
    cost_per_unit : float
        annual costs for pipeline per kilogram of ammonia transported.
     : string
        size of pipeline to build
    '''
    # Converts kg to t
    quantity = quantity / 1000
    all_parameters = pd.read_excel(pipeline_params_filepath,
                                   sheet_name='All',
                                    index_col = 'Parameter'
                                    ).squeeze('columns')
    opex = all_parameters['Opex (% of capex)']
    availability = all_parameters['Availability']
    lifetime_pipeline = all_parameters['Pipeline lifetime (a)']
    electricity_demand = all_parameters['Electricity demand (kWh/kg*km)']
    # t to kg
    large_max_flow = all_parameters['Large pipeline max capacity (t NH3/a)']*availability
    large_min_flow = all_parameters['Large pipeline min capacity (t NH3/a)']*availability
    med_min_flow = all_parameters['Medium pipeline min capacity (t NH3/a)']*availability
    small_min_flow = all_parameters['Small pipeline min capcity (t NH3/a)']*availability

    # If demand is large enough, then split flow into multiple pipelines
    if quantity > large_max_flow:
        n_pipelines = math.ceil(quantity/large_max_flow)
        quantity_per_pipeline = quantity / n_pipelines
    else:
        n_pipelines = 1
        quantity_per_pipeline = quantity
        
    if quantity_per_pipeline >= small_min_flow and quantity_per_pipeline < med_min_flow:
        pipeline_type = 'Small'
    
    elif quantity_per_pipeline >= med_min_flow and quantity_per_pipeline < large_min_flow:
        pipeline_type = 'Medium'
    
    elif quantity_per_pipeline >= large_min_flow and quantity_per_pipeline <= large_max_flow:
        pipeline_type = 'Large'
    
    elif quantity_per_pipeline < small_min_flow:
        return np.nan,'Flow too small for pipeline'
    
    pipeline_parameters = pd.read_excel(pipeline_params_filepath,
                                   sheet_name=pipeline_type,
                                    index_col = 'Parameter'
                                    ).squeeze('columns')
    y_int = pipeline_parameters[f'Capex y-intercept ({currency}/t/yr/100km)']
    print(y_int)
    slope = pipeline_parameters[f'Capex flow coefficient ({currency}/t^2/yr^2/100km)']
    print(slope)
    capex_coeff = (y_int + slope*quantity_per_pipeline)
    print(capex_coeff)
    # Distance coefficients are per 100 km
    capex_annual = (n_pipelines*(capex_coeff*distance/100*quantity_per_pipeline)*CRF(interest,lifetime_pipeline))
    opex_annual = opex*n_pipelines*(capex_coeff*distance/100*quantity_per_pipeline)
    electricity_costs = electricity_demand * distance * quantity * elec_cost

    annual_costs = capex_annual + opex_annual + electricity_costs
    # Converts back to kg
    cost_per_unit = annual_costs/(quantity*1000)

    return cost_per_unit, f"{pipeline_type} Pipeline"
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
import geopy

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
    rounded_up_deliveries_needed = np.ceil(amount_deliveries_needed)
    # Calculate how many deliveries each truck can do per day
    deliveries_per_truck = working_hours/(loading_unloading_time +
                                          (2 * distance/average_truck_speed))
    
    # Deliveries per day / Deliveries per truck = Trucks per day
    trailors_needed = np.ceil(amount_deliveries_needed/
                             deliveries_per_truck)
    total_drives_day = rounded_up_deliveries_needed
    if transport_state == 'NH3':
        trucks_needed = trailors_needed
    else:
        trucks_needed = max(np.ceil(total_drives_day * 2 * distance *
                                        working_days/max_driving_dist),
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
        fuel_costs = (rounded_up_deliveries_needed * 2 * distance * 
                      365/100) * diesel_consumption * diesel_price
        wages = rounded_up_deliveries_needed * (
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
    pipeline_size : string
        size of pipeline to build.
    '''
    pipeline_cost, pipeline_size = pipeline_costs(distance, quantity, 
                                                   elec_cost_grid, 
                                                   pipeline_params_filepath, 
                                                   interest, currency)
    if final_state == 'NH3':
        dist_costs_pipeline = pipeline_cost +\
                h2_conversion_stand(final_state+'_load', quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]  
    else:
        dist_costs_pipeline = pipeline_cost +\
                h2_conversion_stand(final_state, quantity, elec_costs, heat_costs, interest, conversion_params_filepath, currency)[2]

    costs_per_unit = dist_costs_pipeline/quantity
    
    return costs_per_unit, pipeline_size


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
    slope = pipeline_parameters[f'Capex flow coefficient ({currency}/t^2/yr^2/100km)']
    capex_coeff = (y_int + slope*quantity_per_pipeline)

    # Distance coefficients are per 100 km
    capex_annual = (n_pipelines*(capex_coeff*distance/100*quantity_per_pipeline)*CRF(interest,lifetime_pipeline))
    opex_annual = opex*n_pipelines*(capex_coeff*distance/100*quantity_per_pipeline)
    electricity_costs = electricity_demand * distance * quantity * elec_cost

    annual_costs = capex_annual + opex_annual + electricity_costs
    # Converts back to kg
    cost_per_unit = annual_costs/(quantity*1000)

    return cost_per_unit, f"{pipeline_type} Pipeline"

def mineral_conversion_facility(final_state, quantity, interest,
                                conversion_params_filepath):
    """
    Calculate the CAPEX cost and OPEX cost of a new mineral processing facility
    as annual equivalent value (euros/kg)

    Parameters
    ----------
    final_state : TYPE
        DESCRIPTION.
    quantity : TYPE
        DESCRIPTION.
    interest : TYPE
        DESCRIPTION.
    
    Returns
    -------
    None.

    """
    conversion_params = pd.read_excel(conversion_params_filepath,
                                             sheet_name = final_state,
                                             index_col = 'Parameter'
                                             )
    plant_lifetime = conversion_params.loc['Plant lifetime (a)']
    capex_exp_coef_A = conversion_params.loc["CapEx exp coef A (euros/tonne)"]
    capex_exp_coef_B = conversion_params.loc["CapEx exp coef B (euros/tonne)"]
    capex_exp_coef_k = conversion_params.loc["CapEx exp coef k (1/tonne)"]
    
    # Proposed exponential function as used by Baptiste in cost curve analysis
    # tot_capex = (capex_exp_coef_A*np.exp(-capex_exp_coef_k*(quantity/1000)) + capex_exp_coef_B) * (quantity/1000) * plant_lifetime
    
    # the cost curve analysis showed a flat relationship between cost and facility capacities greater than 20,000 tonnes
    # for version 1 a flat rate for capex based on CapEx exp coef B is used until a better function/data can be identified
    
    # to include the affect of interest rate, a fixed 10% is assumed to have been used in analysis of coef_B
    tot_capex = (capex_exp_coef_B * (quantity / 1000)) / CRF(0.1, plant_lifetime)
    
    # equivalent annual levilised capex
    annual_capex = tot_capex * CRF(interest, plant_lifetime)
    capex_per_unit = annual_capex / quantity
        
    # opex cost - based on baptiste breakdown    
    opex_df = conversion_params.loc[["Opex Labor (euros/tonne)",
                                         "Opex Reagents (euros/tonne)",
                                         "Opex Corporate overhead (euros/tonne)",
                                         "Opex Other (euros/tonne)",
                                         "Opex Royalties (euros/tonne)"]]
    annual_opex = opex_df.sum() * (quantity/1000)
    opex_per_unit = annual_opex / quantity
    
    
    return capex_per_unit, opex_per_unit

def calculate_construction_costs(hexagon, infrastructure_interest_rate,
                                infrastructure_lifetime, long_capex,
                                short_capex, short_opex, long_opex,
                                construction_type):
    dist = hexagon[f'{construction_type}_dist']
    if dist==0:
        construction_costs = 0.
    elif dist<10:
        construction_costs = (dist
                                * short_capex
                                * CRF(infrastructure_interest_rate,
                                        infrastructure_lifetime)
                                + dist * short_opex)
    else:
        construction_costs = (dist
                                * long_capex 
                                * CRF(infrastructure_interest_rate,
                                        infrastructure_lifetime)
                                + dist * long_opex)
    return construction_costs
    
def geodesic_matrix(gdf1, gdf2):
    distances = np.empty((gdf1.shape[0], gdf2.shape[0]))
    # the use of centroid throws a UserWarning when a geographic CRS (e.g. EPSG=4326) is used.
    # ideally you should convert to a projected crs (ideally with equal area cylindrical?) then back again.
    for i, p1 in enumerate(gdf1.centroid):
        for j, p2 in enumerate(gdf2.centroid):
            distances[i,j] = geopy.distance.geodesic((p1.y, p1.x),
                                                        (p2.y, p2.x)).km
  
    return pd.DataFrame(distances, index=gdf1.index, columns=gdf2.index)

def find_nearest_hex(idx, hex_to_X_distance):
    hix = hex_to_X_distance.loc[:, idx].sort_values().index[0] # index of the nearest hexagon to the demand centre (contain might not work at low hexagon resolution)
    return hix

def determine_feedstock_sources(feedstock_points_gdf,
                                hexagon_to_feedstock_distance_matrix,
                                hix, feedstock_quantity,file):
    feedstock_ranked = feedstock_points_gdf.merge(hexagon_to_feedstock_distance_matrix.loc[hix,:], left_index=True, right_index=True).sort_values(by=hix)[['Annual capacity [kg/a]']]
    feedstock_ranked["Cumulative [kg/year]"] = feedstock_ranked.cumsum()
    feedstock_ranked["Feedstock used [kg/year]"] = 0.0
    remaining_quantity = feedstock_quantity
    for f, feedstock in feedstock_ranked.iterrows():
        if remaining_quantity <= 0:
            break
        feedstock_used = min(feedstock_points_gdf.loc[f, 'Annual capacity [kg/a]'], remaining_quantity)
        feedstock_ranked.loc[f,"Feedstock used [kg/year]"] = feedstock_used
        remaining_quantity -= feedstock_used
        
    
    feedstock = feedstock_ranked.loc[:feedstock_ranked[feedstock_ranked["Cumulative [kg/year]"] >= feedstock_quantity].index[0], :]
    file.write(f"FEEDSTOCK: {feedstock}\n\n")
    feedstock_ranked_idxs = feedstock.index
    return feedstock, feedstock_ranked_idxs

def calculate_train_costs(transport_state, distance, quantity, interest, 
                             transport_params_filepath, currency):
    '''
    Estimates the annual cost of transporting resource by train.
    Method:
        1. 
        2.
    Assumptions:
        1.
        2.

    Parameters
    ----------
    transport_state : string
        state resource is transported in.
    distance : float
        distance between production site and demand site.
    quantity : float
        annual amount of resource to transport.
    interest : float
        interest rate on capital investments.
    excel_path : string
        path to transport_parameters.xlsx file
        
    Returns
    -------
    annual_costs : float
        annual cost of hydrogen transport with specified method.
    '''
    transport_parameters = pd.read_excel(transport_params_filepath,
                                         sheet_name = transport_state,
                                         index_col = 'Parameter'
                                         ).squeeze('columns')

    average_train_speed = transport_parameters['Average train speed (km/h)']
    working_hours = transport_parameters['Working hours (h/day)']
    diesel_price = transport_parameters[f'Diesel price ({currency}/L)']
    costs_for_driver = transport_parameters[f'Costs for driver ({currency}/h)']
    working_days = transport_parameters['Working days (per year)']

    spec_capex_loco = transport_parameters[f'Spec capex loco ({currency})']
    spec_opex_loco = transport_parameters['Spec opex loco (% of capex/a)']
    diesel_consumption = transport_parameters['Diesel consumption (L/100 km)']
    loco_lifetime = transport_parameters['Loco lifetime (a)']

    spec_capex_wagon = transport_parameters[f'Spec capex wagon ({currency})']
    spec_opex_wagon =transport_parameters['Spec opex wagon (% of capex/a)']
    net_capacity = transport_parameters['Net capacity (kg of commodity)']
    wagon_lifetime = transport_parameters['Wagon lifetime (a)']
    loading_unloading_time = transport_parameters['Loading unloading time (h)']
    max_wagons = transport_parameters['Max wagons per loco']

    # Maximum journeys a single loco can make per year
    max_journeys_per_loco = (working_hours * working_days) / (loading_unloading_time + (2 * distance / average_train_speed)) 
    # Maximum quantity a single train (loco + max wagons) can transfer per year
    max_quantity_per_train = max_journeys_per_loco * net_capacity * max_wagons 
    # Minimum number of locos required assuming 1 loco per train
    min_locos = np.ceil(quantity/max_quantity_per_train)
    min_wagons_per_train = np.ceil((quantity/min_locos)/(net_capacity))  
    
    total_train_journeys = quantity / (min_wagons_per_train * net_capacity)

    capex_loco = min_locos * spec_capex_loco
    capex_wagon = min_wagons_per_train * spec_capex_wagon
    
    fuel_costs = ((total_train_journeys * 2 * distance) / 100) * diesel_consumption * diesel_price
    wages = total_train_journeys * ((distance / average_train_speed) * 2 + loading_unloading_time) * costs_for_driver

    annual_costs = ((capex_loco * CRF(interest, loco_lifetime) + capex_wagon * CRF(interest, wagon_lifetime))
                    + capex_loco * spec_opex_loco 
                    + capex_wagon * spec_opex_wagon 
                    + fuel_costs 
                    + wages)
    
    return annual_costs
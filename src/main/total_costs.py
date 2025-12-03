"""
@authors:
 - Claire Halloran
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk

Bring together all previous data to calculate lowest-cost of commodity.
"""
import geopandas as gpd
import numpy as np
import pandas as pd

from utils import check_folder_exists

def main():
    hexagons = gpd.read_file(snakemake.input.hexagons)
    len_hexagons = len(hexagons)
    demand_params_filepath = snakemake.input.demand_parameters
    demand_center_list = pd.read_excel(demand_params_filepath,
                                    index_col='Demand center')
    demand_centers = demand_center_list.index
    plant_type = snakemake.wildcards.plant_type
    currency = snakemake.config["currency"]

    check_folder_exists("results")

    # Get lowest cost for each transport type
    for demand_center in demand_centers:
        print(f"Calculating total costs for {demand_center} begins...\n")
        if plant_type == 'minx':
            # Calculating total grid LCOP
            hexagons[f'{demand_center} total grid LCOP'] = (hexagons[f'{demand_center}  grid total energy cost ({currency}/kg/year)']
                                    + hexagons[f'{demand_center} annual facility costs ({currency}/kg/year)']
                                    + hexagons[f'{demand_center} total road transport costs ({currency}/kg/year)']
                                    + hexagons[f'{demand_center} feedstock cost ({currency}/kg/year)']
                                    + hexagons[f'{demand_center} lowest water cost ({currency}/kg/year)'])
            
            # Calculating total hybrid LCOP
            hexagons[f'{demand_center} total hybrid LCOP'] = (hexagons[f'{demand_center} hybrid total energy cost ({currency}/kg/year)']
                                + hexagons[f'{demand_center} annual facility costs ({currency}/kg/year)']
                                + hexagons[f'{demand_center} total road transport costs ({currency}/kg/year)']
                                + hexagons[f'{demand_center} feedstock cost ({currency}/kg/year)']
                                + hexagons[f'{demand_center} lowest water cost ({currency}/kg/year)'])
        else:
            if plant_type == "hydrogen":
                trucking_tranport_costs = hexagons[f'{demand_center} trucking transport and conversion costs']
                pipeline_transport_costs = hexagons[f'{demand_center} pipeline transport and conversion costs']
            elif plant_type == "ammonia":
                trucking_tranport_costs = hexagons[f'{demand_center} trucking transport costs']
                pipeline_transport_costs = hexagons[f'{demand_center} pipeline transport costs']

            # Calculating trucking total cost
            hexagons[f'{demand_center} trucking total cost'] =\
                hexagons[f'{demand_center} road construction costs'] +\
                    trucking_tranport_costs +\
                        hexagons[f'{demand_center} production cost'] +\
                            hexagons['Lowest water cost']
            # Calculating trucking total cost
            hexagons[f'{demand_center} pipeline total cost'] =\
                    pipeline_transport_costs +\
                        hexagons[f'{demand_center} production cost'] +\
                            hexagons['Lowest water cost']

        for i in range(len_hexagons):
            if plant_type == 'minx':
                # Calculating total offgrid LCOP
                if hexagons[f'{demand_center} offgrid total energy cost ({currency}/kg/year)'][i] == np.nan:
                    hexagons.loc[i, f'{demand_center} total offgrid LCOP'] = np.nan
                else:
                    hexagons.loc[i, f'{demand_center} total offgrid LCOP'] = (hexagons.loc[i, f'{demand_center} offgrid total energy cost ({currency}/kg/year)']
                                            + hexagons.loc[i, f'{demand_center} annual facility costs ({currency}/kg/year)']
                                            + hexagons.loc[i, f'{demand_center} total road transport costs ({currency}/kg/year)']
                                            + hexagons.loc[i, f'{demand_center} feedstock cost ({currency}/kg/year)']
                                            + hexagons.loc[i, f'{demand_center} lowest water cost ({currency}/kg/year)'])
            else:
                # Get the lowest between the trucking and the pipeline options
                print(f"Determining lowest cost for {i+1} of {len_hexagons}...")
                hexagons.loc[i, f'{demand_center} lowest cost'] = np.nanmin(
                                        [hexagons.loc[i, f'{demand_center} trucking total cost'],
                                        hexagons.loc[i, f'{demand_center} pipeline total cost']])
            
    print("\nCalculations complete.\n")
    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')

if __name__ == "__main__":
    main()
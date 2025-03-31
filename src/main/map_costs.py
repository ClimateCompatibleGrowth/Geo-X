"""
@authors: 
 - Claire Halloran
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk

This script visualizes the spatial cost of the commodity for each demand center.
"""
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

from utils import check_folder_exists

def plot_and_save(crs, hexagons, name, legend_kwds, output_folder, figsize=(10,5), legend=True, cmap='viridis_r', 
             missing_kwds={"color": "lightgrey", "label": "Missing values",},
             bbox_inches="tight",
             ):
    """
    Plots into a figure and saves figure into a file.

    Parameters
    ----------
    crs : cartopy Orthographic Object
        determines the coordinate reference systems (crs) of a given 
        central latitude and longitude.
    hexagons : geodataframe
        hexagon GeoJSON file.
    name : string
        used for column title, axes title and file naming.
    legend_kwds : dictionary
        keyword arguments to pass to matplotlib.pyplot.legend() - 'label' will
        overwrite the auto-generated label.
    output_folder : string
        path to the output folder.
    figsize : tuple
        size of figure. Default is (10,5).
    legend : boolean
        whether to plot a legend or not. Default is True.
    cmap : string
        the name of a colormap. Default is 'viridis_r'.
    missing_kwds : dictionary
        keyword arguments specifying the colour for missing values. If None is 
        specified, missing values are not plotted. Default is 
        {"color": "lightgrey", "label": "Missing values",}.
    bbox_inches : string
        bounding box in inches - will only save the given portion of the 
        figure. Default is "tight".
    """
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    # Only plots if there is at least one non-null value
    if hexagons[name].isnull().all()==False:
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column = name,
            legend = legend,
            cmap = cmap,
            legend_kwds = legend_kwds,
            missing_kwds = missing_kwds,
        )

    ax.set_title(name)
    fig.savefig(output_folder + f"/{name}.png", bbox_inches=bbox_inches)
    plt.close()

if __name__ == "__main__":
    currency = snakemake.config["currency"]
    plant_type = str(snakemake.wildcards.plant_type)
    hexagons = gpd.read_file(str(snakemake.input.hexagons))
    demand_excel_path = str(snakemake.input.demand_parameters)
    demand_parameters = pd.read_excel(demand_excel_path,index_col='Demand center')
    demand_centers = demand_parameters.index
    transport_methods = ["trucking", "pipeline"]

    # Update central coordinates for area considered
    hexagon_bounds = hexagons.geometry.bounds    
    min_lon, min_lat = hexagon_bounds[['minx','miny']].min()
    max_lon, max_lat = hexagon_bounds[['maxx','maxy']].max()

    central_lon = (min_lon + max_lon)/2
    central_lat = (min_lat + max_lat)/2

    crs = ccrs.Orthographic(central_longitude = central_lon, central_latitude= central_lat)
    generators = dict(snakemake.config['generators_dict'])

    output_folder = str(snakemake.output)
    check_folder_exists(output_folder)

    for demand_center in demand_centers:
        print(f"\nPlotting for {demand_center} begins...")

        # Lowest LC in each location
        plot_and_save(crs, hexagons, f'{demand_center} lowest cost', 
                    {'label' : f'LC [{currency}/kg]'}, output_folder)

        for transport_method in transport_methods:
            # Trucking/ pipeline production cost
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} production cost',
                        {'label' : f'Production LC [{currency}/kg]'}, output_folder)

            # Trucking/ pipeline costs
            if plant_type == "hydrogen":
                if transport_method == "trucking":
                    hexagons[f'{demand_center} total {transport_method} costs'] =\
                        hexagons[f'{demand_center} {transport_method} transport and conversion costs']+\
                            hexagons[f'{demand_center} road construction costs']
                    
                    plot_and_save(crs, hexagons, f'{demand_center} total {transport_method} costs', 
                        {'label' : f'{transport_method} cost [{currency}/kg]'}, output_folder)
                elif transport_method == "pipeline":
                    plot_and_save(crs, hexagons, f'{demand_center} {transport_method} transport and conversion costs', 
                                {'label' : f'{transport_method} cost [{currency}/kg]'}, output_folder)
            elif plant_type == "ammonia":
                if transport_method == "trucking":
                    hexagons[f'{demand_center} total {transport_method} costs'] =\
                        hexagons[f'{demand_center} {transport_method} transport costs']+\
                            hexagons[f'{demand_center} road construction costs']
            
                    plot_and_save(crs, hexagons, f'{demand_center} total {transport_method} costs', 
                                {'label' : f'{transport_method} cost [{currency}/kg]'}, output_folder)
                elif transport_method == "pipeline":
                    plot_and_save(crs, hexagons, f'{demand_center} {transport_method} transport costs', 
                                {'label': f'{transport_method} cost [{currency}/kg]'}, output_folder)

            # Total cost
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} total cost',
                        {'label': f'LC [{currency}/kg]'}, output_folder)

            # Electrolyzer capacity
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} electrolyzer capacity',
                        {'label': 'Size (MW)'}, output_folder)
        
            # Electrolyzer costs
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} electrolyzer costs',
                        {'label': f'{currency}'}, output_folder)

            # H2 storage capacity
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} H2 storage capacity',
                        {'label': 'Size (MW)'}, output_folder)
        
            # H2 storage costs
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} H2 storage costs',
                        {'label': f'{currency}'}, output_folder)
            
            # Battery capcity
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} battery capacity', 
                        {'label': 'Size (MW)'}, output_folder)
        
            # Battery costs
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} battery costs',
                        {'label': f'{currency}'}, output_folder)

            for generator in generators:
                # Generator capacity
                plot_and_save(crs, hexagons, f'{demand_center} {transport_method} {generator} capacity',
                            {'label': 'Capacity (MW)'}, output_folder)
        
                # Generator costs
                plot_and_save(crs, hexagons, f'{demand_center} {transport_method} {generator} costs',
                            {'label': f'{currency}'}, output_folder)
        
    # Ocean water costs
    plot_and_save(crs, hexagons, 'Ocean water costs',
                {'label': f'Water cost [{currency}/kg H2]'}, output_folder)

    plt.ticklabel_format(style='plain')

    # Freshwater costs
    plot_and_save(crs, hexagons, 'Freshwater costs',
                {'label': f'Water cost [{currency}/kg H2]'}, output_folder)

    plt.ticklabel_format(style='plain')

    print("\nPlotting complete\n")
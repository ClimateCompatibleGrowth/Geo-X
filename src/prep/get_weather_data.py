"""
@authors: 
 - Claire Halloran
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk
 - Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk

This script fetches historical weather data to calculate wind and solar potential
from the ERA-5 reanalysis dataset using Atlite.

Create cutouts with `atlite <https://atlite.readthedocs.io/en/latest/>`_.
"""
import atlite
import geopandas as gpd
import logging

from utils import check_folder_exists

def calculate_coords(hexagons):
    """
    Calculates mininum and maximum coordinates using bounds from the hexagon file.

    ...
    Parameters
    ----------
    hexagons : geodataframe
        hexagon file from data folder.

    Returns
    -------
    min_lon : integer
        mininum longitude.
    min_lat : integer
        minimum latitude.
    max_lon : integer
        maximum longitude.
    max_lat : integer
        maximum latitude.
    """
    hexagon_bounds = hexagons.geometry.bounds
    min_lon, min_lat = hexagon_bounds[['minx','miny']].min()
    max_lon, max_lat = hexagon_bounds[['maxx','maxy']].max()
    
    return min_lon, min_lat, max_lon, max_lat

def prepare_cutout(min_lon, min_lat, max_lon, max_lat, start_date, end_date):
    """
    Creates and prepares the cutout, using Atlite.

    ...
    Parameters
    ----------
    min_lon : integer
        mininum longitude.
    min_lat : integer
        minimum latitude.
    max_lon : integer
        maximum longitude.
    max_lat : integer
        maximum latitude.
    start_date : string
        start date for weather collection in 'YYYY-MM-DD' format.
    end_date : string
        end date for weather collection in 'YYYY-MM-DD' format.
    """
    cutout = atlite.Cutout(
        path=str(snakemake.output),
        module="era5",
        x=slice(min_lon, max_lon),
        y=slice(min_lat, max_lat),
        time=slice(start_date, end_date),
    )
    
    # Without monthly_requests=True, CDS can reject requests due to being 
    # "too large"
    cutout.prepare(tmpdir="temp", monthly_requests=True)

def main():
    hexagons = gpd.read_file(str(snakemake.input.hexagons))
    # Displays information on process as it runs.
    logging.basicConfig(level=logging.INFO)

    min_lon, min_lat, max_lon, max_lat = calculate_coords(hexagons)

    start_weather_year = int(snakemake.wildcards.weather_year)
    end_weather_year = int(snakemake.wildcards.weather_year)+int(snakemake.config["years_to_check"])
    start_date = f'{start_weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'

    check_folder_exists("cutouts")
    check_folder_exists("temp")

    prepare_cutout(min_lon, min_lat, max_lon, max_lat, start_date, end_date)

if __name__ == "__main__":
    main()
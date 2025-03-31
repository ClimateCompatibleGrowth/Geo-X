"""
@authors:
 - Samiyha Naqvi, University of Oxford, samiyha.naqvi@eng.ox.ac.uk

This script contains functions needed by the following files:
 - get_weather_data.py
"""
from os import makedirs
from os.path import exists

def check_folder_exists(name):
    """
    Create folders if they do not already exist.

    ...
    Parameters
    ----------
    name : string
        File name to be checked.
    """
    if not exists(name):
        makedirs(name)
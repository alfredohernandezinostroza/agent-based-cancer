import pandas as pd
import ast

def load_simulation_configs(path):
    """
    Loads the config file, reading the csv given in path, and adding their values to the global scope
    Returns the names of the values in a list
    Input:
        string: path
        return: array of all the names of the variables
    """ 
    df_configs = pd.read_csv(path, header=0, converters={"Values": ast.literal_eval})
    # if the configs were loaded from a previous 
    # simulation, drop the last three rows
    # they shoud not be loaded as they change for each simulation according to user input
    if path != "simulations_configs.csv":
        df_configs = df_configs[:-3]
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    globals().update(dict_configs)  
    return(list(df_configs["Names"]))
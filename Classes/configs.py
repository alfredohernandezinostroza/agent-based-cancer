import pandas as pd
import ast

def init_simulation_configs(path):
    """
    Loads the config file, reading the csv given in path, and adding their values to the global scope
    Returns the names of the values in a list
    Input:
        string: path
        return: array of all the names of the variables
    """ 
    df_configs = pd.read_csv(path, header=0, converters={"Values": ast.literal_eval})
    if len(df_configs) < 29:
        raise Exception("Less than 29 configuration options! Are there some missing?")
    if len(df_configs) > 29:
        raise Exception("More than 29 configuration options!")
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    globals().update(dict_configs)  
    return(list(df_configs["Names"]))

def load_simulation_configs_for_data_generation(path):
    """
    Loads the config file, reading the csv given in path, and adding their values to the global scope
    Returns the names of the values in a list
    Input:
        string: path
        return: array of all the names of the variables
    """ 
    df_configs = pd.read_csv(path, header=0, converters={"Values": ast.literal_eval})
    if len(df_configs) < 32:
        raise Exception("Less than 32 configuration options! Are there some missing?")
    if len(df_configs) > 32:
        raise Exception("More than 32 configuration options!")
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    globals().update(dict_configs)  
    return(list(df_configs["Names"]))

def load_simulation_configs_for_reloaded_simulation(path):
    """
    Loads the config file, reading the csv given in path, and adding their values to the global scope
    Returns the names of the values in a list
    Input:
        string: path
        return: array of all the names of the variables
    """ 
    df_configs = pd.read_csv(path, header=0, converters={"Values": ast.literal_eval})
    # if the configs were loaded to continue from a previous 
    # simulation, drop the last three rows
    # they shoud not be loaded as they change for each simulation according to user input
    # maxSteps, dataCollectionPeriod and grids_number)
    df_configs = df_configs[:-3]
    if len(df_configs) < 29:
        raise Exception("Less than 29 configuration options! Are there some missing?")
    if len(df_configs) > 29:
        raise Exception("More than 29 configuration options!")
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    globals().update(dict_configs)  
    return(list(df_configs["Names"]))
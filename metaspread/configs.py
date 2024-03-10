import pandas as pd
import os
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
    if len(df_configs) < 31:
        raise Exception("Less than 31 configuration options! Are there some missing?")
    if len(df_configs) > 31:
        raise Exception("More than 31 configuration options!")
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    globals().update(dict_configs)  
    error_string = ""
    if sum(extravasation_probs) != 1:
        error_string += "Extravasation probabilities must sum 1!\n"
    if len(extravasation_probs) != grids_number - 1:
        error_string += "There must be as many Extravasation probabilities as the value of (grids_number - 1)!\n"
    if len(secondary_sites_vessels) != grids_number-1:
        error_string += "There must be as many secondary site vessels as the value of grids_number - 1!\n"
    if mesenchymal_proportion + epithelial_proportion != 1:
        error_string += "Mesenchymal_proportion + epithelial_proportion must be 1!\n"
    if error_string != "":
        raise ValueError(error_string)
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
    if len(df_configs) < 33:
        raise Exception("Less than 33 configuration options! Are there some missing?")
    if len(df_configs) > 33:
        raise Exception("More than 33 configuration options!")
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
    # simulation, drop the last two rows
    # they shoud not be loaded as they change for each simulation according to user input
    # (max steps and step size)
    df_configs = df_configs[:-2]
    if len(df_configs) < 31:
        raise Exception("Less than 33 configuration options! Are there some missing?")
    if len(df_configs) > 31:
        raise Exception("More than 33 configuration options!")
    dict_configs = dict(zip(df_configs["Names"], df_configs["Values"]))
    globals().update(dict_configs)  
    return(list(df_configs["Names"]))

def generate_default_configs():
    """Creates a default simulation_configs.csv file"""
    names = ["th","tha","xh","xha","dM","dE","phiM","phiE","dmmp","theta","Lambda","gamma1","gamma2","vasculature_time","doubling_time_M","doubling_time_E","single_cell_survival","cluster_survival","extravasation_probs","dissagreggation_prob","carrying_capacity","normal_vessels_primary","ruptured_vessels_primary","secondary_sites_vessels","n_center_points_for_tumor","n_center_points_for_Vessels","gridsize","grids_number","mesenchymal_proportion","epithelial_proportion","number_of_initial_cells"]
    values = [0.001,0.001,0.005,0.005,1e-4,5e-4,0.0005,0.0005,0.001,0.195,0.1,1,1,180,2000,3000,5e-04,0.025,[0.75, 0.25],0.5,4,8,2,[10, 10],97,200,201,3,0.6,0.4,388]
    default_configs = pd.DataFrame({"Names": names, "Values": values})
    default_configs.to_csv("simulations_configs.csv", index=False)

def check_if_configs_are_present():
    """Checks if a configs file is present in the current directory"""
    if not os.path.isfile("simulations_configs.csv"):
        print("Configs file simulations_configs.csv not found! Creating a configs file with default values.")
        generate_default_configs()
    else:
        print("Configs file found!")
import shutil
import re
import ast
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# from Classes.CancerModel import *
import mesa
import os
# from Classes.configs import *
import Classes.configs

# To run this code you must be in the parent folder of the program

def main_Batch(maxSteps, dataCollectionPeriod, loadedSimulationPath=""):

    # load configs file from a previous simulation or loads the general configs file
    print(loadedSimulationPath)
    loadedSimulationPath= loadedSimulationPath.strip('\"')
    if loadedSimulationPath != "":
        configs_path = os.path.join(loadedSimulationPath, "configs.csv")
        # config_var_names = Classes.configs.load_simulation_configs_for_data_generation(configs_path)
        config_var_names = Classes.configs.load_simulation_configs_for_reloaded_simulation(configs_path)
    else:
        configs_path = "simulations_configs.csv"
        config_var_names = Classes.configs.init_simulation_configs(configs_path)
    
    # Parameters for this simulation
    N = Classes.configs.number_of_initial_cells # Number of cancer cells
    gridsize     = Classes.configs.gridsize #PDF: 201
    grids_number = Classes.configs.grids_number #PDF: 201
    width        = gridsize
    height       = gridsize

    # Name of the directories
    simulations_dir = "Simulations"
    os.makedirs(simulations_dir, exist_ok=True)
    if loadedSimulationPath != "":
        new_simulation_folder = loadedSimulationPath.split('\\')[-1]
        new_simulation_folder = "Continuing_" + new_simulation_folder
        new_simulation_path = os.path.join(simulations_dir, new_simulation_folder)
        shutil.copytree(loadedSimulationPath, new_simulation_path)
        cells_path = os.path.join(new_simulation_path, "CellsData.csv")
        df = pd.read_csv(cells_path, index_col = 0, converters={"Position": ast.literal_eval})
        loaded_max_step = max(df["Step"])
        # new_simulation_folder = f"Sim maxSteps-{loaded_max_step}+{maxSteps} stepsize-{dataCollectionPeriod} N-{N} gridsNumber-{grids_number}"
    else:
        df = pd.DataFrame()
        loaded_max_step = 0
        new_simulation_folder = f"Sim_maxSteps-{maxSteps}_stepsize-{dataCollectionPeriod}_N-{N}_gridsNumber-{grids_number}"

        # Creates the path for the new simulation
        path = os.path.join(simulations_dir, new_simulation_folder)
        pathMmp2 = os.path.join(path, "Mmp2")
        pathEcm = os.path.join(path, "Ecm")
        pathVasculature = os.path.join(path, "Vasculature")
        pathTimeOfPopulation = os.path.join(path, "Time when grids were populated")

        # Create folder for all cells analysis, for Mmp2 matrices and Ecm matrices
        if not os.path.exists(path):
            print(f'\t Folder for this simulation: {path}')
            print(f'\t Saving agents data at: {path}')
            print(f'\t Saving Mmp2 data at: {pathMmp2}')
            print(f'\t Saving Ecm data at: {pathEcm}')
            print(f'\t Saving Vasculature data at: {pathVasculature}')

            os.makedirs(path)
            os.makedirs(pathMmp2)
            os.makedirs(pathEcm)
            os.makedirs(pathVasculature)
            os.makedirs(pathTimeOfPopulation)
        # If there is already a simulation you skip it
        else:
            return print("This simulation already exists!")

        # Run the simulation and saves the data
    run_simulation(Classes.CancerModel, N, width, height, grids_number, maxSteps, loaded_max_step, dataCollectionPeriod, new_simulation_folder, simulations_dir, config_var_names, df, loadedSimulationPath)


    return print('Finished the simulation')


def run_simulation(CancerModel, N, width, height, grids_number, maxSteps, loaded_max_step, dataCollectionPeriod, new_simulation_folder, simulations_dir, config_var_names, loaded_df, loadedSimulationPath=""):

    # Setting parameters for mesa.batch_run
    params = {"N": N, "width": width, "height": height, "grids_number": grids_number, "maxSteps": maxSteps, "dataCollectionPeriod": dataCollectionPeriod, "newSimulationFolder": new_simulation_folder }
    if loadedSimulationPath != "":
        params["loadedSimulationPath"] = loadedSimulationPath

    # Saves the simulation configuration
    print(f"Saving all the simulations parameters at: {os.path.join(simulations_dir, new_simulation_folder, 'configs.csv')}")
    values = [getattr(Classes.configs, i) for i in config_var_names]
    names = config_var_names

    #add configurations that are not in the global variables
    # if loadedSimulationPath == "":
    #     names += ['maxSteps', 'dataCollectionPeriod', 'grids_number']
    #     values += [maxSteps, dataCollectionPeriod, grids_number]
    names += ['maxSteps', 'dataCollectionPeriod']
    values += [maxSteps, dataCollectionPeriod]
    df_vars = pd.DataFrame({"Names": names, "Values": values})
    df_vars = df_vars.set_index("Names")
    path = os.path.join(simulations_dir, new_simulation_folder, 'configs.csv')
    df_vars.to_csv(path)

    # The whole simulation occurs here
    print(f'\n\tStarting the simulation')
    results = mesa.batch_run(
        CancerModel,
        parameters=params,
        iterations=1, # for the code to run and save in the specified folders this must be always 1
        max_steps=maxSteps,
        number_processes=None,
        data_collection_period=dataCollectionPeriod,
        display_progress=False,
    )

    # Create data frames for the cells
    cells_df = pd.DataFrame(results)
    cells_df = cells_df[["AgentID","Step", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]
    cells_df = cells_df.iloc[1:]
    if loaded_max_step != 0:
        cells_df["Step"] = cells_df["Step"] + loaded_max_step
        cells_df = pd.concat([loaded_df, cells_df])
    print(f'Example of data collected:\n{cells_df.head(10)}')

    # Saves data analysis
    nameOfCsv = f'CellsData.csv'
    pathToSave = os.path.join(simulations_dir, new_simulation_folder, nameOfCsv)
    cells_df.to_csv(pathToSave)
    print(f'All data saved')

if __name__ == "__main__":
    maxSteps = 100
    dataCollectionPeriod = 10
    # This runs all the code to create folders, run the simulation and save the data
    main_Batch(maxSteps, dataCollectionPeriod)
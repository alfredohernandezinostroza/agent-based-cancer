import re
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
    if loadedSimulationPath != "":
        configs_path = os.path.join(loadedSimulationPath, "configs.csv")
    else:
        configs_path = "simulations_configs.csv"
    config_var_names = Classes.configs.init_simulation_configs(configs_path)
    # for var in config_var_names:
    #     globals()[var] = getattr(Classes.configs, var)
    # from Classes.CancerModel import CancerModel
    
    # Parameters for this simulation
    N = 388 # Number of cancer cells
    gridsize     = Classes.configs.gridsize_utils #PDF: 201
    width        = gridsize
    height       = gridsize
    grids_number = 3

    # Name of the directories
    simulations_dir = "Simulations"
    if loadedSimulationPath != "":
        pattern = r'maxSteps-(\d+)'
        match = re.search(pattern, loadedSimulationPath)
        if match:
            loadedMaxSteps = int(match.group(1))
        else:
            raise ValueError("Error finding loaded simulation maxSteps!")
        newSimulationFolder = f"Sim maxStepz-{loadedMaxSteps}+{maxSteps} stepsize-{dataCollectionPeriod} N-{N} gridsNumber-{grids_number}"
    else:
        newSimulationFolder = f"Sim maxSteps-{maxSteps} stepsize-{dataCollectionPeriod} N-{N} gridsNumber-{grids_number}"

    # Create parent directory if it doesn't exist
    if not os.path.exists(simulations_dir):
        print(f'Creating the folder for this simulation at: {simulations_dir}')
        os.makedirs(simulations_dir)

    # Creates the path for the new simulation
    path = os.path.join(simulations_dir, newSimulationFolder)
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

        # Run the simulation and saves the data
        run_simulation(Classes.CancerModel, N, width, height, grids_number, maxSteps, dataCollectionPeriod, newSimulationFolder, simulations_dir, config_var_names, loadedSimulationPath)

    # If there is already a simulation you skip it
    else:
        return print("This simulation already exists!")

    return print('Finished the simulation')


def run_simulation(CancerModel, N, width, height, grids_number, maxSteps, dataCollectionPeriod, newSimulationFolder, simulations_dir, config_var_names, loadedSimulationPath=""):

    # Setting parameters for mesa.batch_run
    params = {"N": N, "width": width, "height": height, "grids_number": grids_number, "maxSteps": maxSteps, "dataCollectionPeriod": dataCollectionPeriod, "newSimulationFolder": newSimulationFolder }
    if loadedSimulationPath != "":
        params["loadedSimulationPath"] = loadedSimulationPath

    # Saves the simulation configuration
    print(f"Saving all the simulations parameters at: {os.path.join(simulations_dir, newSimulationFolder, 'configs.csv')}")
    values = [getattr(Classes.configs, i) for i in config_var_names]
    #add configurations that are not in the global variables
    names = config_var_names + ['maxSteps', 'dataCollectionPeriod', 'grids_number']
    values += [maxSteps, dataCollectionPeriod, grids_number]
    df_vars = pd.DataFrame({"Names": names, "Values": values})
    df_vars = df_vars.set_index("Names")
    path = os.path.join(simulations_dir, newSimulationFolder, 'configs.csv')
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
        display_progress=True,
    )

    # Create data frames for the cells
    cells_df = pd.DataFrame(results)
    cells_df = cells_df[["Step", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]
    print(f'Example of data collected: {cells_df.head(10)}')

    # Saves data analysis
    nameOfCsv = f'CellsData.csv'
    pathToSave = os.path.join(simulations_dir, newSimulationFolder, nameOfCsv)
    cells_df[1:].to_csv(pathToSave)
    print(f'All data saved')

if __name__ == "__main__":
    maxSteps = 100
    dataCollectionPeriod = 10
    # This runs all the code to create folders, run the simulation and save the data
    main_Batch(maxSteps, dataCollectionPeriod)
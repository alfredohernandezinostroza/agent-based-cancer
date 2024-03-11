import metaspread.configs
import pandas as pd
import shutil
import mesa
import ast
import os

# To run this code you must be in the parent folder of the program

def main(max_steps, data_collection_period, loaded_simulation_path=""):

    # load configs file from a previous simulation or loads the general configs file
    print(loaded_simulation_path)
    loaded_simulation_path= loaded_simulation_path.strip('\"')
    if loaded_simulation_path != "":
        configs_path = os.path.join(loaded_simulation_path, "configs.csv")
        config_var_names = metaspread.configs.load_simulation_configs_for_reloaded_simulation(configs_path)
    else:
        configs_path = "simulations_configs.csv"
        config_var_names = metaspread.configs.init_simulation_configs(configs_path)
    
    # Parameters for this simulation
    number_of_initial_cells = metaspread.configs.number_of_initial_cells # Number of cancer cells
    gridsize     = metaspread.configs.gridsize #PDF: 201
    grids_number = metaspread.configs.grids_number #PDF: 201
    width        = gridsize
    height       = gridsize

    # Name of the directories
    simulations_dir = "Simulations"
    os.makedirs(simulations_dir, exist_ok=True)
    if loaded_simulation_path != "":
        cells_path = os.path.join(loaded_simulation_path, "CellsData.csv")
        df = pd.read_csv(cells_path)
        loaded_max_step = max(df["Step"])
        new_simulation_folder = os.path.normpath(loaded_simulation_path)
        new_simulation_folder = os.path.basename(new_simulation_folder)
        new_simulation_folder = f"Continue-from-{loaded_max_step}-until-{loaded_max_step+max_steps}-" + new_simulation_folder
        new_simulation_path = os.path.join(simulations_dir, new_simulation_folder)
        print(f"Copying simulation in {loaded_simulation_path} to {new_simulation_path}")
        shutil.copytree(loaded_simulation_path, new_simulation_path)
    else:
        df = pd.DataFrame()
        loaded_max_step = 0
        new_simulation_folder = f"Sim-max_steps-{max_steps}-collection_period-{data_collection_period}-cells-{number_of_initial_cells}-grids_number-{grids_number}"

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
    run_simulation(metaspread.CancerModel, number_of_initial_cells, width, height, grids_number, max_steps, loaded_max_step, data_collection_period, new_simulation_folder, simulations_dir, config_var_names, df, loaded_simulation_path)


    return print('Finished the simulation')


def run_simulation(CancerModel, number_of_initial_cells, width, height, grids_number, max_steps, loaded_max_step, data_collection_period, new_simulation_folder, simulations_dir, config_var_names, loaded_df, loaded_simulation_path=""):

    # Setting parameters for mesa.batch_run
    params = {"number_of_initial_cells": number_of_initial_cells, "width": width, "height": height, "grids_number": grids_number, "max_steps": max_steps, "data_collection_period": data_collection_period, "new_simulation_folder": new_simulation_folder }
    if loaded_simulation_path != "":
        params["loaded_simulation_path"] = loaded_simulation_path

    # Saves the simulation configuration
    print(f"\t Saving all the simulations parameters at: {os.path.join(simulations_dir, new_simulation_folder, 'configs.csv')}")
    values = [getattr(metaspread.configs, i) for i in config_var_names]
    names = config_var_names

    #add configurations that are not in the global variables
    names += ['max_steps', 'data_collection_period']
    values += [max_steps, data_collection_period]
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
        max_steps=max_steps,
        number_processes=None,
        data_collection_period=data_collection_period,
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
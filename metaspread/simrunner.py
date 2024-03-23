import metaspread.configs
import pandas as pd
import shutil
import mesa
import ast
import os

# To run this code you must be in the parent folder of the program

def save_configs(simulations_dir, new_simulation_folder, config_var_names, max_steps, data_collection_period):
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


def run_simulation(max_steps, data_collection_period, loaded_simulation_path=""):

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
        # new_simulation_folder = f"Continue-from-{loaded_max_step}-until-{loaded_max_step+max_steps}-" + new_simulation_folder
        new_simulation_path = os.path.join(simulations_dir, new_simulation_folder)
        # print(f"Copying simulation in {loaded_simulation_path} to {new_simulation_path}")
        # shutil.copytree(os.path.join(loaded_simulation_path,"Ecm"),os.path.join(new_simulation_path,"Ecm"))
        # shutil.copytree(os.path.join(loaded_simulation_path,"Mmp2"),os.path.join(new_simulation_path,"Mmp2"))
        # shutil.copytree(os.path.join(loaded_simulation_path,"Vasculature"),os.path.join(new_simulation_path,"Vasculature"))
        # shutil.copytree(os.path.join(loaded_simulation_path,"Time when grids were populated"),os.path.join(new_simulation_path,"Time when grids were populated"))
        # shutil.copy2(os.path.join(loaded_simulation_path,"CellsData.csv"),os.path.join(new_simulation_path,"CellsData.csv"))
        # shutil.copy2(os.path.join(loaded_simulation_path,"CellsData.csv"),os.path.join(new_simulation_path,"configs.csv"))
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
    save_configs(simulations_dir, new_simulation_folder, config_var_names, max_steps, data_collection_period)
    model = metaspread.CancerModel(number_of_initial_cells, width, height, grids_number, max_steps, data_collection_period, new_simulation_folder, loaded_simulation_path)
    for i in range(max_steps):
        model.step()
    print(f'Finished the simulation at time step {model.schedule.time}!')
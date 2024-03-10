import ast
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
import re
import os
import json
import metaspread.configs

# To run this code you must be in the parent folder of agent-based-cancer

def read_coords_for_plot(step, all_cells_dataframe, grid_id):
    df = all_cells_dataframe[["Step", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]

    # Selected step to plot
    df_selected_step_and_grid = df.loc[(df["Step"] == step) & (df["Grid"] == grid_id)]

    # Save the position data
    mPoints = list(df_selected_step_and_grid.loc[df_selected_step_and_grid["Phenotype"] == "mesenchymal"]["Position"])
    ePoints = list(df_selected_step_and_grid.loc[df_selected_step_and_grid["Phenotype"] == "epithelial"]["Position"])
    vPoints = list(df_selected_step_and_grid.loc[(df_selected_step_and_grid["Agent Type"] == "vessel") & (df_selected_step_and_grid["Ruptured"] == False)]["Position"])
    vRupturedPoints = list(df_selected_step_and_grid.loc[(df_selected_step_and_grid["Agent Type"] == "vessel") & (df_selected_step_and_grid["Ruptured"] == True)]["Position"])

    X_position_m, Y_position_m = [i[0] for i in mPoints], [i[1] for i in mPoints]
    X_position_e, Y_position_e = [i[0] for i in ePoints], [i[1] for i in ePoints]
    X_position_v, Y_position_v = [i[0] for i in vPoints], [i[1] for i in vPoints]
    X_position_vr, Y_position_vr = [i[0] for i in vRupturedPoints], [i[1] for i in vRupturedPoints]

    return [X_position_m, Y_position_m, X_position_e, Y_position_e, X_position_v, Y_position_v, X_position_vr, Y_position_vr]

def save_cancer(coords_list, grid_id, step, real_time_at_step, TumorDataPath):
    path = os.path.join(TumorDataPath, f'Cells-grid{grid_id}-step{step} - Tumor size at {real_time_at_step/(3600*24):.2f} days.csv')
    if not os.path.isfile(path):
        Xm, Ym, Xe, Ye, Xv, Yv, Xvr, Yvr = coords_list[0], coords_list[1], coords_list[2], coords_list[3], coords_list[4], coords_list[5], coords_list[6], coords_list[7]

        # save the data
        df_export = pd.DataFrame([Xm, Ym, Xe, Ye, Xv, Yv, Xvr, Yvr]) 
        df_export.to_csv(path)

    #create histogram of positions
    path = os.path.join(TumorDataPath, f'Cells-grid{grid_id}-step{step} - Histogram at {real_time_at_step/(3600*24):.2f} days.csv')
    if not os.path.isfile(path):
        Xm, Ym, Xe, Ye, Xv, Yv, Xvr, Yvr = coords_list[0], coords_list[1], coords_list[2], coords_list[3], coords_list[4], coords_list[5], coords_list[6], coords_list[7]
        df_positions = pd.DataFrame({'Position': zip(Xm + Xe, Ym + Ye)})
        position_repetition_count = df_positions['Position'].value_counts()
        histogram = position_repetition_count.value_counts()
        histogram = pd.DataFrame({'Bins': histogram.index, 'Frequency': histogram.values})
        number_of_empty_positions = metaspread.configs.gridsize * metaspread.configs.gridsize - len(position_repetition_count)
        new_row = pd.DataFrame({'Bins': [0], 'Frequency': [number_of_empty_positions]})
        histogram = pd.concat([histogram, new_row])
        histogram.to_csv(path)
        #this will return the radius and diameter of the processed tumor
        return get_cluster_centroid_radius_and_diameter(df_positions, grid_id)
    return ([np.nan, np.nan], np.nan, np.nan)
    
def save_growth_data_n(all_cells_dataframe, grid_id, CellsDataPath, step_number, real_time_at_step, real_delta_time, df_csv_last_step):
    path_to_save = os.path.join(CellsDataPath, f'CellsGrowth-grid{grid_id}-step{step_number} - {real_time_at_step/(3600*24):.2f} days.csv')
    if os.path.isfile(path_to_save):
        return
    df = all_cells_dataframe.loc[(all_cells_dataframe["Grid"] == grid_id) & (all_cells_dataframe["Step"] == step_number)]
    amount_of_mesenchymal = len(df[df["Phenotype"] == "mesenchymal"])
    amount_of_epithelial = len(df[df["Phenotype"] == "epithelial"])
    df = pd.DataFrame({"Number of Epithelial Cells": [amount_of_epithelial], "Number of Mesenchymal Cells": [amount_of_mesenchymal], "Steps": [step_number], "Days": [real_delta_time*step_number/(3600*24)]})
    if not df_csv_last_step.empty:
        df_csv = pd.concat([df_csv_last_step, df], ignore_index=True)
    else:
        df_csv = df
    df_csv.to_csv(path_to_save)
    return df_csv

def save_growth_data(all_cells_dataframe, stepsize, grid_id, CellsDataPath, step_number, real_time_at_step, real_delta_time):
    path_to_save = os.path.join(CellsDataPath, f'CellsGrowth-grid{grid_id}-step{step_number} - {real_time_at_step/(3600*24):.2f} days.csv')
    if os.path.isfile(path_to_save):
        return
    df = all_cells_dataframe[["Step", "Phenotype", "Grid"]]
    df = df.loc[(df["Grid"] == grid_id) & (df["Step"] <= step_number)]

    df_m = df.loc[(df["Step"] % stepsize == 0)]
    if df_m.empty:
        return
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_m["Step"]) + 1, stepsize)
    numberMesenchymalEachStep = [df_m.loc[(df_m["Step"] == step) & (df_m["Phenotype"] == "mesenchymal")].shape[0] for step in steps]

    df_e = df.loc[(df["Step"] % stepsize == 0)]
    if df_e.empty:
        return
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_e["Step"]) + 1, stepsize)
    numberEpithelialEachStep = [df_e.loc[(df_e["Step"] == step) & (df_e["Phenotype"] == "epithelial")].shape[0] for step in steps]
    #save data from plot into csv
    df_csv = pd.DataFrame({"Number of Epithelial Cells": numberEpithelialEachStep, "Number of Mesenchymal Cells": numberMesenchymalEachStep, "Steps": steps, "Days": real_delta_time* steps/(3600*24)})
    df_csv.to_csv(path_to_save)

def get_vasculature_state_at_step(pathToSave, vasculature_json_path, step):
    # Reads the dict in the json file
    with open(vasculature_json_path, 'r') as f:
        vasculature_dict = json.load(f)

    # Change keys to int
    vasculature_dict = {int(k): v for k, v in vasculature_dict.items()}

    # Prepare the data for the bar chart
    mesenchymal_count = 0
    epithelial_count = 0
    multicellular_cluster_count = 0
    total_cluster_count = 0
    time_steps = sorted(vasculature_dict.keys())
    for time_step in time_steps:
        clusters = vasculature_dict[time_step]
        for cluster in clusters:
            total_cluster_count += 1
            mesenchymal_count += cluster[0]
            epithelial_count += cluster[1]
            if cluster[0] + cluster[1] > 1:
                multicellular_cluster_count += 1
    df_export = pd.DataFrame({ "Time": [step], "Mesenchymal cells": [mesenchymal_count], "Epithelial cells" :[epithelial_count], "Multicellular clusters": [multicellular_cluster_count], "Total clusters": [total_cluster_count]})
    return df_export

def saveVasculatureData(pathToSave, vasculature_json_path, max_step):
    # Reads the dict in the json file
    with open(vasculature_json_path, 'r') as f:
        vasculature_dict = json.load(f)

    # Change keys to int and only add the key-value pair that is before the given max_step
    vasculature_dict = {int(k): v for k, v in vasculature_dict.items() if int(k) <= max_step}

    # Prepare the data for the bar chart
    mesenchymal_data = []
    epithelial_data = []
    multicellular_cluster_data = []
    total_cluster_data = []
    time_steps = sorted(vasculature_dict.keys())
    for time_step in time_steps:
        clusters = vasculature_dict[time_step]
        mesenchymal_count = 0
        epithelial_count = 0
        multicellular_cluster_count = 0
        total_cluster_count = 0
        for cluster in clusters:
            total_cluster_count += 1
            mesenchymal_count += cluster[0]
            epithelial_count += cluster[1]
            if cluster[0] + cluster[1] > 1:
                multicellular_cluster_count += 1
        mesenchymal_data.append(mesenchymal_count)
        epithelial_data.append(epithelial_count)
        multicellular_cluster_data.append(multicellular_cluster_count)
        total_cluster_data.append(total_cluster_count)
    df_export = pd.DataFrame({ "Time": time_steps, "Mesenchymal cells": mesenchymal_data, "Epithelial cells" :epithelial_data, "Multicellular clusters": multicellular_cluster_data, "Total clusters": total_cluster_data})
    return df_export

def generate_data(nameOfTheSimulation):
    SimulationPath = os.path.join("Simulations", nameOfTheSimulation)
    
    configs_path = os.path.join(SimulationPath, "configs.csv")
    metaspread.configs.load_simulation_configs_for_data_generation(configs_path)
    
    print(f'\tAnalyzing data in the folder {SimulationPath}\n')

    # Get the Ecm and Mmp2 data filenames 
    EcmPath = os.path.join(SimulationPath, "Ecm")
    Mmp2Path = os.path.join(SimulationPath, "Mmp2")
    ecm_files_name = [f for f in sorted(os.listdir(EcmPath), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(EcmPath, f)) and f.endswith(".csv")]
    mmp2_files_name = [f for f in sorted(os.listdir(Mmp2Path), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(Mmp2Path, f)) and f.endswith(".csv")]

    # Get the vasculature data filename
    VasculaturePath = os.path.join(SimulationPath, "Vasculature")
    vasculature_files_name = [f for f in os.listdir(VasculaturePath) if os.path.isfile(os.path.join(VasculaturePath, f)) and f.endswith(".json")]

    all_cells_filename = os.path.join(SimulationPath, "CellsData.csv")
    if os.path.isfile(all_cells_filename):
        print("Using cells data at:", all_cells_filename)
    else:
        print("No CellsData.csv cell data found in directory:", SimulationPath)
        return

    if ecm_files_name:
        ecm_files_path = [os.path.join(EcmPath, p) for p in ecm_files_name]
        print("Using Ecm data in the folder:", EcmPath)
    else:
        print("No .csv Ecm data found in directory:", EcmPath)
        return

    if mmp2_files_name:
        mmp2_files_path = [os.path.join(Mmp2Path, p) for p in mmp2_files_name]
        print("Using Mmp2 data in the folder:", Mmp2Path)
    else:
        print("No .csv Mmp2 data found in directory:", Mmp2Path)
        return

    if vasculature_files_name:
        # first_vasculature_name = vasculature_files_name[0]
        # first_vasculature_path = os.path.join(VasculaturePath, first_vasculature_name)
        print("Using vasculature data at:", VasculaturePath)
    else:
        print("No .json vasculature data found in directory:", SimulationPath)
        return

    step_size = metaspread.configs.data_collection_period
    real_delta_time = 40 * metaspread.configs.th/0.001 #in seconds (the original ratio is 40 seconds/0.001 non-dimensional time)
    grids_number = metaspread.configs.grids_number
    configs_max_step = metaspread.configs.max_steps
    all_cells_dataframe = pd.read_csv(all_cells_filename, converters={"Position": ast.literal_eval})
    all_cells_dataframe = all_cells_dataframe[["Step", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]
    max_step = max(all_cells_dataframe["Step"])
    if configs_max_step != max_step:
        print(f"Warning: the run for this simulation terminated early")
        print(f"Max step reached is {max_step} while {configs_max_step} was expected.")

    # Path to save the data:
    dataFolder = "Data analysis"
    dataPath = os.path.join(SimulationPath, dataFolder)

    TumorDataPath = os.path.join(dataPath, "Tumor growth")
    CellsDataPath = os.path.join(dataPath, "Cells growth")
    VasculatureDataPath = os.path.join(dataPath, "Vasculature evolution")

    # Create folder for all the data analysis
    os.makedirs(dataPath, exist_ok = True)
    os.makedirs(TumorDataPath, exist_ok = True)
    os.makedirs(CellsDataPath, exist_ok = True)
    os.makedirs(VasculatureDataPath, exist_ok = True)

    print(f"\nSaving tumor data in the folder:", TumorDataPath)
    print("Saving cells numbers data in the folder:", CellsDataPath)
    print("Saving vasculature data in the folder:", VasculatureDataPath)
    
    for grid_id in range(1, grids_number+1):
        print(f'\nGrid: {grid_id}')

        # Plot the cells graphs
        print(f'\tSaving tumor data...')
        df_radius_diameter_history = pd.DataFrame(columns=['Centroid x', 'Centroid y', 'Radius', 'Diameter', 'Step', 'Grid Id'])
        for id, step in enumerate(range(step_size,max_step+1,step_size)):
            ccells_coords = read_coords_for_plot(step, all_cells_dataframe, grid_id)
            real_time_at_step = real_delta_time * step
            if grid_id == 1:
                (centroid, radius, diameter) = save_cancer(ccells_coords, grid_id, step, real_time_at_step, TumorDataPath)
                new_row = pd.DataFrame({'Centroid x': [centroid[0]], 'Centroid y': [centroid[1]],'Radius': [radius], 'Diameter': [diameter], 'Step': [step], 'Grid Id': [grid_id]})
                df_radius_diameter_history = pd.concat([df_radius_diameter_history, new_row])
            else:
                save_cancer(ccells_coords, grid_id, step, real_time_at_step, TumorDataPath)
        if grid_id == 1:
            path = os.path.join(TumorDataPath, f'Tumor radius and diameter history in grid {grid_id}.csv')
            if not os.path.isfile(path):
                df_radius_diameter_history.to_csv(path)


        # Plot the growth of ephitelial and mesenchymal cells
        print(f'\tSaving cells numbers graph data...')
        df_csv_last_step = pd.DataFrame()
        for id, step in enumerate(range(step_size,max_step+1,step_size)):
            real_time_at_step = real_delta_time * step
            # save_growth_data(all_cells_dataframe, step_size, grid_id , CellsDataPath, step, real_time_at_step, real_delta_time)
            df_csv_last_step = save_growth_data_n(all_cells_dataframe, grid_id , CellsDataPath, step, real_time_at_step, real_delta_time, df_csv_last_step)

    # Plot the vasculature data
    print(f'Saving vasculature...')
    df_export = pd.DataFrame(columns=["Time", "Mesenchymal cells", "Epithelial cells", "Multicellular clusters", "Total clusters"])
    for step in range(step_size,max_step+1,step_size):
        path = os.path.join(VasculatureDataPath, f'Vasculature-step{step}.csv')
        if os.path.isfile(path):
            continue
        file_vasculature_on_step = os.path.join(VasculaturePath, f"Vasculature-{step}step.json")
        row = get_vasculature_state_at_step(VasculatureDataPath, file_vasculature_on_step, step)
        df_export = pd.concat([df_export, row])
        df_export.to_csv(path)

def generate_data_vasculature_only(nameOfTheSimulation):
    SimulationPath = os.path.join("Simulations", nameOfTheSimulation)
    
    configs_path = os.path.join(SimulationPath, "configs.csv")
    metaspread.configs.load_simulation_configs_for_data_generation

    # Get the vasculature data filename
    VasculaturePath = os.path.join(SimulationPath, "Vasculature")
    vasculature_files_name = [f for f in os.listdir(VasculaturePath) if os.path.isfile(os.path.join(VasculaturePath, f)) and f.endswith(".json")]

    all_cells_filename = os.path.join(SimulationPath, "CellsData.csv")
    if os.path.isfile(all_cells_filename):
        print("Using cells data at:", all_cells_filename)
    else:
        print("No CellsData.csv cell data found in directory:", SimulationPath)
        return
    if vasculature_files_name:
        # first_vasculature_name = vasculature_files_name[0]
        # first_vasculature_path = os.path.join(VasculaturePath, first_vasculature_name)
        # print("Using vasculature data at:", first_vasculature_path)
        print("Using vasculature data at:", VasculaturePath)
    else:
        print("No .json vasculature data found in directory:", SimulationPath)
        return

    print("Loading CellsData.csv. This might take a minute...")
    step_size = metaspread.configs.data_collection_period
    configs_max_step = metaspread.configs.max_steps
    all_cells_dataframe = pd.read_csv(all_cells_filename, converters={"Position": ast.literal_eval})
    all_cells_dataframe = all_cells_dataframe[["Step", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]
    max_step = max(all_cells_dataframe["Step"])
    if configs_max_step != max_step:
        print(f"Warning: the run for this simulation terminated early")
        print(f"Max step reached is {max_step} while {configs_max_step} was expected.")

    # Path to save the vasculature data:
    dataFolder = "Data analysis"
    dataPath = os.path.join(SimulationPath, dataFolder)
    VasculatureDataPath = os.path.join(dataPath, "Vasculature evolution")

    print(f'Saving vasculature...')
    df_export = pd.DataFrame(columns=["Time", "Mesenchymal cells", "Epithelial cells", "Multicellular clusters", "Total clusters"])
    for step in range(step_size,max_step+1,step_size):
        file_vasculature_on_step = os.path.join(VasculaturePath, f"Vasculature-{step}step.json")
        row = get_vasculature_state_at_step(VasculatureDataPath, file_vasculature_on_step, step)
        df_export = pd.concat([df_export, row])
        path = os.path.join(VasculatureDataPath, f'Vasculature-step{step}.csv')
        df_export.to_csv(path)

def get_cluster_centroid_radius_and_diameter(ccells_positions, grid_id):
    """"
    Calculates the radius and diameter of the cancer cells in a given site.

    Input:
        model_dataframe: the dataframe containing al information about the model
        grid_id: the grid number for which the radius and diameter will be calculated.
    Returns:
        radius: the maximum of all the distances from each cell to the cell's centroid.
        diameter: the maximum of all the cell-cell distances.
    """
    if ccells_positions.empty or grid_id != 1:
        return ([np.nan, np.nan], np.nan, np.nan)
    #centroid = ccells_positions.mean(axis=0)
    ccells_positions= list(ccells_positions['Position'].unique())
    # ccells_positions= list(set(ccells_positions))#delete repeated entries to spped up computation
    centroid = np.average(ccells_positions, axis=0)
    #calculating radius
    radii  = np.linalg.norm(ccells_positions - centroid, axis=1)
    radius = radii.max()
    #calculating diameter
    dist_matrix = get_distance_matrix(ccells_positions)
    diameter = dist_matrix.max()
    return (centroid, radius, diameter)

def get_distance_matrix(vectors):
    """
    Calculates the distance matrix between all the vectors in a list

    Input:
        vectors (list of numpy 2 by 1 arrays): 
    Returns:
        distance matrix: len(vectors) by len(Vectors) numpy array with the distances
        for all vectors.
    """
    vectors = np.array(vectors)
    x = np.sum(vectors**2,axis=1)
    xx = np.matmul(vectors,vectors.T)
    x2 = x.reshape(-1,1) #transposing the vector
    return np.sqrt(x2-2*xx+x)


if __name__ == "__main__":

    # CHANGE THIS LINE according to the simulation you want to plot the graphs  
    name_of_the_simulation = "Sim max_steps-48000 stepsize-150 N-388 gridsNumber-3"

    # This runs all the code to generate the graphs in the folder
    generate_data(name_of_the_simulation)
    # generate_data_vasculature_only(name_of_the_simulation)
import ast
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
# from Classes.utils import gridsize_utils as gridsize
import re
import os
import json
import Classes.configs

# To run this code you must be in the parent folder of agent-based-cancer

def getCoordsForPlot(step, allCellsCsvPath, grid_id):
    df = pd.read_csv(allCellsCsvPath, converters={"Position": ast.literal_eval})
    df = df[["Step", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]

    # Select the step you want to plot, from 0 to 24000 (~11 days)
    df_step0 = df.loc[(df["Step"] == step) & (df["Grid"] == grid_id)]

    # Save the position data
    mPoints = list(df_step0.loc[df_step0["Phenotype"] == "mesenchymal"]["Position"])
    ePoints = list(df_step0.loc[df_step0["Phenotype"] == "epithelial"]["Position"])
    vPoints = list(df_step0.loc[(df_step0["Agent Type"] == "vessel") & (df_step0["Ruptured"] == False)]["Position"])
    vRupturedPoints = list(df_step0.loc[(df_step0["Agent Type"] == "vessel") & (df_step0["Ruptured"] == True)]["Position"])

    Xm, Ym = [i[0] for i in mPoints], [i[1] for i in mPoints]
    Xe, Ye = [i[0] for i in ePoints], [i[1] for i in ePoints]
    Xv, Yv = [i[0] for i in vPoints], [i[1] for i in vPoints]
    Xvr, Yvr = [i[0] for i in vRupturedPoints], [i[1] for i in vRupturedPoints]

    return [Xm, Ym, Xe, Ye, Xv, Yv, Xvr, Yvr]

def save_cancer(coordsList, grid_id, step, TumorDataPath):

    Xm, Ym, Xe, Ye, Xv, Yv, Xvr, Yvr = coordsList[0], coordsList[1], coordsList[2], coordsList[3], coordsList[4], coordsList[5], coordsList[6], coordsList[7]

    # save the data
    df_export = pd.DataFrame([Xm, Ym, Xe, Ye, Xv, Yv, Xvr, Yvr])
    path = os.path.join(TumorDataPath, f'Cells-grid{grid_id}-step{step} - Tumor size at {11/24000 * step:.2f} days.csv')
    df_export.to_csv(path)

    #create histogram of positions
    df_positions = pd.DataFrame({'Position': zip(Xm + Xe, Ym + Ye)})
    position_repetition_count = df_positions['Position'].value_counts()
    histogram = position_repetition_count.value_counts()
    histogram = pd.DataFrame({'Bins': histogram.index, 'Frequency': histogram.values})
    number_of_empty_positions = Classes.configs.gridsize_utils * Classes.configs.gridsize_utils - len(position_repetition_count)
    new_row = pd.DataFrame({'Bins': [0], 'Frequency': [number_of_empty_positions]})
    histogram = pd.concat([histogram, new_row])
    path = os.path.join(TumorDataPath, f'Cells-grid{grid_id}-step{step} - Histogram at {11/24000 * step:.2f} days.csv')
    histogram.to_csv(path)

    #this will return the radius and diameter of the processed tumor
    return get_cluster_radius_and_diameter(df_positions, grid_id)

def saveGrowthData(allCellsCsvPath, stepsize, grid_id, CellsDataPath, step_number):
    # get data at each step
    df = pd.read_csv(allCellsCsvPath, converters={"Position": ast.literal_eval})
    df = df[["Step", "Phenotype", "Grid"]]
    df = df.loc[(df["Grid"] == grid_id) & (df["Step"] <= step_number)]

    # For mesenchymal
    #stepsize = 3000 doubling rate, since it is different it will stay horizontal sometimes
    df_m = df.loc[(df["Step"] % stepsize == 0)]
    if df_m.empty:
        return
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_m["Step"]) + 1, stepsize)
    numberMesenchymalEachStep = [df_m.loc[(df_m["Step"] == step) & (df_m["Phenotype"] == "mesenchymal")].shape[0] for step in steps]

    # For epithelial
    #stepsize = 2000 doubling rate
    df_e = df.loc[(df["Step"] % stepsize == 0)]
    if df_e.empty:
        return
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_e["Step"]) + 1, stepsize)
    numberEpithelialEachStep = [df_e.loc[(df_e["Step"] == step) & (df_e["Phenotype"] == "epithelial")].shape[0] for step in steps]

    path_to_save = os.path.join(CellsDataPath, f'CellsGrowth-grid{grid_id}-step{step_number} - {11/24000 * step_number:.2f} days')

    #save data from plot into csv
    df_csv = pd.DataFrame({"Number of Epithelial Cells": numberEpithelialEachStep, "Number of Mesenchymal Cells": numberMesenchymalEachStep, "Days": steps*11/24000})
    df_csv.to_csv(path_to_save + ".csv")

def get_vasculature_state_at_step(pathToSave, vasculature_json_path, step):
    # Reads the dict in the json file
    with open(vasculature_json_path, 'r') as f:
        vasculature_dict = json.load(f)

    # Change keys to int and only add the key-value pair that is before the given max_step
    vasculature_dict = {int(k): v for k, v in vasculature_dict.items()}
    # vasculature_dict = {int(k): v for k, v in vasculature_dict.items() if int(k) <= max_step}

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
    # path = os.path.join(pathToSave, f'Vasculature-step{max_step}.csv')
    # df_export.to_csv(path)

def generate_data(nameOfTheSimulation):
    SimulationPath = os.path.join("Simulations", nameOfTheSimulation)
    
    configs_path = os.path.join(SimulationPath, "configs.csv")
    Classes.configs.load_simulation_configs_for_data_generation(configs_path)
    
    print(f'\tAnalyzing data in the folder {SimulationPath}\n')

    # Get the agents' data filename 
    csv_files_name = [f for f in os.listdir(SimulationPath) if os.path.isfile(os.path.join(SimulationPath, f)) and f.endswith(".csv")]
    
    # Get the Ecm and Mmp2 data filenames 
    EcmPath = os.path.join(SimulationPath, "Ecm")
    Mmp2Path = os.path.join(SimulationPath, "Mmp2")
    ecm_files_name = [f for f in sorted(os.listdir(EcmPath), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(EcmPath, f)) and f.endswith(".csv")]
    mmp2_files_name = [f for f in sorted(os.listdir(Mmp2Path), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(Mmp2Path, f)) and f.endswith(".csv")]

    # Get the vasculature data filename
    VasculaturePath = os.path.join(SimulationPath, "Vasculature")
    vasculature_files_name = [f for f in os.listdir(VasculaturePath) if os.path.isfile(os.path.join(VasculaturePath, f)) and f.endswith(".json")]

    if csv_files_name:
        first_csv_name = csv_files_name[0]
        first_csv_path = os.path.join(SimulationPath, first_csv_name)
        print("Using cells data at:", first_csv_path)
    else:
        print("No .csv cell data found in directory:", SimulationPath)
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
        first_vasculature_name = vasculature_files_name[0]
        first_vasculature_path = os.path.join(VasculaturePath, first_vasculature_name)
        print("Using vasculature data at:", first_vasculature_path)
    else:
        print("No .json vasculature data found in directory:", SimulationPath)
        return

    step_size = Classes.configs.dataCollectionPeriod
    grids_number = Classes.configs.grids_number
    configs_max_step = Classes.configs.maxSteps
    df = pd.read_csv(first_csv_path, converters={"Position": ast.literal_eval})
    df = df[["Step", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]
    max_step = max(df["Step"])
    if configs_max_step != max_step:
        print(f"Warning: the run for this simulation terminated early")
        print("Max step reached is {max_step} while {configs_max_step} was expected.")

    # Path to save the data:
    dataFolder = "Data analysis"
    dataPath = os.path.join(SimulationPath, dataFolder)

    TumorDataPath = os.path.join(dataPath, "Tumor growth")
    CellsDataPath = os.path.join(dataPath, "Cells growth")
    # EcmDataPath = os.path.join(dataPath, "Ecm evolution")
    # Mmp2DataPath = os.path.join(dataPath, "Mmp2 evolution")
    VasculatureDataPath = os.path.join(dataPath, "Vasculature evolution")

    # Create folder for all the data analysis
    if not os.path.exists(dataPath):
        os.makedirs(dataPath)
        os.makedirs(TumorDataPath)
        # os.makedirs(Mmp2DataPath)
        # os.makedirs(EcmDataPath)
        os.makedirs(CellsDataPath)
        os.makedirs(VasculatureDataPath)

        print(f"\nSaving tumor data in the folder:", TumorDataPath)
        # print("Saving Mmp2 data in the folder:", Mmp2DataPath)
        # print("Saving Ecm data in the folder:", EcmDataPath)
        print("Saving cells numbers data in the folder:", CellsDataPath)
        print("Saving vasculature data in the folder:", VasculatureDataPath)

    # If the data analysis is already done, tell the user
    else:
        print("This data has already been saved!")
        return
    
    for grid_id in range(1, grids_number+1):
        print(f'\nGrid: {grid_id}')

        # Plot the cells graphs
        print(f'\tSaving tumor data...')
        df_radius_diameter_history = pd.DataFrame(columns=['Centroid x', 'Centroid y', 'Radius', 'Diameter', 'Step', 'Grid Id'])
        for id, step in enumerate(range(step_size,max_step+1,step_size)):
            ccells_coords = getCoordsForPlot(step, first_csv_path, grid_id)
            # save_cancer(ccells_coords, grid_id, step, TumorDataPath)
            (centroid, radius, diameter) = save_cancer(ccells_coords, grid_id, step, TumorDataPath)
            new_row = pd.DataFrame({'Centroid x': [centroid[0]], 'Centroid y': [centroid[1]],'Radius': [radius], 'Diameter': [diameter], 'Step': [step], 'Grid Id': [grid_id]})
            df_radius_diameter_history = pd.concat([df_radius_diameter_history, new_row])
        path = os.path.join(TumorDataPath, f'Tumor radius and diameter history in grid {grid_id} at {11/24000 * step:.2f} days.csv')
        df_radius_diameter_history.to_csv(path)


        # Plot the growth of ephitelial and mesenchymal cells
        print(f'\tSaving cells numbers graph data...')
        for id, step in enumerate(range(step_size,max_step+1,step_size)):
            saveGrowthData(first_csv_path, step_size, grid_id , CellsDataPath, step)

    # Plot the vasculature data
    print(f'Saving vasculature...')
    df_export = pd.DataFrame(columns=["Time", "Mesenchymal cells", "Epithelial cells", "Multicellular clusters", "Total clusters"])
    for step in range(step_size,max_step+1,step_size):
        file_vasculature_on_step = os.path.join(VasculaturePath, f"Vasculature-{step}step.json")
        row = get_vasculature_state_at_step(VasculatureDataPath, file_vasculature_on_step, step)
        df_export = pd.concat([df_export, row])
        path = os.path.join(VasculatureDataPath, f'Vasculature-step{step}.csv')
        df_export.to_csv(path)

def generate_data_vasculature_only(nameOfTheSimulation):
    SimulationPath = os.path.join("Simulations", nameOfTheSimulation)
    
    configs_path = os.path.join(SimulationPath, "configs.csv")
    Classes.configs.load_simulation_configs_for_data_generation(configs_path)
    
    # Get the agents' data filename 
    csv_files_name = [f for f in os.listdir(SimulationPath) if os.path.isfile(os.path.join(SimulationPath, f)) and f.endswith(".csv")]

    # Get the vasculature data filename
    VasculaturePath = os.path.join(SimulationPath, "Vasculature")
    vasculature_files_name = [f for f in os.listdir(VasculaturePath) if os.path.isfile(os.path.join(VasculaturePath, f)) and f.endswith(".json")]

    if csv_files_name:
        first_csv_name = csv_files_name[0]
        first_csv_path = os.path.join(SimulationPath, first_csv_name)
        print("Using cells data at:", first_csv_path)
    else:
        print("No .csv cell data found in directory:", SimulationPath)
        return
    if vasculature_files_name:
        first_vasculature_name = vasculature_files_name[0]
        first_vasculature_path = os.path.join(VasculaturePath, first_vasculature_name)
        print("Using vasculature data at:", first_vasculature_path)
    else:
        print("No .json vasculature data found in directory:", SimulationPath)
        return

    step_size = Classes.configs.dataCollectionPeriod
    configs_max_step = Classes.configs.maxSteps
    df = pd.read_csv(first_csv_path, converters={"Position": ast.literal_eval})
    df = df[["Step", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]
    max_step = max(df["Step"])
    if configs_max_step != max_step:
        print(f"Warning: the run for this simulation terminated early")
        print("Max step reached is {max_step} while {configs_max_step} was expected.")

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

def get_cluster_radius_and_diameter(ccells_positions, grid_id):
    """"
    Calculates the radius and diameter of the cancer cells in a given site.

    Input:
        model_dataframe: the dataframe containing al information about the model
        grid_id: the grid number for which the radius and diameter will be calculated.
    Returns:
        radius: the maximum of all the distances from each cell to the cell's centroid.
        diameter: the maximum of all the cell-cell distances.
    """
    if ccells_positions.empty:
        return ([np.nan, np.nan], np.nan, np.nan)
    #centroid = ccells_positions.mean(axis=0)
    ccells_positions= list(ccells_positions['Position'])
    ccells_positions= list(set(ccells_positions))#delete repeated entries to spped up computation
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
    name_of_the_simulation = "Sim maxSteps-22 stepsize-2 N-388 gridsNumber-3"

    # This runs all the code to generate the graphs in the folder
    # generate_data(name_of_the_simulation)
    generate_data_vasculature_only(name_of_the_simulation)
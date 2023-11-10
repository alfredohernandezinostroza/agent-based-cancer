import ast
import pandas as pd
import numpy as np
from Classes.utils import parent_dir
from Classes.utils import gridsize_utils as gridsize
import re
import os
import json

# To run this code you must be in the parent folder of agent-based-cancer

def getCoordsForPlot(step, allCellsCsvPath, grid_id):
    df = pd.read_csv(allCellsCsvPath, converters={"Position": ast.literal_eval})
    df = df[["Step", "Total cells", "Position", "Phenotype", "Grid", "Agent Type", "Ruptured"]]

    # Select the step you want to plot, from 0 to 24000 (~11 days)
    df_step0 = df.loc[(df["Step"] == step) & (df["Grid"] == grid_id)]

    # Save the position data
    mPoints = df_step0.loc[df_step0["Phenotype"] == "mesenchymal"]["Position"] # Series object
    ePoints = df_step0.loc[df_step0["Phenotype"] == "epithelial"]["Position"]
    vPoints = df_step0.loc[(df_step0["Agent Type"] == "vessel") & (df_step0["Ruptured"] == False)]["Position"]
    vRupturedPoints = df_step0.loc[(df_step0["Agent Type"] == "vessel") & (df_step0["Ruptured"] == True)]["Position"]

    mPoints = list(mPoints.map(eval)) # [(104, 101), (101, 97), (101, 95)]
    ePoints = list(ePoints.map(eval))
    vPoints = list(vPoints.map(eval))
    vRupturedPoints = list(vRupturedPoints.map(eval))

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
    histogram = pd.DataFrame({'Bins': histogram.values, 'Frequency': histogram.index})
    number_of_empty_positions = gridsize * gridsize - len(position_repetition_count)
    new_row = pd.DataFrame({'Bins': [0], 'Frequency': [number_of_empty_positions]})
    histogram = pd.concat([histogram, new_row])
    path = os.path.join(TumorDataPath, f'Cells-grid{grid_id}-step{step} - Histogram at {11/24000 * step:.2f} days.csv')
    histogram.to_csv(path)

    #this will return the radius and diameter of the processed tumor
    return get_cluster_radius_and_diameter(df_positions, grid_id)

def saveGrowthData(allCellsCsvPath, stepsize, grid_id, CellsDataPath, step_number):
    # get data at each step
    df = pd.read_csv(allCellsCsvPath, converters={"Position": ast.literal_eval})
    df = df[["Step", "Total cells", "Phenotype", "Grid"]]
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

def saveVasculatureData(pathToSave, vasculature_json_path, max_step):
    # Reads the dict in the json file
    with open(vasculature_json_path, 'r') as f:
        vasculature_dict = json.load(f)

    # Change keys to int and only add the key-value pair that is before the given max_step
    vasculature_dict = {int(k): v for k, v in vasculature_dict.items() if int(k) <= max_step}

    # Prepare the data for the bar chart
    mesenchymal_data = []
    epithelial_data = []
    cluster_data = []
    time_steps = sorted(vasculature_dict.keys())
    for time_step in time_steps:
        clusters = vasculature_dict[time_step]
        mesenchymal_count = 0
        epithelial_count = 0
        cluster_count = 0
        for cluster in clusters:
            mesenchymal_count += cluster[0]
            epithelial_count += cluster[1]
            if cluster[0] + cluster[1] > 1:
                cluster_count += 1

        mesenchymal_data.append(mesenchymal_count)
        epithelial_data.append(epithelial_count)
        cluster_data.append(cluster_count)
    df_export = pd.DataFrame({"Mesenchymal cells": mesenchymal_data, "Epithelial cells" :epithelial_data, "Multicellular clusters": cluster_data, "Time": time_steps})
    path = os.path.join(pathToSave, f'Vasculature-step{max_step}.csv')
    df_export.to_csv(path)

def generate_data(nameOfTheSimulation):
    SimulationPath = os.path.join(parent_dir, nameOfTheSimulation)
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
        #all_cells_analysis_path = os.path.join(parent_dir, SimulationPath, first_csv_path)
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


    # use regex to find the values before 'steps', 'stepsize' and 'grids'. Ex: 50steps-10stepsize-cells
    max_step = int(re.search(r"(\d+)steps", first_csv_name).group(1))
    step_size = int(re.search(r"(\d+)stepsize", first_csv_name).group(1))
    grids_number = int(re.search(r"(\d+)grids", first_csv_name).group(1))

    # Path to save the data:
    dataFolder = "Data analysis"
    dataPath = os.path.join(SimulationPath, dataFolder)

    TumorDataPath = os.path.join(dataPath, "Tumor growth")
    CellsDataPath = os.path.join(dataPath, "Cells growth")
    EcmDataPath = os.path.join(dataPath, "Ecm evolution")
    Mmp2DataPath = os.path.join(dataPath, "Mmp2 evolution")
    VasculatureDataPath = os.path.join(dataPath, "Vasculature evolution")
    histogram_data_path = os.path.join(dataPath, "Positions histogram")

    # Create folder for all the data analysis
    if not os.path.exists(dataPath):
        os.makedirs(dataPath)
        os.makedirs(TumorDataPath)
        os.makedirs(Mmp2DataPath)
        os.makedirs(EcmDataPath)
        os.makedirs(CellsDataPath)
        os.makedirs(VasculatureDataPath)

        print(f"\nSaving tumor data in the folder:", TumorDataPath)
        print("Saving Mmp2 data in the folder:", Mmp2DataPath)
        print("Saving Ecm data in the folder:", EcmDataPath)
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
        df_radius_diameter_history = pd.DataFrame(columns=['Radius', 'Diameter', 'Step', 'Grid Id'])
        for id, step in enumerate(range(1,max_step+1,step_size)):
            (radius, diameter) = save_cancer(getCoordsForPlot(step, first_csv_path, grid_id), grid_id, step, TumorDataPath)
            new_row = pd.DataFrame({'Radius': [radius], 'Diameter': [diameter], 'Step': [step], 'Grid Id': [grid_id]})
            df_radius_diameter_history = pd.concat([df_radius_diameter_history, new_row])
        path = os.path.join(TumorDataPath, f'Tumor radius and diameter history at {11/24000 * step:.2f} days.csv')
        df_radius_diameter_history.to_csv(path)


        # Plot the growth of ephitelial and mesenchymal cells
        print(f'\tPlotting cells numbers graph...')
        for id, step in enumerate(range(1,max_step+1,step_size)):
            saveGrowthData(first_csv_path, step_size, grid_id , CellsDataPath, step)

    # Plot the vasculature data
    print(f'Plotting vasculature...')
    for id, step in enumerate(range(0,max_step+1,step_size)):
        saveVasculatureData(VasculatureDataPath, first_vasculature_path, step)

def get_cluster_radius_and_diameter(ccells_positions, grid_id):
    """"
    Calculates the radius and diameter of the cancer cells in a given site.

    Input:
        model_dataframe: the dataframe containing al information about the model
        grid_id: the grid number for which the radius and diameter will be calculated.
    Returns:
        radius: the maximum of all the distances from each cell to the cell's centroid.
        diameter: the maximum of all the cell-cell distanes.
    """
    if ccells_positions.empty:
        return (float("NaN"), float("NaN"))
    centroid = ccells_positions.mean(axis=0)
    #calculating radius
    radii  = np.linalg.norm(ccells_positions - centroid, axis=1)
    radius = radii.max()
    #calculating diameter
    dist_matrix = get_distance_matrix(ccells_positions)
    diameter = dist_matrix.max()
    return (radius, diameter)

def get_distance_matrix(vectors):
    """
    Calculates the distance matrix between all the vectors in a list

    Input:
        vectors (list of numpy 2 by 1 arrays): 
    Returns:
        distance matrix: len(vectors) by len(Vectors) numpy array with the distances
        for all vectors.
    """
    x = np.sum(vectors**2,axis=1)
    xx = np.matmul(vectors,vectors.T)
    x2 = x.reshape(-1,1) #transposing the vector
    return np.sqrt(x2-2*xx+x)

if __name__ == "__main__":

    # CHANGE THIS LINE according to the simulation you want to plot the graphs  
    nameOfTheSimulation = "Sim maxSteps-1100 stepsize-100 N-388 gridsNumber-3"

    # This runs all the code to generate the graphs in the folder
    generate_data(nameOfTheSimulation)







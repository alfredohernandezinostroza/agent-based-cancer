import ast
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Classes.utils import gridsize_utils as gridsize
from Classes.utils import carrying_capacity
import re
import os
import sys

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

def plotCancer(coordsList, figCounter, imagesFolder, grid_id, step, TumorImagesPath):

    Xm, Ym, Xe, Ye, Xv, Yv, Xvr, Yvr = coordsList[0], coordsList[1], coordsList[2], coordsList[3], coordsList[4], coordsList[5], coordsList[6], coordsList[7]
    plt.figure(figCounter, figsize=(6, 5))
    
    plt.scatter(Xm, Ym, marker='o', color='blue', alpha=0.10, label="Mesenchymal cells")
    plt.scatter(Xe, Ye, marker='h', color='orange', alpha=0.05, label="Epithelial cells")
    plt.scatter(Xv, Yv, marker='.', color='red', alpha=0.8, label="Vasculature points")
    plt.scatter(Xvr, Yvr, marker='+', color='darkred', alpha=0.8, label="Ruptured vasculature points")
    plt.xlim(0, gridsize)
    plt.ylim(0, gridsize)

    xticks = np.arange(0, gridsize, step=int(gridsize/6)) # 6 ticks
    xticklabels = [str(round(j,1)) for j in np.arange(0, 2.1, step = 2/201*(gridsize/6))]
    plt.xticks(xticks, xticklabels)
    plt.yticks(xticks, xticklabels)
    plt.xlabel("mm")
    plt.ylabel("mm")
    plt.grid(False)

    
    if step == 0:
        plt.title(f'Initial Tumor at {11/24000 * step:.2f} days ({step} steps) - grid {grid_id}', fontsize = 13)
        plt.legend(loc="upper left")
    else:
        plt.title(f'Tumor size at {11/24000 * step:.2f} days ({step} steps) - grid {grid_id}', fontsize = 13)

    # save the figure
    figure_path = os.path.join(TumorImagesPath, f'Cells-grid{grid_id}-step{step} - Tumor size at {11/24000 * step:.2f} days.png')
    plt.savefig(figure_path)
    # df_export = pd.DataFrame([Xm, Ym, Xe, Ye, Xv, Yv, Xvr, Yvr])
    # df_export.save_to_csv("dadsa")


def plotGrowthData(fig_index, allCellsCsvPath, stepsize, grid_id, CellsImagesPath, step_number):
    # get data at each step
    df = pd.read_csv(allCellsCsvPath, converters={"Position": ast.literal_eval})
    df = df[["Step", "Total cells", "Phenotype", "Grid"]]
    df = df.loc[(df["Grid"] == grid_id) & (df["Step"] <= step_number)]
    plt.figure(fig_index)

    # For mesenchymal
    #stepsize = 3000 doubling rate, since it is different it will stay horizontal sometimes
    df_m = df.loc[(df["Step"] % stepsize == 0)]
    if df_m.empty:
        return
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_m["Step"]) + 1, stepsize)
    numberMesenchymalEachStep = [df_m.loc[(df_m["Step"] == step) & (df_m["Phenotype"] == "mesenchymal")].shape[0] for step in steps]
    plt.plot(steps*11/24000, numberMesenchymalEachStep, label="Mesenchymal cells", color='tab:blue')

    # For epithelial
    #stepsize = 2000 doubling rate
    df_e = df.loc[(df["Step"] % stepsize == 0)]
    if df_e.empty:
        return
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_e["Step"]) + 1, stepsize)
    numberEpithelialEachStep = [df_e.loc[(df_e["Step"] == step) & (df_e["Phenotype"] == "epithelial")].shape[0] for step in steps]
    
    plt.plot(steps*11/24000, numberEpithelialEachStep, label="Epithelial cells", color='tab:orange')

    plt.title(f'Cells growth - grid {grid_id}', fontsize = 13)
    plt.xlabel('Days')
    plt.ylabel('Number of cells')
    plt.legend(loc="upper left")
    plt.ylim(0)


    # save the figure
    path_to_save = os.path.join(CellsImagesPath, f'CellsGrowth-grid{grid_id}-step{step_number} - {11/24000 * step_number:.2f} days')
    plt.savefig(path_to_save + ".png")

    # #save data from plot into csv
    # df_csv = pd.DataFrame({"Number of Epithelial Cells": numberEpithelialEachStep, "Number of Mesenchymal Cells": numberMesenchymalEachStep, "Days": steps*11/24000})
    # df_csv.to_csv(path_to_save + ".csv")


def plotMMP2orECM(i, step, files_path, figCounter, grid_id, pathToSave, type="Mmp2"):
    df = pd.read_csv(files_path[i], index_col=0)
    plt.figure(figCounter, figsize=(6, 5))

    if type=="Mmp2":
        plt.imshow(df, vmin=0, vmax=2) #vmax = 3 in PDF
        #print(f'mmp2 max at step {step} is {df.values.max()}')
    elif type == "Ecm":
        plt.imshow(df, vmin=0, vmax=1)
        #print(f'ecm min at step {step} is {df.values.min()}')
        
    plt.colorbar()

    plt.xlim(0, gridsize)
    plt.ylim(0, gridsize)

    xticks = np.arange(0, gridsize, step=int(gridsize/6)) # 6 ticks
    xticklabels = [str(round(j,1)) for j in np.arange(0, 2.1, step = 2/201*(gridsize/6))]
    plt.xticks(xticks, xticklabels)
    plt.yticks(xticks, xticklabels)
    plt.xlabel("mm")
    plt.ylabel("mm")
    plt.grid(visible=None)

    plt.title(f'{type} at {11/24000 * step:.2f} days ({step} steps) - grid {grid_id}', fontsize = 13)
    
    # save the figure
    if type=="Mmp2":
        figure_path = os.path.join(pathToSave, f'{type}-grid{grid_id}-step{step} - {11/24000 * step:.2f} days.png')
    elif type=="Ecm":
        figure_path = os.path.join(pathToSave, f'{type}-grid{grid_id}-step{step} - {11/24000 * step:.2f} days.png')
    plt.savefig(figure_path)

def plot_histogram(data_folder_path, histogram_images_path):
    csv_files_names = [f for f in os.listdir(data_folder_path) if os.path.isfile(os.path.join(data_folder_path, f)) and f.endswith(".csv") and "Histogram" in f]
    if csv_files_names:
        for csv_file in csv_files_names:
            histogram = pd.read_csv(os.path.join(data_folder_path, csv_file), header=0, index_col=0)
            if histogram.empty:
                continue
            plt.bar(histogram['Bins'], histogram['Frequency'])
            plt.xticks(range(carrying_capacity + 1))
            plt.xlim([-1, carrying_capacity + 1])
            path_to_save = os.path.join(histogram_images_path, csv_file[0:-4]+".png")
            plt.savefig(path_to_save)
    else:
        raise Exception(f'No CSV files found at {data_folder_path}. Did you run the data analysis first?')

def plotVasculatureGraphs(vasculature_df, pathToSave, max_step):
    # Prepare the data for the bar chart
    mesenchymal_data = vasculature_df["Mesenchymal cells"]
    epithelial_data = vasculature_df["Epithelial cells"]
    total_cluster_data = vasculature_df["Total clusters"]
    multicellular_cluster_data = vasculature_df["Multicellular clusters"]
    time_steps = vasculature_df["Time"]
    plt.style.use("Solarize_Light2")

    #second plot, clusters
    plt.figure()
    
    # Plot the data for each category
    plt.step(time_steps, mesenchymal_data, where='mid', color='tab:blue', label='Mesenchymal Cells', linewidth=1.5)
    plt.step(time_steps, epithelial_data, where='mid', color='tab:orange', label='Epithelial Cells', linewidth=1.5)

    # Set the chart labels, title, and legend
    plt.xlabel('Day')
    plt.ylabel('Cell Count')
    plt.title('Vasculature Cell Counts')
    plt.xticks([0, max_step/2, max_step], [0, f"{(max_step/2)*11/24000:.2f}", f"{max_step*11/24000:.2f}"])
    
    plt.xlim([-0.5, max_step + 0.5])
    plt.legend(loc="upper left")

    # save the figure
    figure_path = os.path.join(pathToSave, f'Vasculature-cells-step{max_step}.png')
    plt.savefig(figure_path)
    plt.close()
    #second plot, clusters
    plt.figure()
    
    # Plot the data for each category
    plt.step(time_steps, total_cluster_data, where='mid', color='darkred', label='Total clusters', linewidth=1.5)
    plt.step(time_steps, multicellular_cluster_data, where='mid', color='green', label='Multicellular clusters', linewidth=1.5)

    # Set the chart labels, title, and legend
    plt.xlabel('Day')
    plt.ylabel('Cluster Count')
    plt.title('Vasculature Cluster Counts')
    plt.xticks([0, max_step/2, max_step], [0, f"{(max_step/2)*11/24000:.2f}", f"{max_step*11/24000:.2f}"])
    

    plt.xlim([-0.5, max_step + 0.5])
    plt.legend(loc="upper left")

    # save the figure
    figure_path = os.path.join(pathToSave, f'Vasculature-clusters-step{max_step}.png')
    plt.savefig(figure_path)
    plt.close()
    
def generate_graphs(nameOfTheSimulation):
    simulations_dir = "Simulations"
    simulation_path = os.path.join(simulations_dir, nameOfTheSimulation)
    print(f'\tAnalyzing data in the folder {simulation_path}\n')

    # Get the agents' data filename 
    csv_files_name = [f for f in os.listdir(simulation_path) if os.path.isfile(os.path.join(simulation_path, f)) and f.endswith(".csv")]
    
    # Get the Ecm and Mmp2 data filenames 
    EcmPath = os.path.join(simulation_path, "Ecm")
    Mmp2Path = os.path.join(simulation_path, "Mmp2")
    ecm_files_name = [f for f in sorted(os.listdir(EcmPath), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(EcmPath, f)) and f.endswith(".csv")]
    mmp2_files_name = [f for f in sorted(os.listdir(Mmp2Path), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(Mmp2Path, f)) and f.endswith(".csv")]

    # Get the vasculature data filename
    VasculaturePath = os.path.join(simulation_path, "Vasculature")
    vasculature_files_name = [f for f in os.listdir(VasculaturePath) if os.path.isfile(os.path.join(VasculaturePath, f)) and f.endswith(".json")]

    if csv_files_name:
        first_csv_name = csv_files_name[0]
        first_csv_path = os.path.join(simulation_path, first_csv_name)
        print("Using cells data at:", first_csv_path)
    else:
        print("No .csv cell data found in directory:", simulation_path)
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
        print("No .json vasculature data found in directory:", simulation_path)
        return


    #Path of the data
    dataFolder = "Data analysis"
    dataPath = os.path.join(simulation_path, dataFolder)

    TumorDataPath = os.path.join(dataPath, "Tumor growth")
    CellsDataPath = os.path.join(dataPath, "Cells growth")
    EcmDataPath = os.path.join(dataPath, "Ecm evolution")
    Mmp2DataPath = os.path.join(dataPath, "Mmp2 evolution")
    VasculatureDataPath = os.path.join(dataPath, "Vasculature evolution")


    # use regex to find the values before 'steps', 'stepsize' and 'grids'. Ex: 50steps-10stepsize-cells
    max_step = int(re.search(r"(\d+)steps", first_csv_name).group(1))
    step_size = int(re.search(r"(\d+)stepsize", first_csv_name).group(1))
    grids_number = int(re.search(r"(\d+)grids", first_csv_name).group(1))

    # Path to save all the images:
    imagesFolder = "Visual analysis"
    imagesPath = os.path.join(simulation_path, imagesFolder)

    TumorImagesPath = os.path.join(imagesPath, "Tumor growth")
    CellsImagesPath = os.path.join(imagesPath, "Cells growth")
    EcmImagesPath = os.path.join(imagesPath, "Ecm evolution")
    Mmp2ImagesPath = os.path.join(imagesPath, "Mmp2 evolution")
    vasculature_images_path = os.path.join(imagesPath, "Vasculature evolution")
    histogram_images_path = os.path.join(imagesPath, "Positions histogram")

    # Create folder for all the visual analysis
    if not os.path.exists(imagesPath):
        os.makedirs(imagesPath)
        os.makedirs(TumorImagesPath)
        os.makedirs(Mmp2ImagesPath)
        os.makedirs(EcmImagesPath)
        os.makedirs(CellsImagesPath)
        os.makedirs(vasculature_images_path)
        os.makedirs(histogram_images_path)

        print(f"\nSaving tumor images in the folder:", TumorImagesPath)
        print("Saving Mmp2 images in the folder:", Mmp2ImagesPath)
        print("Saving Ecm images in the folder:", EcmImagesPath)
        print("Saving cells numbers images in the folder:", CellsImagesPath)
        print("Saving vasculature images in the folder:", vasculature_images_path)
        print("Saving histogram images in the folder:", histogram_images_path)

    # If there the visual analysis is already done, tell the user
    else:
        print("This visual analysis already exists!")
        return
    
    plt.style.use("Solarize_Light2")
    figCounter = 1
    for grid_id in range(1, grids_number+1):
        print(f'\nGrid: {grid_id}')

        # Plot the cells graphs
        print(f'\tPlotting tumor graphs...')
        for id, step in enumerate(range(1,max_step+1,step_size)):
            plotCancer(getCoordsForPlot(step, first_csv_path, grid_id), figCounter, imagesFolder, grid_id, step, TumorImagesPath)
            plt.close()
            figCounter += 1

        # Plot the histogram graphs
        print(f'\tPlotting histogram graphs...')
        for id, step in enumerate(range(1,max_step+1,step_size)):
            plot_histogram(TumorDataPath, histogram_images_path)
            plt.close()
            figCounter += 1

        # Plot the Mmp2 graphs
        print(f'\tPlotting Mmp2 graphs...')
        for id, step in enumerate(range(1,max_step+1,step_size)):
            mmp2_files_path_this_grid = [path for path in mmp2_files_path if f"Mmp2-{grid_id}grid-" in path]
            plotMMP2orECM(id, step, mmp2_files_path_this_grid, figCounter, grid_id, Mmp2ImagesPath, type="Mmp2")
            plt.close()
            figCounter += 1

        # Plot the Ecm graphs
        print(f'\tPlotting Ecm graphs...')
        for id, step in enumerate(range(1,max_step+1,step_size)):
            ecm_files_path_this_grid = [path for path in ecm_files_path if f"Ecm-{grid_id}grid-" in path]
            plotMMP2orECM(id, step, ecm_files_path_this_grid, figCounter, grid_id, EcmImagesPath, type="Ecm")
            plt.close()
            figCounter += 1

        # Plot the growth of ephitelial and mesenchymal cells
        print(f'\tPlotting cells numbers graph...')
        for id, step in enumerate(range(1,max_step+1,step_size)):
            plotGrowthData(figCounter, first_csv_path, step_size, grid_id , CellsImagesPath, step)
            plt.close()
            figCounter += 1

    # Plot the vasculature data
    print(f'Plotting vasculature...')
    for id, step in enumerate(range(0,max_step+1,step_size)):
        folder_path = os.path.join(simulation_path, "Data analysis", "Vasculature evolution")
        try:
            file_path = os.path.join(folder_path, f"Vasculature-step{step}.csv")
            vasculature_data = pd.read_csv(file_path, header=0)
        except:
            print(f"Error while reading the vasculature data for step {step} in grid {grid_id}", file=sys.stderr)
            print("Did you run the 'Data analysis' in the postprocessing menu first?", file=sys.stderr)
            os._exit(1)
        plotVasculatureGraphs(vasculature_data, vasculature_images_path, step)
        plt.close()
        figCounter += 1

    plt.show()
    plt.close()


###########################################################################################################

if __name__ == "__main__":

    # CHANGE THIS LINE according to the simulation you want to plot the graphs
    nameOfTheSimulation = "Sim maxSteps-1100 stepsize-100 N-388 gridsNumber-3"

    # This runs all the code to generate the graphs in the folder
    generate_graphs(nameOfTheSimulation)







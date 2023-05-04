import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Classes.utils import parent_dir, imagesFolder_utils
from Classes.utils import gridsize_utils as gridsize
import re
import os

# To run this code you must be in the parent folder of agent-based-cancer

def getCoordsForPlot(step, allCellsCsvPath, grid_id):
    df = pd.read_csv(allCellsCsvPath)
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

    
    if step == 0:
        plt.title(f'Initial Tumor at {11/24000 * step:.2f} days ({step} steps) - grid {grid_id}')
        plt.legend(loc="upper left")
    else:
        plt.title(f'Tumor size at {11/24000 * step:.2f} days ({step} steps) - grid {grid_id}')

    # save the figure
    figure_path = os.path.join(TumorImagesPath, f'Cells-grid{grid_id}-step{step} - Tumor size at {11/24000 * step:.2f} days.png')
    plt.savefig(figure_path)


def plotGrowthData(fig_index, allCellsCsvPath, imagesFolder, stepsize, grid_id, CellsImagesPath):
    # get data at each step
    df = pd.read_csv(allCellsCsvPath)
    df = df[["Step", "Total cells", "Phenotype", "Grid"]]
    df = df.loc[df["Grid"] == grid_id]
    plt.figure(fig_index)

    # For mesenchymal
    #stepsize = 3000 doubling rate, since it is different it will stay horizontal sometimes
    df_m = df.loc[(df["Step"] % stepsize == 0)]
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_m["Step"]) + 1, stepsize)
    numberMesenchymalEachStep = [df_m.loc[(df_m["Step"] == step) & (df_m["Phenotype"] == "mesenchymal")].shape[0] for step in steps]
    plt.plot(steps*11/24000, numberMesenchymalEachStep, label="Mesenchymal cells")

    # For epithelial
    #stepsize = 2000 doubling rate
    df_e = df.loc[(df["Step"] % stepsize == 0)]
    # create arrays with the step number and number of cells
    steps = np.arange(0, max(df_e["Step"]) + 1, stepsize)
    numberEpithelialEachStep = [df_e.loc[(df_e["Step"] == step) & (df_e["Phenotype"] == "epithelial")].shape[0] for step in steps]
    
    plt.plot(steps*11/24000, numberEpithelialEachStep, label="Epithelial cells")

    plt.title(f'Cells growth - grid {grid_id}')
    plt.xlabel('Days')
    plt.ylabel('Number of cells')
    plt.legend(loc="upper left")

    # save the figure
    figure_path = os.path.join(CellsImagesPath, f'CellsGrowth-grid{grid_id}.png')
    plt.savefig(figure_path)
    


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

    plt.title(f'{type} at {11/24000 * step:.2f} days ({step} steps) - grid {grid_id}')
    
    # save the figure
    if type=="Mmp2":
        figure_path = os.path.join(pathToSave, f'{type}-grid{grid_id}-step{step} - {11/24000 * step:.2f} days.png')
    elif type=="Ecm":
        figure_path = os.path.join(pathToSave, f'{type}-grid{grid_id}-step{step} - {11/24000 * step:.2f} days.png')
    plt.savefig(figure_path)



def main_graphs():
    SimulationPath = os.path.join(parent_dir, nameOfTheSimulation)
    print(f'\tAnalyzing data in the folder {SimulationPath}')

    csv_files_name = [f for f in os.listdir(SimulationPath) if os.path.isfile(os.path.join(SimulationPath, f)) and f.endswith(".csv")]
    
    EcmPath = os.path.join(SimulationPath, "Ecm")
    Mmp2Path = os.path.join(SimulationPath, "Mmp2")
    ecm_files_name = [f for f in sorted(os.listdir(EcmPath), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(EcmPath, f)) and f.endswith(".csv")]
    mmp2_files_name = [f for f in sorted(os.listdir(Mmp2Path), key=lambda x: int(re.findall(r'\d+(?=step)', x)[0])) if os.path.isfile(os.path.join(Mmp2Path, f)) and f.endswith(".csv")]


    if csv_files_name:
        first_csv_name = csv_files_name[0]
        first_csv_path = os.path.join(SimulationPath, first_csv_name)
        all_cells_analysis_path = os.path.join(parent_dir, SimulationPath, first_csv_path)
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


    # use regex to find the values before 'steps', 'stepsize' and 'grids'. Ex: 50steps-10stepsize-cells
    max_step = int(re.search(r"(\d+)steps", first_csv_name).group(1))
    step_size = int(re.search(r"(\d+)stepsize", first_csv_name).group(1))
    grids_number = int(re.search(r"(\d+)grids", first_csv_name).group(1))

    # Path to save all the images:
    imagesFolder = imagesFolder_utils
    imagesPath = os.path.join(SimulationPath, imagesFolder)

    TumorImagesPath = os.path.join(imagesPath, "Tumor growth")
    CellsImagesPath = os.path.join(imagesPath, "Cells growth")
    EcmImagesPath = os.path.join(imagesPath, "Ecm evolution")
    Mmp2ImagesPath = os.path.join(imagesPath, "Mmp2 evolution")

    # Create folder for all the visual analysis
    if not os.path.exists(imagesPath):
        os.makedirs(imagesPath)
        os.makedirs(TumorImagesPath)
        os.makedirs(Mmp2ImagesPath)
        os.makedirs(EcmImagesPath)
        os.makedirs(CellsImagesPath)

        print("Saving tumor images in the folder:", TumorImagesPath)
        print("Saving Mmp2 images in the folder:", Mmp2ImagesPath)
        print("Saving Ecm images in the folder:", EcmImagesPath)
        print("Saving cells numbers images in the folder:", CellsImagesPath)

    # If there the visual analysis is already done, tell the user
    else:
        print("This visual analysis already exists!")
        return
    
    plt.style.use("Solarize_Light2")
    figCounter = 1
    for grid_id in range(1, grids_number+1):
        print(f'\nGrid: {grid_id}')

        # Plot the cells graphs
        print(f'\tPlotting tumor graphs')
        for id, step in enumerate(range(0,max_step+1,step_size)):
            plotCancer(getCoordsForPlot(step, first_csv_path, grid_id), figCounter, imagesFolder, grid_id, step, TumorImagesPath)
            plt.close()
            figCounter += 1

        # Plot the Mmp2 graphs
        print(f'\tPlotting Mmp2 graphs')
        for id, step in enumerate(range(0,max_step+1,step_size)):
            mmp2_files_path_this_grid = [path for path in mmp2_files_path if f"Mmp2-{grid_id}grid-" in path]
            plotMMP2orECM(id, step, mmp2_files_path_this_grid, figCounter, grid_id, Mmp2ImagesPath, type="Mmp2")
            plt.close()
            figCounter += 1

        # Plot the Ecm graphs
        print(f'\tPlotting Ecm graphs')
        for id, step in enumerate(range(0,max_step+1,step_size)):
            ecm_files_path_this_grid = [path for path in ecm_files_path if f"Ecm-{grid_id}grid-" in path]
            plotMMP2orECM(id, step, ecm_files_path_this_grid, figCounter, grid_id, EcmImagesPath, type="Ecm")
            plt.close()
            figCounter += 1

        # Plot the growth of ephitelial and mesenchymal cells
        print(f'\tPlotting cells number')
        plotGrowthData(figCounter, first_csv_path, imagesFolder, step_size, grid_id , CellsImagesPath)
        plt.close()
        figCounter += 1

    plt.show()
    plt.close()


###########################################################################################################

if __name__ == "__main__":

    # CHANGE THIS LINE according to the simulation you want to plot the graphs
    nameOfTheSimulation = "Sim maxSteps-30 stepsize-10 N-388 gridsNumber-3"

    # This runs all the code to generate the graphs in the folder
    main_graphs()







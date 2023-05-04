import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Classes.CancerModel import *
import os
from Classes.utils import parent_dir, gridsize_utils
from Classes import utils

# To run this code you must be in the parent folder of agent-based-cancer
# Before running this code always check if isBatchRun = True is in Utils
# otherwise it won't save the Mmp2 and Ecm data

# Parameters for this simulation
maxSteps = 2000
dataCollectionPeriod = 10
N = 388 # Number of cancer cells
gridsize     = gridsize_utils #PDF: 201
width        = gridsize
height       = gridsize
grids_number = 3

# Name of the directories
simulations_dir = parent_dir
newSimulationFolder = f"Sim maxSteps-{maxSteps} stepsize-{dataCollectionPeriod} N-{N} gridsNumber-{grids_number}"


def run_simulation():

    # Setting parameters for mesa.batch_run
    params = {"N": N, "width": width, "height": height, "grids_number": grids_number}


    # Saves the simulation configuration
    var_names = dir(utils) + ['maxSteps', 'dataCollectionPeriod', 'N', 'grids_number', 'simulations_dir']

    # Open the file for writing
    print(f"Saving all the simulations parameters at: {os.path.join(simulations_dir, newSimulationFolder, 'configs.txt')}")
    with open(os.path.join(simulations_dir, newSimulationFolder, 'configs.txt'), 'w') as f:
        # Write the values of each variable to the file
        for name, value in globals().items():
            if name in var_names:
                f.write(f'{name}: {value}\n')


    # The whole simulation occurs here
    print(f'\n\tStarting the simulation')
    results = mesa.batch_run(
        CancerModel,
        parameters=params,
        iterations=1, # I think for the code to run and save in the specified folders this must be always 1
        max_steps=maxSteps,
        number_processes=1,
        data_collection_period=dataCollectionPeriod,
        display_progress=True,
    )

    # Create data frames for the cells
    cells_df = pd.DataFrame(results)
    print(f'Example of data collected: {cells_df.head(10)}')

    # Saves data analysis
    nameOfCsv = f'CellsData-{maxSteps}steps-{dataCollectionPeriod}stepsize-{grids_number}grids.csv'
    pathToSave = os.path.join(simulations_dir, newSimulationFolder, nameOfCsv)
    cells_df.to_csv(pathToSave)
    print(f'All data saved')



def main_Batch():

    # Create parent directory if it doesn't exist
    if not os.path.exists(simulations_dir):
        print(f'Creating the folder for this simulation at: {simulations_dir}')
        os.makedirs(simulations_dir)

    # Creates the path for the new simulation
    path = os.path.join(simulations_dir, newSimulationFolder)
    pathMmp2 = os.path.join(path, "Mmp2")
    pathEcm = os.path.join(path, "Ecm")

    # Create folder for all cells analysis, for Mmp2 matrices and Ecm matrices
    if not os.path.exists(path):
        print(f'\t Folder for this simulation: {path}')
        print(f'\t Saving agents data at: {path}')
        print(f'\t Saving Mmp2 data at: {pathMmp2}')
        print(f'\t Folder Ecm data at: {pathEcm}')

        os.makedirs(path)
        os.makedirs(pathMmp2)
        os.makedirs(pathEcm)

        # Run the simulation and saves the data
        run_simulation()

    # If there is already a simulation you skip it
    else:
        return print("This simulation already exists!")

    return print('Finished the simulation')


if __name__ == "__main__":

    # This runs all the code to create folders, run the simulation and save the data
    main_Batch()

    




# Error with the size of the array beeing too big -->

#numpy.core._exceptions._ArrayMemoryError: Unable to allocate 3.61 GiB for an array
#with shape (24000, 201, 201) and data type int32


# Error in proliferation function --> problem with Id
# Erro:   File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 75, in step
#    self.proliferate("mesenchymal")
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 89, in proliferate
#    self.schedule.add(new_cell)
#  File "C:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\.venv\lib\site-packages\mesa\time.py", line 68, in add
#    raise Exception(
#Exception: Agent with unique id 4809 already


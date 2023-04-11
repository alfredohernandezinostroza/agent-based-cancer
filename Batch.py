import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Classes.CancerModel import *
import os
from Classes.utils import parent_dir, gridsize_utils

# To run this code you must be in the parent folder of agent-based-cancer
# Before running this code always check if isBatchRun = True is in Utils
# otherwise it won't save the Mmp2 and Ecm data

# Parameters for this simulation
maxSteps = 12001
dataCollectionPeriod = 2000

N = 388 # Number of cancer cells
gridsize     = gridsize_utils #PDF: 201
width        = gridsize
height       = gridsize
grids_number = 2

# Name of the directories
simulations_dir = parent_dir
newSimulationFolder = f"Sim maxSteps-{maxSteps} stepsize-{dataCollectionPeriod} N-{N} gridsNumber-{grids_number}"

def main():

    params = {"N": N, "width": width, "height": height, "grids_number": grids_number} # N=388 pdf

    # The whole simulation occurs here
    results = mesa.batch_run(
        CancerModel,
        parameters=params,
        iterations=1, # I think for the code to run and save in the specified folders this must be always 1
        max_steps=maxSteps,
        number_processes=1,
        data_collection_period=dataCollectionPeriod,
        display_progress=True,
    )

    cells_df = pd.DataFrame(results)
    # Create data frames for the cells
    print(cells_df)

    nameOfCsv = f'{maxSteps}steps-{dataCollectionPeriod}stepsize-cells.csv'
    pathToSave = os.path.join(simulations_dir, newSimulationFolder, nameOfCsv)
    cells_df.to_csv(pathToSave)
    
    #print('cells_df df')
    #print(cells_df.head())
    #print(cells_df.shape)


if __name__ == "__main__":

    # Create parent directory if it doesn't exist
    if not os.path.exists(simulations_dir):
        os.makedirs(simulations_dir)

    # Creates the path for the new simulation
    path = os.path.join(simulations_dir, newSimulationFolder)
    pathMmp2 = os.path.join(path, "Mmp2")
    pathEcm = os.path.join(path, "Ecm")

    # Create folder for all cells analysis, for Mmp2 matrices and Ecm matrices
    if not os.path.exists(path):
        os.makedirs(path)
        os.makedirs(pathMmp2)
        os.makedirs(pathEcm)
        main()

    # If there is already a simulation you skip it
    else:
        print("This simulation already exists!")




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


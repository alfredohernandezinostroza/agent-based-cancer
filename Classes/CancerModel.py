import mesa
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import os
import json
import ast
from Classes.CancerCell import CancerCell
from Classes.Vessel import Vessel
from Classes.QuasiCircle import find_quasi_circle
from matplotlib import pyplot as plt
from matplotlib import cm
# from Classes.configs import *
import Classes.configs



def get_cluster_survival_probability(cluster):
    """
    Takes in a tuple representing a cluster, returns the survival probabiltiy.
    
    Input:
        Cluster: a tuple representing the cancer cells cluster, where the first 
        element corresponds to the amount of Mesenchymal cells, and the second
        corresponds to the amount of Epithelial cells.
    Returns:
        The corresponding survival probability, according to if its a songle-cell
        cluster or a multi-cellular one.
    """
    if cluster[0] < 0:
        raise Exception(f"Error! Mesenchymal cells are negative: {cluster[0]}")
    if cluster[1] < 0:
        raise Exception(f"Error! Epithelial cells are negative: {cluster[1]}")
    if sum(cluster) == 1:
        return (single_cell_survival)
    elif sum(cluster) > 1:
        return (cluster_survival)
    elif sum(cluster) == 0:
        raise Exception(f"Error, no cells in cluster!")
    else:
        raise Exception(f"Error, nothing returned for cluster survival probability" )
    

def count_total_cells(model):
    """Counts all the cells present in the model, in all sites.

    Input:
        model: CancerModel object.
    Returns:
        amount_of_cells (int): the total amount of cells in every site,
        NOT considering the vasculature    
    """
    amount_of_cells = len([1 for agent in model.schedule.agents if agent.agent_type == "cell"])
    # print(amount_of_cells)
    return amount_of_cells

def count_vasculature_cells(model):
    """"
    Counts the total amount of cells in the vasculature

    Input: 
        model: CancerModel object.
    Returns:
        amount_of_cells (int): the total amount of cells in the vasculature
    """
    amount_of_cells = sum([len(value) for value in model.vasculature.values()])
    return amount_of_cells

class CancerModel(mesa.Model):
    """
    Class for the model.

    Attributes:
    ---------------
    N: int
        number of initial cancer cells
    width: int
        the width of each grid
    height: int
        the height of each grid
    grids_number: int
        the amount of sites, consideiring the intial tumour site + the secondary sites
    seed: int
        the seed used for the random number generation for all the simulation.
        If None, the random one will be selected by default.

    Methods:
    ---------------
    _initialize_grids__()
        initialize the grid with the initial bessel and cancer cell population
    proliferate(cellType)
        Duplicates every cancer cell in the model of the cellType phenotype
    calculateEnvironments(mmp2, ecm)
        Calculates the next step for the given arrays of mmp2 and ecm concentrations
    disaggregate_clusters(time)
        For a given time, it will dissagregate single cells from clusters
    """

    def __init__(self, N, width, height, grids_number, maxSteps, dataCollectionPeriod, newSimulationFolder, loadedSimulationPath="", seed=None):
        super().__init__()  
        self.simulations_dir = "Simulations"
        self.vasculature = {}
        self.num_agents = N
        self.width = width
        self.height = height
        self.phenotypes = ["mesenchymal", "epithelial"]
        self.grid_vessels_positions = [[],[],[]]
        self.current_agent_id = 0
        self.maxSteps = maxSteps
        self.dataCollectionPeriod = dataCollectionPeriod
        self.newSimulationFolder  = newSimulationFolder 
        self.mesenchymalCount = [np.zeros((width, height), dtype=float) for _ in range(grids_number)]
        self.epithelialCount = [np.zeros((width, height), dtype=float) for _ in range(grids_number)]
        self.grids_number = grids_number
        self.grids = [mesa.space.MultiGrid(width, height, False) for _ in range(self.grids_number)]
        self.grid_ids = [i+1 for i in range(self.grids_number)] # need a number to appear in the data analysis (.csv)
        self.time_grid_got_populated = [-1 for _ in range(self.grids_number)]
        self.schedule = mesa.time.RandomActivation(self)
        #list of numpy arrays, representing mmp2 and ecm concentration in each grid
        self.mmp2 = [np.zeros((2, width, height), dtype=float) for _ in range(grids_number)]
        self.ecm = [np.ones((2, width, height), dtype=float) for _ in range(grids_number)]

        if loadedSimulationPath != "":
            # configs_path = os.path.join(loadedSimulationPath, "configs.csv")
            # config_var_names = configs.load_simulation_configs(configs_path)
            # for var in config_var_names:
            #     globals()[var] = getattr(configs, var)
            #     self.load_previous_simulation(loadedSimulationPath)
            self.load_previous_simulation(loadedSimulationPath)
            configs_path = os.path.join(loadedSimulationPath, "configs.csv")
            config_var_names = Classes.configs.load_simulation_configs_for_reloaded_simulation(configs_path)
            for var in config_var_names:
                globals()[var] = getattr(Classes.configs, var)
        else:
            # load_simulation_configs("simulations_configs.csv")

            configs_path = "simulations_configs.csv"
            config_var_names = Classes.configs.init_simulation_configs(configs_path)
            for var in config_var_names:
                globals()[var] = getattr(Classes.configs, var)
            self._initialize_grids()
        self.datacollector = mesa.DataCollector(
            model_reporters={"Total cells": count_total_cells}, agent_reporters={"Position": "pos", "Agent Type": "agent_type", "Phenotype": "phenotype", "Ruptured": "ruptured", "Grid": "grid_id"})
            # model_reporters={"Total cells": count_total_cells, "Cluster radius and diameter": get_cluster_radius_and_diameter,
            #                   "Amount of vasculature cells": count_vasculature_cells,}, agent_reporters={"Position": "pos", "Agent Type": "agent_type", "Phenotype": "phenotype", "Ruptured": "ruptured", "Grid": "grid_id"})
        #model_reporters={"Mmp2": "mmp2", "Grid": "grid"},

    def step(self):
        """Advance the model by one step.
        
        The step function of the model will be called on each of the siulations steps.
        It will correctly allocate the cells coming from the vaculature, if any.
        Then, it will calculate the ECM and MMP2 concentration changes for this step,
        proliferate the cells if its due, and collect the data if it is the appropiate time.

        Input: none
        Returns: none
        """        
        # print("=========================================")
        if self.schedule.time in self.vasculature: # Add keys
            self.disaggregate_clusters(self.schedule.time)
            surviving_clusters = [cluster for cluster in self.vasculature[self.schedule.time] if self.random.random() < get_cluster_survival_probability(cluster)]
            del self.vasculature[self.schedule.time]
            for cluster in surviving_clusters:
                selected_site = self.random.choices(range(1,self.grids_number), weights=extravasation_probs[0:self.grids_number-1])[0]
                arriving_point = self.random.choice(self.grid_vessels_positions[selected_site])
                x,y = arriving_point
                onLeftBorder    = self.grids[selected_site].out_of_bounds((x-1,y))
                onRightBorder   = self.grids[selected_site].out_of_bounds((x+1,y))
                onTopBorder     = self.grids[selected_site].out_of_bounds((x,y+1))
                onBottomBorder  = self.grids[selected_site].out_of_bounds((x,y-1))
                possible_places = self.grids[selected_site].get_neighborhood(arriving_point, moore=False, include_center=False)
                number_of_ccells_in_arriving_point ={}
                for x2,y2 in possible_places:
                    number_of_ccells_in_arriving_point[x2,y2] = len([agent for agent in self.grids[selected_site].get_cell_list_contents([(x2,y2)]) if agent.agent_type == "cell"])
                for tuple_index, ccells_amount in enumerate(cluster):
                    cell_type = "mesenchymal" if tuple_index == 0 else "epithelial"
                    while ccells_amount > 0:
                        if not onLeftBorder and carrying_capacity > number_of_ccells_in_arriving_point[x-1,y]:
                            ccell = CancerCell(self.current_agent_id, self, self.grids[selected_site], self.grid_ids[selected_site], cell_type, self.ecm[selected_site], self.mmp2[selected_site])
                            self.current_agent_id += 1
                            self.grids[selected_site].place_agent(ccell, (x-1,y)) 
                            number_of_ccells_in_arriving_point[x-1,y] += 1
                            self.schedule.add(ccell)
                        elif not onRightBorder and carrying_capacity > number_of_ccells_in_arriving_point[x+1,y]:
                            ccell = CancerCell(self.current_agent_id, self, self.grids[selected_site], self.grid_ids[selected_site], cell_type, self.ecm[selected_site], self.mmp2[selected_site])
                            self.current_agent_id += 1
                            self.grids[selected_site].place_agent(ccell, (x+1,y))
                            number_of_ccells_in_arriving_point[x+1,y] += 1
                            self.schedule.add(ccell)
                        elif not onBottomBorder and carrying_capacity > number_of_ccells_in_arriving_point[x,y-1]:
                            ccell = CancerCell(self.current_agent_id, self, self.grids[selected_site], self.grid_ids[selected_site], cell_type, self.ecm[selected_site], self.mmp2[selected_site])
                            self.current_agent_id += 1
                            self.grids[selected_site].place_agent(ccell, (x,y-1))
                            number_of_ccells_in_arriving_point[x,y-1] += 1
                            self.schedule.add(ccell)
                        elif not onTopBorder and carrying_capacity > number_of_ccells_in_arriving_point[x,y+1]:
                            ccell = CancerCell(self.current_agent_id, self, self.grids[selected_site], self.grid_ids[selected_site], cell_type, self.ecm[selected_site], self.mmp2[selected_site])
                            self.current_agent_id += 1
                            self.grids[selected_site].place_agent(ccell, (x,y+1))
                            number_of_ccells_in_arriving_point[x,y+1] += 1
                            self.schedule.add(ccell)
                        ccells_amount -= 1
                    
        #Calculo do quimico que fomenta haptotaxis e da matriz extracelular
        self.calculateEnvironment(self.mmp2, self.ecm)
        
        # Reprodução
        if (self.schedule.time % doubling_time_M == 0 and self.schedule.time != 0):
            # print(f"Before proliferating: {count_total_cells(self)}")
            self.proliferate("mesenchymal")
            # print(f"After proliferating: {count_total_cells(self)}")

        if (self.schedule.time % doubling_time_E == 0 and self.schedule.time != 0):
            self.proliferate("epithelial")

        print(f'step number: {self.schedule.time}')
        self.schedule.step()
        self.datacollector.collect(self)

        #at the end of each step, check if the grid has been populated, and if it happened, store the time step when it did
        for index, time in enumerate(self.time_grid_got_populated):
            if time == -1: #if it has not been populated already, we check:
                cell_count = len([1 for agent in self.schedule.agents if agent.agent_type == "cell" and (agent.grid_id - 1) == index])
                if cell_count > 0:
                    self.time_grid_got_populated[index] = self.schedule.time

        # Saving of non agents data
        if (self.schedule.time != 0 and (self.schedule.time % self.dataCollectionPeriod == 0)) \
            or self.schedule.time == self.maxSteps:
            df_time_grids_got_populated = pd.DataFrame()
            for grid_id in self.grid_ids:
                new_mmp2_df = pd.DataFrame(self.mmp2[grid_id-1][0,:,:])
                mmp2CsvName = f"Mmp2-{grid_id}grid-{self.schedule.time}step.csv"
                pathToSave = os.path.join(self.simulations_dir, self.newSimulationFolder, "Mmp2", mmp2CsvName)
                new_mmp2_df.to_csv(pathToSave)

                new_ecm_df = pd.DataFrame(self.ecm[grid_id-1][0,:,:])
                EcmCsvName = f"Ecm-{grid_id}grid-{self.schedule.time}step.csv"
                pathToSave = os.path.join(self.simulations_dir, self.newSimulationFolder, "Ecm", EcmCsvName)
                new_ecm_df.to_csv(pathToSave)

                df_time_grids_got_populated[f"Time when grid {grid_id} was first populated"] = [self.time_grid_got_populated[grid_id-1]]
                df_time_grids_got_populated_csv_name = f"Cells-are-present-grid-{grid_id}-{self.schedule.time}step.csv"
            pathToSave = os.path.join(self.simulations_dir, self.newSimulationFolder, "Time when grids were populated", df_time_grids_got_populated_csv_name)
            df_time_grids_got_populated.to_csv(pathToSave)

            # Saves vasculature data
            vasculature_json = json.dumps(self.vasculature)
            # {key: list of clusters} -> {timestep: [(number of Mcells, number of Ecells), ..., (..., ...)]}
            
            vasculatureJsonName = f"Vasculature-{self.schedule.time}step.json"
            pathToSave = os.path.join(self.simulations_dir, self.newSimulationFolder, "Vasculature", vasculatureJsonName)
            
            with open(pathToSave, 'w') as f:
                f.write(vasculature_json)
                
            # Saves cancer cells data as a backup in case the simulation fails
            _, current_model_data = mesa.batchrunner._collect_data(self, self.dataCollectionPeriod-1)
            df_current_model_data = pd.DataFrame(current_model_data)
            df_current_model_data["Step"] = self.dataCollectionPeriod
            for step in range(self.dataCollectionPeriod * 2, self.schedule.time, self.dataCollectionPeriod):
                _, step_model_data = mesa.batchrunner._collect_data(self, step-1)
                df_step_model_data = pd.DataFrame(step_model_data)
                df_step_model_data["Step"] = step
                df_current_model_data = pd.concat([df_current_model_data, df_step_model_data])
            pathToSave = os.path.join(self.simulations_dir, self.newSimulationFolder, f'CellsData.csv')
            df_current_model_data.to_csv(pathToSave)



    def proliferate(self, cellType):
        """"
        Duplicates every cell of cellType phenotype in every site of the model

        Input: none
        Returns: none
        """
        all_agents = [agent for agent in self.schedule.agents]
        total_amount_of_agents = len(all_agents)
        for agent in all_agents:
            if agent.agent_type == "cell":
                x, y = agent.pos
                amount_of_cells = len([cell for cell in agent.grid.get_cell_list_contents([(x, y)]) if cell.agent_type == "cell"])
                if carrying_capacity > amount_of_cells and agent.phenotype == cellType:
                    # print("Created new cell!!")
                    new_cell = CancerCell(self.current_agent_id, self, agent.grid, agent.grid_id, agent.phenotype, agent.ecm, agent.mmp2)
                    self.current_agent_id += 1
                    self.schedule.add(new_cell)
                    agent.grid.place_agent(new_cell, (x,y))
                    total_amount_of_agents +=1
        


    def load_previous_simulation(self, pathToSimulation):
        """
        Loads the last step of a previously computed simulation as the initial condition of this model

        Input: The simulation's path to be loaded
        Returns: none
        """
        
        #load mmp2 and ecm
        for grid_number in range(self.grids_number):
            mmp2_files_path = os.path.join(pathToSimulation, "Mmp2")
            ecm_files_path  = os.path.join(pathToSimulation, "Ecm")
            mmp2_files = os.listdir(mmp2_files_path)
            ecm_files  = os.listdir(ecm_files_path)
            mmp2_files.sort(key = lambda file_name: int(file_name.split('step')[0][11:]))
            ecm_files.sort(key  = lambda file_name: int(file_name.split('step')[0][10:]))
            last_state_of_mmp2_filepath = os.path.join(mmp2_files_path,mmp2_files[-1])
            last_state_of_ecm_filepath  = os.path.join(ecm_files_path,ecm_files[-1])
            self.ecm[grid_number][0,:,:]  = pd.read_csv(last_state_of_ecm_filepath, index_col=0).to_numpy(dtype=float)
            self.mmp2[grid_number][0,:,:] = pd.read_csv(last_state_of_mmp2_filepath, index_col=0).to_numpy(dtype=float)

        path = os.path.join(pathToSimulation, "CellsData.csv")
        previous_sim_df = pd.read_csv(path, converters={"Position": ast.literal_eval})
        last_step = previous_sim_df["Step"].max()
        previous_sim_df = previous_sim_df[previous_sim_df["Step"] == last_step]
        last_step_cells = previous_sim_df[previous_sim_df["Agent Type"] == "cell"]
        last_step_vessels = previous_sim_df[previous_sim_df["Agent Type"] == "vessel"]
        for index, row in last_step_cells.iterrows():
            grid = int(row["Grid"]) - 1
            ccell = CancerCell(self.current_agent_id, self, self.grids[grid], self.grid_ids[grid], row["Phenotype"], self.ecm[grid], self.mmp2[grid])
            self.current_agent_id += 1
            self.schedule.add(ccell)
            self.grids[grid].place_agent(ccell, row["Position"])
        for index, row in last_step_vessels.iterrows():
            grid = int(row["Grid"]) - 1
            ruptured_state = bool(row["Ruptured"])
            vessel = Vessel(self.current_agent_id, self, ruptured_state, self.grids[grid], self.grid_ids[grid])
            self.current_agent_id += 1
            self.schedule.add(vessel)
            self.grids[grid].place_agent(vessel, row["Position"])

        #load vasculature
        vasculature_path = os.path.join(pathToSimulation, "Vasculature")
        vasculature_files = os.listdir(vasculature_path)
        vasculature_files.sort(key = lambda file_name: int(file_name.split('step')[0][12:]))
        last_state_of_vasculature_filepath = os.path.join(vasculature_path,vasculature_files[-1])
        with open(last_state_of_vasculature_filepath, 'r') as f:
            last_state_of_vasculature = json.load(f)
        # Change keys to int
        last_state_of_vasculature = {int(k): v for k, v in last_state_of_vasculature.items()}
        self.vasculature = last_state_of_vasculature



    def _initialize_grids(self):
        """
        Places the initial cancer cell and vessel in the initial grid in a circle

        Input: none
        Returns: none
        """
        mesenchymal_number = round(self.num_agents * mesenchymal_proportion)
        possible_places = find_quasi_circle(n_center_points_for_tumor, self.width, self.height)[1]
        # Place all the agents in the quasi-circle area in the center of the grid
        for i in range(self.num_agents):
            if mesenchymal_number > 0:
                cell_type = "mesenchymal"
                mesenchymal_number -= 1
            elif mesenchymal_number == 0:
                cell_type = "epithelial"

            a = CancerCell(self.current_agent_id, self, self.grids[0], self.grid_ids[0], cell_type, self.ecm[0], self.mmp2[0])
            self.current_agent_id += 1
            j = self.random.randrange(len(possible_places))
            x = int(possible_places[j][0])
            y = int(possible_places[j][1])

            self.schedule.add(a)
            self.grids[0].place_agent(a, (x, y))

            # Remove the point after it has 4 cells
            possible_places[j][2] += 1
            if possible_places[j][2] == carrying_capacity:
                possible_places.pop(j)


        # Create agents at second grid
        amount_of_second_grid_CAcells=0
        for i in range(amount_of_second_grid_CAcells):
            a = CancerCell(self.current_agent_id, self, self.grids[1], self.grid_ids[1], "mesenchymal", self.ecm[1], self.mmp2[1])
            self.current_agent_id += 1
            self.schedule.add(a)
        
            # Add the agent to a random grid cell
            x = self.random.randrange(3,7)
            y = self.random.randrange(3,7)
            self.grids[1].place_agent(a, (x, y))

        # Create vessels
        num_normal_vessels = normal_vessels_primary
        num_ruptured_vessels = ruptured_vessels_primary
        numVesselsSecondary = 10
        numVesselsThird = 10 # just to test it, final code will not have 1 var to each grid

        # bad code, reduce number of for and make a counter to save the index to de put in each vessel
        # creates grid with 1 where vessels must not be placed
        not_possible_array = find_quasi_circle(n_center_points_for_Vessels, self.width, self.height)[0]
        not_possible_array[:2,:] = 1
        not_possible_array[-2:,:] = 1
        not_possible_array[:,:2] = 1
        not_possible_array[:,-2:] = 1
        possible_places = np.where(not_possible_array == 0)
        pos_coords = [list(tup) for tup in zip(possible_places[0], possible_places[1])]

        # range(2) for 1 secondary site
        #for i in range(2):
        for i in range(len(self.grids)):
            if i == 0: # primary grid
                temp = num_ruptured_vessels
                while temp > 0:
                    j = num_ruptured_vessels - temp
                    cell_to_place = [self.random.randrange(self.width), self.random.randrange(self.height)]
                    if cell_to_place in pos_coords:
                        a = Vessel(self.current_agent_id, self, True, self.grids[0], self.grid_ids[0])
                        self.current_agent_id += 1
                        self.schedule.add(a)
                        self.grids[0].place_agent(a, (int(cell_to_place[0]), int(cell_to_place[1])))
                        # tenho que adicionar a cruz de ruptured e remover 5 cells de pos coords
                        not_possible_array[cell_to_place[0], cell_to_place[1]] = 1
                        pos_coords.remove(cell_to_place)
                        temp -= 1

                temp = num_normal_vessels
                while temp > 0:
                    j = num_normal_vessels - temp
                    cell_to_place = [self.random.randrange(self.width), self.random.randrange(self.height)]
                    if cell_to_place in pos_coords:
                        a = Vessel(self.current_agent_id, self, False, self.grids[0], self.grid_ids[0])
                        self.current_agent_id += 1
                        self.schedule.add(a)
                        self.grids[0].place_agent(a, (int(cell_to_place[0]), int(cell_to_place[1])))

                        not_possible_array[cell_to_place[0], cell_to_place[1]] = 1
                        pos_coords.remove(cell_to_place)
                        temp -= 1
            elif i > 0: # secondary and third grid
                    for m in range(secondary_sites_vessels[i-1]):
                        # make if to only create a vessel if given random value of x and y doesnt already has a vessel
                        a = Vessel(self.current_agent_id, self, False, self.grids[i], self.grid_ids[i])
                        self.current_agent_id += 1
                        self.schedule.add(a)
                        x = self.random.randrange(self.width)
                        y = self.random.randrange(self.height)
                        self.grids[i].place_agent(a, (x,y))
                        self.grid_vessels_positions[i] += [(x,y)]
                
    def calculateEnvironment(self, mmp2, ecm):
        global th
        for i in range(len(mmp2)):
            for cell in self.grids[i].coord_iter():
                cell_contents, (x, y) = cell
                diff = 0
                self.mesenchymalCount[i][x,y] = 0
                self.epithelialCount[i][x,y] = 0
                for cancerCell in cell_contents:
                    if isinstance(cancerCell, CancerCell):
                        if cancerCell.phenotype == "mesenchymal":
                            self.mesenchymalCount[i][x,y] += 1
                            diff = dM
                        elif cancerCell.phenotype == "epithelial":
                            self.epithelialCount[i][x,y] += 1
                            diff = dE
                        else:
                            raise Exception("Unknown phenotype")
                onLeftBorder = self.grids[i].out_of_bounds((x-1,y))
                onRightBorder = self.grids[i].out_of_bounds((x+1,y))
                onTopBorder = self.grids[i].out_of_bounds((x,y-1))
                onBottomBorder = self.grids[i].out_of_bounds((x,y+1))
                mmp2[i][1,x,y]=dmmp*tha/xha**2*((mmp2[i][0,x+1,y] if not onRightBorder else mmp2[i][0,x-1,y])\
                        +(mmp2[i][0,x-1,y] if not onLeftBorder else mmp2[i][0,x+1,y])\
                        +(mmp2[i][0,x,y+1] if not onBottomBorder else mmp2[i][0,x,y-1])\
                        +(mmp2[i][0,x,y-1] if not onTopBorder else mmp2[i][0,x,y+1])\
                        )\
                        +mmp2[i][0,x,y]*(1-4*dmmp*tha/xha**2-th*Lambda)+tha*theta*self.mesenchymalCount[i][x,y]
                ecm[i][1,x,y] = ecm[i][0,x,y]*(1-tha*(gamma1*self.mesenchymalCount[i][x,y]+gamma2*mmp2[i][1,x,y]))
                if ecm[i][1,x,y] < 0:
                    print(f"<0 ecm in [i][1,{x},{y}] is {ecm[i][1,x,y]}")
                    print(".")
                if ecm[i][1,x,y] > 1:
                    print(f">1 ecm in [i][1,{x},{y}] is {ecm[i][1,x,y]}")
                    print(".")
            mmp2[i][0,:,:] = mmp2[i][1,:,:]
            ecm[i][0,:,:] = ecm[i][1,:,:]

    def disaggregate_clusters(self, time):
        """
        Dissagregates cells from clusters into single-cell clusters, according to
        the dissagregation probability

        Input:
            time: given time for which the dissagreggation will occur
        Returns: 
            None
        """
        big_clusters = [cluster for cluster in self.vasculature[time] if sum(cluster) > 1]
        new_vasculature = [cluster for cluster in self.vasculature[time] if sum(cluster) == 1]
        for cluster in big_clusters:
            new_mesenchymal, new_epithelial = cluster
            for ccell_type, ccells_amount in enumerate(cluster):
                for i in range(ccells_amount):
                    if self.random.random() > dissagreggation_prob:
                        if ccell_type == 0:
                            new_vasculature += [(1, 0)]
                            new_mesenchymal -= 1
                        if ccell_type == 1:
                            new_vasculature += [(0, 1)]
                            new_epithelial -= 1
            if new_mesenchymal + new_epithelial > 0:
                new_vasculature += [(new_mesenchymal,new_epithelial)]
        self.vasculature[time] = new_vasculature



#
    #def graph_ecm_mmp2(self, time):
    #    if self.schedule.time == time:
    #        print("SADSADDSA")
    #        fig = plt.figure(figsize=plt.figaspect(0.5))
    #        ax = fig.add_subplot(1, 2, 1, projection='3d')
    #        X = np.arange(0, self.width, 1)
    #        Y = np.arange(0, self.height,1)
    #        X, Y = np.meshgrid(X, Y)
    #        Z = self.mmp2[0][0, :, :]
    #        print(Z)
    #        # ax.scatter(X, Y, Z, marker='o')
    #        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
    #                    linewidth=0, antialiased=False)
    #        ax.set_zlim(-1.01, 1.01)
    #        fig.colorbar(surf, shrink=0.5, aspect=10)
    #        
    #                    
    #        ax = fig.add_subplot(1, 2, 2, projection='3d')
    #        Z = self.ecm[0][0, :, :]
    #        print(Z)
    #        # ax.scatter(X, Y, Z, marker=m)
    #        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
    #                    linewidth=0, antialiased=False)
    #        ax.set_zlim(-1.01, 2.01)
    #        fig.colorbar(surf, shrink=0.5, aspect=10)
    #        
    #        plt.show()
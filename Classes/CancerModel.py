import mesa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import json
# from Classes import *
from Classes.CancerCell import CancerCell
from Classes.Vessel import Vessel
from Classes.utils import *
from Batch import maxSteps, dataCollectionPeriod, newSimulationFolder # Used to save mmp2 and ecm during runtime

from Classes.QuasiCircle import find_quasi_circle
from matplotlib import pyplot as plt
from matplotlib import cm

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
        raise Exception(f"Error, nothing returned for cluster survival probability, time {self.schedule.time}")
    

def count_total_cells(model):
    """Counts all the cells present in the model, in all sites.

    Input:
        model: CancerModel object.
    Returns:
        amount_of_cells (int): the total amount of cells in every site,
        NOT considering the vasculature    
    """
    amount_of_cells = len([1 for agent in model.schedule.agents if agent.agent_type == "cell"])
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

def get_cluster_radius_and_diameter(model,grid_id):
    """"
    Calculates the radius and diameter of the cancer cells in a given site.

    Input:
        model: CancerModel object.
        grid_id: the grid number for which the radius and diameter will be calculated.
    Returns:
        radius: the maximum of all the distances from each cell to the cell's centroid.
        diameter: the maximum of all the cell-cell distanes.
    """
    ccells_positions = np.array([[agent.pos[0],agent.pos.y[1]] for agent in model.schedule.agents if agent.agent_type == "cell" and agent.grid == grid_id])
    centroid = ccells_positions.mean(axis=0)
    #calculating radius
    radii  = np.linalg.norm(a - centroid, axis=1)
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
    return numpy.sqrt(x2-2*xx+x)

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
    def __init__(self, N, width, height, grids_number, seed=None):
        super().__init__()  
        self.vasculature = {}
        self.num_agents = N
        self.width = width
        self.height = height
        self.phenotypes = ["mesenchymal", "epithelial"]
        self.grid_vessels_positions = [[],[],[]]
        self.current_agent_id = 0

        self.mesenchymalCount = [np.zeros((width, height), dtype=np.float) for _ in range(grids_number)]
        self.epithelialCount = [np.zeros((width, height), dtype=np.float) for _ in range(grids_number)]

        self.grids_number = grids_number
        
        self.grids = [mesa.space.MultiGrid(width, height, False) for _ in range(self.grids_number)]
        self.grid_ids = [i+1 for i in range(self.grids_number)] # need a number to appear in the data analysis (.csv)
        
        self.schedule = mesa.time.RandomActivation(self)
        #list of numpy arrays, representing mmp2 and ecm concentration in each grid
        self.mmp2 = [np.zeros((2, width, height), dtype=float) for _ in range(grids_number)]
        self.ecm = [np.ones((2, width, height), dtype=float) for _ in range(grids_number)]

        self._initialize_grids()

        self.datacollector = mesa.DataCollector(
            model_reporters={"Total cells": count_total_cells}, agent_reporters={"Position": "pos", "Agent Type": "agent_type", "Phenotype": "phenotype", "Ruptured": "ruptured", "Grid": "grid_id"})

        #model_reporters={"Mmp2": "mmp2", "Grid": "grid"},

        self.mmp2Data = np.zeros((1, self.width, self.height), dtype=float)
        self.ecmData = np.ones((1, self.width, self.height), dtype=float)

    def step(self):
        """Advance the model by one step.
        
        The step function of the model will be called on each of the siulations steps.
        It will correctly allocate the cells coming from the vaculature, if any.
        Then, it will calculate the ECM and MMP2 concentration changes for this step,
        proliferate the cells if its due, and collect the data if it is the appropiate time.

        Input: none
        Returns: none
        """        
        if self.schedule.time in self.vasculature: # Add keys
            self.disaggregate_clusters(self.schedule.time)
            surviving_clusters = [cluster for cluster in self.vasculature[self.schedule.time] if self.random.random() < get_cluster_survival_probability(cluster)]
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
        if (self.schedule.time % doublingTimeM == 0 and self.schedule.time != 0):
            self.proliferate("mesenchymal")

        if (self.schedule.time % doublingTimeE == 0 and self.schedule.time != 0):
            self.proliferate("epithelial")

        # Saving of non agents data
        if isBatchRun and (self.schedule.time % dataCollectionPeriod == 0):
            # Saving Mmp2 and Ecm data
            for grid_id in self.grid_ids:
                new_mmp2_df = pd.DataFrame(self.mmp2[grid_id-1][0,:,:])
                mmp2CsvName = f"Mmp2-{grid_id}grid-{self.schedule.time}step.csv"
                pathToSave = os.path.join(parent_dir, newSimulationFolder, "Mmp2", mmp2CsvName)
                new_mmp2_df.to_csv(pathToSave)

                new_ecm_df = pd.DataFrame(self.ecm[grid_id-1][0,:,:])
                EcmCsvName = f"Ecm-{grid_id}grid-{self.schedule.time}step.csv"
                pathToSave = os.path.join(parent_dir, newSimulationFolder, "Ecm", EcmCsvName)
                new_ecm_df.to_csv(pathToSave)

            # Saves vasculature data
            if self.schedule.time == maxSteps:
                vasculature_json = json.dumps(self.vasculature)
                # {key: list of clusters} -> {timestep: [(number of Mcells, number of Ecells), ..., (..., ...)]}
                
                vasculatureJsonName = f"Vasculature-{self.schedule.time}step.json"
                pathToSave = os.path.join(parent_dir, newSimulationFolder, "Vasculature", vasculatureJsonName)
                
                with open(pathToSave, 'w') as f:
                    f.write(vasculature_json)

        print(f'step number: {self.schedule.time}')
        self.schedule.step()
        self.datacollector.collect(self)


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
                    new_cell = CancerCell(self.current_agent_id, self, agent.grid, agent.grid_id, agent.phenotype, agent.ecm, agent.mmp2)
                    self.current_agent_id += 1
                    self.schedule.add(new_cell)
                    agent.grid.place_agent(new_cell, (x,y))
                    total_amount_of_agents +=1
        


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

        # Testing code - creates a ruptured vessels at a given point
        a = Vessel(self.current_agent_id, self, True, self.grids[0], self.grid_ids[0])
        self.current_agent_id += 1
        self.schedule.add(a)
        self.grids[0].place_agent(a, (90,100))

        # Create vessels
        numNormalVessels = 0
        numRupturedVessels = 10
        numVesselsSecondary = 10
        numVesselsThird = 2 # just to test it, final code will not have 1 var to each grid

        # bad code, reduce number of for and make a counter to save the index to de put in each vessel
        # creates grid with 1 where vessels must not be placed
        # n_center_points_for_Vessels PDF = 200 
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
                temp = numRupturedVessels
                while temp > 0:
                    j = numRupturedVessels - temp
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

                temp = numNormalVessels
                while temp > 0:
                    j = numNormalVessels - temp
                    cell_to_place = [self.random.randrange(self.width), self.random.randrange(self.height)]
                    if cell_to_place in pos_coords:
                        a = Vessel(self.current_agent_id, self, False, self.grids[0], self.grid_ids[0])
                        self.current_agent_id += 1
                        self.schedule.add(a)
                        self.grids[0].place_agent(a, (int(cell_to_place[0]), int(cell_to_place[1])))

                        not_possible_array[cell_to_place[0], cell_to_place[1]] = 1
                        pos_coords.remove(cell_to_place)
                        temp -= 1


            if i > 0: # secondary and third grid
                if i == 1:
                    for m in range(numVesselsSecondary):
                        # make if to only create a vessel if given random value of x and y doesnt already has a vessel
                        a = Vessel(self.current_agent_id, self, False, self.grids[i], self.grid_ids[i])
                        self.current_agent_id += 1
                        self.schedule.add(a)
                        x = self.random.randrange(self.width)
                        y = self.random.randrange(self.height)
                        self.grids[i].place_agent(a, (x,y))
                        self.grid_vessels_positions[i] += [(x,y)]

                if i == 2:  
                    for m in range(numVesselsThird):
                        # make if to only create a vessel if given random value of x and y doesnt already has a vessel
                        a = Vessel(self.current_agent_id, self, False, self.grids[i], self.grid_ids[i])
                        self.current_agent_id += 1
                        self.schedule.add(a)
                        x = self.random.randrange(self.width)
                        y = self.random.randrange(self.height)
                        self.grids[i].place_agent(a, (x,y))
                        self.grid_vessels_positions[i] += [(x,y)]
                



    def calculateEnvironment(self, mmp2, ecm):
        for i in range(len(mmp2)):
            for cell in self.grids[i].coord_iter():
                cell_contents, x, y = cell
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



                        #ahora hay que mover la celula de acuerdo a las posibilidades

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



#step number: 488 and vasculature: {488: [(5, 0)]}
#step number: 488
#  0%|                                                                                                    | 0/1 [06:26<?, ?it/s] 
#Traceback (most recent call last):
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Batch.py", line 82, in <module>
#    main()
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Batch.py", line 32, in main
#    results = mesa.batch_run(
#  File "C:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\.venv\lib\site-packages\mesa\batchrunner.py", line 87, in batch_run
#    data = process_func(run)
#  File "C:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\.venv\lib\site-packages\mesa\batchrunner.py", line 157, in _model_run_func
#    model.step()
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 89, in step
#    possible_places = self.grid.get_neighborhood(arriving_point, moore=False, include_center=False)
#AttributeError: 'CancerModel' object has no attribute 'grid'. Did you mean: 'grids'?


#step number: 1998 and vasculature: {437: [(1, 0), (1, 0), (3, 0)], 461: [(1, 0), (0, 1)], 529: [(1, 0), (1, 0), (1, 0), (2, 0)], 539: [(1, 0), (1, 0), (1, 0), (2, 0)], 547: [(1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (0, 1), (1, 1)], 573: [(1, 0), (1, 0), (1, 0), (7, 0)], 575: [(1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (5, 0)],
import mesa
import matplotlib.pyplot as plt
import numpy as np
from Classes import *
from Classes.Vessel import Vessel
from Classes.utils import *
from QuasiCircle import find_quasi_circle

def count_total_cells(model):
    amount_of_cells = len([1 for agent in model.schedule.agents if agent.agent_type == "cell"])
    return amount_of_cells

def count_vasculature_cells(model):
    amount_of_cells = sum([len(value) for value in model.vasculature.values()])
    return amount_of_cells

class CancerModel(mesa.Model):

    def __init__(self, N, width, height, grids_number):
        super().__init__()  
        self.vasculature = {}
        self.num_agents = N
        self.width = width
        self.height = height
        self.phenotypes = ["mesenchymal", "epithelial"]
        self.mesenchymalCount = [np.zeros((totalTime, width, height), dtype=int) for _ in range(grids_number)]
        self.epithelialCount = [np.zeros((totalTime, width, height), dtype=int) for _ in range(grids_number)]
        self.grids_number = grids_number
        
        self.grids = [mesa.space.MultiGrid(width, height, False) for _ in range(self.grids_number)]
        
        self.schedule = mesa.time.RandomActivation(self)
        #list of numpy arrays, representing mmp2 and ecm concentration in each grid
        self.mmp2 = [np.zeros((totalTime, width, height), dtype=int) for _ in range(grids_number)]
        self.ecm = [np.zeros((totalTime, width, height), dtype=int) for _ in range(grids_number)]

        self._initialize_grids()

        self.datacollector = mesa.DataCollector(
            model_reporters={"Total cells": count_total_cells}#, agent_reporters={"Wealth": "wealth"}
            # model_reporters={"Cells in vasculature": count_vasculature_cells}#, agent_reporters={"Wealth": "wealth"}
        )

    def step(self):
        """Advance the model by one step."""
        self.datacollector.collect(self)
        if self.schedule.time in self.vasculature: # Add keys
            surviving_cells = [ccell for ccell in self.vasculature[self.schedule.time] if self.random.random() < single_cell_survival]
            n_cells_in_arriving_point = len([agent for agent in self.grids[1].get_cell_list_contents([(30,30)]) if agent.agent_type == "cell"])
            for ccell in surviving_cells:
                if carrying_capacity > n_cells_in_arriving_point:
                    self.grids[1].place_agent(ccell, (30,30))
                    self.schedule.add(ccell)
                elif carrying_capacity > len([agent for agent in self.grids[1].get_cell_list_contents([(30-1,30)]) if agent.agent_type == "cell"]):
                    self.grids[1].place_agent(ccell, (30-1,30))
                    self.schedule.add(ccell)
                elif carrying_capacity > len([agent for agent in self.grids[1].get_cell_list_contents([(30+1,30)]) if agent.agent_type == "cell"]):
                    self.grids[1].place_agent(ccell, (30+1,30))
                    self.schedule.add(ccell)
                elif carrying_capacity > len([agent for agent in self.grids[1].get_cell_list_contents([(30,30-1)]) if agent.agent_type == "cell"]):
                    self.grids[1].place_agent(ccell, (30,30-1))
                    self.schedule.add(ccell)
                elif carrying_capacity > len([agent for agent in self.grids[1].get_cell_list_contents([(30,30+1)]) if agent.agent_type == "cell"]):
                    self.grids[1].place_agent(ccell, (30,30+1))
                    self.schedule.add(ccell)
                    
        #Calculo do quimico que fomenta haptotaxis e da matriz extracelular
        self.calculateEnvironment(self.mmp2, self.ecm, self.schedule.time)
        self.schedule.step()

        # Reprodução
        if (self.schedule.time % doublingTimeM == 0 and self.schedule.time != 0):
            self.proliferate("mesenchymal")

        if (self.schedule.time % doublingTimeE == 0 and self.schedule.time != 0):
            self.proliferate("epithelial")

    def proliferate(self, cellType):
        all_agents = [agent for agent in self.schedule.agents]
        total_amount_of_agents = len(all_agents)
        for agent in all_agents:
            if agent.agent_type == "cell":
                x, y = agent.pos
                amount_of_cells = len([cell for cell in agent.grid.get_cell_list_contents([(x, y)]) if cell.agent_type == "cell"])
                if carrying_capacity > amount_of_cells and agent.phenotype == cellType:
                    new_cell = CancerCell(total_amount_of_agents + 1, self, agent.grid, agent.phenotype, agent.ecm, agent.mmp2)
                    self.schedule.add(new_cell)
                    agent.grid.place_agent(new_cell, (x,y))
                    total_amount_of_agents +=1
        


    def _initialize_grids(self):

        

        # I think we can put these constants in utils.py
        mesenchymal_proportion = 0.6
        epithelial_proportion = 0.4

        mesenchymal_number = round(self.num_agents * mesenchymal_proportion)
        n_center_points = 20
        possible_places = find_quasi_circle(n_center_points, self.width, self.height)[1]

        # Place all the agents in the quasi-circle area in the center of the grid
        for i in range(self.num_agents):

            if mesenchymal_number > 0:
                cell_type = "mesenchymal"
                mesenchymal_number -= 1
            elif mesenchymal_number == 0:
                cell_type = "epithelial"

            a = CancerCell(i, self, self.grids[0], cell_type, self.ecm[0], self.mmp2[0])
            
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
            a = CancerCell(i+self.num_agents+1, self, self.grids[1], "mesenchymal", self.ecm[1], self.mmp2[1])
            self.schedule.add(a)
        
            # Add the agent to a random grid cell
            x = self.random.randrange(3,7)
            y = self.random.randrange(3,7)
            self.grids[1].place_agent(a, (x, y))
        
        for i in range(1):
            a = Vessel(i+self.num_agents+amount_of_second_grid_CAcells+1, self, True, self.grids[0])
            self.schedule.add(a)
            self.grids[0].place_agent(a, (5, 15))

    def calculateEnvironment(self, mmp2, ecm, time):
        for i in range(len(mmp2)):
            for cell in self.grids[i].coord_iter():
                cell_contents, x, y = cell
                diff = 0
                for cancerCell in cell_contents:
                    if isinstance(cancerCell, CancerCell):
                        if cancerCell.phenotype == "mesenchymal":
                            self.mesenchymalCount[i][x][y] += 1
                            diff = dM
                        elif cancerCell.phenotype == "epithelial":
                            self.epithelialCount[i][x][y] += 1
                            diff = dE
                        else:
                            raise Exception("Unknown phenotype")
                        onLeftBorder = self.grids[i].out_of_bounds((x-1,y))
                        onRightBorder = self.grids[i].out_of_bounds((x+1,y))
                        onTopBorder = self.grids[i].out_of_bounds((x,y-1))
                        onBottomBorder = self.grids[i].out_of_bounds((x,y+1))
                        mmp2[i][time+1,x,y]=dmmp*tha/xha**2*(\
                                (mmp2[i][time,x+1,y] if not onRightBorder else mmp2[i][time,x-1,y])\
                                +(mmp2[i][time,x-1,y] if not onLeftBorder else mmp2[i][time,x+1,y])\
                                +(mmp2[i][time,x,y+1] if not onBottomBorder else mmp2[i][time,x,y-1])\
                                +(mmp2[i][time,x,y-1] if not onTopBorder else mmp2[i][time,x,y+1])\
                                )\
                                +mmp2[i][time,x,y]*(1-4*dmmp*tha/xha**2-th*Lambda)+tha*theta*self.mesenchymalCount[i][time,x,y]
                        ecm[i][time+1,x,y] = ecm[i][time,x,y]*(1-tha*(gamma1*self.mesenchymalCount[i][time,x,y]+gamma2*mmp2[i][time,x,y]))

                    #ahora hay que mover la celula de acuerdo a las posibilidades

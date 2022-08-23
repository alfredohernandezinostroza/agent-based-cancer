#%%
import mesa
import matplotlib.pyplot as plt
import numpy as np
from CustomCanvasGridVisualization import CanvasGridPrimary,CanvasGridSecondary

th    = 1e-3
tha   = th
xh    = 5e-3
xha   = xh
yh    = 1e-4
dM    = 1e-4
dE    = 5e-5
phiM  = 5e-4
phiE  = 5e-4
dmmp    = 1e-3
theta = 0.195
Lambda= 0.1
gamma1= 1
gamma2= 1
doublingTimeE=30
doublingTimeM=20
gridsize = 10
totalTime = 200
patchsize=3
middlePoint = round(gridsize/2)
E1 = 0.5461
E2 = 0.2553
E3 = 0.1986
Ps = 5e-4 #probability of survival for single cells in the vasculature
Pc = 2.5e-2 #probability of survival for clusters in the vasculature
Pd = 0.5 #probability of dissagreggation in vasculature
carrying_capacity = 4

class CancerCell(mesa.Agent):

    def __init__(self, unique_id, model, grid, phenotype, ecm, mmp2): #constructor
        super().__init__(unique_id, model)
        self.grid = grid
        self.phenotype = phenotype
        self.ecm = ecm
        self.mmp2 = mmp2
        self.agent_type = "cell"
    def step(self): #what will the agent do every time a step is made
        self.move()

    def move(self):
        time = self.model.schedule.time
        possible_steps = self.grid.get_neighborhood(
            self.pos,
            moore=False,
            include_center=True)
        x, y = self.pos
        onLeftBorder = self.grid.out_of_bounds((x-1,y))
        onRightBorder = self.grid.out_of_bounds((x+1,y))
        onTopBorder = self.grid.out_of_bounds((x,y-1))
        onBottomBorder = self.grid.out_of_bounds((x,y+1))
        Pleft = 0 if onLeftBorder else (th/xh**2*(dE-phiE/4*(0 if onRightBorder else self.ecm[time,x+1,y]-self.ecm[time,x-1,y])))
        Pright = 0 if onRightBorder else (th/xh**2*(dE+phiE/4*(0 if onLeftBorder else self.ecm[time,x+1,y]-self.ecm[time,x-1,y])))
        Ptop = 0 if onTopBorder else (th/xh**2*(dE+phiE/4*(0 if onBottomBorder else self.ecm[time,x,y+1]-self.ecm[time,x,y-1])))
        Pbottom = 0 if onBottomBorder else (th/xh**2*(dE-phiE/4*(0 if onTopBorder else self.ecm[time,x,y+1]-self.ecm[time,x,y-1])))
        Pstay = 1-(Pleft+Pright+Ptop+Pbottom)

        weights=[]
        for x2,y2 in possible_steps:
            if x2 < x:
                weights.append(Pleft)
            elif x2>x:
                weights.append(Pright)
            elif y2<y:
                weights.append(Pbottom)
            elif y2>y:
                weights.append(Ptop)
            else:
                weights.append(Pstay)


        new_position = self.random.choices(possible_steps,weights,k=1)[0]
        isTravelPoint = False
        for agent in self.grid.get_cell_list_contents([(new_position)]):
            if isinstance(agent, TravelPoint):
                isTravelPoint=True
        if isTravelPoint:
            print(self.grid.get_cell_list_contents([(new_position)]))
            print("Travelled!")
            # self.model.grid[1].place_agent(self,(0,5))
        self.grid.move_agent(self, new_position)

class TravelPoint(mesa.Agent):
    def __init__(self, unique_id, model, ruptured, grid):
        super().__init__(unique_id, model)
        self.model = model
        self.ruptured = ruptured
        self.grid = grid
        self.agent_type = "vessel"

    def step(self):
        pass


def count_total_cells(model):
    amount_of_cells = len([1 for agent in model.schedule.agents if agent.agent_type == "cell"])
    return amount_of_cells

class CancerModel(mesa.Model):

    def __init__(self, N, width, height):
        super().__init__()  
        self.num_agents = N
        self.width = width
        self.height = height
        self.phenotypes = ["mesenchymal", "epithelial"]
        self.mesenchymalCount = [np.zeros((totalTime, width, height), dtype=int), np.zeros((totalTime, width, height), dtype=int)]
        self.epithelialCount = [np.zeros((totalTime, width, height), dtype=int), np.zeros((totalTime, width, height), dtype=int)]
        
        #we have two grids, the primary site and the matastasis site
        self.grids = [mesa.space.MultiGrid(width, height, False), mesa.space.MultiGrid(width, height, False)]
        self.grid = self.grids[0]
        
        self.schedule = mesa.time.RandomActivation(self)
        #list of numpy arrays, representing mmp2 and ecm concentration in each grid
        self.mmp2 = [np.zeros((totalTime, width, height), dtype=float), np.zeros((totalTime, width, height), dtype=float)]
        self.ecm = [np.ones((totalTime, width, height), dtype=float), np.ones((totalTime, width, height), dtype=float)]
        # Create agents
        for i in range(self.num_agents):
            a = CancerCell(i, self, self.grids[0], "mesenchymal", self.ecm[0], self.mmp2[0])
            self.schedule.add(a)
        
            # Add the agent to a random grid cell
            x = self.random.randrange(3,7)
            y = self.random.randrange(3,7)
            self.grids[0].place_agent(a, (x, y))
        # Create agents at second grid
        amount=20
        for i in range(amount):
            a = CancerCell(i+N+1, self, self.grids[1], "mesenchymal", self.ecm[1], self.mmp2[1])
            self.schedule.add(a)
        
            # Add the agent to a random grid cell
            x = self.random.randrange(3,7)
            y = self.random.randrange(3,7)
            self.grids[1].place_agent(a, (x, y))
        
        for i in range(1):
            a = TravelPoint(i+N+amount+1, self, True, self.grids[0])
            self.schedule.add(a)

            x = self.random.randrange(self.width)
            y = self.random.randrange(self.width)
            self.grids[0].place_agent(a, (x, y))

        self.datacollector = mesa.DataCollector(
            model_reporters={"Total cells": count_total_cells}#, agent_reporters={"Wealth": "wealth"}
        )

    def step(self):
        """Advance the model by one step."""
        self.datacollector.collect(self)
        self.calculateEnvironment(self.mmp2, self.ecm, self.schedule.time)
        self.schedule.step()
        if (self.schedule.time % doublingTimeM == 0 and self.schedule.time != 0):
            all_agents = [agent for agent in self.schedule.agents]
            total_amount_of_agents = len(all_agents)
            for agent in all_agents:
                if  agent.agent_type == "cell" and agent.phenotype == "mesenchymal":
                    x, y = agent.pos
                    amount_of_cells = len([cell for cell in agent.grid.get_cell_list_contents([(x, y)]) if cell.agent_type == "cell"])
                    print(amount_of_cells)
                    if carrying_capacity > amount_of_cells:
                        new_cell = CancerCell(total_amount_of_agents + 1, self, agent.grid, "mesenchymal", agent.ecm, agent.mmp2)
                        self.schedule.add(new_cell)
                        agent.grid.place_agent(new_cell, (x, y))
                        total_amount_of_agents += 1
        # if (self.schedule.time % doublingTimeE == 0 and self.schedule.time != 0):
        #     all_agents = [agent for agent in self.schedule.agents]
        #     amount_of_agents = len(all_agents)
        #     for agent in all_agents:
        #         if agent.agent_type == "cell" and agent.phenotype == "mesenchymal":
        #             x, y = agent.pos
        #             new_cell = CancerCell(amount_of_agents, self, agent.grid, "mesenchymal", agent.ecm, agent.mmp2)
        #             self.schedule.add(new_cell)
        #             self.grids[1].place_agent(new_cell, (x, y))
        #             amount_of_agents += 1

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



# model = CancerModel(10,10,10)
# for i in range(10):
#     model.step() 

# plt.show()
# agent_wealth = [a.wealth for a in model.schedule.agents]
# plt.hist(agent_wealth)
#%%
# model = MoneyModel(50, 10, 10)
# for i in range(100):
#     model.step()

# agent_counts = np.zeros((model.grid.width, model.grid.height))
# for cell in model.grid.coord_iter():
#     cell_content, x, y = cell
#     agent_count = len(cell_content)
#     agent_counts[x][y] = agent_count
# plt.imshow(agent_counts, interpolation="nearest")
# plt.colorbar()
# plt.show()

# gini = model.datacollector.get_model_vars_dataframe()
# gini.plot()

# plt.show()
# %%
# all_wealth = []
# # This runs the model 100 times, each model executing 10 steps.
# for j in range(100):
#     # Run the model
#     model = MoneyModel(10)
#     for i in range(10):
#         model.step()

#     # Store the results
#     for agent in model.schedule.agents:
#         all_wealth.append(agent.wealth)

# plt.hist(all_wealth, bins=range(max(all_wealth) + 1))

# # %%
# agent_wealth = model.datacollector.get_agent_vars_dataframe()
# agent_wealth.head()
# # %%
# end_wealth = agent_wealth.xs(99, level="Step")["Wealth"]
# end_wealth.hist(bins=range(agent_wealth.Wealth.max() + 1))
# # %%
# one_agent_wealth = agent_wealth.xs(14, level="AgentID")
# one_agent_wealth.Wealth.plot()
# # %%

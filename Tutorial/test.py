#%%
import mesa
import matplotlib.pyplot as plt
import numpy as np

def compute_gini(model):
    agent_wealths = [agent.wealth for agent in model.schedule.agents]
    x = sorted(agent_wealths)
    N = model.num_agents
    B = sum(xi * (N - i) for i, xi in enumerate(x)) / (N * sum(x))
    return 1 + (1 / N) - 2 * B

class MoneyAgent(mesa.Agent):
    """An agent with fixed initial wealth."""

    def __init__(self, unique_id, model, grid):
        super().__init__(unique_id, model)
        self.wealth = 1
        self.grid = grid

    def step(self):
        self.move()
        if self.wealth > 0:
            self.give_money()

    
    def move(self):
        possible_steps = self.grid.get_neighborhood(
            self.pos,
            moore=False,
            include_center=True)
        print("grid2" if self.grid is self.model.grid2 else "grid1")
        print(possible_steps)

        new_position = self.random.choice(possible_steps)
        self.grid.move_agent(self, new_position)
    
    def give_money(self):
        cellmates = self.model.grid.get_cell_list_contents([self.pos])
        if len(cellmates) > 1:
            other = self.random.choice(cellmates)
            other.wealth += 1
            self.wealth -= 1

class MoneyModel(mesa.Model):
    """A model with some number of agents."""

    def __init__(self, N, width, height):
        super().__init__()   # this fixes the problem 
        self.num_agents = N
        self.grid = mesa.space.MultiGrid(width, height, False)
        self.grid2 = mesa.space.MultiGrid(width-5, height-5, False)
        self.schedule = mesa.time.RandomActivation(self)
        # Create agents
        for i in range(self.num_agents):
            a = MoneyAgent(i, self, self.grid)
            self.schedule.add(a)
        
            # Add the agent to a random grid cell
            x = self.random.randrange(self.grid.width)
            y = self.random.randrange(self.grid.height)
            self.grid.place_agent(a, (x, y))
        # Create agents at second grid
        for i in range(5):
            a = MoneyAgent(i+N+1, self, self.grid2)
            self.schedule.add(a)
        
            # Add the agent to a random grid cell
            x = self.random.randrange(self.grid2.width)
            y = self.random.randrange(self.grid2.height)
            self.grid2.place_agent(a, (x, y))
        
        self.datacollector = mesa.DataCollector(
            model_reporters={"Gini": compute_gini}, agent_reporters={"Wealth": "wealth"}
        )

    def step(self):
        """Advance the model by one step."""
        self.datacollector.collect(self)
        self.schedule.step()
        # print("Time: " + str(self.schedule.time))
        # print("Step: " + str(self.schedule.step))
        # print("thing")
        # print(self.grid.out_of_bounds(0,0))
        # print(self.grid.out_of_bounds(-1,0))
        # print(self.grid.out_of_bounds(self.grid.width,0))
        # print("thing")
        # for cell in self.grid.coord_iter():
        #     cell_content, x, y = cell
        #     for agent in cell_content:
        #         print(agent.wealth)


# model = MoneyModel(10)
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


def agent_portrayal(agent):
    portrayal = {"Shape": "circle",
                 "Filled": "true",
                 "r": 0.5}

    if agent.wealth > 0:
        portrayal["Color"] = "red"
        portrayal["Layer"] = 0
    else:
        portrayal["Color"] = "grey"
        portrayal["Layer"] = 1
        portrayal["r"] = 0.2
    return portrayal

def agent_portrayal2(agent):
    portrayal = {"Shape": "circle",
                 "Filled": "true",
                 "r": 0.5}

    if agent.wealth > 0:
        portrayal["Color"] = "black"
        portrayal["Layer"] = 0
    else:
        portrayal["Color"] = "grey"
        portrayal["Layer"] = 1
        portrayal["r"] = 0.2
    return portrayal

chart = mesa.visualization.ChartModule([{"Label": "Gini",
                      "Color": "Black"}],
                    data_collector_name='datacollector')

grid2 = mesa.visualization.CanvasGrid(agent_portrayal2, 10, 10, 500, 500,)
grid = mesa.visualization.CanvasGrid(agent_portrayal, 10, 10, 500, 500)
server = mesa.visualization.ModularServer(
    MoneyModel, [grid, grid2], "Money Modelssss", {"N": 100, "width": 10, "height": 10}
)
server.port = 8521 # The default
server.launch()

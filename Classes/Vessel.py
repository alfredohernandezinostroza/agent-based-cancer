import mesa
from Classes.utils import *

class Vessel(mesa.Agent):
    def __init__(self, unique_id, model, ruptured, grid, grid_id):
        super().__init__(unique_id, model)
        self.model = model
        self.ruptured = ruptured
        self.grid = grid
        self.grid_id = grid_id
        self.agent_type = "vessel"
        self.phenotype = False #need to be able do use data collector on agents

    def step(self):
        pass
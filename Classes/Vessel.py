import mesa
from Classes.utils import *

class Vessel(mesa.Agent):
    def __init__(self, unique_id, model, ruptured, grid):
        super().__init__(unique_id, model)
        self.model = model
        self.ruptured = ruptured
        self.grid = grid
        self.agent_type = "vessel"

    def step(self):
        pass
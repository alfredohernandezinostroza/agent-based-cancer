import mesa
from metaspread.vessel import Vessel
from metaspread.configs import *

class CancerCell(mesa.Agent):

    def __init__(self, unique_id, model, grid, grid_id, phenotype, ecm, mmp2):
        super().__init__(unique_id, model)
        self.grid = grid
        self.grid_id = grid_id
        self.phenotype = phenotype
        if self.phenotype == "mesenchymal":
            self.diff_coeff = dM
            self.phi = phiM
        else:
            self.diff_coeff = dE
            self.phi = phiE
        self.ecm = ecm
        self.mmp2 = mmp2
        self.agent_type = "cell"
        self.ruptured = False #need to be able do use data collector on agents
        
    def step(self): #what will the agent do every time a step is made
        self.move()

    def move(self):
        time = self.model.schedule.time
        possible_steps = self.grid.get_neighborhood(
            self.pos,
            moore=False,
            include_center=True)
        x, y = self.pos
        on_left_border    = self.grid.out_of_bounds((x-1,y))
        on_right_border   = self.grid.out_of_bounds((x+1,y))
        on_top_border     = self.grid.out_of_bounds((x,y+1))
        on_bottom_border  = self.grid.out_of_bounds((x,y-1))
        p_left = 0 if on_left_border else (th/xh**2*(self.diff_coeff-self.phi/4*(0 if on_right_border else self.ecm[0,x+1,y]-self.ecm[0,x-1,y])))
        p_right = 0 if on_right_border else (th/xh**2*(self.diff_coeff+self.phi/4*(0 if on_left_border else self.ecm[0,x+1,y]-self.ecm[0,x-1,y])))
        p_top = 0 if on_top_border else (th/xh**2*(self.diff_coeff+self.phi/4*(0 if on_bottom_border else self.ecm[0,x,y+1]-self.ecm[0,x,y-1])))
        p_bottom = 0 if on_bottom_border else (th/xh**2*(self.diff_coeff-self.phi/4*(0 if on_top_border else self.ecm[0,x,y+1]-self.ecm[0,x,y-1])))
        p_stay = 1-(p_left+p_right+p_top+p_bottom)

        weights=[]
        for x2,y2 in possible_steps:
            if x2 < x:
                weights.append(p_left)
            elif x2>x:
                weights.append(p_right)
            elif y2<y:
                weights.append(p_bottom)
            elif y2>y:
                weights.append(p_top)
            else:
                weights.append(p_stay)


        # new_position = (x,y+1)
        new_position = self.random.choices(possible_steps,weights,k=1)[0]
        is_vessel = False
        is_ruptured = False
        for agent in self.grid.get_cell_list_contents([(new_position)]):
            if isinstance(agent, Vessel):
                is_ruptured = agent.ruptured
                is_vessel=True
                break
        if is_vessel and self.grid_id == 1 and (is_ruptured or self.phenotype == "mesenchymal"): 
                x, y = new_position
                on_left_border    = self.grid.out_of_bounds((x-1,y))
                on_right_border   = self.grid.out_of_bounds((x+1,y))
                on_top_border     = self.grid.out_of_bounds((x,y+1))
                on_bottom_border  = self.grid.out_of_bounds((x,y-1))
                mesenchymal_ccells_to_travel  = [agent for agent in self.grid.get_cell_list_contents([(x,y)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                mesenchymal_ccells_to_travel += [] if on_left_border   else [agent for agent in self.grid.get_cell_list_contents([(x-1,y)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                mesenchymal_ccells_to_travel += [] if on_right_border  else [agent for agent in self.grid.get_cell_list_contents([(x+1,y)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                mesenchymal_ccells_to_travel += [] if on_top_border    else [agent for agent in self.grid.get_cell_list_contents([(x,y-1)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                mesenchymal_ccells_to_travel += [] if on_bottom_border else [agent for agent in self.grid.get_cell_list_contents([(x,y+1)]) if agent.agent_type == 'cell' and agent.phenotype == "mesenchymal"]
                epithelial_ccells_to_travel  = [agent for agent in self.grid.get_cell_list_contents([(x,y)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
                epithelial_ccells_to_travel += [] if on_left_border   else [agent for agent in self.grid.get_cell_list_contents([(x-1,y)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
                epithelial_ccells_to_travel += [] if on_right_border  else [agent for agent in self.grid.get_cell_list_contents([(x+1,y)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
                epithelial_ccells_to_travel += [] if on_top_border    else [agent for agent in self.grid.get_cell_list_contents([(x,y-1)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
                epithelial_ccells_to_travel += [] if on_bottom_border else [agent for agent in self.grid.get_cell_list_contents([(x,y+1)]) if agent.agent_type == 'cell' and agent.phenotype == "epithelial"]
        
                #if there are not clusters at that time in the vasculature dict, create a new key for that time
                #and add the tuple

                if self.model.vasculature.get(time + vasculature_time,False):
                    self.model.vasculature[time + vasculature_time] += [(len(mesenchymal_ccells_to_travel), len(epithelial_ccells_to_travel))]
                # if there are clusters, add the tuple to that key
                else:
                    self.model.vasculature[time + vasculature_time] = [(len(mesenchymal_ccells_to_travel), len(epithelial_ccells_to_travel))]
                for ccell in mesenchymal_ccells_to_travel + epithelial_ccells_to_travel:
                    ccell.grid.remove_agent(ccell)
                    ccell.model.schedule.remove(ccell)
        else:
            if carrying_capacity > len([cell for cell in self.grid.get_cell_list_contents([new_position]) if agent.agent_type == 'cell']):
                self.grid.move_agent(self, new_position)

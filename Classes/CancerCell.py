import mesa
from Classes.utils import *
from Classes.Vessel import Vessel

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


        # new_position = (x,y+1)
        new_position = self.random.choices(possible_steps,weights,k=1)[0]
        isVessel = False
        isRuptured = False
        for agent in self.grid.get_cell_list_contents([(new_position)]):
            if isinstance(agent, Vessel):
                isRuptured = agent.ruptured
                isVessel=True
        if isVessel:
            print("Begin travel!")
            x, y = new_position
            onLeftBorder = self.grid.out_of_bounds((x-1,y))
            onRightBorder = self.grid.out_of_bounds((x+1,y))
            onTopBorder = self.grid.out_of_bounds((x,y-1))
            onBottomBorder = self.grid.out_of_bounds((x,y+1))
            ccells_to_travel = [agent for agent in self.grid.get_cell_list_contents([(x,y)]) if agent.agent_type == 'cell']
            ccells_to_travel += [] if onLeftBorder else [agent for agent in self.grid.get_cell_list_contents([(x-1,y)]) if agent.agent_type == 'cell']
            ccells_to_travel += [] if onRightBorder else [agent for agent in self.grid.get_cell_list_contents([(x+1,y)]) if agent.agent_type == 'cell']
            ccells_to_travel += [] if onTopBorder else [agent for agent in self.grid.get_cell_list_contents([(x,y-1)]) if agent.agent_type == 'cell']
            ccells_to_travel += [] if onBottomBorder else [agent for agent in self.grid.get_cell_list_contents([(x,y+1)]) if agent.agent_type == 'cell']
            # print(ccells_to_travel)
            if self.model.vasculature.get(time + vasculature_time,False):
                self.model.vasculature[time + vasculature_time] += ccells_to_travel
            else:
                self.model.vasculature[time + vasculature_time] = ccells_to_travel
            for ccell in ccells_to_travel:
                ccell.grid.remove_agent(ccell)
                ccell.model.schedule.remove(ccell)
            print(self.model.vasculature)
            # print(self.grid.get_cell_list_contents([(new_position)]))
            # print("Travelled!")
            # self.model.grid[1].place_agent(self,(0,5))
        else:
            self.grid.move_agent(self, new_position)

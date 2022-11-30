import mesa
# import *
from Classes.CancerCell import CancerCell
from Classes.Vessel import Vessel
from Classes.utils import *
from Classes.QuasiCircle import find_quasi_circle
from Classes.MultipleCanvasGrid import MultipleCanvasGrid
from Classes.CancerModel import CancerModel

def agent_portrayal(agent):
    if False or isinstance(agent, Vessel):
        portrayal = {"Shape": "rect",
                 "Filled": "true",
                 "w": 0.8,
                 "h":0.8,
                 "r":0.6,

                 "Layer": 0}
        if agent.ruptured:
            portrayal["Color"] = "rgba(255, 0, 0, 1)"
        else:
            portrayal["Color"] = "rgba(255, 0, 0, 0.5)"

    else:
        portrayal = {"Shape": "circle",
                    "Filled": "true",
                    "r": 0.5}
        if agent.phenotype == "mesenchymal":
            portrayal["Color"] = "rgba(0, 0, 200, 0.25)"
            portrayal["Layer"] = 1
        else:
            portrayal["Color"] = "rgba(100, 100, 100, 0.25)"
            portrayal["Layer"] = 2
            portrayal["r"] = 0.2
    return portrayal


def main():
    gridsize     = 201 #PDF: 201
    width        = gridsize
    height       = gridsize
    grids_number = 2

    grids = [MultipleCanvasGrid(agent_portrayal, width, height, 402, 402, site=grid_number) for grid_number in range(grids_number)]

    chart = mesa.visualization.ChartModule([{"Label": "Total cells", "Color": "Black", 'w': 100}],
                        canvas_height=50, canvas_width=100,
                        data_collector_name='datacollector')

    visual_elements = grids + [chart]

    server = mesa.visualization.ModularServer(
        CancerModel, visual_elements, "Cancer model", {"N": 388, "width": width, "height": height, "grids_number": grids_number}
    )

    server.port = 8521
    server.launch()

if __name__ == "__main__":
    main()
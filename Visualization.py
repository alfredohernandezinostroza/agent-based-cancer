from ClassDefinitions import *

def agent_portrayal(agent):
    if False or isinstance(agent, TravelPoint):
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

def agent_portrayal2(agent):
    portrayal = {"Shape": "circle",
                 "Filled": "true",
                 "r": 0.5}

    if agent.phenotype == "mesenchymal":
        portrayal["Color"] = "red"
        portrayal["Layer"] = 0
    else:
        portrayal["Color"] = "grey"
        portrayal["Layer"] = 1
        portrayal["r"] = 0.2
    return portrayal

chart = mesa.visualization.ChartModule([{"Label": "Gini",
                      "Color": "Black"}],
                    data_collector_name='datacollector')
gridsize=50
width=gridsize
height=gridsize
grid = CanvasGridPrimary(agent_portrayal, width, height, 800, 800)
grid2 = CanvasGridSecondary(agent_portrayal, width, height, 800, 800)
server = mesa.visualization.ModularServer(
    CancerModel, [grid, grid2], "Cancer model", {"N": 50, "width": width, "height": height}
)
server.port = 8521 # The default
server.launch()

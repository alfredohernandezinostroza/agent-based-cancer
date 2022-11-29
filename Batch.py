import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Classes.CancerModel import *



gridsize     = 201 #PDF: 201
width        = gridsize
height       = gridsize
grids_number = 2

params = {"N": 388, "width": width, "height": height, "grids_number": grids_number}

results = mesa.batch_run(
    CancerModel,
    parameters=params,
    iterations=1,
    max_steps=995,
    number_processes=1,
    data_collection_period=1,
    display_progress=True
)

results_df = pd.DataFrame(results)
print(results_df.head())

# Error with the size of the array beeing too big -->
''' numpy.core._exceptions._ArrayMemoryError: Unable to allocate 3.61 GiB for an array
with shape (24000, 201, 201) and data type int32
'''

# Error in proliferation function --> problem with Id
''' Erro:   File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 75, in step
    self.proliferate("mesenchymal")
  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 89, in proliferate
    self.schedule.add(new_cell)
  File "C:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\.venv\lib\site-packages\mesa\time.py", line 68, in add
    raise Exception(
Exception: Agent with unique id 4809 already added to scheduler
'''
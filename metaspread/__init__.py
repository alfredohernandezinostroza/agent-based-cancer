from metaspread.configs import init_simulation_configs, check_if_configs_are_present
check_if_configs_are_present()
init_simulation_configs("simulations_configs.csv")
from metaspread.cancercell import CancerCell
from metaspread.cancermodel import CancerModel
from metaspread.vessel import Vessel
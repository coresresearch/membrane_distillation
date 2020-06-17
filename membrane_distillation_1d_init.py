""" 
Initialize variables, parameters, and Cantera objects for the simulation.

- Reads user inputs from membrane_distillation_inputs.yaml
- Creates necessary Cantera objects.
- Reads and stores necessary parameters.
- Initializes the solution vector


Returns: 
- Initial solution vector SV_0
- Dictionary of cantera objects obj
- Dictionary of paramters params
"""
import numpy as np
from ruamel.yaml import YAML
from pathlib import Path

path = Path('membrane_distillation_inputs.yaml')
yaml = YAML(typ='safe')
input = yaml.load(path)
SV_0 = np.zeros(1)

params = {}

obj = {}



# If tolerances are not specified, load in the defaults:
if 'rtol' not in params:
    params['rtol'] = 1e-3
if 'atol' not in params:
    params['atol'] = 1e-6

# If a final simulation time is not specified, load the default:
if 't final' not in params:
    params['t final'] = 1000 # seconds

# If an integrator method is not specified, load BDF:
if 'method' not in params:
    params['method'] = 'BDF'
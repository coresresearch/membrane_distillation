""" 
One-dimensional direct contact membrane distillation simulation.

This is the main file that runs the simulation code.  
- It initializes the solution vector using membrane_distillation_init.py
- The initialization routine includes reading in inputs from 
    membrane_distillation_input.yaml
- It runs the model for a range of conditions and possible microstructure 
    parameter sets
- Output is saved in the `output` folder, with a date-time string appended  
    to whatever other filename info it provided/specified.
- Currently, this function accepts the following inputs:
    * transport: 'Fick' or 'DGM'
"""

from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
from membrane_distillation_1d_output import *

def membrane_distillation_1d_model(save = None, membrane=None, tau_g=None, 
    transport=None, feed_temp=None, permeate_temp=None, SV_init=None):

    # Initialize the model, including reading of inputs, creating Cantera 
    #   objects, and reading and storing model parameters:
    from membrane_distillation_1d_init import initialize as init
    SV_0, obj, params = init(membrane, tau_g, transport, 
        feed_temp, permeate_temp)
    
    if SV_init is not None:
        SV_0 = SV_init

    # Read in the residual function, to use below:
    from membrane_distillation_functions import residual_1d

    # Run the model:
    solution = solve_ivp(lambda t, y: residual_1d(t, y, obj, params), 
        [0, params['t final']], SV_0, method = params['method'],
        rtol = params['rtol'], atol = params['atol'])

    if save:
        membrane_distillation_1d_output(solution, params, obj, save)
    else:
        data = membrane_distillation_1d_output(solution, params, obj, save)
        return data

# If this file is being executed as the "Main" python script, run the model, 
#   with any provided keyword arguments:
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--transport')
    parser.add_argument('--feed_temp')
    parser.add_argument('--permeate_temp')
    parser.add_argument('--membrane')
    parser.add_argument('--tau_g')
    parser.add_argument('--save')
    parser.add_argument('--init')
    args = parser.parse_args()
    
    membrane_distillation_1d_model(args.save, args.membrane, args.tau_g, 
        args.transport, args.feed_temp, args.permeate_temp, args.SV_init)

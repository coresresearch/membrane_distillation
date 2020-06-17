def initialize(membrane, transport, feed_temp, permeate_temp):
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
    import cantera as ct

    input_file = 'membrane_distillation_inputs.yaml'
    path = Path(input_file)
    yaml = YAML(typ='safe')
    inputs = yaml.load(path)
    SV_0 = np.zeros(1)

    # Create Cantera objects and store in 'obj' dict:
    gas = ct.Solution(input_file, inputs['phase-names']['gas'])
    liquid = ct.Solution(input_file, inputs['phase-names']['liquid'])
    interface = ct.Interface(input_file, inputs['phase-names']['interface'],\
        [gas, liquid])
    obj = {'gas': gas, 'liquid': liquid, 'interface': interface}

    # Read out the selected membrane.  Specification via command-line overrides #   all other input.  If no command-line input is provided, look in the 
    #   input file.  Otherwise, raise an exception:
    if membrane==None:
        if 'membrane' in inputs['simulation-params']:
            membrane = inputs['simulation-params']['membrane']
        else:
            raise ValueError('Please specify a membrane for your simulation'+\
                'either via command line (--membrane=[name string]) or as a'+\
                ' field in membrane_distillation_inputs.yaml:simulation-params')

    # This tracks whether or not the named membrane was found in the list of parameters in the input file.  Initialize to "not found":
    found = 0 
    for mem in inputs['membrane-data']:
        if mem['name']==membrane:
            if found: 
                # More than one set of membrane params with that name. Throw error.
                raise ValueError('There were two or more membranes with the name '\
                    +membrane+' in your input file.  Pleease edit.')

            membrane_params = mem
            found = 1



    # Create dictionary of parameters, 'params'
    params = {}
    params['n_y'] = inputs['simulation-params']['n_y']
    params['dyInv'] = params['n_y']/membrane_params['thickness']
    params['dy'] = 1./params['dyInv']

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

    return SV_0, obj, params
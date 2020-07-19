def initialize(membrane, tau_g, transport, feed_temp, permeate_temp):
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

    #===========================================================================
    #   READ IN INPUTS
    #===========================================================================
    input_file = 'membrane_distillation_inputs.yaml'
    path = Path(input_file)
    yaml = YAML(typ='safe')
    inputs = yaml.load(path)
    
    #===========================================================================
    #   READ OUT MEMBRANE PROPERTIES
    #===========================================================================
    # Read out the selected membrane.  Specification via command-line overrides #   all other input.  If no command-line input is provided, look in the 
    #   input file.  Otherwise, raise an exception:
    if membrane is None:
        if 'membrane' in inputs['simulation-params']:
            membrane = inputs['simulation-params']['membrane']
        else:
            raise ValueError('Please specify a membrane for your simulation'
                'either via command line (--membrane=[name string]) or as a'
                ' field in membrane_distillation_inputs.yaml:simulation-params')
    # print(membrane)
    # This tracks whether or not the named membrane was found in the list of 
    #   parameters in the input file.  Initialize to "not found":
    found = 0 
    for mem in inputs['membrane-data']:
        if mem['name']==str(membrane):
            if found: 
                # More than one set of membrane params has that name. 
                #   Throw an error.
                raise ValueError('There were two or more membranes with the'
                    ' name ' + membrane + ' in your input file. Pleease edit.')

            membrane_params = mem
            found = 1
    if not found:
        raise ValueError('Specified membrane ' + membrane + ' not found in'+
            ' membrane_distillation_inputs.yaml. Pleease edit.')
    #===========================================================================
    #   LOAD PARAMETERS INTO 'param' DICT
    #===========================================================================
    # Create dictionary of parameters, 'params'
    params = {}
    params['n_y'] = inputs['simulation-params']['n_y']
    params['dyInv'] = params['n_y']/membrane_params['thickness']
    params['dy'] = 1./params['dyInv']
    params['eps_g'] = membrane_params['porosity']
    params['r_pore'] = 0.5*membrane_params['pore-diameter']
    params['d_part'] = membrane_params['particle-diameter']
    params['kappa_membrane'] = membrane_params['thermal-conductivity']
    params['c_v_membrane'] = membrane_params['const-vol-specific-heat']
    params['rho_membrane'] = membrane_params['solid-phase-density']
    params['h_fg_feed'] = inputs['temp-data']['h_fg_feed']
    params['h_fg_permeate'] = inputs['temp-data']['h_fg_permeate']

    #===========================================================================
    #   LOAD PARAMETERS THAT CAN COME FROM EITHER INPUT FILE OR COMMAND LINE:
    #===========================================================================
    # Load the transport model:
    if transport:
        params['transport'] = transport
    elif 'transport-model' in inputs['simulation-params']:
        params['transport'] = inputs['simulation-params']['transport-model']
    else:
        raise ValueError('Please select a transport model: DGM or Fick.')

    # Load the tortuosity:
    if tau_g:
        params['tau_g'] = float(tau_g)
    elif 'tau_g' in membrane_params:
        params['tau_g'] = membrane_params['tau_g']
    else:
        raise ValueError('Please specify a tortuosity tau_g, either via the '
        'command line (--tau_g=[value]) or as a field in '
        'membrane_distillation_inputs.yaml:membrane-params')

    # For the Fick model, we manually store the membrane permeability:
    if params['transport']=='Fick':
        params['K_g'] = (4. * params['d_part']**2 * params['eps_g']**3
            / (72.*params['tau_g']**2*(1.-params['eps_g'])**2))
  
    Press = 101325.0     # Initial gas pressure [Pa]
    # Read out the feed flow temperature and convert to K:
    if feed_temp is not None:
        params['T_feed'] = 273.15 + feed_temp
    elif 'T_feed' in inputs['temp-data']:
        params['T_feed'] = 273.15 + inputs['temp-data']['T_feed']
    else:
        raise ValueError('Please specify a feed temperature, either via the'
            ' command line (--feed_temp=[temp]) or as a field in '
            'membrane_distillation_inputs.yaml:temp-data')      
    
    # Read out the permeate flow temperature and convert to K:
    if permeate_temp is not None:
        params['T_permeate'] = 273.15 + permeate_temp
    elif 'T_permeate' in inputs['temp-data']:
        params['T_permeate'] = 273.15 + inputs['temp-data']['T_permeate']
    else:
        raise ValueError('Please specify a permeate temperature, either via '
            'the command line (--permeate_temp=[temp]) or as a field in '
            'membrane_distillation_inputs.yaml:temp-data')
    
    #===========================================================================
    #   CHECK FOR AND LOAD INTEGRATOR PARAMETERS:
    #===========================================================================
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

    #===========================================================================
    #   Create Cantera objects and store in 'obj' dict:
    #===========================================================================
    if params['transport']=='Fick':
        gas = ct.Solution(input_file, inputs['phase-names']['gas'])
    elif params['transport']=='DGM':
        gas = ct.DustyGas(input_file, inputs['phase-names']['gas'])
        gas.tortuosity = params['tau_g']        
        gas.porosity = params['eps_g']
        gas.mean_pore_radius = params['r_pore']
        gas.mean_particle_diameter = params['d_part']

    liquid = ct.Solution(input_file, inputs['phase-names']['liquid'])
    interface = ct.Interface(input_file, inputs['phase-names']['interface'],
        [gas, liquid])
    obj = {'gas': gas, 'liquid': liquid, 'interface': interface}


    #===========================================================================
    #   POINTERS TO SOLUTION VECTOR LOCATIONS
    #===========================================================================
    # Each volume has n_species + 1 variables (the extra being temperature)
    params['n_vars'] = gas.n_species+1
    size_SV = params['n_y']*params['n_vars']
    SV_0 = np.zeros(size_SV)

    params['ptr_temp']=range(0, size_SV, params['n_vars'])
    params['ptr_rho_k'] = np.ndarray(shape=(params['n_y'], gas.n_species),
        dtype='int')
    for j in range(params['n_y']):
        params['ptr_rho_k'][j,:] = range(1+j*params['n_vars'],
                                         1 + j*params['n_vars'] + gas.n_species)

    #===========================================================================
    #   INITIALIZE THE SOLUTION VECTOR
    #===========================================================================
    # Calculate the mole fraction of water vapor at the feed side, assuming 
    #   saturation.
    # Generic string for gas-phase mole fractions:
    comp_string = 'N2:{},H2O:{}'
    X_h2o = inputs['temp-data']['P_H2O_feed']/Press
    X_k_feed = comp_string.format(1.-X_h2o,X_h2o)

    # Calculate the mole fraction of water vapor at the permeate side,
    #   assuming saturation:
    X_h2o = inputs['temp-data']['P_H2O_permeate']/Press
    X_k_permeate = comp_string.format(1.-X_h2o,X_h2o)

    # Use cantera gas object to calcualte species mass densities:
    # Feed side:
    gas.TPX = params['T_feed'], Press, X_k_feed
    rho_k_feed = gas.density*gas.Y

    # Permeate side:
    gas.TPX = params['T_permeate'], Press, X_k_permeate
    rho_k_permeate = gas.density*gas.Y
    
    # Starting guess assumes linear gradients:
    SV_0[params['ptr_temp']] = np.linspace(params['T_feed'],
        params['T_permeate'], params['n_y'])
    for k in np.arange(gas.n_species):
        SV_0[params['ptr_rho_k'][:,k]] = np.linspace(rho_k_feed[k],
            rho_k_permeate[k],params['n_y'])

    return SV_0, obj, params

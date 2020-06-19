"""
Defines functions for running membrane distillation simulation.

- Residual_1d: Calcualtes time-based derivative of solution vector, dSV/dt, for 
    a 1-D simulation.
"""
import numpy as np

def residual_1d(t, SV, obj, params):
    """ 
    This function selects the appropriate flux function (DGM or Fickean),
    Defines the time-derivative governing equations, and then calls the
    integrator function.  It returns the vector of times and an array
    of solution variables, one for each time step.
    """

    # Select the appropriate flux function:
    if params['transport']=='DGM':
        from membrane_distillation_1d_fluxes import DGM_fluxes as flux_calc
    elif params['transport']=='Fick':
        from membrane_distillation_1d_fluxes import Fick_fluxes as flux_calc
    else:
        raise ValueError('Please choose an available transport model. Options'\
            +' are DGM and Fick.')

    # Initialize the residual vector:
    dSV_dt = np.zeros_like(SV)

    # Save a local copy of the cantera objects, for readability:
    gas = obj['gas']

    n_y = params['n_y']
    dyInv = params['dyInv']
    eps_g_Inv = 1./params['eps_g']
    n_vars = params['n_vars']

    # Call the flux calculation, which returns the following arrays:
    #   - mass fluxes J_k,
    #   - gas-phase species enthalpies h_k,
    #   - Inverse of density x specific heat (weighted by phase vol. fractions)
    #   - Volume-weighted average themal conductivities 
    J_k, h_k, rhoCvInv, kappa = flux_calc(SV, obj, params)

    dSV_dt = (J_k[:n_vars*n_y] - J_k[n_vars:])*dyInv*eps_g_Inv

    q_conv = np.zeros(n_y+1)
    for k in np.arange(gas.n_species):
        q_conv = q_conv + J_k[k::n_vars]*h_k[k::n_vars]

    #T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = temp_data

    q_chem = np.zeros(n_y+1)
    q_chem[0] = -params['h_fg_feed']*J_k[n_vars+1]
    q_chem[-1] = -params['h_fg_permeate']*J_k[-2]

    Temps = np.append(np.append(params['T_feed'],SV[0::n_vars]), \
        params['T_permeate'])
    q_cond = kappa*(Temps[:-1]-Temps[1:])*dyInv
    q_cond[-1] = q_cond[-1]*2.

    q_tot =  q_cond + q_conv + q_chem
    dSV_dt[::n_vars] = (q_tot[:-1] - q_tot[1:])*dyInv*rhoCvInv

    return dSV_dt
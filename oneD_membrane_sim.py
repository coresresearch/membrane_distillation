"""
1D Direct Contact Membrane Distilation model: 1D_membrane_sim.py

Author: Spencer Gilleon and Steven C. DeCaluwe
Colorado School of Mines


"""
from string import Template
import numpy as np
from scipy.integrate import solve_ivp

def load_model(gas, membrane_data, params, temp_data, n_points):
    "This function sets up the initial guess for the solution vector"

    # Read in the experimental data and temperature data.
    #   Not all of this is used:
    T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = temp_data
    print(T_h,T_c)
    
    Press = 101325.0     # Initial gas pressure [Pa]
    T1 = 273.15+T_h;     # Hot feed flow temperature
    T2 = 273.15+T_c;     # Cold permeate flow

    # Calculate the mole fraction of water vapor at the feed side,
    #   assuming saturation:
    X_h2o = P_H2O_h/Press
    comp_string = 'N2:{},H2O:{}'
    X_k_h = comp_string.format(1.-X_h2o,X_h2o)

    # Calculate the mole fraction of water vapor at the permeate side,
    #   assuming saturation:
    X_h2o = P_H2O_c/Press
    X_k_c = comp_string.format(1.-X_h2o,X_h2o)
    dy = membrane_data['H']/n_points

    params['ptr_temp']=0
    params['ptr_rho_k']=1

    # Use cantera gas object to calcualte species mass densities:
    # Feed side:
    gas.TPX = T1, Press, X_k_h
    rho_k_h = gas.density*gas.Y

    # Permeate side:
    gas.TPX = T2, Press, X_k_c
    rho_k_c = gas.density*gas.Y

    # Initialize solution vector:
    # Each volume has n_species + 1 variables (the extra being temperature)
    n_vars = gas.n_species+1
    SV_0 = np.zeros(n_points*n_vars)

    # Starting guess assumes linear gradients:
    SV_0[0:n_vars*n_points:n_vars] = np.linspace(T1,T2,n_points)[None,:]
    for k in np.arange(gas.n_species):
        SV_0[k+1:n_vars*n_points+1:n_vars] = \
            np.linspace(rho_k_h[k],rho_k_c[k],n_points)[None,:]

    return SV_0, params


def run_model(t_sim, SV_0, obj, membrane_params, temp_data, params, trans_flag, method):
    """ This function selects the appropriate flux function (DGM or Fickean),
        Defines the time-derivative governing equations, and then calls the
        integrator function.  It returns the vector of times and an array
        of solution variables, one for each time step."""

    # Select the appropriate flux function:
    if trans_flag==0:
        from flux_functions import DGM_fluxes as flux_calc
    elif trans_flag==1:
        from flux_functions import Fick_fluxes as flux_calc

    # Define the system of ODEs to solve
    def dSVdt_func(t,SV,obj,membrane_params,temp_data,params):
        # Save a local copy of the cantera gas object, for readability:
        gas = obj['gas']

        # Initialize the vector of time-based derivatives:
        dSVdt = np.zeros(SV.shape)

        n_points = params['n_points']
        dyInv = params['dyInv']
        eps_g_Inv = 1./membrane_params['eps_g']
        n_vars = gas.n_species+1

        J_k, h_k, rhoCvInv, kappa = flux_calc(SV,obj,membrane_params,temp_data,params)

        dSVdt = (J_k[:n_vars*n_points] - J_k[n_vars:])*dyInv*eps_g_Inv

        q_conv = np.zeros(n_points+1)
        for k in np.arange(gas.n_species):
            q_conv = q_conv + J_k[k::n_vars]*h_k[k::n_vars]

        T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = temp_data

        q_chem = np.zeros(n_points+1)
        q_chem[0] = -h_fg_h*J_k[n_vars+1]
        q_chem[-1] = -h_fg_c*J_k[-2]

        Temps = np.append(np.append(T_h+273.15,SV[0::n_vars]), T_c+273.15)
        q_cond = kappa*(Temps[:-1]-Temps[1:])*dyInv
        q_cond[-1] = q_cond[-1]*2.

        q_tot =  q_cond + q_conv + q_chem
        dSVdt[::n_vars] = (q_tot[:-1] - q_tot[1:])*dyInv*rhoCvInv

        return dSVdt

    sol = solve_ivp(lambda t, y: dSVdt_func(t, y, obj, membrane_params, temp_data, params), \
        [0, t_sim], SV_0, method=method)

    return sol

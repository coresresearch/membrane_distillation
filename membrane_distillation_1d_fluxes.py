"""
Calculate terms related to mass and energy fluxes.

Returns:
- Mass fluxes 'J_k' (column vector of size (n_species+1)*(n_y+1)) [kg-k/m2/s]
- Species enthalpies 'h_k' per unit mass (column vector of size 
    (n_species+1)*(n_y+1)) [kJ/kg-k]
- Volume-averaged thermal conductivity 'kappa' for each interface between   
    volumes (column vector of size (n_y+1)) [W/m2/K]
- Inverse of density times constant volume specific heat 'rhoCvInv1 for each 
    volume (column vector of size n_y) [m3-K/kJ]
"""
import cantera as ct
import numpy as np
from string import Template
from copy import copy

def Fick_fluxes(SV, obj, params):
    # Read out local copies of Cantera objects:
    gas, liquid, interface = read_objs(obj)

    # Read out local copies of other parameters:
    dyInv, _, eps_g, D_eff_mult, T_feed, T_permeate = read_params(params)

    # Initialize return vectors:
    J_k = np.zeros((params['n_y']+1)*params['n_vars'])
    h_k = np.zeros((params['n_y']+1)*params['n_vars'])
    kappa = np.zeros(params['n_y']+1)
    rhoCvInv = np.zeros(params['n_y'])

    # Read out gas properties in first finite volume:
    j = 0
    T_2, rho_2, Y_k2 = read_properties(SV, params, j)

    # Set the cantera object states:
    gas.TDY = T_2, rho_2, Y_k2
    liquid.TP = T_feed, ct.one_atm
    interface.TP = T_feed, ct.one_atm
    
    # Read out gas-phase properties:
    D_k2, D_T2, P_2, mu_2, X_k2 = transport_properties(gas, D_eff_mult)

    # Flux into the first volume equals the reation rate at the gas-liquid 
    #   interface at the feed side.  We assume the interface temperature equals 
    #   the feed temperature
    sdot = interface.get_net_production_rates(gas)
    J_k[params['ptr_rho_k'][j]] = gas.molecular_weights*sdot
    # Species enthalpies:
    h_k[params['ptr_rho_k'][j]] = gas.partial_molar_enthalpies\
        /gas.molecular_weights
    # Constant-volume specific heats (per unit volume)
    rhoCvInv[j] = inv_spec_heat(gas, rho_2, params, eps_g)
    # Volume-averaged thermal conductivity:
    kappa[j] = eps_g*gas.thermal_conductivity + \
        (1.-eps_g)*params['kappa_membrane']

    # Loop through the other membrane volumes:
    for j in range(1,params['n_y']):
        # Re-assign previous "volume 2" variables to volume 1:
        T_1, Y_k1, rho_1, D_k1, D_T1, P_1, mu_1, X_k1 = \
            T_2, Y_k2, rho_2, D_k2, D_T2, P_2, mu_2, X_k2

        # Read out conditions from "next" volume:
        T_2, rho_2, Y_k2 = read_properties(SV, params, j)

        # Set gas object properties and read out other properties:
        gas.TDY = T_2, rho_2, Y_k2
        D_k2, D_T2, P_2, mu_2, X_k2 = transport_properties(gas, D_eff_mult)

        # Inverse specific heat at volume center:
        rhoCvInv[j] = inv_spec_heat(gas, rho_2, params, eps_g)

        # Properties at interface between volumes:
        mu_int, rho_int, D_k_int, Y_k_int, X_k_int, P_int, D_T_int = \
            0.5*(mu_1 + mu_2), 0.5*(rho_1+rho_2), 0.5*(D_k1+D_k2), \
            0.5*(Y_k1+Y_k2), 0.5*(X_k1+X_k2), 0.5*(P_1 + P_2), 0.5*(D_T1 + D_T2)

        #Diffusion operator:
        d_k_int = (X_k2 - X_k1 + (X_k_int - Y_k_int)*(P_2 - P_1)/P_int)*dyInv
        #Diffusive velocity, relative to mass-averaged velocity:
        V_k_int_o = -D_k_int*d_k_int/X_k_int \
            - D_T_int*(T_2 - T_1)*dyInv/rho_int/Y_k_int
        #Correction term, so that mass-weighted diffusive fluxes sum to zero:
        V_corr_int = -np.dot(Y_k_int, V_k_int_o)
        #Corrected diffusive velocities:
        V_diff_int = V_k_int_o + V_corr_int

        #Convective velocity, based on Darcy flow:
        V_conv_int = -params['K_g']*(P_2 - P_1)*dyInv/mu_int

        # Mass flux equals mass fraction times mass density times net velocity:
        J_k[params['ptr_rho_k'][j]] =  Y_k_int*rho_int*(V_diff_int + V_conv_int)

        # Set gas to interface properties to calculate thermal conductivity:
        gas.TDY = 0.5*(T_1+T_2), rho_int, Y_k_int
        h_k[params['ptr_rho_k'][j]] = gas.partial_molar_enthalpies / \
            gas.molecular_weights
        # Volume-weighted average 
        kappa[j] = eps_g*gas.thermal_conductivity + (1.-eps_g) \
            * params['kappa_membrane']

    # Arrays for last interface can't use SV pointers (size(array)>size(SV)).  #    We'll have to manually create an 'offset':
    offset = params['n_y']*params['n_vars']

    # Set Cantera object states:
    liquid.TP = T_permeate, ct.one_atm
    interface.TP = T_permeate, ct.one_atm
    # Read out chemical production rates at the permeate liquid interface:
    sdot = interface.get_net_production_rates(gas)

    # If gas phase species are created at this interface, this represents a 
    #   negative flux (flux in the negative direction):
    J_k[offset+1:offset+gas.n_species+1] = -gas.molecular_weights * sdot

    rhoCvInv[-1] = inv_spec_heat(gas, rho_2, params, eps_g)
    # Species enthalpies (per unit mass):
    h_k[offset+1:offset+gas.n_species+1] = gas.partial_molar_enthalpies / \
        gas.molecular_weights
    # Volume-weighted average thermal conductivity:
    kappa[-1] = eps_g*gas.thermal_conductivity + \
        (1.-eps_g)*params['kappa_membrane']

    return J_k, h_k, rhoCvInv, kappa

def DGM_fluxes(SV,obj,params):
    # Read out local copies of Cantera objects:
    gas, liquid, interface = read_objs(obj)

    # Read out local copies of other parameters:
    _, dY, eps_g, _, T_feed, T_permeate = read_params(params)

    # Initialize return vectors:
    J_k = np.zeros((params['n_y']+1)*params['n_vars'])
    h_k = np.zeros((params['n_y']+1)*params['n_vars'])
    kappa = np.zeros(params['n_y']+1)
    rhoCvInv = np.zeros(params['n_y'])  

    # Read out gas properties in first finite volume:
    j = 0
    T_2, rho_2, Y_k2 = read_properties(SV, params, j)
    
    # Set the cantera object states:
    gas.TDY = T_2, rho_2, Y_k2
    liquid.TP = T_feed, ct.one_atm
    interface.TP = T_feed, ct.one_atm

    # Flux into the first volume equals the reation rate at the gas-liquid 
    #   interface at the feed side.  We assume the interface temperature equals 
    #   the feed temperature
    sdot = interface.get_net_production_rates(gas)
    J_k[params['ptr_rho_k'][j]] = gas.molecular_weights*sdot
    
    # Species enthalpies:
    h_k[params['ptr_rho_k'][j]] = gas.partial_molar_enthalpies\
        /gas.molecular_weights
    # Inverse of volume-averaged specific heat:
    rhoCvInv[j] = inv_spec_heat(gas, rho_2, params, eps_g)

    # Volume-averaged thermal conductivity:
    # NOTE: at present the Cantera 'DustyGas' object has no method to calculate #   the gas-phase thermal conductivity. This requires creation of a 
    #   parallel Ideal Gas object:
    gas2 = ct.Solution('membrane_distillation_inputs.yaml','humid-air')  
    gas2.TDY = T_2, rho_2, Y_k2
    kappa[j] = eps_g*gas2.thermal_conductivity + \
        (1.-eps_g)*params['kappa_membrane']

    # Loop through the other membrane volumes:
    for j in range(1, params['n_y']):
        # Re-assign previous "volume 2" variables to volume 1:
        T_1, Y_k1, rho_1 = T_2, Y_k2, rho_2

        # Read out conditions from "next" volume:
        T_2, rho_2, Y_k2 = read_properties(SV, params, j)

        # Send properties to Cantera molar flux calculator and convert to 
        #   species mass fluxes:
        J_k[params['ptr_rho_k'][j]] = gas.molecular_weights * \
            gas.molar_fluxes(T_1, T_2, rho_1, rho_2, Y_k1, Y_k2, dY)

        # Properties at volume center:
        rhoCvInv[j] = inv_spec_heat(gas, rho_2, params, eps_g)

        # Set gas properties to interface conditions, and read out species 
        #   enthalpies, and volume-averaged thermal conductivity:
        gas.TDY = 0.5*(T_1+T_2),0.5*(rho_1+rho_2),0.5*(Y_k1+Y_k2)
        h_k[params['ptr_rho_k'][j]] = gas.partial_molar_enthalpies / \
            gas.molecular_weights
        # TODO #2
        gas2.TDY = 0.5*(T_1+T_2),0.5*(rho_1+rho_2),0.5*(Y_k1+Y_k2)
        kappa[j] = eps_g*gas2.thermal_conductivity + \
            (1.-eps_g)*params['kappa_membrane']

    
    # Arrays for last interface can't use SV pointers (size(array)>size(SV)).  #    We'll have to manually create an 'offset':
    offset = params['n_y']*params['n_vars']

    # Set Cantera object states:
    liquid.TP = T_permeate, ct.one_atm
    interface.TP = T_permeate, ct.one_atm
    # Read out chemical production rates at the permeate liquid interface:
    sdot = interface.get_net_production_rates(gas)

    # If gas phase species are created at this interface, this represents a 
    #   negative flux (flux in the negative direction):
    J_k[offset+1:offset+gas.n_species+1] = -gas.molecular_weights * sdot

    # Inverse of constant-volume specific heats (per unit volume)
    rhoCvInv[-1] = inv_spec_heat(gas, rho_2, params, eps_g)
    # Species enthalpies (per unit mass):
    h_k[offset+1:offset+gas.n_species+1] = gas.partial_molar_enthalpies / \
        gas.molecular_weights
    # Set 'gas2' properties back to those of final volume:
    gas2.TDY = T_2, rho_2, Y_k2
    # Volume-weighted average thermal conductivity:
    kappa[-1] = eps_g*gas2.thermal_conductivity + \
        (1.-eps_g)*params['kappa_membrane']

    return J_k, h_k, rhoCvInv, kappa

#===============================================================================#   HELPER FUNCTIONS
#===============================================================================

def read_objs(obj):
    gas = obj['gas']
    liquid = obj['liquid']
    interface = obj['interface']

    return gas, liquid, interface

def read_params(params):
    # Inverse volume thickness [1/m]
    dyInv = params['dyInv']
    # Volume thickness [m]
    dY = params['dy']
    # Porosity:
    eps_g = params['eps_g']
    # Scaling factor to obtain effective diffusion coefficients.  Note that 
    #   tortuosity is the tortuosity factor, not the geometric tortuosity:
    D_eff_mult = eps_g/params['tau_g']
    # Feed temperature [K]:
    T_feed = params['T_feed']
    # Permeate temperature [K]:
    T_permeate = params['T_permeate']

    return dyInv, dY, eps_g, D_eff_mult, T_feed, T_permeate

def read_properties(SV, params, j):
    T = SV[params['ptr_temp'][j]]
    rho_k = SV[params['ptr_rho_k'][j]]
    rho = sum(rho_k)
    Y_k = rho_k/rho

    return T, rho, Y_k

def transport_properties(gas, D_eff_mult):
    
    D_k = D_eff_mult*gas.mix_diff_coeffs
    #   Thermal diffusion coefficients:
    D_T = D_eff_mult*gas.thermal_diff_coeffs
    #   Pressure:
    P = gas.P
    #   Viscosity:
    mu = gas.viscosity
    #   Mole fractions:
    X_k = gas.X

    return D_k, D_T, P, mu, X_k

def inv_spec_heat(gas, rho_2, params, eps_g):
    
    cv_vol_gas = gas.cv_mass*rho_2
    cv_vol_membrane = params['c_v_membrane']*params['rho_membrane']
    # Inverse of volume-averaged specific heat:
    rhoCvInv = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)

    return rhoCvInv
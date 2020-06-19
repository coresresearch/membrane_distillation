"""
Calculate terms related to mass and energy fluxes.
"""
import cantera as ct
import numpy as np
from string import Template
from copy import copy

def Fick_fluxes(SV, obj, params):
    gas = obj['gas']
    liquid = obj['liquid']
    interface = obj['interface']

    dyInv = params['dyInv']

    eps_g = params['eps_g']
    D_eff_mult = eps_g/params['tau_g']

    J_k = np.zeros((params['n_y']+1)*params['n_vars'])
    h_k = np.zeros((params['n_y']+1)*params['n_vars'])
    kappa = np.zeros(params['n_y']+1)
    rhoCvInv = np.zeros(params['n_y'])

    T_feed = params['T_feed']
    T_permeate = params['T_permeate']

    # Read out gas properties in first finite volume:
    j = 0
    T_2 = SV[params['ptr_temp'][j]]
    rho_k2 = SV[params['ptr_rho_k'][j]]
    rho_2 = sum(rho_k2)
    Y_k2 = rho_k2/rho_2

    # Set the cantera object states:
    gas.TDY = T_2, rho_2, Y_k2
    liquid.TP = T_feed, ct.one_atm
    interface.TP = T_feed, ct.one_atm
    
    # Read out gas-phase properties:
    #   Mixture-averaged diffusin coefficients:
    D_k2 = D_eff_mult*gas.mix_diff_coeffs
    #   Thermal diffusion coefficients:
    D_T2 = D_eff_mult*gas.thermal_diff_coeffs
    #   Pressure:
    P_2 = gas.P
    #   Viscosity:
    mu_2 = gas.viscosity
    #   Mole fractions:
    X_k2 = gas.X

    # Flux into the first volume equals the reation rate at the gas-liquid 
    #   interface at the feed side.  We assume the interface temperature equals 
    #   the feed temperature
    sdot = interface.get_net_production_rates(gas)
    J_k[params['ptr_rho_k'][j]] = gas.molecular_weights*sdot
    # Species enthalpies:
    h_k[params['ptr_rho_k'][j]] = gas.partial_molar_enthalpies\
        /gas.molecular_weights
    # Constant-volume specific heats (per unit volume)
    cv_vol_gas = gas.cv_mass*rho_2
    cv_vol_membrane = params['c_v_membrane']*params['rho_membrane']
    # Volume-averaged specific heat:
    rhoCvInv[j] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)
    # Volume-averaged thermal conductivity:
    kappa[j] = eps_g*gas.thermal_conductivity + \
        (1.-eps_g)*params['kappa_membrane']

    # Loop through the other membrane volumes:
    for j in range(1,params['n_y']):
        # Re-assign previous "volume 2" variables to volume 1:
        T_1 = T_2
        Y_k1 = Y_k2
        rho_1 = rho_2
        D_k1 = D_k2
        D_T1 = D_T2
        P_1 = P_2
        mu1 = mu_2
        X_k1 = X_k2

        # Read out conditions from "next" volume:
        T_2 = SV[params['ptr_temp'][j]]
        rho_k2 = SV[params['ptr_rho_k'][j]]
        rho_2 = sum(rho_k2)
        Y_k2 = rho_k2/rho_2

        # Set gas object properties and read out other properties:
        gas.TDY = T_2, rho_2, Y_k2
        D_k2 = D_eff_mult*gas.mix_diff_coeffs
        D_T2 = D_eff_mult*gas.thermal_diff_coeffs
        P_2 = gas.P
        mu_2 = gas.viscosity
        X_k2 = gas.X

        # Properties at volume center:
        cv_vol_gas = gas.cv_mass*rho_2
        rhoCvInv[j] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)

        # Properties at interface between volumes:
        mu_int = 0.5*(mu1 + mu_2)
        rho_int = 0.5*(rho_1+rho_2)
        D_k_int = 0.5*(D_k1+D_k2)
        Y_k_int = 0.5*(Y_k1+Y_k2)
        X_k_int = 0.5*(X_k1+X_k2)
        P_int = 0.5*(P_1 + P_2)
        D_T_int = 0.5*(D_T1 + D_T2)

        #Diffusion operator:
        d_k_int = (X_k2 - X_k1 + (X_k_int - Y_k_int)*(P_2 - P_1)/P_int)*dyInv
        #Diffusive velocity, relative to mass-averaged velocity:
        V_k_int_o = -D_k_int*d_k_int/X_k_int \
            - D_T_int*(T_2 - T_1)*dyInv/rho_int/Y_k_int
        #Correction term, so that mass-weighted diffusive fluxes sum to zero:
        V_corr_int = -np.dot(Y_k_int, V_k_int_o)
        #Corrected diffusive velocities:
        V_diff_int = V_k_int_o + V_corr_int

        #Convective velocity:
        V_conv_int = -params['K_g']*(P_2 - P_1)*dyInv/mu_int

        # Mass flux equals mass fraction times mass density times net velocity:
        J_k[params['ptr_rho_k'][j]] =  Y_k_int*rho_int*(V_diff_int + V_conv_int)

        # Set gas to interface properties to calculate thermal conductivity:
        gas.TDY = 0.5*(T_1+T_2),rho_int,Y_k_int
        h_k[params['ptr_rho_k'][j]] = gas.partial_molar_enthalpies / \
            gas.molecular_weights

        kappa[j] = eps_g*gas.thermal_conductivity + (1.-eps_g) \
            * params['kappa_membrane']

    # Arrays for last interface can't use SV pointers (size(array)>size(SV)).  #    We'll have to manually create an 'offset':
    offset = (params['n_y'])*(gas.n_species+1)

    # Set Cantera object states:
    liquid.TP = T_permeate, ct.one_atm
    interface.TP = T_permeate, ct.one_atm
    # Read out chemical production rates at the permeate liquid interface:
    sdot = interface.get_net_production_rates(gas)

    # If gas phase species are created at this interface, this represents a 
    #   negative flux (flux in the negative direction):
    J_k[offset+1:offset+gas.n_species+1] = -gas.molecular_weights * sdot

    # Constant-volume specific heats (per unit volume)
    cv_vol_gas = gas.cv_mass*rho_2
    # Volume-weighted average specific heat:
    rhoCvInv[-1] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)
    # Species enthalpies (per unit mass):
    h_k[offset+1:offset+gas.n_species+1] = gas.partial_molar_enthalpies / \
        gas.molecular_weights
    # Volume-weighted average thermal conductivity:
    kappa[-1] = eps_g*gas.thermal_conductivity + \
        (1.-eps_g)*params['kappa_membrane']

    return J_k, h_k, rhoCvInv, kappa

def DGM_fluxes(SV,obj,params):

    gas = obj['gas']
    liquid = obj['liquid']
    interface = obj['interface']

    dY = params['dy']

    J_k = np.zeros((params['n_y']+1)*(gas.n_species+1))
    h_k = np.zeros((params['n_y']+1)*(gas.n_species+1))
    kappa = np.zeros(params['n_y']+1)
    rhoCvInv = np.zeros(params['n_y'])

    eps_g = params['eps_g']

    gas2 = ct.Solution('air_h2o.cti','air')    

    T_feed = params['T_feed']
    T_permeate = params['T_permeate']

    # Read out gas properties in first finite volume:
    j = 0
    T_2 = SV[params['ptr_temp'][j]]
    rho_k2 = SV[params['ptr_rho_k'][j]]
    rho_2 = sum(rho_k2)
    Y_k2 = rho_k2/rho_2
    gas.TDY = T_2, rho_2, Y_k2
    liquid.TP = T_feed, ct.one_atm
    interface.TP = T_feed, ct.one_atm

    sdot = interface.get_net_production_rates(gas)
    J_k[1:1+gas.n_species] = gas.molecular_weights*sdot

    h_k[1:gas.n_species+1] = gas.partial_molar_enthalpies/gas.molecular_weights
    cv_vol_gas = gas.cv_mass*rho_2
    cv_vol_membrane = params['c_v_membrane']*params['rho_membrane']
    rhoCvInv[0] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)

    gas2.TDY = T_2,rho_2,Y_k2
    kappa[j] = eps_g*gas2.thermal_conductivity + \
        (1.-eps_g)*params['kappa_membrane']

    for j in range(1, params['n_y']):
        T_1 = T_2
        Y_k1 = Y_k2
        rho_1 = rho_2
        offset = (j)*(gas.n_species+1)

        T_2 = SV[params['ptr_temp'][j]]
        rho_k2 = SV[params['ptr_rho_k'][j]]
        rho_2 = sum(rho_k2)
        Y_k2 = rho_k2/rho_2
        J_k[offset+1:offset+gas.n_species+1] = gas.molecular_weights * \
            gas.molar_fluxes(T_1,T_2,rho_1,rho_2,Y_k1,Y_k2,dY)

        gas2.TDY = 0.5*(T_1+T_2),0.5*(rho_1+rho_2),0.5*(Y_k1+Y_k2)
        h_k[offset+1:offset+gas.n_species+1] = gas2.partial_molar_enthalpies / \
            gas.molecular_weights

        kappa[j] = eps_g*gas2.thermal_conductivity + \
            (1.-eps_g)*params['kappa_membrane']

        gas.TDY = T_2,rho_2,Y_k2
        cv_vol_gas = gas.cv_mass*rho_2
        rhoCvInv[j] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)

    offset = (params['n_y'])*(gas.n_species+1)

    liquid.TP = T_permeate, ct.one_atm
    interface.TP = T_permeate, ct.one_atm
    sdot = interface.get_net_production_rates(gas)

    J_k[offset+1:offset+gas.n_species+1] = -gas.molecular_weights * sdot

    gas2.TDY = T_2,rho_2,Y_k2
    h_k[offset+1:offset+gas.n_species+1] = gas.partial_molar_enthalpies / \
        gas.molecular_weights
    kappa[-1] = eps_g*gas2.thermal_conductivity + (1.-eps_g)*params['kappa_membrane']
    cv_vol_gas = gas.cv_mass*rho_2
    rhoCvInv[-1] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)
    return J_k, h_k, rhoCvInv, kappa

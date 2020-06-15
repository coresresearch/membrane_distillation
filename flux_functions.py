"""
1D Direct Contact Membrane Distilation model: 1D_membrane_sim.py

Author: Spencer Gilleon and Steven C. DeCaluwe
Colorado School of Mines


"""
import cantera as ct
import numpy as np
from string import Template
from copy import copy

def Fick_fluxes(SV,obj,membrane_data,temp_data,params):
    gas = obj['gas']
    water = obj['liq']
    int = obj['int']

    dyInv = params['dyInv']

    eps_g = membrane_data['eps_g']
    tau_g = params['tau_g']

    J_k = np.zeros((params['n_points']+1)*(gas.n_species+1))
    h_k = np.zeros((params['n_points']+1)*(gas.n_species+1))
    kappa = np.zeros(params['n_points']+1)
    rhoCvInv = np.zeros(params['n_points'])

    T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = temp_data

    T2 = SV[params['ptr_temp']]
    rho_k2 = SV[params['ptr_rho_k']:params['ptr_rho_k']+gas.n_species]
    rho2 = sum(rho_k2)
    Y_k2 = rho_k2/rho2

    gas.TDY = T2, rho2, Y_k2
    water.TP = T_h+273.15, ct.one_atm
    int.TP = T_h+273.15, ct.one_atm

    # Read out gas properties in first finite volume:
    D_k2 = eps_g*gas.mix_diff_coeffs/tau_g
    D_T2 = eps_g*gas.thermal_diff_coeffs/tau_g
    C2 = gas.density_mole
    P2 = gas.P
    mu2 = gas.viscosity
    X_k2 = gas.X

    sdot = int.net_production_rates
    J_k[params['ptr_rho_k']:params['ptr_rho_k']+gas.n_species] = \
        gas.molecular_weights*sdot[params['sdot_gas_ptr']]

    h_k[params['ptr_rho_k']:params['ptr_rho_k']+gas.n_species] = \
        gas.partial_molar_enthalpies/gas.molecular_weights
    cv_vol_gas = gas.cv_mass*rho2
    cv_vol_membrane = membrane_data['c_v']*membrane_data['rho']
    rhoCvInv[0] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)
    kappa[0] = eps_g*gas.thermal_conductivity + \
        (1.-eps_g)*membrane_data['kappa']

    for j in range(params['n_points']-1):
        T1 = T2
        Y_k1 = Y_k2
        rho1 = rho2
        D_k1 = D_k2
        D_T1 = D_T2
        C1 = C2
        P1 = P2
        mu1 = mu2
        X_k1 = X_k2

        offset = (j+1)*(gas.n_species+1)

        T2 = SV[offset+params['ptr_temp']]
        rho_k2 = SV[offset+params['ptr_rho_k']:offset+params['ptr_rho_k']+gas.n_species]
        rho2 = sum(rho_k2)
        Y_k2 = rho_k2/rho2

        # Set gas propoerties and read out other transport-relevant properties:
        gas.TDY = T2, rho2, Y_k2
        D_k2 = eps_g*gas.mix_diff_coeffs/tau_g
        D_T2 = eps_g*gas.thermal_diff_coeffs/tau_g
        C2 = gas.density_mole
        P2 = gas.P
        mu2 = gas.viscosity
        X_k2 = gas.X

        # Properties at volume center:
        cv_vol_gas = gas.cv_mass*rho2
        rhoCvInv[j] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)

        # Properties at interface:
        mu_p = 0.5*(mu1 + mu2)
        rho_p = 0.5*(rho1+rho2)
        D_kp = 0.5*(D_k1+D_k2)
        Y_kp = 0.5*(Y_k1+Y_k2)
        X_kp = 0.5*(X_k1+X_k2)
        P_p = 0.5*(P1 + P2)
        D_Tp = 0.5*(D_T1 + D_T2)

        #Diffusion operator:
        #d_kp = (X_k2 - X_k1)*dyInv
        d_kp = (X_k2 - X_k1 + (X_kp - Y_kp)*(P2 - P1)/P_p)*dyInv
        #Diffusive velocity, relative to mass-averaged velocity:
        #V_kp_o = -D_kp*d_kp/X_kp
        V_kp_o  =-D_kp*d_kp/X_kp - D_Tp*(T2-T1)*dyInv/rho_p/Y_kp
        #Correction term, so that mass-weighted fluxes sum to zero:
        V_corr_p = -np.dot(Y_kp,V_kp_o)
        #Corrected diffusive velocities:
        V_diff_p = V_kp_o + V_corr_p

        #Convective velocity:
        V_conv_p = -params['K_g']*(P2 - P1)*dyInv/mu_p

        J_k[offset+1:offset+gas.n_species+1] =  Y_kp*rho_p*(V_diff_p + V_conv_p)

        # Calculate c_v:
        cv_vol_gas = gas.cv_mass*rho2
        rhoCvInv[j+1] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)

        # Set gas to interface properties to calculate thermal conductivity:
        gas.TDY = 0.5*(T1+T2),rho_p,Y_kp
        h_k[offset+1:offset+gas.n_species+1] = gas.partial_molar_enthalpies / \
            gas.molecular_weights

        kappa[j+1] = eps_g*gas.thermal_conductivity + (1.-eps_g)*membrane_data['kappa']


    offset = (params['n_points'])*(gas.n_species+1)

    water.TP = T_c+273.15, ct.one_atm
    int.TP = T_c+273.15, ct.one_atm
    sdot = int.net_production_rates

    J_k[offset+1:offset+gas.n_species+1] = -gas.molecular_weights * \
        sdot[params['sdot_gas_ptr']]


    cv_vol_gas = gas.cv_mass*rho2
    rhoCvInv[-1] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)
    h_k[offset+1:offset+gas.n_species+1] = gas.partial_molar_enthalpies / \
        gas.molecular_weights
    kappa[-1] = eps_g*gas.thermal_conductivity + (1.-eps_g)*membrane_data['kappa']

    return J_k, h_k, rhoCvInv, kappa

def DGM_fluxes(SV,obj,membrane_data,temp_data,params):

    gas = obj['gas']
    water = obj['liq']
    int = obj['int']

    dY = params['dY']

    J_k = np.zeros((params['n_points']+1)*(gas.n_species+1))
    h_k = np.zeros((params['n_points']+1)*(gas.n_species+1))
    kappa = np.zeros(params['n_points']+1)
    rhoCvInv = np.zeros(params['n_points'])

    eps_g = membrane_data['eps_g']

    gas2 = ct.Solution('air_h2o.cti','air')

    T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = temp_data

    T2 = SV[0]
    rho_k2 = SV[1:gas.n_species+1]
    rho2 = sum(rho_k2)
    Y_k2 = rho_k2/rho2
    gas.TDY = T2, rho2, Y_k2
    water.TP = T_h+273.15, ct.one_atm
    int.TP = T_h+273.15, ct.one_atm

    sdot = int.net_production_rates
    J_k[1:1+gas.n_species] = gas.molecular_weights*sdot[params['sdot_gas_ptr']]

    h_k[1:gas.n_species+1] = gas.partial_molar_enthalpies/gas.molecular_weights
    cv_vol_gas = gas.cv_mass*rho2
    cv_vol_membrane = membrane_data['c_v']*membrane_data['rho']
    rhoCvInv[0] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)


    gas2.TDY = T2,rho2,Y_k2
    kappa[0] = eps_g*gas2.thermal_conductivity + (1.-eps_g)*membrane_data['kappa']

    for j in range(params['n_points']-1):
        T1 = T2
        Y_k1 = Y_k2
        rho1 = rho2
        offset = (j+1)*(gas.n_species+1)

        T2 = SV[offset]
        rho_k2 = SV[offset+1:offset+gas.n_species+1]
        rho2 = sum(rho_k2)
        Y_k2 = rho_k2/rho2

        J_k[offset+1:offset+gas.n_species+1] = gas.molecular_weights * \
            gas.molar_fluxes(T1,T2,rho1,rho2,Y_k1,Y_k2,dY)

        gas2.TDY = 0.5*(T1+T2),0.5*(rho1+rho2),0.5*(Y_k1+Y_k2)
        h_k[offset+1:offset+gas.n_species+1] = gas2.partial_molar_enthalpies / \
            gas.molecular_weights

        kappa[j+1] = eps_g*gas2.thermal_conductivity + (1.-eps_g)*membrane_data['kappa']

        gas.TDY = T2,rho2,Y_k2
        cv_vol_gas = gas.cv_mass*rho2
        rhoCvInv[j+1] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)

    offset = (params['n_points'])*(gas.n_species+1)

    water.TP = T_c+273.15, ct.one_atm
    int.TP = T_c+273.15, ct.one_atm
    sdot = int.net_production_rates

    J_k[offset+1:offset+gas.n_species+1] = -gas.molecular_weights * \
        sdot[params['sdot_gas_ptr']]

    gas2.TDY = T2,rho2,Y_k2
    h_k[offset+1:offset+gas.n_species+1] = gas.partial_molar_enthalpies / \
        gas.molecular_weights
    kappa[-1] = eps_g*gas2.thermal_conductivity + (1.-eps_g)*membrane_data['kappa']
    cv_vol_gas = gas.cv_mass*rho2
    rhoCvInv[-1] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)

    return J_k, h_k, rhoCvInv, kappa

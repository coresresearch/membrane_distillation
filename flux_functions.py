"""
1D Direct Contact Membrane Distilation model: 1D_membrane_sim.py

Author: Spencer Gilleon and Steven C. DeCaluwe
Colorado School of Mines


"""
import cantera as ct
import numpy as np
from string import Template
from copy import copy

def DGM_fluxes(SV,obj,membrane_data,temp_data,params):

    gas = obj['gas']
    water = obj['liq']
    int = obj['int']

    dY = params['dY']

    J_k = np.zeros((params['n_points']+1)*(gas.n_species+1))
    wdot = np.zeros((params['n_points']+1)*(gas.n_species+1))
    h_k = np.zeros((params['n_points']+1)*(gas.n_species+1))
    kappa = np.zeros(params['n_points']+1)
    rhoCvInv = np.zeros(params['n_points'])

    eps_g = membrane_data['eps_g']

    gas2 = ct.Solution('air_h2o.cti','air')

    T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = temp_data

    T3 = SV[3]
    rhok3 = copy(SV[4:6])
    #rhok3[0] = rhok3[0]*
    rho3= sum(rhok3)
    Y3 = rhok3/rho3

    water.TP = T_h, ct.one_atm
    int.TP = T_h, ct.one_atm
    gas.TDY = T3, rho3, Y3

    sdot = int.net_production_rates
    print(sdot)

    T2 = SV[0]
    rho_k2 = SV[1:gas.n_species+1]
    rho2 = sum(rho_k2)
    Y2 = rho_k2/rho2
    gas.TDY = T2, rho2, Y2

    h_k[1:gas.n_species+1] = gas.partial_molar_enthalpies/gas.molecular_weights
    cv_vol_gas = gas.cv_mass*rho2
    #print(cv_vol_gas)
    cv_vol_membrane = membrane_data['c_v']*membrane_data['rho']
    rhoCvInv[0] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)


    gas2.TDY = T2,rho2,Y2
    kappa[0] = eps_g*gas2.thermal_conductivity + (1.-eps_g)*membrane_data['kappa']

    for j in range(params['n_points']-1):
        T1 = T2
        Y1 = Y2
        rho1 = rho2
        offset = (j+1)*(gas.n_species+1)

        T2 = SV[offset]
        rho_k2 = SV[offset+1:offset+gas.n_species+1]
        rho2 = sum(rho_k2)
        Y2 = rho_k2/rho2

        J_k[offset+1:offset+gas.n_species+1] = gas.molecular_weights * \
            gas.molar_fluxes(T1,T2,rho1,rho2,Y1,Y2,dY)

        gas2.TDY = 0.5*(T1+T2),0.5*(rho1+rho2),0.5*(Y1+Y2)
        h_k[offset+1:offset+gas.n_species+1] = gas2.partial_molar_enthalpies / \
            gas.molecular_weights

        kappa[j+1] = eps_g*gas2.thermal_conductivity + (1.-eps_g)*membrane_data['kappa']

        gas.TDY = T2,rho2,Y2
        cv_vol_gas = gas.cv_mass*rho2
        rhoCvInv[j+1] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)



    # Routine to extrapolate interface pressure:
    gas.TDY = T1,rho1,Y1
    P_1 = gas.P

    gas.TDY = T2,rho2,Y2
    P_2 = gas.P

    P2 = P_2 + 0.5*(P_2 - P_1)
    comp_string = Template('N2:$XN2,H2O:$XH2O')
    X2 = comp_string.substitute(XH2O = P_H2O_c/P2, XN2 =1-P_H2O_c/P2 )
    gas.TPX = T_c+273.15,P2,X2


    T1 = T2
    Y1 = Y2
    rho1 = rho2

    T2 = T_c+273.15
    Y2 = gas.Y
    rho2 = gas.density

    offset = (params['n_points'])*(gas.n_species+1)
    dY = 0.5*params['dY']
    J_k[offset+1:offset+gas.n_species+1] = gas.molecular_weights *\
        gas.molar_fluxes(T1,T2,rho1,rho2,Y1,Y2,dY)

    gas2.TDY = T2,rho2,Y2
    h_k[offset+1:offset+gas.n_species+1] = gas.partial_molar_enthalpies / \
        gas.molecular_weights
    kappa[-1] = eps_g*gas2.thermal_conductivity + (1.-eps_g)*membrane_data['kappa']
    cv_vol_gas = gas.cv_mass*rho2
    rhoCvInv[-1] = 1./(eps_g*cv_vol_gas + (1.-eps_g)*cv_vol_membrane)

    #print(J_k)
    return J_k, wdot, h_k, rhoCvInv, kappa

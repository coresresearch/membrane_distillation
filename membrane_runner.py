# -*- coding: utf-8 -*-
"""
1D Direct Contact Membrane Distilation model: membrane_runner.py

Author: Spencer Gilleon and Steven C. DeCaluwe
Colorado School of Mines


"""
print('hello1')

# Import necessary modules:
import numpy as np
import cantera as ct
import time
import pandas as pdpython
import timeit

from datetime import datetime

from flux_functions import DGM_fluxes, Fick_fluxes

print('hello2')
"""
This file defines a series of models and runs them, looping over some
combination of temperature boudnary conditions and/or tortuosity values, and
then enables plotting of results, to compare to experimental data.

The simulation computes and plots the predicted flux curve for a given membrane
and set of temperature boundary conditions. Temperature boundary conditions and
experimentally-measured fluxes are loaded from the 'Exp_data.csv' spreadsheet.
Individual membrane properties are input manually by the user, either in this
file or in the simulation file 'MD_sim.py'
"""


" --- BEGIN USER INPUTS ---"
# Read in the data:
temp_data = pd.read_csv('Exp_data.csv')

" Membrane Data: 3M 0.2um, 3M 0.45um "
# names: one string per membrane
# eps_g: Porosity [-]
# r_p: Pore radius [m]
# H: membrane thickness [m]
# d_p: average 'particle' diameter [m]
# kappa: Thermal conductivity of membrane solid [W/m-K]
# c_v: specific heat of membrane solid [J/kg-K]
# rho: mass density of membrane solid [kg/m3]
membrane_data = pd.DataFrame(data={'name':['200nm', '450nm'], \
    'eps_g':[0.85,0.85], 'r_p':[0.5*0.59e-6,0.5*0.79e-6], 'H':[110e-6, 110e-6],\
     'd_p':[0.19e-6, 0.19e-6], 'kappa':[0.16, 0.16], 'c_v':[1920,1920], \
     'rho':[946,946] })

# Vector of tortuosities to simulate:
tau_g_vec_0 = np.linspace(1.0,2.1,num=23)

# Specify which membranes to model.  This should be an numpy array,
#   corresponding to the proprerties in membrane_eps_g, membrane_r_p, and
#   membrane_H :
membranes = np.array([2])

# These flags specify whether or not to run Fickean and/or DGM transport models:
DGM = 1
Fick = 1

# Specify cti file and phase name(s):
ctifile = 'air_h2o.cti'
gasphase = 'air'
liqphase = 'water'
intphase = 'water_surf'

# In the cti file, specify the ordinal location of the gas phase in the list of
#   interface phases (0 = first)
int_ord_gas = 0

# Provide the number of volumes:
n_points = 20

# Simtulation time [s]
t_sim = 1000.

# Solver method.  Use RK45 for non-stiff problems, BDF for stiff problems
method = 'BDF'#'RK45'

" --- END USER INPUTS --- "

# Predicted tortuosity from correlation:
eps_g = membrane_data['eps_g']
membrane_data['tau_corr'] = (2. - eps_g)**2/eps_g

gas = ct.Solution(ctifile)
size = (gas.n_species+1)*n_points + 2

if DGM:
    Data_save_DGM = np.zeros((membranes.shape[0],temp_data.shape[0],tau_g_vec_0.shape[0],size))
if Fick:
    Data_save_Fick = np.zeros((membranes.shape[0],temp_data.shape[0],tau_g_vec_0.shape[0],size))

i_membrane = -1
for m in membranes:
    print(m)
    i_membrane += 1
    # Load/re-load initial vector of tortuosity factors:
    tau_g_vec = tau_g_vec_0
    membrane_params = membrane_data.iloc[m]
    print(membrane_params)

    params = {}
    params['dyInv'] = n_points/membrane_params['H']
    params['dY'] = membrane_params['H']/n_points
    params['n_points'] = n_points

    # Add correlation value:
    #tau_g_vec = np.append(tau_g_vec,membrane_params['tau_corr'])
    #tau_g_vec.sort()

    for i_temp, row_temp in temp_data.iterrows():
        print('i_temp = ',i_temp)

        T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = temp_data.loc[i_temp,:]

        int_temps_string = (str(int(row_temp['Feed Temp [C]']))+'_'+str(int(row_temp['Perm Temp [C]']))+'_profiles.csv')
        profiles = np.loadtxt(int_temps_string,delimiter=',')
        T_f = profiles[:,0]
        T_d = profiles[:,1]
        x = profiles[:,2]

        T_f_avg = np.sum(T_f[:-2]*(x[1:-1]-x[:-2]))/x[-1]
        T_d_avg = np.sum(T_d[:-2]*(x[1:-1]-x[:-2]))/x[-1]

        print(T_f_avg,T_d_avg)
        temp_data.at[i_temp,'Feed Temp [C]']=T_f_avg
        temp_data.at[i_temp,'Perm Temp [C]']=T_d_avg


        t = datetime.now()
        dtStr=t.strftime('%m/%d/%Y')

        save_string = (membrane_params['name']+'_'+str(int(row_temp['Feed Temp [C]']))+'_'+str(int(row_temp['Perm Temp [C]']))+'_'+dtStr)
        print(save_string)
        i_tau=-1
        for tau_g in tau_g_vec:
            i_tau += 1
            print('Tortusity factor = ',tau_g)
            tic = timeit.default_timer()
            if DGM:

                # Transport flag: 0 for DGM, 1 for Fickean
                trans_flag = 0

                gas = ct.DustyGas(ctifile)
                gas.porosity = membrane_params['eps_g']
                gas.tortuosity = tau_g
                gas.mean_pore_radius = 2.*membrane_params['r_p']
                gas.mean_particle_diameter = membrane_params['d_p']


                liq = ct.Solution(ctifile,liqphase)
                liq_int = ct.Interface(ctifile,intphase,[gas,liq])

                params['sdot_gas_ptr'] = np.arange(int_ord_gas,int_ord_gas+gas.n_species)

                from oneD_membrane_sim import load_model, run_model

                SV_0, params = load_model(gas, membrane_params, params, temp_data.loc[i_temp,:], n_points)

                obj = {'gas':gas, 'liq':liq, 'int':liq_int}

                solution = run_model(t_sim,SV_0,obj,membrane_params, temp_data.loc[i_temp,:],params,trans_flag,method)
                toc = timeit.default_timer()
                print(toc-tic,' seconds elapsed for DGM.')
                tic = toc

                n_vars = 1+gas.n_species
                SV = solution.y[:,-1].T


                J_k, _ , _ , _ = DGM_fluxes(SV,obj,membrane_params,temp_data.loc[i_temp,:],params)
                dataP=np.array([[tau_g], [J_k[4]]])
                data = np.vstack((dataP, SV[:,np.newaxis]))
                Data_save_DGM[i_membrane,i_temp,i_tau,:] = data[:,0]

                T_DGM = SV[0::n_vars]
                rho_k_h2o_DGM = SV[1::n_vars]
                rho_k_n2_DGM = SV[2::n_vars]
                rho = rho_k_h2o_DGM+rho_k_n2_DGM
                Y_k_h2o_DGM = rho_k_h2o_DGM/rho

                P_DGM = np.zeros(n_points)
                for j in np.arange(n_points):
                    temp = T_DGM[j]
                    Yk = [rho_k_h2o_DGM[j],rho_k_n2_DGM[j]]
                    rho = np.sum(Yk)
                    gas.TDY = temp,rho,Yk
                    P_DGM[j] = gas.P

            if Fick:
                params['tau_g'] = tau_g
                params['K_g'] = 4*membrane_params['d_p']**2 \
                    *membrane_params['eps_g']**3 \
                    /(72*tau_g**2*(1-membrane_params['eps_g'])**2)
                # Transport flag: 0 for DGM, 1 for Fickean
                trans_flag = 1

                gas = ct.Solution(ctifile)

                liq = ct.Solution(ctifile,liqphase)
                liq_int = ct.Interface(ctifile,intphase,[gas,liq])

                params['sdot_gas_ptr'] = np.arange(int_ord_gas,int_ord_gas+gas.n_species)

                from oneD_membrane_sim import load_model, run_model

                SV_0, params = load_model(gas, membrane_params, params, temp_data.loc[i_temp,:], n_points)

                obj = {'gas':gas, 'liq':liq, 'int':liq_int}

                solution = run_model(t_sim,SV_0,obj,membrane_params, temp_data.loc[i_temp,:],params,trans_flag,method)
                toc = timeit.default_timer()
                print(toc-tic,' seconds elapsed for Fick.')
                tic = toc

                n_vars = 1+gas.n_species
                SV = solution.y[:,-1].T

                J_k, _ , _ , _ = Fick_fluxes(SV,obj,membrane_params,temp_data.loc[i_temp,:],params)
                dataP=np.array([[tau_g], [J_k[4]]])
                data = np.vstack((dataP, SV[:,np.newaxis]))
                Data_save_Fick[i_membrane,i_temp,i_tau,:] = data[:,0]

        if DGM:
            np.savetxt(save_string+'_DGM.csv', Data_save_DGM[i_membrane,i_temp,:,:], delimiter=',')
        if Fick:
            np.savetxt(save_string+'_Fick.csv', Data_save_Fick[i_membrane,i_temp,:,:], delimiter=',')



        """T = SV[0::n_vars]
        rho_k_h2o = SV[1::n_vars]
        rho_k_n2 = SV[2::n_vars]
        rho = rho_k_h2o+rho_k_n2
        Y_k_h2o = rho_k_h2o/rho

        P = np.zeros(n_points)
        for j in np.arange(n_points):
            temp = T[j]
            Yk = [rho_k_h2o[j],rho_k_n2[j]]
            rho = np.sum(Yk)
            gas.TDY = temp,rho,Yk
            P[j] = gas.P"""

        """plt.figure
        plt.plot(np.arange(1,n_points+1),Y_k_h2o,np.arange(1,n_points+1),rho_k_h2o)
        plt.plot(np.arange(1,n_points+1),Y_k_h2o_DGM,np.arange(1,n_points+1),rho_k_h2o_DGM)

        plt.figure
        plt.plot(np.arange(1,n_points+1),P)
        plt.plot(np.arange(1,n_points+1),P_DGM)

        plt.figure
        plt.plot(np.arange(1,n_points+1),T)
        plt.plot(np.arange(1,n_points+1),T_DGM)
        plt.show()"""

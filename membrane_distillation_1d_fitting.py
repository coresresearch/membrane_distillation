"""
This script runs over a user-specified set of inputs to fit 1d membrane distillation model results against corresponding experimental data.
"""
# Import necessary modules:
import numpy as np
import timeit
from membrane_distillation_1d_model import membrane_distillation_1d_model

from datetime import datetime

import cantera as ct
import time
import pandas as pdpython

import pandas as pd

#===============================================================================
#   USER INPUTS
#===============================================================================

# Vector of tortuosities to simulate:
tau_g_vec_0 = np.linspace(0.5,2.1,num=3)#23)

# Specify which membranes to model.  The named membrane should exist in the 
#   file 'membrane_distillation_inputs.yaml'
membrane = '200nm'

# These flags specify whether or not to run Fickean and/or DGM transport models:
DGM = 1
Fick = 1

# Folder in which to save results. Date-time will be appended to this
savename = 'fitting_demo'

#===============================================================================
#   SET UP OUTPUTS
#===============================================================================

outputs = list()
t = datetime.now()
date_string = t.strftime('%m_%d_%Y_%H%M')
save_folder = savename + '_' + date_string

if DGM:
    DGM_fluxes = np.zeros((tau_g_vec_0.shape[0],2))
if Fick:
    Fick_fluxes = np.zeros((tau_g_vec_0.shape[0],2))

for i_tau, tau_g in enumerate(tau_g_vec_0):
    print('Tortusity factor = ',tau_g)
    tic = timeit.default_timer()

    if DGM:
        output = 'outputs/' + save_folder + '/DGM_' + str(tau_g) +'.csv'
        outputs.append(output)
        membrane_distillation_1d_model(membrane=membrane, tau_g=tau_g, 
            transport='DGM', save = output)
        toc = timeit.default_timer()
        print(toc-tic,' seconds elapsed for DGM.')
        tic = toc

    if Fick:
        output = 'outputs/' + save_folder + '/Fick_' + str(tau_g) +'.csv'
        outputs.append(output)
        membrane_distillation_1d_model(membrane=membrane, tau_g=tau_g, 
            transport='DGM', save = output)
        toc = timeit.default_timer()
        print(toc-tic,' seconds elapsed for Fick.')
        tic = toc


"""  J_k, _ , _ , _ = DGM_fluxes(SV,obj,membrane_params,temp_data.loc[i_temp,:],params)
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

            """
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

# -*- coding: utf-8 -*-
"""
1D Direct Contact Membrane Distilation model: membrane_runner.py

Author: Spencer Gilleon and Steven C. DeCaluwe
Colorado School of Mines


"""

# Import necessary modules:
import pylab as plt
import numpy as np
import cantera as ct
import time
import pandas as pd
import timeit




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
membrane_data = pd.DataFrame(data={'name':['0.2um', '0.45um'], \
    'eps_g':[0.85,0.85], 'r_p':[0.5*0.59e-6,0.5*0.79e-6], 'H':[110e-6, 110e-6],\
     'd_p':[0.19e-6, 0.19e-6], 'kappa':[0.16, 0.16], 'c_v':[1920,1920], \
     'rho':[946,946] })

# Vector of tortuosities to simulate:
tau_g_vec_0 = np.linspace(1.0,3.0,num=41)

# Specify which membranes to model.  This should be an numpy array,
#   corresponding to the proprerties in membrane_eps_g, membrane_r_p, and
#   membrane_H :
membranes = np.array([0]);

# These flags specify whether or not to run Fickean and/or DGM transport models:
DGM = 1
Fick = 1

# Specify cti file and phase name(s):
ctifile = 'air_h2o.cti'
gasphase = 'air'
liqphase = 'water'
intphase = 'water_surf'

# Provide the number of volumes:
n_points = 20

# Simtulation time [s]
t_sim = 1000.

# Solver method.  Use RK45 for non-stiff problems, LSODA for stiff problems
method = 'BDF'#'RK45'

" --- END USER INPUTS --- "

# Predicted tortuosity from correlation:
eps_g = membrane_data['eps_g']
membrane_data['tau_corr'] = (2. - eps_g)**2/eps_g

if DGM:
    model_fluxes_DGM = np.zeros((membranes.shape[0],tau_g_vec_0.shape[0]+1))

if Fick:
    model_fluxes_Fick = np.zeros((membranes.shape[0],tau_g_vec_0.shape[0]+1))


for m in membranes:

    # Load/re-load initial vector of tortuosity factors:
    tau_g_vec = tau_g_vec_0
    membrane_params = membrane_data.iloc[m]
    print(membrane_params)

    params = {}
    params['dyInv'] = n_points/membrane_params['H']
    params['dY'] = membrane_params['H']/n_points
    params['n_points'] = n_points

    # Add correlation value:
    tau_g_vec = np.append(tau_g_vec,membrane_params['tau_corr'])
    tau_g_vec.sort()

    for i_temp, row_temp in temp_data.iterrows():
        #T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = \
        #    temp_data.loc[index,:]


        for tau_g in tau_g_vec:
            print('Tortusity factor = ',tau_g)
            tic = timeit.default_timer()
            if DGM:

                # Transport flag: 0 for DGM, 1 for Fickean
                trans_flag = 0

                gas = ct.DustyGas(ctifile)
                gas.porosity = membrane_params['eps_g']
                gas.tortuosity = tau_g
                gas.mean_pore_radius = membrane_params['r_p']
                gas.mean_particle_diameter = membrane_params['d_p']

                liq = ct.Solution(ctifile,liqphase)
                int = ct.Interface(ctifile,intphase,[gas,liq])

                from oneD_membrane_sim import load_model, run_model

                SV_0 = load_model(gas, membrane_params, temp_data.loc[i_temp,:], n_points)

                obj = {'gas':gas, 'liq':liq, 'int':int}

                solution = run_model(t_sim,SV_0,obj,membrane_params, temp_data.loc[i_temp,:],params,trans_flag,method)
                toc = timeit.default_timer()
                print(toc-tic,' seconds elapsed.')
                print(solution)
                #plt.plot(solution.t,solution.y[4:-2:3,:].T)
                n_vars = 1+gas.n_species
                SV = solution.y[:,-1].T
                T = SV[0::n_vars]
                rho_k_h2o = SV[1::n_vars]
                rho_k_n2 = SV[2::n_vars]
                rho = rho_k_h2o+rho_k_n2
                Y_k_h2o = rho_k_h2o/rho
                plt.figure
                plt.plot(np.arange(1,n_points+1),Y_k_h2o,np.arange(1,n_points+1),rho_k_h2o)
                plt.show()

                P = np.zeros(n_points)
                for j in np.arange(n_points):
                    temp = T[j]
                    Yk = [rho_k_h2o[j],rho_k_n2[j]]
                    rho = np.sum(Yk)
                    gas.TDY = temp,rho,Yk
                    P[j] = gas.P


                plt.figure
                plt.plot(np.arange(1,n_points+1),P)
                plt.show()

                plt.figure
                plt.plot(np.arange(1,n_points+1),T)
                plt.show()
                print(Xk)



"""
class ep():

    # The Crate is the rate of charge/discharge - how many charges/discharges can
    #   be carried out in 1 hour? This sets the current density:
    C_rate = 0.05

    # For these half-cell simulations, let us assume a constant cathode voltage for
    #   our energy storage and efficiency calcs:
    V_cathode = 3.0  # [V]

    # Simulation temperature
    T = 300  # [K]

    # Number of nodes in the y-direction:
    npoints = 1

    # Number of "shells" in an anode particle for the Li intercalation diffusion:
    nshells = 5

    # Initial conditions:
    X_an_init = 'LiC6:0.0001, C6:0.9999'
    X_elyte_init = 'Li(e):0.98, solvent:0.02'
    V_an_init = 1  # Initial anode voltage

    # Cti file info:
    ctifile = 'graphite_anode.cti'
    anodephase = 'graphite'
    elytephase = 'electrolyte'
    interfacephase = 'anode_surf'

    # Cutoff values for lithiation and delithiation of anode:
    X_Li_an_max = 1 - 1e-6
    X_Li_an_min = 1e-6

    # Microstructure:
    phi_an = 0.6   # Graphite volume fraction [-]
    tau_an = 1.6   # Tortuosity - assume equal values for carbon and elyte [-]
    r_p = 5e-6     # Average pore radius [m]
    d_part = 5e-6  # Average particle diameter for graphite [m]
    H = 50e-6      # Anode thickness [m]

    # Other parameters:
    C_dl = 1.5e-2       # Double-layer capacitance [F/m^2]
    sigma_an = 75.0     # Bulk anode electrical conductivity [S/m]
    D_Li_elyte = 1e-10  # Bulk diffusion coefficient for Li+ in elyte [m^2/s]
    D_Li_an = 7.5e-16   # Bulk diffusion coefficient for Li in graphite [m^2/s]

    nVars = nshells + 3

    # Pointers
    ptr = {}

    # Solution vector at each point:
    # kmol of Li per kmol C6 [%] - one for each 'shell'
    # kmol of Li+ per kmol elyte [%]
    # anode electric potential [V]
    # electrolyte electric potential [V]

    ptr['iFar'] = 2    # Location of LiC6 in 'sdot' vector

    ptr['X_an'] = np.arange(0, nshells)
    ptr['X_elyte'] = nshells
    ptr['V_an'] = nshells + 1
    ptr['V_elyte'] = nshells + 2



    # For spherical particles, the total surface area per unit volume can be
    #   calculated from the geometry. Since some particles will overlap, we take
    #   a percentage (60%, below) of this theoretical value:
    A_surf = 0.6*6*phi_an/d_part  # Anode/elyte interface area [m^2/m_3]

    # Create Cantera objects:

    # Graphite anode:
    anode = ct.Solution(ctifile, anodephase)
    anode.X = (X_an_init)

    # Read out init conditions in vector form:
    X_an_0 = anode.X

    # Electrolyte:
    elyte = ct.Solution(ctifile, elytephase)
    elyte.X = X_elyte_init

    # Read out init conditions in vector form:
    X_elyte_0 = elyte.X

    # Graphite/electrolyte interface:
    surf = ct.Interface(ctifile, interfacephase, [anode, elyte])

    # Set up solution vector:
    nSV = npoints*nVars
    SV_0 = np.zeros([nSV])

    for j in range(npoints):
        offset = (j)*nVars
        SV_0[offset + ptr['X_an']] = np.ones([nshells])*X_an_0[0]
        SV_0[offset + ptr['X_elyte']] = X_elyte_0[0]
        SV_0[offset + ptr['V_an']] = V_an_init
        # Electrolyte electric potential = 0, which is already the initialized
        #   value

    # Create dict to hold geometrical parameters:
    geom = {}
    geom['phi_an'] = phi_an
    geom['phi_elyte'] = 1 - phi_an
    geom['tau_an'] = tau_an
    geom['r_p'] = r_p
    geom['d_part'] = d_part
    geom['A_surf'] = A_surf
    geom['dyInv'] = npoints/H

    # Calculate the current density [A/m^2] corresponding to a C_rate of 1:
    oneC = geom['phi_an']*anode.density_mole*H*96485e3/3600

    # Calculate the actual current density:
    # The minus sign is because we begin with the charging reaction, which
    #   delivers negative charge to the anode:
    i_ext = -C_rate*oneC

    # Create dict to hold 'other' parameters

    params = {}
    params['npoints'] = npoints
    params['nshells'] = nshells
    params['nVars'] = nVars
    params['T'] = T
    params['C_dl'] = C_dl
    params['X_Li_max'] = X_Li_an_max
    params['X_Li_min'] = X_Li_an_min
    params['D_Li_an'] = D_Li_an
    params['i_ext'] = i_ext

    # Calculate the percent volume of a single graphite particle that exists in
    #   each 'shell'. I.e. for shell j, what is the volume of that shell,
    #   relative to the total particle volume? The radius of the volume is
    #   currently discretized evenly (i.e. 'dr' the differential radius is
    #   constant). Certainly other discretizations (such as constant
    #   differential volume) are possible (and perhaps better).
    #
    #   Because the volume is 4/3 pi*r^3, the volume of the shell relative to
    #   the total volume is (r_shell/r_particle)^3, and the differential volume
    #   relative to the total, for shell 'j' is:
    #       (r_shell(j)^3 - r_shell(j-1)^3)/r_particle^3
    #   Because the radius is discretized evenly, the radius of shell j, r_j,
    #   relative to the total radius r_particle, is:
    #       r_j/r_particle = j/nshells

    params['V_shell'] = np.ones([nshells, 1])
    params['V_shell'][0] = 1/nshells/nshells/nshells
    for j in np.arange(1, nshells, 1):
        params['V_shell'][j] = ((j + 1)**3 - (j)**3)/nshells/nshells/nshells

    params['sigma_eff_an'] = sigma_an*geom['phi_an']/geom['tau_an']**3
    params['u_Li_elyte'] = D_Li_elyte*geom['phi_elyte']/ct.gas_constant/T/geom['tau_an']**3

    # Set up mass matrix:
    algvar = np.zeros([nSV])
    for j in range(npoints):
        offset = (j)*nVars
        # X_Li_an has a differential equation:
        for i in range(nshells):
            algvar[offset + ptr['X_an'][i]] = 1

        # X_Li_elyte has a differential equation:
        algvar[offset + ptr['X_elyte']] = 1

        # Anode elec potential has a differential equation involving phi_an and
        #   phi_elyte:
        algvar[offset + ptr['V_an']] = 1
        algvar[offset + ptr['V_elyte']] = 1

        # For points 1 to npoints-1, the electrolyte electric potential is governed
        #   by an algebraic equation. For the final point, the electrolyte electric
        #   potential is our reference value, and is held constant at zero. Hence
        #   it is governed by a differential equation, with dSVdt = 0
        if j < npoints - 1:
            algvar[offset + ptr['V_elyte']] = 0
        else:
            algvar[offset + ptr['V_elyte']] = 1

    # Store Cantera objects
    obj = {}
    obj['anode'] = anode
    obj['elyte'] = elyte
    obj['surf'] = surf



def main():

    t_count = time.time()

    # Calculate the time span, which should be enough to charge or discharge fully
    #   (i.e. 3600 seconds divided by the C-rate):
    t_f = 3600/ep.C_rate
    t_0 = 0
    SV_0 = ep.SV_0
    SV_dot_0 = ep.SV_0*0

        # Create problem instance
    Battery_equil = Implicit_Problem(Battery_Func, SV_0, SV_dot_0, t_0)

    Battery_equil.algvar = ep.algvar

    # Simulation parameters
    equil_sim = IDA(Battery_equil)           # Create simulation instance
    equil_sim.atol = 1e-8                  # Solver absolute tolerance
    equil_sim.rtol = 1e-6                  # Solver relative tolerance
    equil_sim.verbosity = 50
    equil_sim.make_consistent('IDA_YA_YDP_INIT')

    # Equilibrate by integrating at zero current:
    ep.params['i_ext'] = 0

    print('\nEquilibrating\n')

    t, SV, SV_dot = equil_sim.simulate(t_f)
#    t_df = pd.DataFrame(t)
    SV_df = pd.DataFrame(SV)
    SV_plot = SV_df.plot()
    SV_plot.legend(loc = 'upper left')
#    equil_sim.plot()


    print('\nDone equilibrating\n')
#    input("Done equilibrating, Press Enter to continue...")
    # New initial conditions are the final equilibrium conditions
    SV_0 = SV[-1, :]
    SV_dot_0 = SV_0*0

    # Charge the battery

    # Create problem instance
    Battery_charge = Implicit_Problem(Battery_Func, SV_0, SV_dot_0, t_0)
    Battery_charge.algvar = ep.algvar

    # Simulation parameters
    charge_sim = IDA(Battery_charge)
    charge_sim.atol = 1e-8
    charge_sim.rtol = 1e-6
    charge_sim.verbosity = 50
    charge_sim.make_consistent('IDA_YA_YDP_INIT')

    print('Charging\n')
    ep.params['i_ext'] = ep.i_ext

    t_charge, SV_charge, SV_dot_charge = charge_sim.simulate(t_f)
    SV_charge_df = pd.DataFrame(SV_charge)
    SV_charge_plot = SV_charge_df.plot()
    SV_charge_plot.legend(loc = 'upper left')
#    charge_sim.plot()
#    input("Done charging, press enter to continue...")
    print('\nDone charging\n')

    data = np.zeros([len(t_charge), ep.nSV + 1])
    data[:, 0] = t_charge
    data[:, 1:] = SV_charge
    np.savetxt('p5_charge.csv', data, delimiter=',')



    # New initial conditions are the final charge conditions
    SV_0 = SV_charge[-1, :]

    # Equilibrate again. Note - this is a specific choice to reflect
    #   equilibration after the charging steps. We may want, at times, to
    #   simulate a situation where the battery is not equilibrated between
    #   charge and discharge, or is equilibrated for a shorter amount of time.

    print('Re-equilibrating\n')

    Battery_re_equil = Implicit_Problem(Battery_Func, SV_0, SV_dot_0, t_0)
    Battery_re_equil.algvar = ep.algvar

    # Simulation parameters
    re_equil_sim = IDA(Battery_charge)
    re_equil_sim.atol = 1e-8
    re_equil_sim.rtol = 1e-6
    re_equil_sim.verbosity = 50
    re_equil_sim.make_consistent('IDA_YA_YDP_INIT')

    ep.params['i_ext'] = 0

    t, SV, SV_dot = re_equil_sim.simulate(t_f)
    SV_df = pd.DataFrame(SV)
    SV_plot = SV_df.plot()
    SV_plot.legend(loc = 'upper left')
#    re_equil_sim.plot()

    print('\nDone equilibrating. Press any key to continue\n')

    print('Discharging...\n')
    ep.params['i_ext'] = -ep.i_ext
    SV_0 = SV[-1, :]

    Battery_discharge = Implicit_Problem(Battery_Func, SV_0, SV_dot_0, t_0)
    Battery_discharge.algvar = ep.algvar

    # Simulation parameters
    Battery_discharge = IDA(Battery_charge)
    Battery_discharge.atol = 1e-8
    Battery_discharge.rtol = 1e-6
    Battery_discharge.verbosity = 50
    Battery_discharge.make_consistent('IDA_YA_YDP_INIT')

    t_discharge, SV_discharge, SV_dot_discharge = Battery_discharge.simulate(t_f)
    SV_discharge_df = pd.DataFrame(SV_discharge)
    SV_discharge_plot = SV_discharge_df.plot()
    SV_discharge_plot.legend(loc = 'upper left')
#    Battery_discharge.plot()

    print('\nDone discharging. Press any key to continue\n')

    # Post processing

    V_charge = SV_charge[:, ep.ptr['V_an']]
    V_discharge = SV_discharge[:, ep.ptr['V_an']]
    t_charge = np.array(t_charge)
    t_discharge = np.array(t_discharge)

    # Plot charge-discharge curve

    Capacity_charge = t_charge*ep.params['i_ext']/3600  # A-h/m^2
    Capacity_discharge = t_discharge*ep.params['i_ext']/3600  # A-h/m^2

#    fig1, ax1 = plt.subplots()
#    ax1.plot(Capacity_charge, V_charge)
#    ax1.plot(Capacity_discharge, V_discharge)

    # Calculate battery energy storage/recovery and calculate round-trip
    #   efficiency. Anode voltage is referenced to its initial equilibrium
    #   value (i.e. in the discharged state).

    # NOTE: This is in W-h/m^2, per unit area of battery. For the specific
    #   capacity, you want W-h/g of anode material.
    E_stored = 0
    E_recovered = 0

    for k in np.arange(1, len(t_charge)):
        E_stored = (E_stored - (ep.V_cathode - 0.5*(V_charge[k] + V_charge[k-1]))
                    *ep.params['i_ext']*(t_charge[k] - t_charge[k-1]))

    for k in np.arange(1, len(t_discharge)):
        E_recovered = (E_recovered - (ep.V_cathode - 0.5*(V_discharge[k] + V_discharge[k-1]))
                        *ep.params['i_ext']*(t_discharge[k] - t_discharge[k-1]))

    Cap_recovered = Capacity_discharge[-1]
    Eta_RT = E_recovered/E_stored

    print('Cap_recovered = ', Cap_recovered, '\n')
    print('Eta_RT = ', Eta_RT, '\n')

    elapsed = time.time() - t_count
    print('t_cpu=', elapsed, '\n')

    return t, SV, SV_dot


if __name__ == "__main__":
    t, SV, SV_dot = main()"""

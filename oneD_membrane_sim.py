"""
1D Direct Contact Membrane Distilation model: 1D_membrane_sim.py

Author: Spencer Gilleon and Steven C. DeCaluwe
Colorado School of Mines


"""
from string import Template
import numpy as np
from scipy.integrate import solve_ivp

""" Set up inputs and givens """
"-----------------------------------------------------------------------------"

def load_model(gas, membrane_data, temp_data, n_points):

    T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = temp_data
    # Geometry and physical constants:
    Press = 101325.0     # Initial gas pressure [Pa]
    T1 = 273.15+T_h;     # Hot feed flow temperature
    T2 = 273.15+T_c;     # Cold permeate flow

    comp_string = Template('N2:1.0,H2O:$PH2O')
    X_k_h = comp_string.substitute(PH2O = P_H2O_h/Press)

    X_k_c = comp_string.substitute(PH2O = P_H2O_c/Press)
    dy = membrane_data['H']/n_points

    # Use gas object to determine gas-phase properties at interfaces:

    # Feed side:
    gas.TPX = T1, Press, X_k_h
    rho_k_h = gas.density*gas.Y

    # Permeate side:
    gas.TPX = T2, Press, X_k_c
    rho_k_c = gas.density*gas.Y

    # Initialize solution vector:
    n_species = gas.n_species
    n_vars = n_species+1
    SV_0 = np.zeros(n_points*n_vars)

    # Starting guess is a linear gradients:
    SV_0[0:n_vars*n_points:n_vars] = np.linspace(T1,T2,n_points)[None,:]
    for k in np.arange(n_species):
        SV_0[k+1:n_vars*n_points+1:n_vars] = np.linspace(rho_k_h[k],rho_k_c[k],n_points)[None,:]

    return SV_0


def run_model(t_sim, SV_0, obj, membrane_params, temp_data, params, trans_flag, method):


    """BOUNDARY CONDITION IS WRONG!!!!  N2 FLUX IS NOT ZERO AT BOUNDARIES!!!!"""
    if trans_flag==0:
        from flux_functions import DGM_fluxes as flux_calc
    elif trans_flag==1:
        from flux_functions import Fick_fluxes as flux_calc

    def dSVdt_func(t,SV,obj,membrane_params,temp_data,params): # System of ODEs to solve
        gas = obj['gas']
        dSVdt = np.zeros(SV.shape)

        n_points = params['n_points']
        dyInv = params['dyInv']
        n_vars = gas.n_species+1

        J_k, wdot, h_k, rhoCvInv, kappa = flux_calc(SV,obj,membrane_params,temp_data,params)

        J_k[1:gas.n_species+1] = J_k[n_vars+1:n_vars+gas.n_species+1]
        #print(J_k)

        q_conv = np.zeros(n_points+1)
        for k in np.arange(gas.n_species):
            dSVdt[n_vars+k+1:n_vars*n_points+1:n_vars] = \
                (J_k[n_vars+k+1:n_vars*n_points+1:n_vars] \
                - J_k[n_vars+k+1+n_vars:n_vars*(n_points+1):n_vars])*dyInv/membrane_params['eps_g']

            q_conv = q_conv + J_k[k::n_vars]*h_k[k::n_vars]

        T_h, T_c, h_fg_h, h_fg_c, P_H2O_h, P_H2O_c, Flux_p45, Flux_p2 = temp_data

        #print(h_fg_h,h_fg_c)
        q_chem = np.zeros(n_points+1)
        q_chem[1] = -h_fg_h*J_k[n_vars+1]
        q_chem[-1] = -h_fg_c*J_k[-2]

        Temps = np.append(np.append(T_h+273.15,SV[0::n_vars]), T_c+273.15)
        q_cond = kappa*(Temps[:-1]-Temps[1:])*dyInv
        q_cond[-1] = q_cond[-1]*2.
        #print(q_cond, q_conv,q_chem)
        q_tot =  q_cond + q_conv + q_chem

        dSVdt[n_vars::n_vars] = (q_tot[1:-1] - q_tot[2:])*dyInv*rhoCvInv[1:]
        #print(dSVdt)
        return dSVdt

    sol = solve_ivp(lambda t, y: dSVdt_func(t, y, obj, membrane_params, temp_data, params), \
        [0, t_sim], SV_0, method=method)

    return sol



def run_model_old():
    #print('hewllo3')
    """ Sub-functions for integrator """
    "-----------------------------------------------------------------------------"
    if Diff_Mod == 0:
        from Ficks_func import Flux_Calc
    elif Diff_Mod == 1:
        from DGM_func import Flux_Calc

    def dSVdt_func(t,SV): # System of ODEs to solve
        dSVdt = np.zeros(Nx*Ny*Nspecies)
        Fluxes_X, Fluxes_Y = Flux_Calc(SV,Nx,dX,Ny,dY,Nspecies,BC_in,inlet_BC,gas,phi_g,tau_g,d_p)
        Fluxes_Y[Ny*Nx*Nspecies:] = J_BC # Constant flux out BC

    #    #print(np.reshape(Fluxes_X[iO2::Nspecies],(Ny,Nx+1)))
    #    #print(np.reshape(Fluxes_Y[iO2::Nspecies],(Ny+1,Nx)))
    #    #print(np.reshape(SV[iO2::Nspecies],(Ny,Nx)))
    #    input('Press Enter to continue')

        # Initialize the fluxes into the first row (y-direction)
        Flux_Y_in = Fluxes_Y[0:Nx*Nspecies]

        # Vector math for inside domain all but last row of cells
        for j in range(Ny):
            ind1 = j*Nx*Nspecies # index for first cell in each row
            ind2 = (j+1)*Nx*Nspecies # index for last cell in each row

            Flux_Y_out = Fluxes_Y[(j+1)*Nx*Nspecies:(j+2)*Nx*Nspecies]

            Flux_X_in = Fluxes_X[j*(Nx+1)*Nspecies:(j+1)*(Nx+1)*Nspecies-Nspecies]
            Flux_X_out = Fluxes_X[j*(Nx+1)*Nspecies+Nspecies:(j+1)*(Nx+1)*Nspecies]

            dSVdt[ind1:ind2] = phi_inv*((Flux_Y_in - Flux_Y_out)*dY_inv \
                             + (Flux_X_in - Flux_X_out)*dX_inv)

            # The fluxes leaving the current row are the inlets to the next row
            Flux_Y_in = Flux_Y_out

        return dSVdt

    """ Call ODE integrator and process results """
    "-----------------------------------------------------------------------------"
    #print('hewllo4')

    sol = solve_ivp(dSVdt_func, [0, t_sim], SV_0, method=method) # Use ODE solver

    plt_species_ind = gas.species_index(plt_species) # Extract index for plot

    # Define variable vectors to be plotted
    x_plt = np.linspace(0, rad+max_space, Nx) * 1e6
    y_plt = -1*np.linspace(0, t_glass, Ny)
    SV = sol.y.T
    sol_t = sol.t
    SV_plt = np.reshape(SV[-1,plt_species_ind::Nspecies], (Ny,Nx))

    # Create contour plot
    plt.imshow(SV_plt,interpolation='nearest')
    #plt.contourf(x_plt,y_plt,SV_plt,
    #             levels=np.linspace(np.min(SV_plt),np.max(SV_plt),Nc))
    plt.title('%s Density [kg/m^3 gas]' %plt_species)
    plt.xlabel('Horizontal distance, x [micron]')
    plt.ylabel('Glass thickness, y [micron]')
    plt.colorbar()
    plt.show()

    # Create movie of entire transient solution over time interval
    if movie == 1:
        from matplotlib.animation import FuncAnimation
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set(xlim=(0, (rad+max_space)*1e6), ylim=(-t_glass, 0))

        div = make_axes_locatable(ax)
        cax = div.append_axes('right', '5%', '5%')

        SV_cb = SV[:,plt_species_ind::Nspecies]

        ax.set_xlabel('Horizontal distance, x [micron]')
        ax.set_ylabel('Glass thickness, y [micron]')

        def animate(i):
            ax.collections = []
            SV_plt = np.reshape(SV[i,plt_species_ind::Nspecies],(Ny,Nx))
            cf = ax.contourf(x_plt,y_plt,SV_plt,
                             levels=np.linspace(np.min(SV_cb),np.max(SV_cb),Nc))
            cax.cla()
            fig.colorbar(cf, cax=cax)
            ax.set_title('Time (s): ' + str(round(sol_t[i])))

        anim = FuncAnimation(fig, animate, interval=10, frames=len(sol_t)-1)

        anim.save('animation.html')


















""" Comments and future steps """
"-----------------------------------------------------------------------------"
# Determine appropriate discretization in time and space and apply for a time
# that allows for a steady-state solution to be reached.

# Exercise model for different geometries.

"""
This script runs over a user-specified set of inputs to fit 1d membrane distillation model results against corresponding experimental data.
"""
# Import necessary modules:
import numpy as np
import timeit
import pandas as pd
import pylab as plt
from datetime import datetime
import time
import cantera as ct

from membrane_distillation_1d_model import membrane_distillation_1d_model
from membrane_distillation_1d_output import ensure_dir


#===============================================================================
#   USER INPUTS
#===============================================================================

# Vector of tortuosities to simulate:
tau_g_vec_0 = np.linspace(0.85,2.0,num=15)

# Specify which membranes to model.  The named membrane should exist in the 
#   file 'membrane_distillation_inputs.yaml'
membranes = ['450nm']

# Specify the feed and permeate temperatures to simulate.  Each tuple should 
# list the feed temperature first, in degrees C.
temps = [(30, 20), (40, 20), (50, 20), (60, 20), (40, 30), (50, 30), (60, 30), 
    (70, 30)]

# These flags specify whether or not to run Fickean and/or DGM transport models:
#   0 for do not run.
#   1 for run.
DGM = 1
Fick = 1

# Plotting parameters
x_UL = 2.0 # upper limit on x-axis for flux vs. tortuosity
fontname = 'Arial'
fontsize = 10

# Folder in which to save results. Date-time will be appended to this
savename = 'tortuosity_fitting'

#===============================================================================
#   VALIDATION DATA
#===============================================================================
temp_data = pd.read_csv('Experimental_and_CFD_data/Exp_data.csv')

#===============================================================================
#   FUNCTION TO LOOP OVER TORTUOSITIES
#===============================================================================
def run_tortuosities(tau_g_vec, membrane, DGM, Fick, SV_0, T_f, T_p):
    data_DGM = np.zeros((SV_0.shape[0]+2, np.size(tau_g_vec_0,0)))
    data_Fick = np.zeros((SV_0.shape[0]+2, np.size(tau_g_vec_0,0)))

    SV_0_DGM = SV_0
    SV_0_Fick = SV_0
    for i_tau, tau_g in enumerate(tau_g_vec):
        print('      Tortusity factor = ',tau_g)
        tic = timeit.default_timer()

        if DGM:
            data = membrane_distillation_1d_model(membrane=membrane,    
                tau_g=tau_g, transport='DGM', save = None, SV_init=SV_0_DGM,feed_temp=T_f, permeate_temp=T_p)
            toc = timeit.default_timer()
            print('         {:.2f} seconds elapsed for DGM.'.format(toc-tic))
            tic = toc
            data_DGM[:,i_tau] = np.insert(data, 0, tau_g)
            SV_0_DGM = data[1:]


        if Fick:
            data = membrane_distillation_1d_model(membrane=membrane, 
                tau_g=tau_g, transport='Fick', save = None, SV_init=SV_0_Fick,feed_temp=T_f, permeate_temp=T_p)
            toc = timeit.default_timer()
            print('         {:.2f} seconds elapsed for Fick.'.format(toc-tic))
            tic = toc
            data_Fick[:,i_tau] = np.insert(data, 0, tau_g)
            SV_0_Fick = data[1:]
    
    return data_Fick, data_DGM

#===============================================================================
#   SET UP AND RUN SIMULATIONS
#===============================================================================
t = datetime.now()
date_string = t.strftime('%m_%d_%Y_%H%M')
save_folder = 'outputs/' + savename + '_' + date_string

# Create a solution vector so we know its size.
from membrane_distillation_1d_init import initialize as init
# The particular parameters don't matter.  We just need SV to be the correct 
# length:
SV_0, _, _ = init(None, 1.0, 'DGM', 300., 275.)

# Set up plotting:
font = plt.matplotlib.font_manager.FontProperties(family=fontname,size=fontsize)

cmap = plt.get_cmap('hsv')
ndata = 11
color_ind = np.linspace(0,1,ndata)
colors = list()

for i in np.arange(ndata):
    colors.append(cmap(color_ind[i]))

colors[2] = '#E5BC3B'#'#8B8B00'
colors[4] = '#215E21'
colors[6] = '#1E90FF'
colors = colors[6:]

# Line width:
lw=1.5

# Loop over the membranes and temperature inputs, and simulate MD for a range 
# of tortuosities at each condition:
for membrane in membranes:
    print(membrane)
    membrane_folder = save_folder + '/' + membrane + '/'
    tau_fitted_DGM = []
    tau_fitted_Fick = []
    
    if DGM:
        f_DGM = plt.figure()
        ax_DGM = f_DGM.add_axes([0.2,0.2,0.75,0.75])
        f_DGM.set_size_inches((3.0,2.25))
    if Fick:
        f_Fick = plt.figure()
        ax_Fick = f_Fick.add_axes([0.2,0.2,0.75,0.75])
        f_Fick.set_size_inches((3.0,2.25))

    tau_Fitting_Fick = np.zeros((len(temps), 3))
    tau_Fitting_DGM = np.zeros((len(temps),3))
    index = 0
    for T_feed, T_perm in temps:
        print('   T_feed = {}, T_perm = {}'.format(T_feed, T_perm))

        if T_perm == 20:
            style = '-'
        elif T_perm == 30:
            style = '--'

        exp_data = temp_data.loc[(temp_data['Feed Temp [C]']==T_feed) & 
            (temp_data['Perm Temp [C]']==T_perm)]
        if membrane=='200nm':
            Flux_exp = exp_data['Flux_p2 [L/m2/h]']
        elif membrane=='450nm':
            Flux_exp = exp_data['Flux_p45 [L/m2/h]']

        SV_0, _, _ = init(membrane, 1.0, 'DGM', T_feed, T_perm)
        # Depending on transport flags, one of these outputs may be all zeros:
        data_Fick, data_DGM = run_tortuosities(tau_g_vec_0, 
            membrane, DGM, Fick, SV_0, T_feed, T_perm)

        # Check to see if the intended folder exists.  If not, create it:
        ensure_dir(membrane_folder)
        save_string = (membrane_folder+'TortuosityVariation_T_f_'+str(T_feed)
            +'_T_p_'+str(T_perm))

        # Calculate the density of the permeate water:
        liq = ct.Water()
        liq.TP = T_perm + 273.15, ct.one_atm
        rho = liq.density

        # Set up plot color:
        dT = T_feed - T_perm
        i_color = int(dT/10  - 1)

        if DGM:
            np.savetxt(save_string+'_DGM.csv', data_DGM, delimiter=',')
            
            Fluxes = 3600*1000*data_DGM[1,:]/rho
            tau_fit = np.interp(Flux_exp, np.flip(Fluxes,0), 
                np.flip(tau_g_vec_0,0))
            tau_Fitting_DGM[index,:] = tau_fit[0], T_feed, T_perm
                        
            ax_DGM.plot(tau_g_vec_0,Fluxes,color=colors[i_color],linewidth=lw,
                linestyle=style)
            ax_DGM.plot(tau_fit[0],Flux_exp,'o',linewidth=lw,
                color=colors[i_color],label=None)

        if Fick:
            np.savetxt(save_string+'_Fick.csv', data_Fick, delimiter=',')
            
            Fluxes = 3600*1000*data_Fick[1,:]/rho
            tau_fit = np.interp(Flux_exp, np.flip(Fluxes,0), 
                np.flip(tau_g_vec_0,0))
            tau_Fitting_Fick[index,:] = tau_fit[0], T_feed, T_perm
                        
            if T_perm == 20:
                ax_Fick.plot(tau_g_vec_0,Fluxes,color=colors[i_color],
                    linewidth=lw, label=r'$\Delta T = $'+str(dT)+r'$ ^\circ$'+'C',linestyle=style)
            else:
                ax_Fick.plot(tau_g_vec_0,Fluxes,color=colors[i_color],
                    linewidth=lw, linestyle=style)

            ax_Fick.plot(tau_fit[0],Flux_exp,'o',linewidth=lw,
                color=colors[i_color],label=None)

        index+=1

    if DGM: 
        np.savetxt(membrane_folder+'tortuosity_vs_Temps_DGM.csv',   
            tau_Fitting_DGM, delimiter=',')

        plots = list()
        plots.append(plt.matplotlib.lines.Line2D([0, 1],[0,1],linestyle='-',
            linewidth=1.5,color='k'))
        plots.append(plt.matplotlib.lines.Line2D([0, 1],[0,1],linestyle='--',
            linewidth=1.5,color='k'))

        tau_20 = []
        tau_30 = []
        for i, tau in enumerate(tau_Fitting_DGM[:,0]):
            if tau_Fitting_DGM[i,2] == 20:
                tau_20.append(tau)
            else:
                tau_30.append(tau)

        tau_20_mu = np.average(tau_20)
        tau_20_std = np.std(tau_20)
        
        tau_30_mu = np.average(tau_30)
        tau_30_std = np.std(tau_30)

        ax_DGM.axvspan(tau_20_mu-tau_20_std, tau_20_mu+tau_20_std,alpha=0.15, 
            color='0.1')
        ax_DGM.axvspan(tau_30_mu-tau_30_std, tau_30_mu+tau_30_std,alpha=0.15, 
            color='0.1')
        
        ax_DGM.set_xlim((1.0,x_UL))
        ax_DGM.set_ylim(0.0,125.)
        ax_DGM.set_xlabel('Tortuosity Factor',fontsize=10)
        ax_DGM.set_ylabel('Water flux (L m'+r'$^{-2}$'+'h'+r'$^{-1}$'+')',
            fontsize=10,labelpad=0.5)

        for tick in ax_DGM.xaxis.get_major_ticks():
            tick.label1.set_fontsize(10)
            tick.label1.set_fontname('Arial')

        ax_DGM.legend(plots,[r'$T_d = 20^\circ$'+'C',r'$T_d = 30^\circ$'+'C'],
            prop=font, borderaxespad=0.5,loc=1,frameon=False,handletextpad=0.2,handlelength=1.5)
        ax_DGM.set_xticks([1.,1.25,1.5,1.75,2.0])
        ax_DGM.set_yticks([0,25,50,75,100,125])
        #ax_DGM.annotate('b) DGM',xy=[1.05,90],color='k',va='bottom',family='Arial',size=10,annotation_clip=False)
        ax_DGM.text(1.03, 110.5, 'b) DGM', fontsize=10,bbox={'facecolor':'white', 'alpha':1.0, 'edgecolor':'1.0'})
        f_DGM.savefig(membrane_folder+'DGM_fig.pdf',fmt='pdf',dpi=350)
    if Fick:
        np.savetxt(membrane_folder+'tortuosity_vs_Temps_Fick.csv', 
            tau_Fitting_Fick, delimiter=',')

        plots = list()
        plots.append(plt.matplotlib.lines.Line2D([0, 1],[0,1],linestyle='-',
            linewidth=1.5,color='k'))
        plots.append(plt.matplotlib.lines.Line2D([0, 1],[0,1],linestyle='--',
            linewidth=1.5,color='k'))
            

        tau_20 = []
        tau_30 = []
        for i, tau in enumerate(tau_Fitting_Fick[:,0]):
            if tau_Fitting_DGM[i,2] == 20:
                tau_20.append(tau)
            else:
                tau_30.append(tau)

        tau_20_mu = np.average(tau_20)
        tau_20_std = np.std(tau_20)
        
        tau_30_mu = np.average(tau_30)
        tau_30_std = np.std(tau_30)

        ax_Fick.axvspan(tau_20_mu-tau_20_std, tau_20_mu+tau_20_std,alpha=0.15, 
            color='0.1')
        ax_Fick.axvspan(tau_30_mu-tau_30_std, tau_30_mu+tau_30_std,alpha=0.15, 
            color='0.1')
        
        ax_Fick.set_xlim((1.0,x_UL))
        ax_Fick.set_ylim(0.0,125.)
        ax_Fick.set_xlabel('Tortuosity Factor',fontsize=10)
        ax_Fick.set_ylabel('Water flux (L m'+r'$^{-2}$'+'h'+r'$^{-1}$'+')',
            fontsize=10,labelpad=0.5)

        for tick in ax_Fick.xaxis.get_major_ticks():
            tick.label1.set_fontsize(10)
            tick.label1.set_fontname('Arial')

        ax_Fick.legend(frameon=False,fontsize=8,borderaxespad=0.25,
            labelspacing=0.2,handletextpad=0.2,handlelength=1.5,loc=1)
        ax_Fick.set_xticks([1.,1.25,1.5,1.75, 2.0])
        ax_Fick.set_yticks([0,25,50,75,100,125])
        ax_Fick.text(1.03, 110.5, 'a) Fick', fontsize=10,bbox=
            {'facecolor':'white', 'alpha':1.0, 'edgecolor':'1.0'})
        f_Fick.savefig(membrane_folder+'Fick_fig.pdf',fmt='pdf',dpi=350)

    # plt.show()
    
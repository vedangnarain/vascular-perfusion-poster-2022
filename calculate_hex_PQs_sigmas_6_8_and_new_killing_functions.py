#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Apr  6 19:06:30 2022

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)
  
Tested in Python 3.7.4.

Calculate metrics of the hexagonal networks with log-normally distributed radii 
and different pruning thresholds with sigma 6 and 8. 
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
#import matplotlib.colors as colors
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import time

# Starts stopwatch to clock execution time
start_time = time.time()

# Set LaTex-style font
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to get the PQ from a file
def get_pq_df(pruning_parameter):
    
    if pruning_parameter=='threshold':
        max_pruning_parameter = max_threshold
        filename = '/home/narain/Desktop/Results/Link to Hexagonal/Log Normal Distribution/Radius Threshold Pruning in Hexagonal Network with 100 Selections and New Killing Function and Beta 6 and 8/TestHexagonalNetwork/hex_lognormal_radius_threshold_perfusion_quotients.txt'
        pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "alpha", "selection", "threshold", "PQ"])
    
    elif pruning_parameter=='kills':
        max_pruning_parameter = max_kills
        filename = '/home/narain/Desktop/Results/Link to Hexagonal/Log Normal Distribution/Individual Pruning in Hexagonal Network with 100 Selections and New Killing Function and Beta 6 and 8/hex_lognormal_individual_perfusion_quotients.txt'
        pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "alpha", "selection", "kills", "PQ"])
    
    elif pruning_parameter=='beta':
        max_pruning_parameter = max_beta
        filename = '/home/narain/Desktop/Results/Link to Hexagonal/Log Normal Distribution/Stochastic Pruning in Hexagonal Network with 100 Selections and New Killing Function and Beta 6 and 8/hex_lognormal_stochastic_perfusion_quotients.txt'
        pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "alpha", "selection", "beta", "PQ"])

    return pq_df, max_pruning_parameter

# Define a function to compute the averages for the betas in the alpha group
def get_alpha_line(alpha_group, max_pruning_parameter, pruning_parameter):
    
    pq_table = np.array([])
    
    if pruning_parameter=='threshold':
        alpha_threshold_grouped = alpha_group.groupby(alpha_group.threshold)
        for pruning_parameter in range(1, max_pruning_parameter+1):
            pq_table_entry = np.array([])
            alpha_threshold_group = alpha_threshold_grouped.get_group(pruning_parameter)
            pq_table_entry = np.array([alpha_threshold_group["alpha"].mean(), alpha_threshold_group["threshold"].mean(), alpha_threshold_group["PQ"].mean(), alpha_threshold_group["PQ"].std()])
            pq_table = np.vstack([pq_table, pq_table_entry]) if pq_table.size else pq_table_entry
            
    elif pruning_parameter=='kills':
        alpha_kills_grouped = alpha_group.groupby(alpha_group.kills)
        for pruning_parameter in range(1, max_pruning_parameter+1):
            pq_table_entry = np.array([])
            alpha_kills_group = alpha_kills_grouped.get_group(pruning_parameter)
            pq_table_entry = np.array([alpha_kills_group["alpha"].mean(), alpha_kills_group["kills"].mean(), alpha_kills_group["PQ"].mean(), alpha_kills_group["PQ"].std()])
            pq_table = np.vstack([pq_table, pq_table_entry]) if pq_table.size else pq_table_entry        
    
    elif pruning_parameter=='beta':
        alpha_beta_grouped = alpha_group.groupby(alpha_group.beta)
        for pruning_parameter in range(1, max_pruning_parameter+1):
            pq_table_entry = np.array([])
            alpha_beta_group = alpha_beta_grouped.get_group(pruning_parameter)
            pq_table_entry = np.array([alpha_beta_group["alpha"].mean(), alpha_beta_group["beta"].mean(), alpha_beta_group["PQ"].mean(), alpha_beta_group["PQ"].std()])
            pq_table = np.vstack([pq_table, pq_table_entry]) if pq_table.size else pq_table_entry

    # Plot the alpha data
#    ax.plot(pq_table[:,0], pq_table[:,1], pq_table[:,2])
    
    # Return the table for reference
    return pq_table

# =============================================================================
# PARAMETERS
# =============================================================================

# Enter details to allow looping over folders
solver_list = ['Constant']
alpha_list = ['6', '8']
max_threshold = 13
threshold_list = [str(x) for x in range(0, max_threshold + 1)]
max_beta = 35
max_kills = 150
kills_list = [str(x) for x in range(max_kills + 1)]
max_trials = 100
max_layouts = 100

# Drop extra data
#max_beta = 35
#pq_df = pq_df.loc[(pq_df["beta"] <= max_beta)]

# Pick the pruning method
pruning_parameter = 'threshold'
#pruning_parameter = 'kills'
#pruning_parameter = 'beta'

# =============================================================================
# PLOTS
# =============================================================================

# Set solver name
solver_name = solver_list[0]

# Get the PQs
pq_df, max_pruning_parameter = get_pq_df(pruning_parameter) 

# Filter PQ data for multiple solvers
solver_filter = solver_name + 'Haematocrit'
pq_df = pq_df.loc[(pq_df["solver_name"] == solver_filter)]

# Separate by alpha 
alpha_grouped = pq_df.groupby(pq_df.alpha)
alpha_1 = alpha_grouped.get_group(6)
alpha_2 = alpha_grouped.get_group(8)

# Compute average of all selections for PQ
line_1 = get_alpha_line(alpha_1, max_pruning_parameter, pruning_parameter)
line_2 = get_alpha_line(alpha_2, max_pruning_parameter, pruning_parameter)

# Combine the PQs
pq_composite = np.hstack([line_1, line_2])

# Set the figure layout
fig, axs = plt.subplots(1, len(alpha_list), figsize=(10, 4), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)

# Set plot title
if pruning_parameter=='threshold':
    plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with radius threshold pruning')
elif pruning_parameter=='kills':
    plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with individual vessel pruning')
elif pruning_parameter=='beta':
    plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with stochastic pruning')

# Plot the distribution stats for a solver
axs = axs.ravel()
for i in range(len(alpha_list)):
    axs[i].set_ylim([0,1.1])  # set PQ limits
    axs[i].plot(line_1[:,1], pq_composite[:, (4*i)+2])  
    axs[i].fill_between(line_1[:,1], pq_composite[:, (4*i)+2]+pq_composite[:, (4*i)+3], pq_composite[:, (4*i)+2]-pq_composite[:, (4*i)+3], color='grey', alpha=0.5, label='PQ ± SD')
    axs[i].set_xlim(0)
    axs[i].set_ylim(0)
    if pruning_parameter=='threshold':
        axs[i].set_xlabel('radius threshold (μm)')    
    elif pruning_parameter=='kills':
        axs[i].set_xlabel('kills (vessels)')    
    elif pruning_parameter=='beta':
        axs[i].set_xlabel('${γ}$ (μm)')    
    if i==0:
        axs[i].set_ylabel('PQ') 
    axs[i].legend(loc="upper right", prop={'size': 13})
#    axs[i].legend()
    axs[i].grid()

# Show plots
plt.show()

'''
# Save image
file_path = Path('~/Desktop/Final Figures/' + solver_name + '_lognormal_hexagonal_radius_threshold_pruning.svg').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
file_path = Path('~/Desktop/Final Figures/' + solver_name + '_lognormal_hexagonal_radius_threshold_pruning.png').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
'''

# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

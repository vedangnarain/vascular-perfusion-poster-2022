#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Jun 27 15:51:06 2021

Plot Perfusion Quotient for Stochastic Killing (Dichotomous Network)

Tested in Python 3.7.4.

"""

# Initialise libraries
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import time

# Starts stopwatch to clock execution time
start_time = time.time()

# Define a function to generate the data for a single alpha value
def plot_alpha_line(alpha_group, max_beta):
    
    # Compute the averages for the betas in the alpha group
    pq_table = np.array([])
    alpha_beta_grouped = alpha_group.groupby(alpha_group.beta)
    for beta_value in range(1, max_beta+1):
        pq_table_entry = np.array([])
        alpha_beta_group = alpha_beta_grouped.get_group(beta_value)
        pq_table_entry = np.array([alpha_beta_group["alpha"].mean(), alpha_beta_group["beta"].mean(), alpha_beta_group["PQ"].mean()])
        pq_table = np.vstack([pq_table, pq_table_entry]) if pq_table.size else pq_table_entry
    
    # Plot the alpha data
    ax.plot(pq_table[:,0], pq_table[:,1], pq_table[:,2])
    
    # Return the table for reference
    return pq_table

# Read PQ file
filename = '/home/narain/Desktop/Data Dump/TestDichotomousNetworkBeta10/perfusion_quotients.txt'
pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "lambda", "alpha", "beta", "trial", "PQ"])

# Drop extra data
max_beta = 10
solver_name = 'Fung'
solver_filter = solver_name + 'Haematocrit/'
pq_df = pq_df.loc[(pq_df["beta"] <= max_beta)]
pq_df = pq_df.loc[(pq_df["solver_name"] == solver_filter)]

# Separate by alpha 
alpha_grouped = pq_df.groupby(pq_df.alpha)
alpha_0 = alpha_grouped.get_group(1.00)
alpha_1 = alpha_grouped.get_group(1.10)
alpha_2 = alpha_grouped.get_group(1.20)
alpha_3 = alpha_grouped.get_group(1.30)
alpha_4 = alpha_grouped.get_group(1.40)

# Plot the hypoxic fractions
fig = plt.figure()
ax = fig.gca(projection='3d')
plt.gca().invert_xaxis()
plt.suptitle(solver_name + ' haematocrit solver in the dichotomous vessel network')
line_0 = plot_alpha_line(alpha_0, max_beta)
line_1 = plot_alpha_line(alpha_1, max_beta)
line_2 = plot_alpha_line(alpha_2, max_beta)
line_3 = plot_alpha_line(alpha_3, max_beta)
line_4 = plot_alpha_line(alpha_4, max_beta)
ax.set_xlabel('α')
ax.set_ylabel('β')
ax.set_zlabel('perfusion quotient')
ax.view_init(elev=40., azim=60)
plt.show() 

# prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

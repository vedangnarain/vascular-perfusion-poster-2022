#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri May 21 16:04:22 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Plot Perfusion Quotient (Dichotomous Network)

Tested in Python 3.7.4.

Version 3: Use contour plots.

"""

# Initialise libraries
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import time

# Starts stopwatch to clock execution time
start_time = time.time()

# Read PQ file
filename = '/home/narain/Chaste/perfusion_quotients.txt'
pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "lambda", "alpha", "radius_threshold", "PQ"])

# Drop extra data
max_radius = 35
pq_df = pq_df.loc[(pq_df["radius_threshold"] <= max_radius)]
pq_df = pq_df.loc[(pq_df["solver_name"] == 'PriesHaematocrit/')]

# Separate by alpha 
alpha_grouped = pq_df.groupby(pq_df.alpha)
alpha_0 = alpha_grouped.get_group(1.00)
alpha_1 = alpha_grouped.get_group(1.10)
alpha_2 = alpha_grouped.get_group(1.20)
alpha_3 = alpha_grouped.get_group(1.30)
alpha_4 = alpha_grouped.get_group(1.40)

# Define a function to plot a 3D plot for the hypoxic fraction
fig = plt.figure()
ax = fig.gca(projection='3d')
plt.gca().invert_xaxis()
plt.suptitle('Constant haematocrit solver in the dichotomous vessel network')
ax.plot(alpha_0['alpha'], alpha_0['radius_threshold'], alpha_0['PQ'])
ax.plot(alpha_1['alpha'], alpha_1['radius_threshold'], alpha_1['PQ'])
ax.plot(alpha_2['alpha'], alpha_2['radius_threshold'], alpha_2['PQ'])
ax.plot(alpha_3['alpha'], alpha_3['radius_threshold'], alpha_3['PQ'])
ax.plot(alpha_4['alpha'], alpha_4['radius_threshold'], alpha_4['PQ'])
ax.set_xlabel('α')
ax.set_ylabel('radius threshold (μm)')
ax.set_zlabel('perfusion quotient')
ax.view_init(elev=40., azim=60)
plt.show() 

# prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

# =============================================================================
# SNIPPETS
# =============================================================================

'''
fig = plt.figure()
plt.suptitle('Constant haematocrit solver in the dichotomous vessel network (α = 1.4)')
plt.plot(alpha_4['radius_threshold'], alpha_4['PQ'], label='PQ')
plt.plot(alpha_4['radius_threshold'], constant_data[:,14], label='HF')
plt.xlabel('radius threshold (μm)')
plt.ylabel('HF/PQ')
plt.legend()
'''
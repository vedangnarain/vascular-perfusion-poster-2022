#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  1 23:23:54 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)
  
Calculate Hypoxic Fraction (Dichotomous Network)

Tested in Python 3.7.4.

Version 6: Compute for hexagonal networks.

"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
from collections import deque
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import scipy.interpolate
import time
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# Starts stopwatch to clock execution time
start_time = time.time()

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read a .vti file and return the field distribution
def get_vti_data(input_filename):
    
    # Import the file
    extension = input_filename.split('.').pop()
    reader = None
    if extension == 'vtk':
        reader = vtk.vtkDataSetReader() 
    elif extension == 'vti':
        reader = vtk.vtkXMLImageDataReader() 
    else:
        raise RuntimeError('Unknown File Type: %s ' % input_filename)
    reader.SetFileName( "%s" % (input_filename) ) 
    reader.Update()
    image_data = reader.GetOutput()
    
    # Extract the dimensions, spacing, and origin
    spacing = image_data.GetSpacing()
    
    # Extract the point values 
    field_point_data = image_data.GetPointData() 
    field_values = vtk_to_numpy(field_point_data.GetArray(0)) 
    
    # Get the coordinates of each point
    position_list = deque()
    for index in range(len(field_values)):  # do something 
        position = image_data.GetPoint(index)
        position_list.append(position)
    position_array = np.array(position_list)
    
    # Return the field distribution
    distribution_array = np.column_stack((position_array,field_values))
    distribution = pd.DataFrame(data=distribution_array[:,[0,1,3]], columns=["x", "y", "oxygen"])
    return distribution, spacing#, dimensions, origin

# Define a function to convert the .vti data to a plottable field with some basic stats
def get_field(paraview_data):
    
    # Change the data type
    paraview_data = paraview_data.astype('float32')
    
    # Calculate concentration statistics
    O2_stats = paraview_data['oxygen'].describe()
    
    # Downsample data if needed
    paraview_data = paraview_data[::int(1)]
    
    # Convert dataframe into NumPy matrix
    mat = paraview_data.to_numpy()
    
    # Get the x and y axes
    x = np.unique(mat[:,0])
    y = np.unique(mat[:,1])
    
    # Create a mesh from the x and y axes 
    X,Y = np.meshgrid(x, y)
    
    # Interpolate the concentration values over the mesh
    Z = scipy.interpolate.griddata((mat[:,0], mat[:,1]), mat[:,2], (X,Y), method='nearest')

    # Return the oxygen stats
    return Z, O2_stats

# Define a function to extract the middle distribution of the field
def get_hex_domain(field_dataset):
#    reader = vtk.vtkXMLPolyDataReader()
#    reader.SetFileName(node_file_path)
#    reader.Update()
#    polydata = reader.GetOutput()
#    points = polydata.GetPoints()
#    array = points.GetData()
#    point_coordinates = vtk_to_numpy(array)
    x_start = 100
    x_end = 1950
    y_start = 170#3.205
    y_end = 1730#2.05
    field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end)) & (field_dataset['y'] >= int(y_start)) & (field_dataset['y'] <= int(y_end))]
    return field_dataset
#    return vessel_network, oxygen_distribution, middle_x
    
# Define a function to read an O2 distribution .vti file and return the hypoxic fraction and other stats
def get_hypoxic_fraction(solver_name, layout_selection, nhet, kills, hypoxic_threshold, plot=0):    
    
# Set the file path
#    folder_path = '/home/narain/Desktop/Stochastic Pruning with 100 Trials/' + solver_name + 'Haematocrit/Lambda4/Alpha' + alpha_value + '/Beta' + beta + '/Trial' + trial
    folder_path = '/home/narain/Desktop/Deterministic Pruning/TestHexagonalNetwork/' + solver_name + 'Haematocrit/Selection' + layout_selection + '/NHet' + nhet + '/Kills' + kills
    field_path =  folder_path + '/oxygen_solution_0.vti'    
#    network_path = '/home/narain/Desktop/ReferenceNetworks/Alpha' + alpha_value + '/FinalHaematocrit.vtp'

    # Print status update
    print(field_path)

    # Import the data for the network
    field_data, field_spacing = get_vti_data(field_path)
    
    # Get the distribution in the middle of the field (replace with designated file)
    middle_field = get_hex_domain(field_data)

    # Get the basic stats 
    middle_O2_stats = middle_field['oxygen'].describe()
    
    # Plot the O2 distribution
    if plot==1:
        # Get the O2 mesh for plots and stats
        O2_field, _ = get_field(middle_field)
        fig = plt.figure()
        ax = plt.axes()
        colour_map = plt.cm.get_cmap('jet')
        plt.suptitle('O$_2$ distribution generated by the ' + solver_name + ' solver in the hexagonal vessel network with thin vessel selection = ' + layout_selection + ' (α = ' + nhet + ' and kills = ' + kills + '  vessels (1 unit = ' + str(field_spacing [0]) + ' μm)')
        ref_map = ax.imshow(O2_field, cmap=colour_map, origin='lower')
        fig.colorbar(ref_map, ax=ax, label='nM')
        plt.show()
    
    # Calculate the number of points below hypoxic threshold    
    hypoxic_points = (middle_field['oxygen'] < hypoxic_threshold).sum()
    
    # Get the number of total points
    number_of_points = 1
    for dim in np.shape(middle_field['oxygen']): number_of_points *= dim
    
    # Calculate the hypoxic fraction
    hypoxic_fraction = hypoxic_points/number_of_points
    
    # Return the hypoxic fraction
    return hypoxic_fraction, middle_O2_stats['mean'], middle_O2_stats['min'], middle_O2_stats['50%'], middle_O2_stats['max'], middle_O2_stats['std'],

# Define function to compute the average HT of all trials in a kill selection
def compute_average_trial(solver_name, nhet, kills, hypoxic_threshold, max_layouts):
        
    # Create table to store all the kill data trials in an alpha group
    kill_table = np.array([])

    # Extract metrics from all trials and store in beta table
    for layout_selection in range(1, max_layouts+1):    
#        print(layout_selection)
        hypoxic_fraction, mean_value, min_value, half_value, max_value, std_value = get_hypoxic_fraction(solver_name, str(layout_selection), nhet, str(kills), hypoxic_threshold, plot=0)
        table_entry = np.array([hypoxic_fraction, mean_value, min_value, half_value, max_value, std_value])
        kill_table = np.vstack([kill_table, table_entry]) if kill_table.size else table_entry

    # Return the hypoxic fraction, mean, min, 50%, max, and std averaged across all trials
    return np.average(kill_table[:,0]), np.average(kill_table[:,1]), np.average(kill_table[:,2]), np.average(kill_table[:,3]), np.average(kill_table[:,4]), np.average(kill_table[:,5])
#    return np.average(kill_table[0]), np.average(kill_table[1]), np.average(kill_table[2]), np.average(kill_table[3]), np.average(kill_table[4])

# Define a function to return the hypoxic fraction data for all the simulations in a solver
def get_hypoxic_fraction_stats(solver_name, alpha_list, kill_list, max_layouts, hypoxic_threshold):
#    table = np.array([])
    alpha_table = np.array([])
    for nhet in alpha_list:
        for kill in kill_list: 
            hypoxic_fraction, mean_value, min_value, half_value, max_value, std_value = compute_average_trial(solver_name, nhet, kill, hypoxic_threshold, max_layouts)
            table_entry = np.array([float(nhet), float(kill), hypoxic_fraction, mean_value, min_value, half_value, max_value, std_value])
            alpha_table = np.vstack([alpha_table, table_entry]) if alpha_table.size else table_entry
#        table = np.vstack([table, alpha_table]) if table.size else alpha_table
    return alpha_table

# Define a function to get the data for a contour plot for the hypoxic fraction
def get_xyz(solver_data, solver_name, hypoxic_threshold):

    # Prepare arrays
    alphas = solver_data[:,0]
    betas = solver_data[:,1]
    oxygen_points = solver_data[:,2]  
    
    # Get the x and y axes
    x = np.unique(alphas)
    y = np.unique(betas)
    
    # Create a mesh from the x and y axes 
    X,Y = np.meshgrid(x, y)
    
    # Interpolate the concentration values over the mesh
    Z = scipy.interpolate.griddata((alphas, betas), oxygen_points, (X,Y), method='nearest')
    
    # Return the arrays
    return X, Y, Z

'''
# Define a function to plot a comparison of the PQ and HF
def plot_simulation_stats(solver_name, alpha_df, ht_data, fixed_alpha_value, fixed_hypoxic_threshold):
    
    # Group the PQ data by alpha
    alpha_grouped = alpha_df.groupby(alpha_df.alpha)
    alpha_group = alpha_grouped.get_group(fixed_alpha_value)
    
    # Filter the other stats by alpha
    ht_array = ht_data[(ht_data[:,0]==fixed_alpha_value)]
    hf_data = ht_array[:,2]
    mean_data = ht_array[:,3]
    min_data = ht_array[:,4]
    half_data = ht_array[:,5]
    max_data = ht_array[:,6]
    sd_data = ht_array[:,7]
    
    # Plot the comparison
    fig, ax1 = plt.subplots()
    plt.suptitle(solver_name + ' haematocrit solver in the dichotomous vessel network (α = ' + str(fixed_alpha_value) + ', hypoxic threshold = ' + str(int(fixed_hypoxic_threshold)) + ' nM)')

    color = 'tab:red'
    ax1.set_xlabel('radius threshold (μm)')
    ax1.set_ylabel('HF/PQ', color=color)
    pq, = ax1.plot(alpha_group['radius_threshold'], alpha_group['PQ'], color=color, ls='--', label='PQ')
    hf, = ax1.plot(alpha_group['radius_threshold'], hf_data, color=color, label='HF')
    ax1.tick_params(axis='y', labelcolor=color)
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    color = 'tab:blue'
    ax2.set_ylabel('O$_2$ Concentration (nM)', color=color)  # we already handled the x-label with ax1
    ax2.plot(alpha_group['radius_threshold'], mean_data, color=color, ls='solid', label='mean')
    ax2.plot(alpha_group['radius_threshold'], min_data, color=color, ls='dotted', label='min')
    ax2.plot(alpha_group['radius_threshold'], half_data, color=color, ls='dashed', label='50%')
    ax2.plot(alpha_group['radius_threshold'], max_data, color=color, ls='dashdot', label='max')
    ax2.plot(alpha_group['radius_threshold'], sd_data, color=color, ls=(0, (3, 5, 1, 5)), label='SD')
    ax2.tick_params(axis='y', labelcolor=color)
    
#    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    ax1.legend(loc='center left')
    ax2.legend(loc='center right')
    plt.show()
'''

# =============================================================================
# HYPOXIC FRACTION
# =============================================================================

# Enter details to allow looping over folders
solver_list = ['Constant', 'Pries', 'Memory', 'Fung']
alpha_list = ['0', '1', '2', '3', '4']
max_kills = 77
max_layouts = 100
kill_list = [str(x) for x in range(0, max_kills + 1)]
hypoxic_threshold_list = [0, 1000, 5000, 10000, 20000, 28000] 

# Set solver name
solver_name = solver_list[0]

# Normalise the colourbar
lower_bound = 0  # get the overall min.
upper_bound = 1  # get the overall max.
resolution = 10  # set the resolution
bounds = np.linspace(lower_bound, upper_bound, 10)  # set a range and resolution
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)  # normalise    

# Set the figure layout
fig, axs = plt.subplots(2, 3, figsize=(20, 10))
fig.subplots_adjust(hspace = .5, wspace=.25)
plt.suptitle(solver_name + ' haematocrit solver in the hexagonal vessel network')

# Plot the contour maps
axs = axs.ravel()
for i in range(len(hypoxic_threshold_list)):
    solver_data = get_hypoxic_fraction_stats(solver_name, alpha_list, kill_list, max_layouts, hypoxic_threshold_list[i])
    X, Y, Z = get_xyz(solver_data, solver_name, hypoxic_threshold_list[i])
    if i == int(len(hypoxic_threshold_list)/1):  # use the last plot to map the colourbar
        ref_map = axs[i].contourf(X, Y, Z, resolution, cmap='RdGy', norm=norm, origin='lower') 
    else:
        axs[i].contourf(X, Y, Z, resolution, cmap='RdGy', norm=norm, origin='lower')
    axs[i].set_title('hypoxic threshold = ' + str(int(hypoxic_threshold_list[i])) + ' nM')
    axs[i].set_xlabel('α')
    axs[i].set_ylabel('vessels killed')
fig.colorbar(ref_map, ax=axs, label='hypoxic fraction')
plt.show()

# =============================================================================
# PERFUSION QUOTIENT
# =============================================================================
'''
# Read PQ file
filename = '/home/narain/Chaste/perfusion_quotients.txt'
pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "lambda", "alpha", "radius_threshold", "PQ"])

# Drop extra data
max_radius = 35
solver_filter = solver_name + 'Haematocrit/'
pq_df = pq_df.loc[(pq_df["radius_threshold"] <= max_radius)]
pq_df = pq_df.loc[(pq_df["solver_name"] == solver_filter)]

# Separate by alpha 
alpha_grouped = pq_df.groupby(pq_df.alpha)
alpha_0 = alpha_grouped.get_group(1.00)
alpha_1 = alpha_grouped.get_group(1.10)
alpha_2 = alpha_grouped.get_group(1.20)
alpha_3 = alpha_grouped.get_group(1.30)
alpha_4 = alpha_grouped.get_group(1.40)
alpha_5 = alpha_grouped.get_group(1.50)
'''

# =============================================================================
# COMPARE HF & PQ
# =============================================================================
'''
# Enter choice of simulation
fixed_hypoxic_threshold = 5000  #O2_stats['25%'] # set hypoxia threshold (<0.1% mmHg is radioresistance, <2% is hypoxia)
fixed_alpha_value = 1.5

# Get data for specific hypoxic threshold
ht_data = get_hypoxic_fraction_stats(solver_name, alpha_list, radius_threshold_list, fixed_hypoxic_threshold)

# Plot the PQ and HF for a simulation
plot_simulation_stats(solver_name, pq_df, ht_data, fixed_alpha_value, fixed_hypoxic_threshold)
'''

# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Jun  6 21:05:03 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Created on Wed May 26 11:04:33 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)
  
Calculate Additional Metrics (Dichotomous Network)

Tested in Python 3.7.4.

Version 1: Get some metrics for different solvers.
"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
from collections import deque
import matplotlib.colors as colors
import matplotlib.pyplot as plt
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

# Define a function to read a .vti file
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
def get_middle_portion(node_file_path, field_dataset, middle_generation_number):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(node_file_path)
    reader.Update()
    polydata = reader.GetOutput()
    points = polydata.GetPoints()
    array = points.GetData()
    point_coordinates = vtk_to_numpy(array)
    generation_coordinates = np.unique(point_coordinates[:,0])
    x_start = generation_coordinates[middle_generation_number-1]
    x_end = generation_coordinates[-middle_generation_number]
    field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end))]
    return field_dataset
#    return vessel_network, oxygen_distribution, middle_x

# Define a function to read an O2 distribution .vti file and return the basic stats
def get_distribution_stats(solver_name, alpha_value, radius_threshold, plot=0):    
    
    # Set the file path
    folder_path = '/Users/vedang/Simulations/Haematocrit Solver Comparisons/TestDichotomousNetwork4SolversAlpha1.4/' + solver_name + 'Haematocrit/Alpha' + alpha_value + '/LambdaEquals4/RadiusThresholdEquals'
    field_path =  folder_path + radius_threshold + 'um/oxygen_solution_0.vti'    
    network_path = folder_path + '0um/FinalHaematocrit.vtp'

    # Import the data for the network
    field_data, field_spacing = get_vti_data(field_path)
    
    # Get the distribution in the middle of the field (replace with designated file)
    middle_field = get_middle_portion(network_path, field_data, 6)

    # Get the basic stats 
    middle_O2_stats = middle_field['oxygen'].describe()
    
    # Plot the O2 distribution
    if plot==1:
        # Get the O2 mesh for plots and stats
        O2_field, _ = get_field(middle_field)
        fig = plt.figure()
        ax = plt.axes()
        colour_map = plt.cm.get_cmap('jet')
        plt.suptitle('O$_2$ distribution generated by the ' + solver_name + ' solver in the dichotomous vessel network with α = ' + alpha_value + ' and radius threshold = ' + radius_threshold + '  μm (1 unit = ' + str(field_spacing [0]) + ' μm)')
        ref_map = ax.imshow(O2_field, cmap=colour_map, origin='lower')
        fig.colorbar(ref_map, ax=ax, label='nM')
        plt.show()

    # Return the stats
    return middle_O2_stats['mean'], middle_O2_stats['min'], middle_O2_stats['50%'], middle_O2_stats['max'], middle_O2_stats['std'],

# Define a function to return statistics for all the simulations in a solver
def get_solver_stats(solver_name, alpha_list, radius_threshold_list):
    table = np.array([])
    for alpha_value in alpha_list:
        alpha_table = np.array([])
        for radius_threshold in radius_threshold_list:
            mean_value, min_value, half_value, max_value, std_value = get_distribution_stats(solver_name, alpha_value, radius_threshold)
            table_entry = np.array([float(alpha_value), float(radius_threshold), mean_value, min_value, half_value, max_value, std_value])
            alpha_table = np.vstack([alpha_table, table_entry]) if alpha_table.size else table_entry
        table = np.vstack([table, alpha_table]) if table.size else alpha_table
    return table

# =============================================================================
# BASIC STATS
# =============================================================================

# Enter details to allow looping over folders
solver_list = ['Constant', 'Pries', 'Memory', 'Betteridge']
alpha_list = ['1', '1.1', '1.2', '1.3', '1.4']
max_radius = 35
radius_threshold_list = [str(x) for x in range(max_radius + 1)]

# Set solver name
solver_name = solver_list[2]

# Get the stats for all solvers
solver_stats = get_solver_stats(solver_name, alpha_list, radius_threshold_list)

# Filter by alpha
mean_composite = np.array([])
min_composite = np.array([])
half_composite = np.array([])
max_composite = np.array([])
sd_composite = np.array([])
for alpha_value in alpha_list:
    alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
    mean_data = alpha_array[:,2]
    min_data = alpha_array[:,3]
    half_data = alpha_array[:,4]
    max_data = alpha_array[:,5]
    sd_data = alpha_array[:,6]
    mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
    min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
    half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
    max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
    sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data

# Set the figure layout
fig, axs = plt.subplots(2, len(alpha_list), figsize=(20, 10))
fig.subplots_adjust(hspace = .5, wspace=.25)
plt.suptitle(solver_name + ' haematocrit solver in the dichotomous vessel network')

# Plot the stats for a solver
axs = axs.ravel()
for i in range(len(alpha_list)):
    axs[i].plot(alpha_array[:,1], mean_composite[i], ls='dashed', label='mean')
    axs[i].plot(alpha_array[:,1], min_composite[i], ls='dotted', label='min')
    axs[i].plot(alpha_array[:,1], half_composite[i], ls=(0, (3, 5, 1, 5)), label='50%')
    axs[i].plot(alpha_array[:,1], max_composite[i], ls='dashdot', label='max')
    axs[i].plot(alpha_array[:,1], sd_composite[i], ls='solid', label='SD')
    axs[i].set_xlabel('radius threshold (μm)')    
    if i==0:
        axs[i].set_ylabel('O$_2$ Concentration (nM)') 
    axs[i].legend()
    axs[i].grid()
    axs[i].title.set_text('α = ' + alpha_list[i])
#plt.show()

# =============================================================================
# PERFUSION QUOTIENT
# =============================================================================

# Read PQ file
filename = '/Users/vedang/Simulations/Haematocrit Solver Comparisons/perfusion_quotients.txt'
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
#alpha_5 = alpha_grouped.get_group(1.50)
alpha_composite = np.vstack([alpha_0['PQ'], alpha_1['PQ'], alpha_2['PQ'], alpha_3['PQ'], alpha_4['PQ']])
#'''

# Plot the PQ for a solver
for i in range(len(alpha_list),len(alpha_list)*2):
    axs[i].plot(alpha_0['radius_threshold'], alpha_composite[i-len(alpha_list)], label='PQ')
    axs[i].set_xlabel('radius threshold (μm)')    
    if i==len(alpha_list):
        axs[i].set_ylabel('perfusion quotient') 
    axs[i].legend()
    axs[i].grid()
#    axs[i].title.set_text('α = ' + alpha_list[i])
plt.show()

# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

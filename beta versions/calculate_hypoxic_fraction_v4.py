#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 11:04:33 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)
  
Calculate Hypoxic Fraction (Dichotomous Network)

Tested in Python 3.7.4.

Version 4: Use a single solver for a greater resolution of thresholds. 

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

# Define a function to read an O2 distribution .vti file and return the hypoxic fraction and other stats
def get_hypoxic_fraction(solver_name, alpha_value, radius_threshold, hypoxic_threshold, plot=0):    
    
    # Set the file path
#    folder_path = '/home/narain/Desktop/Data Dump/TestDichotomousNetwork4SolversAlpha1.4/' + solver_name + 'Haematocrit/Alpha' + alpha_value + '/LambdaEquals4/RadiusThresholdEquals'
    folder_path = '/tmp/narain/testoutput/TestDichotomousNetwork/' + solver_name + 'Haematocrit/Lambda4/Alpha' + alpha_value + '/RadiusThreshold'
    field_path =  folder_path + radius_threshold + '/oxygen_solution_0.vti'    
    network_path = folder_path + '0/FinalHaematocrit.vtp'

    # Import the data for the network
    field_data, field_spacing = get_vti_data(field_path)
    
    # Get the distribution in the middle of the field (replace with designated file)
    middle_field = get_middle_portion(network_path, field_data, 6)

    # Get some basic stats 
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
    
    # Calculate the number of points below hypoxic threshold    
    hypoxic_points = (middle_field['oxygen'] < hypoxic_threshold).sum()
    
    # Get the number of total points
    number_of_points = 1
    for dim in np.shape(middle_field['oxygen']): number_of_points *= dim
    
    # Calculate the hypoxic fraction
    hypoxic_fraction = hypoxic_points/number_of_points
    
    # Return the hypoxic fraction
    return hypoxic_fraction, middle_O2_stats['mean'], middle_O2_stats['min'], middle_O2_stats['50%'], middle_O2_stats['max'], middle_O2_stats['std'],

# Define a function to return the hypoxic fraction data for all the simulations in a solver
def get_hypoxic_fraction_stats(solver_name, alpha_list, radius_threshold_list, hypoxic_threshold):
    table = np.array([])
    for alpha_value in alpha_list:
        alpha_table = np.array([])
        for radius_threshold in radius_threshold_list:
            hypoxic_fraction, mean_value, min_value, half_value, max_value, std_value = get_hypoxic_fraction(solver_name, alpha_value, radius_threshold, hypoxic_threshold)
            table_entry = np.array([float(alpha_value), float(radius_threshold), hypoxic_fraction, mean_value, min_value, half_value, max_value, std_value])
            alpha_table = np.vstack([alpha_table, table_entry]) if alpha_table.size else table_entry
        table = np.vstack([table, alpha_table]) if table.size else alpha_table
    return table

# Define a function to get the data for a contour plot for the hypoxic fraction
def get_xyz(solver_data, solver_name, hypoxic_threshold):

    # Prepare arrays
    alphas = solver_data[:,0]
    thresholds = solver_data[:,1]
    oxygen_points = solver_data[:,2]  
    
    # Get the x and y axes
    x = np.unique(alphas)
    y = np.unique(thresholds)
    
    # Create a mesh from the x and y axes 
    X,Y = np.meshgrid(x, y)
    
    # Interpolate the concentration values over the mesh
    Z = scipy.interpolate.griddata((alphas, thresholds), oxygen_points, (X,Y), method='nearest')
    
    # Return the arrays
    return X, Y, Z

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

# =============================================================================
# HYPOXIC FRACTION
# =============================================================================

# Enter details to allow looping over folders
solver_name = 'Constant'
max_alpha = 1.5
n_alphas = 6
max_radius = 35  #50
max_hypoxic_threshold = 30000

# Convert alpha list to format of file system
#alpha_names = [int(num) if float(num).is_integer() else round(num,2) for num in np.linspace(1, max_alpha, num=n_alphas)]
alpha_names = ["{:.2f}".format(num) for num in np.linspace(1, max_alpha, num=n_alphas)]
alpha_list = [str(x) for x in alpha_names]
radius_threshold_list = [str(x) for x in range(max_radius + 1)]
#hypoxic_threshold = 7331.5  #O2_stats['25%'] # set hypoxia threshold (<0.1% mmHg is radioresistance, <2% is hypoxia)
hypoxic_threshold_list = np.linspace(0, max_hypoxic_threshold, num=12)

#alpha_list = ['1', '1.1', '1.2', '1.3', '1.35']
#hypoxic_threshold_list = [0, 1000, 5000, 10000, 20000, 28000] 

# Normalise the colourbar
lower_bound = 0  # get the overall min.
upper_bound = 1  # get the overall max.
resolution = 10  # set the resolution
bounds = np.linspace(lower_bound, upper_bound, 10)  # set a range and resolution
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)  # normalise    
'''
# Set the figure layout
fig, axs = plt.subplots(3, 4, figsize=(20, 10))
fig.subplots_adjust(hspace = .5, wspace=.25)
plt.suptitle(solver_name + ' haematocrit solver in the dichotomous vessel network')

# Plot the contour maps
axs = axs.ravel()
for i in range(len(hypoxic_threshold_list)):
    solver_data = get_hypoxic_fraction_stats(solver_name, alpha_list, radius_threshold_list, hypoxic_threshold_list[i])
    X, Y, Z = get_xyz(solver_data, solver_name, hypoxic_threshold_list[i])
    if i == int(len(hypoxic_threshold_list)/2):  # use the middle plot to map the colourbar
        ref_map = axs[i].contourf(X, Y, Z, resolution, cmap='RdGy', norm=norm, origin='lower') 
    else:
        axs[i].contourf(X, Y, Z, resolution, cmap='RdGy', norm=norm, origin='lower')
    axs[i].set_title('hypoxic threshold = ' + str(int(hypoxic_threshold_list[i])) + ' nM')
    axs[i].set_xlabel('α')
    axs[i].set_ylabel('radius threshold (μm)')
fig.colorbar(ref_map, ax=axs, label='hypoxic fraction')
plt.show()
'''
# =============================================================================
# PERFUSION QUOTIENT
# =============================================================================

# Read PQ file
filename = '/home/narain/Chaste/perfusion_quotients.txt'
pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "lambda", "alpha", "radius_threshold", "PQ"])

# Drop extra data
max_radius = 35
solver_filter = solver_name + 'Haematocrit/'
pq_df = pq_df.loc[(pq_df["radius_threshold"] <= max_radius)]
pq_df = pq_df.loc[(pq_df["solver_name"] == solver_filter)]

'''
# Separate by alpha 
alpha_grouped = pq_df.groupby(pq_df.alpha)
alpha_0 = alpha_grouped.get_group(1.00)
alpha_1 = alpha_grouped.get_group(1.10)
alpha_2 = alpha_grouped.get_group(1.20)
alpha_3 = alpha_grouped.get_group(1.30)
alpha_4 = alpha_grouped.get_group(1.40)
alpha_5 = alpha_grouped.get_group(1.50)

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
ax.plot(alpha_5['alpha'], alpha_5['radius_threshold'], alpha_5['PQ'])
ax.set_xlabel('α')
ax.set_ylabel('radius threshold (μm)')
ax.set_zlabel('perfusion quotient')
ax.view_init(elev=40., azim=60)
plt.show() 
'''

# =============================================================================
# COMPARE HF & PQ
# =============================================================================

# Enter choice of simulation
fixed_hypoxic_threshold = 5000  #O2_stats['25%'] # set hypoxia threshold (<0.1% mmHg is radioresistance, <2% is hypoxia)
fixed_alpha_value = 1.5

# Get data for specific hypoxic threshold
ht_data = get_hypoxic_fraction_stats(solver_name, alpha_list, radius_threshold_list, fixed_hypoxic_threshold)

# Plot the PQ and HF for a simulation
plot_simulation_stats(solver_name, pq_df, ht_data, fixed_alpha_value, fixed_hypoxic_threshold)




























# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

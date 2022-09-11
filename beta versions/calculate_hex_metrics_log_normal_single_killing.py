#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 16:13:11 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)
  
Calculate Metrics of Hexagonal Network with thick-thin ordered killing.

Tested in Python 3.7.4.

Calculate Metrics of Hexagonal Network with log-normally distributed radii and
single vessel killing. The metrics are: basic distribution stats, the perfusion 
quotient, and the hypoxic fraction.

Used to generate figure in Transfer Report.

"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
from collections import deque
#import matplotlib.colors as colors
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import scipy.interpolate
import time
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# Starts stopwatch to clock execution time
start_time = time.time()

# Set LaTex-style font
from pathlib import Path
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 22})

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read a .vti file and return the data
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

## Define a function to extract the middle distribution of the field
#def get_middle_portion(node_file_path, field_dataset, middle_generation_number):
#    reader = vtk.vtkXMLPolyDataReader()
#    reader.SetFileName(node_file_path)
#    reader.Update()
#    polydata = reader.GetOutput()
#    points = polydata.GetPoints()
#    array = points.GetData()
#    point_coordinates = vtk_to_numpy(array)
#    generation_coordinates = np.unique(point_coordinates[:,0])
#    x_start = generation_coordinates[middle_generation_number-1]
#    x_end = generation_coordinates[-middle_generation_number]
#    field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end))]
#    return field_dataset
#    return vessel_network, oxygen_distribution, middle_x

# Define a function to extract the densely vascularised portion of the field
def get_hex_domain(field_dataset, field_spacing):
#    reader = vtk.vtkXMLPolyDataReader()
#    reader.SetFileName(node_file_path)
#    reader.Update()
#    polydata = reader.GetOutput()
#    points = polydata.GetPoints()
#    array = points.GetData()
#    point_coordinates = vtk_to_numpy(array)
    x_start = 10*field_spacing[0]
    x_end = 195*field_spacing[0]
    y_start = 17*field_spacing[1]
    y_end = 173*field_spacing[1]
    field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end)) & (field_dataset['y'] >= int(y_start)) & (field_dataset['y'] <= int(y_end))]
    return field_dataset
#    return vessel_network, oxygen_distribution, middle_x

# Define a function to read an O2 distribution .vti file and return the basic stats
def get_distribution_stats(solver_name, layout_selection, sigma, kills, hypoxic_threshold_list, plot=0, read=0):    
    
    # Set the file path
#    folder_path = '/home/narain/Desktop/Stochastic Pruning with 100 Trials/' + solver_name + 'Haematocrit/Lambda4/Alpha' + alpha_value + '/Beta' + beta + '/Trial' + trial
    folder_path = '/home/narain/Desktop/Results/Hexagonal/Log Normal Distribution/Kills Pruning in Hexagonal Network with 100 Selections/TestHexagonalNetwork/' + solver_name + 'Haematocrit/Sigma' + sigma + '/Selection' + layout_selection + '/Kills' + kills
    field_path =  folder_path + '/oxygen_solution_0.vti'    
#    network_path = '/home/narain/Desktop/ReferenceNetworks/Alpha' + alpha_value + '/FinalHaematocrit.vtp'

    # Print status update
    print(field_path)

    # Save file or read from existing file
    if read==0:    
        field_data, field_spacing = get_vti_data(field_path)  # import the data for the network
        middle_field = get_hex_domain(field_data, field_spacing)  # get the distribution in the middle of the field (replace with designated file)
        if plot==1:  # plot the O2 distribution
            O2_field, _ = get_field(middle_field)  # get the O2 mesh for plots and stats
            fig = plt.figure()
            ax = plt.axes()
            colour_map = plt.cm.get_cmap('jet')
            plt.suptitle('O$_2$ distribution generated by the ' + solver_name + ' solver in the hexagonal vessel network with thin vessel selection = ' + layout_selection + ' (σ = ' + sigma + ' and kills = ' + kills + '  vessels (1 unit = ' + str(field_spacing [0]) + ' μm)')
            ref_map = ax.imshow(O2_field, cmap=colour_map, origin='lower')
            fig.colorbar(ref_map, ax=ax, label='nM')
            plt.show()
        flat_field = middle_field['oxygen'].to_numpy().flatten()
#        np.save('/home/narain/Temporary Python Files/' + solver_name + 'Haematocrit_Selection' + layout_selection + '_Sigma' + sigma + '_Kills' + kills + '_distribution.npy', flat_field)
#    else:
#        flat_field = np.load('/home/narain/Temporary Python Files/' + solver_name + 'Haematocrit_Selection' + layout_selection + '_Sigma' + sigma + '_Kills' + kills + '_distribution.npy')

    # Get the basic stats 
#    middle_O2_stats = pd.DataFrame(flat_field).describe()
            
    # Get the number of total points
    number_of_points = 1
    for dim in np.shape(flat_field): number_of_points *= dim

    # Calculate the number of points below the given hypoxic thresholds
    hypoxic_fraction_list = []
    for hypoxic_threshold in hypoxic_threshold_list:
        hypoxic_points = (flat_field < hypoxic_threshold).sum()
        hypoxic_fraction = hypoxic_points/number_of_points  # calculate the hypoxic fraction
#        print(hypoxic_threshold, hypoxic_points)
        hypoxic_fraction_list.append(hypoxic_fraction)

    # Return the stats
    return hypoxic_fraction_list, np.mean(flat_field), np.amin(flat_field), np.percentile(flat_field, 50), np.amax(flat_field), np.std(flat_field)

# Define function to compute the average of all layouts in a kill selection
def compute_average_kill(solver_name, nhet, kills, max_layouts, hypoxic_threshold_list, plot, read):
    
    # Create table to store all the kill data trials in an alpha group
    kill_table = np.array([])

    # Extract metrics from all trials and store in beta table
    for layout_selection in range(1, max_layouts+1):    
#        print(layout_selection)
        hypoxic_fraction_list, mean_value, min_value, half_value, max_value, std_value = get_distribution_stats(solver_name, str(layout_selection), nhet, str(kills), hypoxic_threshold_list, plot, read)
        table_entry = np.hstack([hypoxic_fraction_list, mean_value, min_value, half_value, max_value, std_value])
        kill_table = np.vstack([kill_table, table_entry]) if kill_table.size else table_entry

    # Return the hypoxic fractions, mean, min, 50%, max, and std for a kill count averaged across all layouts
    return np.average(kill_table, axis=0)
#    return np.average(kill_table[0]), np.average(kill_table[1]), np.average(kill_table[2]), np.average(kill_table[3]), np.average(kill_table[4])

# Define a function to return statistics for all the heterogeneities in the data
def get_solver_stats(solver_name, alpha_list, kill_list, max_layouts, hypoxic_threshold_list, plot, read):
#    table = np.array([])
    alpha_table = np.array([])
    for nhet in alpha_list:
        for kill in kill_list:    
            average_kill_data = compute_average_kill(solver_name, nhet, kill, max_layouts, hypoxic_threshold_list, plot, read)
            table_entry = np.hstack([float(nhet), float(kill), average_kill_data])
            alpha_table = np.vstack([alpha_table, table_entry]) if alpha_table.size else table_entry
    return alpha_table

# =============================================================================
# DISTRIBUTION STATS & HYPOXIC FRACTIONS
# =============================================================================

# Enter details to allow looping over folders
solver_list = ['Constant', 'Pries', 'Memory', 'Fung']
alpha_list = ['1', '2', '3', '4']
max_kills = 77
max_layouts = 100
kill_list = [str(x) for x in range(0, max_kills + 1)]
#hypoxic_threshold_list = [2195, 10000, 15000, 20000, 25000, 27441] 
hypoxic_threshold_list = [2195, 5488, 10976, 16465, 21953, 27441] 
#hypoxic_threshold_list_pp = [0.8, 2, 4, 6, 8, 10] 

# Set solver name
solver_name = solver_list[0]

# Get the stats for all solvers (change to read=1 to extract from .vti files directly)
#solver_stats = get_solver_stats(solver_name, alpha_list, kill_list, max_layouts, hypoxic_threshold_list, plot=0, read=0)

# Save array
#np.save('/home/narain/Desktop/Results/Log Normal Distribution/Kills Pruning in Hexagonal Network with 100 Selections/TestHexagonalNetwork/' + solver_name + 'Haematocrit/python_solver_data.npy', solver_stats)
solver_stats = np.load('/home/narain/Desktop/Results/Hexagonal/Log Normal Distribution/Kills Pruning in Hexagonal Network with 100 Selections/TestHexagonalNetwork/' + solver_name + 'Haematocrit/python_solver_data.npy')

# Filter by alpha
mean_composite = np.array([])
min_composite = np.array([])
half_composite = np.array([])
max_composite = np.array([])
sd_composite = np.array([])
hypoxic_fraction_composite = np.array([])
for alpha_value in alpha_list:
    alpha_array = solver_stats[(solver_stats[:,0]==float(alpha_value))]
    hypoxic_fraction_data = alpha_array[:,2:-5]  # extract data between identifiers and basic stats (i.e., hypoxic fractions)
    mean_data = alpha_array[:,-5]
    min_data = alpha_array[:,-4]
    half_data = alpha_array[:,-3]
    max_data = alpha_array[:,-2]
    sd_data = alpha_array[:,-1]
    hypoxic_fraction_composite = np.vstack([hypoxic_fraction_composite, hypoxic_fraction_data]) if hypoxic_fraction_composite.size else hypoxic_fraction_data
    mean_composite = np.vstack([mean_composite, mean_data]) if mean_composite.size else mean_data
    min_composite = np.vstack([min_composite, min_data]) if min_composite.size else min_data
    half_composite = np.vstack([half_composite, half_data]) if half_composite.size else half_data
    max_composite = np.vstack([max_composite, max_data]) if max_composite.size else max_data
    sd_composite = np.vstack([sd_composite, sd_data]) if sd_composite.size else sd_data

# =============================================================================
# PERFUSION QUOTIENTS
# =============================================================================

# Define a function to generate the data for a single alpha value
def get_alpha_line(alpha_group, max_kills):
    
    # Compute the averages for the betas in the alpha group
    pq_table = np.array([])
    alpha_kills_grouped = alpha_group.groupby(alpha_group.kills)
    for kill_count in range(1, max_kills+1):
        pq_table_entry = np.array([])
        alpha_kills_group = alpha_kills_grouped.get_group(kill_count)
        pq_table_entry = np.array([alpha_kills_group["alpha"].mean(), alpha_kills_group["kills"].mean(), alpha_kills_group["PQ"].mean()])
        pq_table = np.vstack([pq_table, pq_table_entry]) if pq_table.size else pq_table_entry
    
    # Plot the alpha data
#    ax.plot(pq_table[:,0], pq_table[:,1], pq_table[:,2])
    
    # Return the table for reference
    return pq_table

# Read PQ file
filename = '/home/narain/Desktop/Results/Hexagonal/Log Normal Distribution/Kills Pruning in Hexagonal Network with 100 Selections/TestHexagonalNetwork/hex_lognormal_kills_perfusion_quotients.txt'
pq_df = pd.read_csv(filename, delim_whitespace=True, names=["network_name", "solver_name", "alpha", "selection", "kills", "PQ"])
#pq_df = pd.read_csv(filename, delim_whitespace=True, names=["alpha", "threshold", "PQ"], skiprows=1)

# Filter PQ data for multiple solvers
solver_filter = solver_name + 'Haematocrit'
pq_df = pq_df.loc[(pq_df["solver_name"] == solver_filter)]

# Drop extra data
#max_beta = 35
#pq_df = pq_df.loc[(pq_df["beta"] <= max_beta)]

# Separate by alpha 
alpha_grouped = pq_df.groupby(pq_df.alpha)
#alpha_0 = alpha_grouped.get_group(0)
alpha_1 = alpha_grouped.get_group(1)
alpha_2 = alpha_grouped.get_group(2)
alpha_3 = alpha_grouped.get_group(3)
alpha_4 = alpha_grouped.get_group(4)

# Compute average of all selections for PQ
line_1 = get_alpha_line(alpha_1, max_kills)
line_2 = get_alpha_line(alpha_2, max_kills)
line_3 = get_alpha_line(alpha_3, max_kills)
line_4 = get_alpha_line(alpha_4, max_kills)
# Combine the PQs
pq_composite = np.hstack([line_1, line_2, line_3 ,line_4])

# =============================================================================
# PLOTS
# =============================================================================

# Set the figure layout
fig, axs = plt.subplots(3, len(alpha_list), figsize=(20, 12), tight_layout = {'pad': 2})
fig.subplots_adjust(hspace = .5, wspace=.25)
linestyles = ['solid', 'dashed', 'dotted', 'dashdot', (0,(5,10)), (0, (3, 1, 1, 1, 1, 1))]
linecolours = ['#1f77b4','b','b','b','b', 'r']
linelegends = ['anoxia', 2, 4, 6, 8, 'hypoxia']
#plt.suptitle(solver_name + ' haematocrit solver in the heterogeneous hexagonal vessel network with individual kills pruning')

# Plot the distribution stats for a solver
axs = axs.ravel()
for i in range(len(alpha_list)):
    axs[i].set_ylim([0,30000])
    axs[i].plot(alpha_array[:,1], mean_composite[i], ls='dashed', label='mean')
    axs[i].plot(alpha_array[:,1], min_composite[i], ls='dotted', label='min')
#    axs[i].plot(alpha_array[:,1], half_composite[i], ls=(0, (3, 5, 1, 5)), label='50%')
    axs[i].plot(alpha_array[:,1], max_composite[i], ls='dashdot', label='max')
    axs[i].plot(alpha_array[:,1], sd_composite[i], ls='solid', label='SD')
    axs[i].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#    axs[i].set_xlabel('vessels killed')    
    if i==0:
        axs[i].set_ylabel('oxygen (nM)') 
        axs[i].legend(loc="best", prop={'size': 15})        
    axs[i].set_xlim(0)
    axs[i].grid()
    axs[i].title.set_text('${σ}$ = ' + alpha_list[i])

# Plot the PQ for a solver
for i in range(len(alpha_list),len(alpha_list)*2):
    axs[i].set_ylim([0,1.1])  # set PQ limits
    axs[i].plot(line_1[:,1], pq_composite[:, (3*i)-10], label='PQ')
#    axs[i].set_xlabel('vessels killed')    
    if i==len(alpha_list):
        axs[i].set_ylabel('perfusion quotient') 
#    axs[i].legend()
    axs[i].set_xlim(0)
    axs[i].grid()

# Plot the HF for a solver
for i in range(len(alpha_list*2),len(alpha_list)*3):
    offset = i-len(alpha_list*2)
    for threshold_index in range(len(hypoxic_threshold_list)):  # plot a line for each threshold
        if threshold_index==0 or threshold_index==5:
            axs[i].plot(alpha_array[:,1], hypoxic_fraction_composite[offset*len(kill_list):(offset+1)*len(kill_list), threshold_index], label=linelegends[threshold_index], ls=linestyles[threshold_index], c=linecolours[threshold_index])
        if i==len(alpha_list*2):
            axs[i].set_ylabel('hypoxic fraction') 
            axs[i].legend(loc="best", prop={'size': 15})
        axs[i].set_xlabel('vessels killed')    
        axs[i].set_xlim(0)
        axs[i].set_ylim([0,1.1])  # set HF limits
        axs[i].grid(b=True)

# Show plots
plt.show()

# Save image
file_path = Path('~/Desktop/Final Figures/' + solver_name + '_lognormal_hexagonal_kills_pruning.svg').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')
file_path = Path('~/Desktop/Final Figures/' + solver_name + '_lognormal_hexagonal_kills_pruning.png').expanduser()
fig.savefig(file_path, dpi=500, bbox_inches = 'tight')

# Prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

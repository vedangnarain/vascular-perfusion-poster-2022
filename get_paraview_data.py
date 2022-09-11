#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Mar 31 16:44:36 2022

@author: narain

Tested in Python 3.7.4.

Read VTK and VTI files to process them in Python.

"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
from collections import deque
#from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import pandas as pd
import scipy.interpolate
import vtk

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to read a .vtk file and return the point coordinates and data and the cell data
def get_vtk_data(vtk_path):
    print(vtk_path)

    # Set up the file reader
#    if 'vtp' in filename:
#        reader = vtk.vtkXMLPolyDataReader()
#    else:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtk_path)
    reader.Update()
    polydata = reader.GetOutput()
#    polydata.ReleaseDataFlagOn()
    
    # Extract the point coordinates
    point_coordinates_data = polydata.GetPoints().GetData()
    point_coordinates_data_array = vtk_to_numpy(point_coordinates_data)
    
    # Extract the point data 
    point_data_array = {}
    point_data = polydata.GetPointData()
    for column_index in range(point_data.GetNumberOfArrays()):
       column_name =  point_data.GetArrayName(column_index)
       point_data_array[column_name] = vtk_to_numpy(point_data.GetArray(column_index))
    
    # Extract the cell data    
    cell_data_array = {}
    cell_data = polydata.GetCellData()
    for column_index in range(cell_data.GetNumberOfArrays()):
       column_name =  cell_data.GetArrayName(column_index)
       cell_data_array[column_name] = vtk_to_numpy(cell_data.GetArray(column_index))
    
    # Return a dictionary with the point coordinates and data and cell data
    return point_coordinates_data_array, point_data_array, cell_data_array, polydata
           
# Define a function to read a .vti file and return the data
def get_vti_data(vti_path):
    
    # Import the file
    extension = vti_path.split('.').pop()
    reader = None
    if extension == 'vtk':
        reader = vtk.vtkDataSetReader() 
    elif extension == 'vti':
        reader = vtk.vtkXMLImageDataReader() 
    else:
        raise RuntimeError('Unknown File Type: %s ' % vti_path)
    reader.SetFileName( "%s" % (vti_path) ) 
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
def get_plottable_field(vti_data):
    
    # Change the data type
    vti_data = vti_data.astype('float32')
    
    # Calculate concentration statistics
    O2_stats = vti_data['oxygen'].describe()
    
    # Downsample data if needed
    vti_data = vti_data[::int(1)]
    
    # Convert dataframe into NumPy matrix
    mat = vti_data.to_numpy()
    
    # Get the x and y axes
    x = np.unique(mat[:,0])
    y = np.unique(mat[:,1])
    
    # Create a mesh from the x and y axes 
    X,Y = np.meshgrid(x, y)
    
    # Interpolate the concentration values over the mesh
    Z = scipy.interpolate.griddata((mat[:,0], mat[:,1]), mat[:,2], (X,Y), method='nearest')

    # Return the oxygen stats
    return Z, O2_stats

# Define a function to extract the region of evaluation of the forking network
def get_forking_domain(field_dataset, middle_generation_number, generation_coordinates):
    x_start = generation_coordinates[middle_generation_number-1]
    x_end = generation_coordinates[-middle_generation_number]
    middle_field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end))]
    return middle_field_dataset

# Define a function to extract the the region of evaluation of the hexagonal network
def get_hex_domain(field_dataset, field_spacing):
    x_start = 10*field_spacing[0]
    x_end = 195*field_spacing[0]
    y_start = 17*field_spacing[1]
    y_end = 173*field_spacing[1]
    field_dataset = field_dataset[(field_dataset['x'] >= int(x_start)) & (field_dataset['x'] <= int(x_end)) & (field_dataset['y'] >= int(y_start)) & (field_dataset['y'] <= int(y_end))]
    return field_dataset
#    return vessel_network, oxygen_distribution, middle_x

# Define a function to extract the predictive metrics from the forking network
def get_forking_predictors(vtk_path, reference_rank_lengths, reference_node_coordinates):
        
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    
    # Convert the radii and lengths in metres into the right units
    segment_diameters_um = cell_data_array['Vessel Radius m']*2*1000000  # convert radii (m) to diameters (um)
    reference_rank_lengths_um = reference_rank_lengths*1000000  # convert length (m to um)
    
    # Get the architectural metrics
    n_vessels = len(cell_data_array['Vessel Radius m'])
    mean_diameter = sum(segment_diameters_um)/n_vessels
    segment_lengths = [reference_rank_lengths_um[int(vessel_rank)] for vessel_rank in cell_data_array['Vessel Owner Rank']]  # convert length (m to um)
    mean_geometric_resistance = sum(segment_lengths/(segment_diameters_um**4))    
    
    # Get a list of vessel segments with node IDs for the adjacency matrices
    cellIds = vtk.vtkIdList()  # cell ids store to
    numberOfCells = polydata.GetNumberOfCells()
    segment_nodes = np.array([])
    for cellIndex in range(numberOfCells):  # for every cell
    #    print('new cell')
        polydata.GetCellPoints(cellIndex, cellIds)  # get IDs of nodes of the given cell
        cell_nodes = np.array([])
        for i in range(0, cellIds.GetNumberOfIds()):  # for every node of the given cell
            coord = polydata.GetPoint(cellIds.GetId(i))  # get coordinates of the node, type: class 'tuple'
            x = np.around(coord[0], 2)  # get x-coordinate of the node, type: class 'float'
    #        print(x)
            y = np.around(coord[1], 2)  # get y-coordinate of the node, type: class 'float'
    #        print(y)
            node_id = np.where((reference_node_coordinates[:,0] == x) & (reference_node_coordinates[:,1] == y))[0]
    #        print(node_id)
            cell_nodes = np.hstack([cell_nodes, node_id])
        segment_nodes = np.vstack([segment_nodes, cell_nodes]) if segment_nodes.size else cell_nodes
    segment_nodes = segment_nodes.astype('int')
    
    # Get the number of nodes
    number_of_nodes = len(reference_node_coordinates)
    
    # Initialise two nxn matrices of zeros, where n is the number of nodes
    diameter_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    length_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))

    # Fill in the adjacency matrices
    for segment_index, segment in enumerate(segment_nodes):
        
        # Get the segment nodes
        row, col = segment
        
        # Fill in the adjacency matrices
        diameter_adjacency_matrix[row,col] = segment_diameters_um[segment_index]
        diameter_adjacency_matrix[col,row] = segment_diameters_um[segment_index]
        length_adjacency_matrix[row,col] = reference_rank_lengths_um[int(cell_data_array['Vessel Owner Rank'][segment_index])]
        length_adjacency_matrix[col,row] = reference_rank_lengths_um[int(cell_data_array['Vessel Owner Rank'][segment_index])]

    # Return the architectural features and the adjacency matrices
    return n_vessels, mean_diameter, mean_geometric_resistance, diameter_adjacency_matrix, length_adjacency_matrix

# Define a function to extract the predictive metrics from the hexagonal network
def get_hex_predictors(vtk_path, vessel_length_m, reference_node_coordinates):
    
    # Get the .vtk data
    point_coordinates_data_array, point_data_array, cell_data_array, polydata = get_vtk_data(vtk_path)
    
    # Convert the radii and lengths in metres into the right units
    segment_diameters_um = cell_data_array['Vessel Radius m']*2*1000000  # convert radii (m) to diameters (um)
    segment_length_um = vessel_length_m*1000000  # convert length (m to um)
       
    # Get the architectural metrics    
    n_vessels = len(cell_data_array['Vessel Radius m'])
    mean_diameter = sum(segment_diameters_um)/n_vessels
    mean_geometric_resistance = sum(segment_length_um/(segment_diameters_um**4))
    
    # Get a list of vessel segments with node IDs for the adjacency matrices
    cellIds = vtk.vtkIdList()  # cell ids store to
    numberOfCells = polydata.GetNumberOfCells()
    segment_nodes = np.array([])
    for cellIndex in range(numberOfCells):  # for every cell
    #    print('new cell')
        polydata.GetCellPoints(cellIndex, cellIds)  # get IDs of nodes of the given cell
        cell_nodes = np.array([])
        for i in range(0, cellIds.GetNumberOfIds()):  # for every node of the given cell
            coord = polydata.GetPoint(cellIds.GetId(i))  # get coordinates of the node, type: class 'tuple'
            x = np.around(coord[0], 2)  # get x-coordinate of the node, type: class 'float'
    #        print(x)
            y = np.around(coord[1], 2)  # get y-coordinate of the node, type: class 'float'
    #        print(y)
            node_id = np.where((reference_node_coordinates[:,0] == x) & (reference_node_coordinates[:,1] == y))[0]
    #        print(node_id)
            cell_nodes = np.hstack([cell_nodes, node_id])
        segment_nodes = np.vstack([segment_nodes, cell_nodes]) if segment_nodes.size else cell_nodes
    segment_nodes = segment_nodes.astype('int')
    
    # Get the number of nodes
    number_of_nodes = len(reference_node_coordinates)
    
    # Initialise two nxn matrices of zeros, where n is the number of nodes
    diameter_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    length_adjacency_matrix = np.zeros((number_of_nodes, number_of_nodes))
    
    # Fill in the adjacency matrices
    for segment_index, segment in enumerate(segment_nodes):
        
        # Get the segment nodes
        row, col = segment
        
        # Fill in the adjacency matrices
        diameter_adjacency_matrix[row,col] = segment_diameters_um[segment_index]
        diameter_adjacency_matrix[col,row] = segment_diameters_um[segment_index]
        length_adjacency_matrix[row,col] = segment_length_um
        length_adjacency_matrix[col,row] = segment_length_um
      
    return n_vessels, mean_diameter, mean_geometric_resistance, diameter_adjacency_matrix, length_adjacency_matrix
    
# Import the generation coordinates, rank lengths, and node coordinates for the forking network 
def get_reference_forking_network(lambda_value=4):
    alpha_value = '1.00'  # only the diameter varies with alpha 
    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/Alpha' + alpha_value + '/FinalHaematocrit.vtp'
    reference_point_coordinates_data_array, point_data_array, reference_cell_data_array, polydata = get_vtk_data(reference_network_path)
    generation_coordinates = np.unique(reference_point_coordinates_data_array[:,0])
    rank_diameters = -np.sort(-np.unique(reference_cell_data_array['Vessel Radius m'])*2)
    reference_rank_lengths = rank_diameters*4
    
    # Round the coordinates to two decimal places
    reference_point_coordinates_data_array = reference_point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    reference_node_coordinates = np.around(reference_point_coordinates_data_array[:,0:2], 2)
 
    return generation_coordinates, reference_rank_lengths, reference_node_coordinates

# Import the node coordinates for the hexagonal network 
def get_reference_hexagonal_network(vessel_length_m=100*(10**-6)):
    sigma = 1
    selection = 1
    radius_threshold = 0
    reference_network_path = '/home/narain/Desktop/Scripts/reference_networks/Sigma' + str(sigma) + '/Selection' + str(selection) +'/RadiusThreshold' + str(radius_threshold) + '/FinalHaematocrit.vtp'
    reference_point_coordinates_data_array, _, _, _ = get_vtk_data(reference_network_path)
        
    # Round the coordinates to two decimal places
    reference_point_coordinates_data_array = reference_point_coordinates_data_array.astype('float64')  # convert to the right data type for later comparison
    reference_node_coordinates = np.around(reference_point_coordinates_data_array[:,0:2], 2)
 
    return reference_node_coordinates

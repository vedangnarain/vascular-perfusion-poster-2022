#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Apr 25 19:17:32 2021

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Modify Hexagonal Network

This script modifies the list of edges that Jakub generated. It also allows us
to generate a subset of edges in a certain spatial range.

Tested in Python 3.7.4.

"""

# =============================================================================
# LIBRARIES & INITIALISATION
# =============================================================================

# Initialise libraries
import pandas as pd
import time

# Starts stopwatch to clock execution time
start_time = time.time()

# =============================================================================
# FUNCTIONS
# =============================================================================

# Define a function to filter vertices based on limits (look at points in Paraview for easy viewing)
def filter_vertices(all_vertices, x_min, x_max, y_min, y_max):
    
    # Filter vertices
    filtered_vertices = edges.loc[(all_vertices['start_x'] >= x_min) & 
                            (all_vertices['start_x'] <= x_max) & 
                            (all_vertices['end_x'] >= x_min) & 
                            (all_vertices['end_x'] <= x_max) & 
                            (all_vertices['start_y'] >= y_min) & 
                            (all_vertices['start_y'] <= y_max) & 
                            (all_vertices['end_y'] >= y_min) & 
                            (all_vertices['end_y'] <= y_max)]
    
    # Set edge numbers
    filtered_vertices['edge_number'] = range(1, 1+len(filtered_vertices))
    
    # Return filtered array
    return filtered_vertices

# =============================================================================
# DATA
# =============================================================================

# Import file
filename = '/Users/vedang/Simulations/scripts/Modify Hexagonal Points/EdgesMatrix.txt'
edges = pd.read_csv(filename, sep='\t', names=['start_x', 'start_y', 'end_x', 'end_y', 'edge_number','NaN'])
edges = edges.drop(columns=['NaN'])

# Offset x-coordinate by 50
edges['start_x'] = edges['start_x'] + 50
edges['end_x'] = edges['end_x'] + 50

# =============================================================================
# FULL NETWORK
# =============================================================================

# Save file with vertices for the full network
edges.to_csv("/Users/vedang/Simulations/scripts/Modify Hexagonal Points/EdgesMatrix_Offset.txt", sep="\t", header=0, index=0)

# =============================================================================
# SINGLE FEED NETWORK
# =============================================================================

# Save file with vertices for single feed network
single_path = filter_vertices(all_vertices = edges, x_min = 0, x_max = 400, y_min = 173, y_max = 347)
single_path.to_csv("/Users/vedang/Simulations/scripts/Modify Hexagonal Points/EdgesMatrix_SingleFeed.txt", sep="\t", header=0, index=0)

# =============================================================================
# MULTI FEED NETWORK
# =============================================================================

# Save file with vertices for multi-feed network
multi_path = filter_vertices(all_vertices = edges, x_min = 0, x_max = 400, y_min = 0, y_max = 436)
multi_path.to_csv("/Users/vedang/Simulations/scripts/Modify Hexagonal Points/EdgesMatrix_MultiFeed.txt", sep="\t", header=0, index=0)

# =============================================================================
# NEIGHBOURHOOD NETWORK
# =============================================================================

# Save file with vertices for network with neighbours
neighbourhood = filter_vertices(all_vertices = edges, x_min = 0, x_max = 700, y_min = 0, y_max = 607)
neighbourhood.to_csv("/Users/vedang/Simulations/scripts/Modify Hexagonal Points/EdgesMatrix_Neighbourhood.txt", sep="\t", header=0, index=0)

# =============================================================================
# VAMPIRE NETWORK
# =============================================================================

# Save file with vertices for network with neighbours
neighbourhood = filter_vertices(all_vertices = edges, x_min = 0, x_max = 700, y_min = 0, y_max = 607)
neighbourhood.to_csv("/Users/vedang/Simulations/scripts/Modify Hexagonal Points/EdgesMatrix_Vampire.txt", sep="\t", header=0, index=0)

# prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

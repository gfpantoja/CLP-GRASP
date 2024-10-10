# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 08:40:24 2024

@author: dellphoto
"""

import plotly.graph_objects as go
import random

# Define the file paths
instance_file_path = 'Instances/MLBR_1_1_1.txt'
solution_file_path = 'Solutions/R_MLBR_1_1_1_MD1_P0_S1_EM2_st_rep1_nt1_d0.txt'

# Function to read instance file and extract box dimensions
def read_instance_file(instance_file_path):
    Lx, Ly, Lz = 0, 0, 0
    with open(instance_file_path, 'r') as f:
        lines = f.readlines()
        line = lines[1]
        data = line.split()
        Lx, Ly, Lz = int(data[0]), int(data[1]), int(data[2])
    return Lx, Ly, Lz

# Function to read the solution file and extract packed box positions
def read_solution_file(solution_file_path):
    packed_boxes = []
    with open(solution_file_path, 'r') as f:
        lines = f.readlines()[1:]  # Skip the first line (General Information)
        for line in lines:
            data = list(map(int, line.split()))
            if len(data) == 9:
                packed_boxes.append(data)
    if len(packed_boxes) == 0:
        print("Solution file was not correctly read")
    else:
        print("Solution file was correctly read")
    return packed_boxes

# Function to plot the packed boxes in a 3D space
def plot_boxes(Lx, Ly, Lz, packed_boxes):
    fig = go.Figure()

    for box in packed_boxes:
        x1, y1, z1, x2, y2, z2 = box[3], box[4], box[5], box[6], box[7], box[8]

        # Generate random color for each box
        color = f'rgb({random.randint(0, 255)}, {random.randint(0, 255)}, {random.randint(0, 255)})'

        # Define the vertices of the box
        vertices = [
            [x1, y1, z1], [x2, y1, z1], [x2, y2, z1], [x1, y2, z1],  # Bottom face
            [x1, y1, z2], [x2, y1, z2], [x2, y2, z2], [x1, y2, z2]   # Top face
        ]
        x, y, z = zip(*vertices)

        # Define the faces of the box (two triangles per face)
        i_faces = [0, 0, 4, 4, 1, 1, 0, 0, 3, 3, 3, 3]
        j_faces = [1, 2, 5, 6, 5, 6, 4, 7, 2, 6, 1, 5]
        k_faces = [2, 3, 6, 7, 6, 2, 7, 3, 6, 7, 5, 4]

        # Add the box to the 3D plot
        fig.add_trace(go.Mesh3d(
            x=x,
            y=y,
            z=z,
            i=i_faces,
            j=j_faces,
            k=k_faces,
            color=color,
            flatshading=True
        ))

    # Configure the layout of the plot
    fig.update_layout(
        scene=dict(
            xaxis_title='X Axis',
            yaxis_title='Y Axis',
            zaxis_title='Z Axis',
            aspectmode='data'  # Automatically scale the axes based on the data
        ),
        title='3D Visualization of Packed Boxes',
        margin=dict(l=0, r=0, b=0, t=30),
    )
    fig.write_html("Image.html")
    fig.show()

# Read instance and solution files
Lx, Ly, Lz = read_instance_file(instance_file_path)
packed_boxes = read_solution_file(solution_file_path)

# Plot the boxes
plot_boxes(Lx, Ly, Lz, packed_boxes)

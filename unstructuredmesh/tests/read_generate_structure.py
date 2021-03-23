#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:08:27 2021

@author: kissami
"""

import os
from unstructuredmesh import ddm

# ... get the mesh directory
try:
    MESH_DIR = os.environ['MESH_DIR']

except:
    BASE_DIR = os.path.dirname(os.path.realpath(__file__))
    BASE_DIR = os.path.join(BASE_DIR ,'..', '..')
    MESH_DIR = os.path.join(BASE_DIR, 'mesh')

dim = 2
#File name
if dim == 2:
    filename = os.path.join(MESH_DIR, "rectangle.msh")
elif dim == 3:
    filename = os.path.join(MESH_DIR, "cube.msh")

#Read gmsh file
ddm.readmesh(filename, dim=dim)

#Create the informations about cells, faces and nodes
grid = ddm.generate_structure(dim=dim)

faces = grid["faces"]
cells = grid["cells"]
nodes = grid["nodes"]

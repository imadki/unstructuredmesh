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

#reading gmsh file and partitioning into size subdomains
filename = os.path.join(MESH_DIR, "rectangle.msh")

#generating local grid for each subdomain
ddm.readmesh(filename)

grid = ddm.generate_structure()

faces = grid["faces"]
cells = grid["cells"]
nodes = grid["nodes"]
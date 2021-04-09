#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:08:27 2021

@author: kissami
"""

import os
from unstructuredmesh import ddm
from mpi4py import MPI
import numpy as np
import timeit

COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

# ... get the mesh directory
try:
    MESH_DIR = os.environ['MESH_DIR']

except:
    BASE_DIR = os.path.dirname(os.path.realpath(__file__))
    BASE_DIR = os.path.join(BASE_DIR ,'..', '..')
    MESH_DIR = os.path.join(BASE_DIR, 'mesh')

dim = 3
#File name
if dim == 2:
    filename = os.path.join(MESH_DIR, "rectangle.msh")
    typeofCells = 'triangle'
elif dim == 3:
    filename = os.path.join(MESH_DIR, "cube.msh")
    typeofCells = 'tetra'

start = timeit.default_timer()

if RANK == 0:
    #Read gmsh file
    ddm.readmesh(filename, dim=dim, size=SIZE)

COMM.Barrier()
#Create the informations about cells, faces and nodes
grid = ddm.generate_structure(dim=dim, size=SIZE)

faces = grid["faces"]
cells = grid["cells"]
nodes = grid["nodes"]
halos = grid["halos"]

stop = timeit.default_timer()

if RANK == 0:
    print("Global CPU time", stop - start , RANK)

nbfaces    = len(faces.nodeid)
nbelements = len(cells.nodeid)
nbnodes    = len(nodes.vertex)

elements = {typeofCells: cells.nodeid}
points = []
for i in nodes.vertex:
    points.append([i[0], i[1], i[2]])

element = np.array(elements)
points = np.array(points)

import meshio
meshio.write_points_cells("visu"+str(RANK)+".vtu",  points, elements, file_format="vtu")

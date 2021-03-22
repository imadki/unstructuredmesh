#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:08:27 2021

@author: kissami
"""


# coding: utf-8
import os
import timeit
import meshio
import numpy as np

__all__ = ['readmesh']

def readmesh(filename, dim):

    if dim == 2 :
        typeOfCells = "triangle"
        typeOfFaces = "line"
    else:
        typeOfCells = "tetra"
        typeOfFaces = "triangle"

    def load_gmsh_mesh(filename):
        #mesh = meshio.gmsh.read(filename)
        mesh = meshio.read(filename)
        return mesh

    def create_cell_nodeid(mesh):
        cell_nodeid = []
        
        if type(mesh.cells) == dict:
            cell_nodeid = mesh.cells[typeOfCells]

        elif type(mesh.cells) == list:
            cell_nodeid = mesh.cells[1].data

        for i in range(len(cell_nodeid)):
            cell_nodeid[i].sort()

        return cell_nodeid

    def define_ghost_node(mesh, nodes):
        ghost_nodes = [0]*len(nodes)

        if type(mesh.cells) == dict:
            for i, j in mesh.cell_data.items():
                if i == typeOfFaces:
                    ghost = j.get('gmsh:physical')

            for i, j in mesh.cells.items():
                if i == typeOfFaces:
                    for k in range(len(j)):
                        for index in range(dim):
                            if ghost[k] > 2:
                                ghost_nodes[j[k][index]] = int(ghost[k])
            for i, j in mesh.cells.items():
                if i == typeOfFaces:
                    for k in range(len(j)):
                        for index in range(dim):
                            if ghost[k] <= 2:
                                ghost_nodes[j[k][index]] = int(ghost[k])

        elif type(mesh.cells) == list:
            ghost = mesh.cell_data['gmsh:physical'][0]
            for i in range(len(mesh.cells[0].data)):
                for j in range(dim):
                    if ghost[i] > 2:
                        ghost_nodes[mesh.cells[0].data[i][j]] = int(ghost[i])

            for i in range(len(mesh.cells[0].data)):
                for j in range(dim):
                    if ghost[i] <= 2:
                        ghost_nodes[mesh.cells[0].data[i][j]] = int(ghost[i])

        return ghost_nodes

    def create_nodes(mesh):
        nodes = []
        nodes = mesh.points
        return nodes

    start = timeit.default_timer()

    #load mesh
    mesh = load_gmsh_mesh(filename)
    #coordinates x, y of each node
    nodes = create_nodes(mesh)
    #nodes of each cell
    cell_nodeid = create_cell_nodeid(mesh)

    ghost_nodes = define_ghost_node(mesh, nodes)
    
    nbelements = len(cell_nodeid)
    nbnodes = len(nodes)

    print("Number of Cells : ", nbelements)
    print("Number of Nodes : ", nbnodes)
    
    # THe mesh file
    filename = "mesh.txt"

    if os.path.exists(filename):
        os.remove(filename)

    with open(filename, "a") as text_file:
        text_file.write("elements\n")
        np.savetxt(text_file, cell_nodeid, fmt='%u')
        text_file.write("endelements\n")
        text_file.write("nodes\n")
        for i in range(len(nodes)):
            for j in range(3):
                text_file.write(str(nodes[i][j])+str(" "))
            text_file.write(str(ghost_nodes[i]))
            text_file.write("\n")
        text_file.write("endnodes\n")

    stop = timeit.default_timer()

    print('Global Execution Time: ', stop - start)
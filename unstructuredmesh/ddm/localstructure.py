#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:08:27 2021

@author: kissami
"""

#from numba.typed import Dict
import timeit
from numpy import zeros, ones, asarray,  double, int64, unique

from unstructuredmesh.ddm.utils import (Compute_2dcentervolumeOfCell, create_info_2dfaces, 
                                        create_2dfaces, create_2doppNodeOfFaces,
                                        create_node_cellid, create_cell_faceid,
                                        create_cellsOfFace, create_NeighborCellByFace, 
                                        create_NormalFacesOfCell, Compute_3dcentervolumeOfCell, 
                                        create_info_3dfaces, create_3doppNodeOfFaces, create_3dfaces)

__all__ = ["generate_structure"]

class Cells:
    nodeid = []
    faceid = []
    center = []
    volume = []
    cellfid = []
    cellnid = []
    nf = [] 

    def __init__(self, nodeid, faceid, center, volume, cellfid, cellnid, nf):
        self.nodeid = nodeid    # instance variable unique to each instance
        self.faceid = faceid
        self.center = center
        self.volume = volume
        self.cellfid = cellfid
        self.cellnid = cellnid
        self.nf = nf
        
class Nodes:
    vertex = []
    name = []
    cellid = []
    ghostcenter = []

    def __init__(self, vertex, name, cellid, ghostcenter):
        self.vertex = vertex    #instance variable unique to each instance
        self.name = name
        self.cellid = cellid
        self.ghostcenter = ghostcenter

class Faces:
    nodeid = []
    cellid = []
    name = []
    normal = []
    mesure = []
    bound = 0
    center = []
    ghostcenter = []
    oppnodeid = []

    def __init__(self, nodeid, cellid, name, normal, mesure, center, bound, ghostcenter, oppnodeid):
        self.nodeid = nodeid    # instance variable unique to each instance
        self.cellid = cellid
        self.name = name
        self.normal = normal
        self.mesure = mesure
        self.bound = bound
        self.center = center
        self.ghostcenter = ghostcenter
        self.oppnodeid = oppnodeid

def CreateStructure(file, dim):

    #Lecture des cellules à partir du fichier mesh..txt
    for line in file:
        #read elements
        if line == "elements\n":
            continue
        if line == "endelements\n":
            continue
        if line == "nodes\n":
            break
        Cells.nodeid.append([int64(x) for x in line.split()])
    #Lecture des coordonnées des noeuds à partir du fichier mesh..txt
    for line in file:
        #read Nodes
        if line == "nodes\n":
            continue
        if line == "endnodes\n":
            continue
        if line == "halosint64\n":
            break
        Nodes.vertex.append([double(x) for x in line.split()])
        
    start = timeit.default_timer()

    Cells.nodeid = asarray(Cells.nodeid, dtype=int64)
    Nodes.vertex = asarray(Nodes.vertex, dtype=double)

    nbelements = len(Cells.nodeid)
    nbnodes    = len(Nodes.vertex)
    
    Nodes.name = zeros(nbnodes, dtype=int64)
    for i in range(nbnodes):
        Nodes.name[i] = int64(Nodes.vertex[i][3])
    
    #Create center and volume for each cell (pycceliser)
    Cells.center = zeros((nbelements, 3), dtype=double)
    Cells.volume = zeros(nbelements, dtype=double)
    if dim == 2:
        Compute_2dcentervolumeOfCell(Cells.nodeid, Nodes.vertex, nbelements, Cells.center, Cells.volume)
    elif dim == 3:
        Compute_3dcentervolumeOfCell(Cells.nodeid, Nodes.vertex, nbelements, Cells.center, Cells.volume)
    #create cells over each node (still numba function)
    Nodes.cellid, Cells.cellnid =  create_node_cellid(Cells.nodeid, Nodes.vertex, nbelements, nbnodes, dim=dim)
    
    #creating faces (pycceliser)
    faces = zeros(((dim+1)*nbelements, dim), dtype=int64)
    cellf = zeros((nbelements, dim+1), dtype=int64)
    if dim == 2:
        create_2dfaces(Cells.nodeid, nbelements, faces, cellf)
    elif dim == 3:
        create_3dfaces(Cells.nodeid, nbelements, faces, cellf)
    
    faces = unique(faces, axis=0, return_inverse="True")
    Faces.nodeid  = faces[0]
    oldTonewIndex = faces[1]
    nbfaces = len(Faces.nodeid)
    
    Cells.faceid = zeros((nbelements, (dim+1)), dtype=int64)
    create_cell_faceid(nbelements, oldTonewIndex, cellf, Cells.faceid, dim=dim)
    
    ############################################################################
    #creater cells left and right of each face (pycceliser)
    Faces.cellid = -1*ones((nbfaces, 2), dtype=int64)
    create_cellsOfFace(Cells.faceid, nbelements, nbfaces, Faces.cellid, dim=dim)
    ############################################################################
    Cells.cellfid = create_NeighborCellByFace(Cells.faceid, Faces.cellid, nbelements, dim=dim)
    
    ############################################################################
    #create info of faces (pycceliser)
    Faces.name   = zeros(nbfaces, dtype=int64)
    Faces.normal = zeros((nbfaces, 3), dtype=double)
    Faces.mesure = zeros(nbfaces, dtype=double)
    Faces.center = zeros((nbfaces, 3), dtype=double)
    if dim == 2:
        create_info_2dfaces(Faces.cellid, Faces.nodeid, Nodes.name, Nodes.vertex, Cells.center, 
                            nbfaces, Faces.normal, Faces.mesure, Faces.center, Faces.name)
    elif dim == 3:
        create_info_3dfaces(Faces.cellid, Faces.nodeid, Nodes.name, Nodes.vertex, Cells.center, 
                            nbfaces, Faces.normal, Faces.mesure, Faces.center, Faces.name)
    #Number of boundary faces
    Faces.bound = len(Faces.name[Faces.name !=0])
    ############################################################################
    #Create outgoing normal vectors (pycceliser)
    Cells.nf = zeros((nbelements, dim+1, 3), dtype=double)
    create_NormalFacesOfCell(Cells.center, Faces.center, Cells.faceid, Faces.normal, nbelements, Cells.nf, dim=dim)
    ###########################################################################   
    #still numba function
    if dim == 2:
        Faces.oppnodeid = create_2doppNodeOfFaces(Cells.nodeid, Cells.faceid, Faces.nodeid, nbelements, nbfaces)
    elif dim == 3:
        Faces.oppnodeid = create_3doppNodeOfFaces(Cells.nodeid, Cells.faceid, Faces.nodeid, nbelements, nbfaces)
    
    stop = timeit.default_timer()
    print("Create "+str(dim)+"d structure ", stop-start)

    cells = Cells(Cells.nodeid,  Cells.faceid, Cells.center,
                  Cells.volume,  Cells.cellfid, Cells.cellnid, Cells.nf)
    
    nodes = Nodes(Nodes.vertex, Nodes.name, Nodes.cellid, Nodes.ghostcenter)

    faces = Faces(Faces.nodeid, Faces.cellid, Faces.name,
                  Faces.normal, Faces.mesure, Faces.center, 
                  Faces.bound, Faces.ghostcenter, Faces.oppnodeid) 

    return cells, nodes, faces


def generate_structure(dim):
    
    filename = 'mesh.txt'
    txt_file = open(filename)

    cells, nodes, faces = CreateStructure(txt_file, dim=dim)

    grid = {}

    grid["cells"] = cells
    grid["nodes"] = nodes
    grid["faces"] = faces

    txt_file.close()

    return grid

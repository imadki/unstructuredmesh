#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:08:27 2021

@author: kissami
"""

from collections import OrderedDict
import numpy as np

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
        self.vertex = vertex    # instance variable unique to each instance
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

def normalvector2d(node_a, node_b, bary):

    norm = np.zeros(2)#[None]*2
    snorm = np.zeros(2)#[None]*2
    center = np.zeros(2)
    normal = np.zeros(2)#[None]*2

    norm[0] = node_a[1] - node_b[1]
    norm[1] = node_b[0] - node_a[0]

    center[0] = 0.5 * (node_a[0] + node_b[0])
    center[1] = 0.5 * (node_a[1] + node_b[1])

    snorm[0] = bary[0] - center[0]
    snorm[1] = bary[1] - center[1]

    if (snorm[0] * norm[0] + snorm[1] * norm[1]) > 0:
        normal[0] = -1*norm[0]
        normal[1] = -1*norm[1]
    else:
        normal[0] = norm[0]
        normal[1] = norm[1]

    normal[0] = normal[0]
    normal[1] = normal[1]
    
    mesure = np.sqrt(normal[0]**2 + normal[1]**2)

    return normal, mesure


def wedge_3d(u, v):
    wedge = np.zeros(3)
    
    wedge[0] = u[1]*v[2] - u[2]*v[1]
    wedge[1] = u[2]*v[0] - u[0]*v[2]
    wedge[2] = u[0]*v[1] - u[1]*v[0]
    
    return wedge

def dot_vec3(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def det_vec_3d(u, v, w):
    return dot_vec3(u, wedge_3d(v, w))

def normalvector3d(node_a, node_b, node_c, bary):
    
    a = node_a
    b = node_b
    c = node_c
    norm = np.zeros(3)
   
    norm = 0.5*wedge_3d(b-a, c-a)

    snorm = np.zeros(3)
    center = np.zeros(3)
    normal = np.zeros(3)

    center[0] = 1./3 * (a[0] + b[0] + c[0])
    center[1] = 1./3 * (a[1] + b[1] + c[1])
    center[2] = 1./3 * (a[2] + b[2] + c[2])

    snorm[0] = bary[0] - center[0]
    snorm[1] = bary[1] - center[1]
    snorm[2] = bary[2] - center[2]

    if (snorm[0] * norm[0] + snorm[1] * norm[1] + snorm[2] * norm[2]) > 0:
        normal[0] = -1*norm[0]
        normal[1] = -1*norm[1]
        normal[2] = -1*norm[2]
    else:
        normal[0] = norm[0]
        normal[1] = norm[1]
        normal[2] = norm[2]

    mesure = np.sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)

    return normal, mesure

def Create2DStructure(file):

    #Lecture des cellules à partir du fichier mesh..txt
    for line in file:
        #read elements
        if line == "elements\n":
            continue
        if line == "endelements\n":
            continue
        if line == "nodes\n":
            break
        Cells.nodeid.append([int(x) for x in line.split()])
    #Lecture des coordonnées des noeuds à partir du fichier mesh..txt
    for line in file:
        #read Nodes
        if line == "nodes\n":
            continue
        if line == "endnodes\n":
            continue
        if line == "halosint\n":
            break
        Nodes.vertex.append([np.double(x) for x in line.split()])
        Nodes.name = [0]*len(Nodes.vertex)
        
        
    nbelements = len(Cells.nodeid)
    nbnodes    = len(Nodes.vertex)

    for i in range(nbnodes):
        Nodes.name[i] = int(Nodes.vertex[i][3])

    #calcul du barycentre
    for i in range(nbelements):
        s_1 = Cells.nodeid[i][0]
        s_2 = Cells.nodeid[i][1]
        s_3 = Cells.nodeid[i][2]

        x_1 = Nodes.vertex[s_1][0]
        y_1 = Nodes.vertex[s_1][1]
        x_2 = Nodes.vertex[s_2][0]
        y_2 = Nodes.vertex[s_2][1]
        x_3 = Nodes.vertex[s_3][0]
        y_3 = Nodes.vertex[s_3][1]

        Cells.center.append((1./3 * (x_1 + x_2 + x_3), 1./3*(y_1 + y_2 + y_3)))
        Cells.volume.append((1./2) * abs((x_1-x_2)*(y_1-y_3)-(x_1-x_3)*(y_1-y_2)))
        var1 = (x_2-x_1)*(y_3-y_1)-(y_2-y_1)*(x_3-x_1)
        if var1 < 0:
            Cells.nodeid[i][0] = s_1;   Cells.nodeid[i][1] = s_3; Cells.nodeid[i][2] = s_2

    tmp = [[] for i in range(nbnodes)]
    longn = [0]*nbnodes
    for i in range(nbelements):
        for j in range(3):
            tmp[Cells.nodeid[i][j]].append(i)
            longn[Cells.nodeid[i][j]] = longn[Cells.nodeid[i][j]] + 1
    
    longc = [0]*nbelements
    tmp2 = [[] for i in range(nbelements)]
    for i in range(nbelements):
        for j in range(3):
            for k in range(len(tmp[Cells.nodeid[i][j]])):
                if (tmp[Cells.nodeid[i][j]][k] not in tmp2[i] and  tmp[Cells.nodeid[i][j]][k] != i):
                    tmp2[i].append(tmp[Cells.nodeid[i][j]][k])
                    longc[i] = longc[i] + 1

    maxlongn = max(longn)
    Nodes.cellid = [[-1]*maxlongn for i in range(nbnodes)]

    for i in range(len(tmp)):
        for j in range(len(tmp[i])):
            Nodes.cellid[i][j] = tmp[i][j]
            
    for i in range(nbnodes):
        Nodes.cellid[i].append(longn[i])
    
    maxlongc = max(longc)
    Cells.cellnid = [[-1]*maxlongc for i in range(nbelements)]
    
    for i in range(len(tmp2)):
        for j in range(len(tmp2[i])):
            Cells.cellnid[i][j] = tmp2[i][j]

    #Création des faces
    cellule = Cells.nodeid
    cellf = []
    faces = []
    k = 0
    
    for i in range(len(cellule)):
        faces.append([cellule[i][0], cellule[i][1]])
        faces.append([cellule[i][1], cellule[i][2]])
        faces.append([cellule[i][2], cellule[i][0]])
        cellf.append([faces[k], faces[k+1], faces[k+2]])
        k = k+3
    for i in range(len(faces)):
        faces[i].sort()
    faces = set(tuple(x) for x in faces)
    faces = list(faces)

    nbfaces = len(faces)
    
    facesdict = OrderedDict()
    for i in range(nbfaces):
        facesdict[faces[i]] = i
        Faces.nodeid.append(faces[i])
        Faces.cellid.append((-1, -1))

    #Création des 3 faces de chaque cellule
    for i in range(len(cellf)):
        Cells.faceid.append([facesdict.get(tuple(cellf[i][0])), facesdict.get(tuple(cellf[i][1])),
                             facesdict.get(tuple(cellf[i][2]))])

    #Faces.cellid = [(-1, -1)]*len(faces)
    for i in range(nbelements):
        for j in range(3):
            if Faces.cellid[Cells.faceid[i][j]] == (-1, -1):
                Faces.cellid[Cells.faceid[i][j]] = (i, -1)

            if Faces.cellid[Cells.faceid[i][j]][0] != i:
                Faces.cellid[Cells.faceid[i][j]] = (Faces.cellid[Cells.faceid[i][j]][0], i)

    Cells.cellfid = [[] for i in range(nbelements)]
    #Création des 3 triangles voisins par face
    for i in range(nbelements):
        for j in range(3):
            f = Cells.faceid[i][j]
            if Faces.cellid[f][1] != -1:
                if i == Faces.cellid[f][0]:
                    Cells.cellfid[i].append(Faces.cellid[f][1])
                else:
                    Cells.cellfid[i].append(Faces.cellid[f][0])
            else:
                Cells.cellfid[i].append(-1)
    
    #Faces aux bords (1,2,3,4), Faces à l'interieur 0    A VOIR !!!!!
    for i in range(nbfaces):
        Faces.name.append(0)
        if (Faces.cellid[i][1] == -1 and Faces.cellid[i][1] != -10):
            if Nodes.name[faces[i][0]] == Nodes.name[faces[i][1]]:
                Faces.name[i] = Nodes.name[faces[i][0]]
            if ((Nodes.name[faces[i][0]] == 3 and Nodes.name[faces[i][1]] != 0) or
                    (Nodes.name[faces[i][0]] != 0 and Nodes.name[faces[i][1]] == 3)):
                Faces.name[i] = 3
            if ((Nodes.name[faces[i][0]] == 4 and Nodes.name[faces[i][1]] != 0) or
                    (Nodes.name[faces[i][0]] != 0 and Nodes.name[faces[i][1]] == 4)):
                Faces.name[i] = 4
        normal, mesure = normalvector2d(np.asarray(Nodes.vertex[faces[i][0]]),
                                        np.asarray(Nodes.vertex[faces[i][1]]),
                                        Cells.center[Faces.cellid[i][0]])
        Faces.normal.append(normal)
        Faces.mesure.append(mesure)
        Faces.center.append([0.5 * (Nodes.vertex[faces[i][0]][0] + Nodes.vertex[faces[i][1]][0]),
                             0.5 * (Nodes.vertex[faces[i][0]][1] + Nodes.vertex[faces[i][1]][1])])
    for i in range(nbfaces):
        if Faces.name[i] != 0:
            Faces.bound += 1
        
    #compute the outgoing normal faces for each cell
    for i in range(nbelements):
        ss = np.zeros((3, 2))

        G = np.asarray(Cells.center[i])

        for j in range(3):
            f = Cells.faceid[i][j]
            c = Faces.center[f]

            if np.dot(G-c, Faces.normal[f]) < 0.0:
                ss[j] = Faces.normal[f]
            else:
                ss[j] = -1.0*Faces.normal[f]
                
        Cells.nf.append(ss)
    
    
    #TODO improve the opp node creation
    Faces.oppnodeid = [[] for i in range(nbfaces)]
    for i in range(nbelements):
        f1 = Cells.faceid[i][0]; f2 = Cells.faceid[i][1]; f3 = Cells.faceid[i][2]
        n1 = Cells.nodeid[i][0]; n2 = Cells.nodeid[i][1]; n3 = Cells.nodeid[i][2] 
        
        if n1 not in Faces.nodeid[f1] :
            Faces.oppnodeid[f1].append(n1)
        if n1 not in Faces.nodeid[f2] :
            Faces.oppnodeid[f2].append(n1)
        if n1 not in Faces.nodeid[f3] :
            Faces.oppnodeid[f3].append(n1)
        
        if n2 not in Faces.nodeid[f1] :
            Faces.oppnodeid[f1].append(n2)
        if n2 not in Faces.nodeid[f2] :
            Faces.oppnodeid[f2].append(n2)
        if n2 not in Faces.nodeid[f3] :
            Faces.oppnodeid[f3].append(n2)
        
        if n3 not in Faces.nodeid[f1] :
            Faces.oppnodeid[f1].append(n3)
        if n3 not in Faces.nodeid[f2] :
            Faces.oppnodeid[f2].append(n3)
        if n3 not in Faces.nodeid[f3] :
            Faces.oppnodeid[f3].append(n3)
            
    for i in range(nbfaces):
        if len(Faces.oppnodeid[i]) < 2:
            Faces.oppnodeid[i].append(-1)
        
        

    cells = Cells(np.asarray(Cells.nodeid), np.asarray(Cells.faceid), np.asarray(Cells.center),
                  np.asarray(Cells.volume), np.asarray(Cells.cellfid), 
                  np.asarray(Cells.cellnid), np.asarray(Cells.nf))
    
    nodes = Nodes(np.asarray(Nodes.vertex), np.asarray(Nodes.name), np.asarray(Nodes.cellid),
                  np.asarray(Nodes.ghostcenter))

    faces = Faces(np.asarray(Faces.nodeid), np.asarray(Faces.cellid), np.asarray(Faces.name),
                  np.asarray(Faces.normal), np.asarray(Faces.mesure), np.asarray(Faces.center), 
                  Faces.bound, np.asarray(Faces.ghostcenter), np.asarray(Faces.oppnodeid) )

    return cells, nodes, faces

def Create3DStructure(file):

    #Lecture des cellules à partir du fichier mesh..txt
    for line in file:
        #read elements
        if line == "elements\n":
            continue
        if line == "endelements\n":
            continue
        if line == "nodes\n":
            break
        Cells.nodeid.append([int(x) for x in line.split()])
    #Lecture des coordonnées des noeuds à partir du fichier mesh..txt
    for line in file:
        #read Nodes
        if line == "nodes\n":
            continue
        if line == "endnodes\n":
            continue
        if line == "halosint\n":
            break
        Nodes.vertex.append([np.double(x) for x in line.split()])
        Nodes.name = [0]*len(Nodes.vertex)

    nbelements = len(Cells.nodeid)
    nbnodes    = len(Nodes.vertex)
    
    
    for i in range(nbnodes):
        Nodes.name[i] = int(Nodes.vertex[i][3])
    
    #calcul du barycentre
    for i in range(nbelements):
        
        s_1 = Cells.nodeid[i][0]
        s_2 = Cells.nodeid[i][1]
        s_3 = Cells.nodeid[i][2]
        s_4 = Cells.nodeid[i][3]
        
        a = np.asarray(Nodes.vertex[s_1])
        b = np.asarray(Nodes.vertex[s_2])
        c = np.asarray(Nodes.vertex[s_3])
        d = np.asarray(Nodes.vertex[s_4])
        
        x_1 = Nodes.vertex[s_1][0]; y_1 = Nodes.vertex[s_1][1]; z_1 = Nodes.vertex[s_1][2]
        x_2 = Nodes.vertex[s_2][0]; y_2 = Nodes.vertex[s_2][1]; z_2 = Nodes.vertex[s_2][2]
        x_3 = Nodes.vertex[s_3][0]; y_3 = Nodes.vertex[s_3][1]; z_3 = Nodes.vertex[s_3][2]
        x_4 = Nodes.vertex[s_4][0]; y_4 = Nodes.vertex[s_4][1]; z_4 = Nodes.vertex[s_4][2]
        
        Cells.center.append((1./4*(x_1 + x_2 + x_3 + x_4), 1./4*(y_1 + y_2 + y_3 + y_4), 1./4*(z_1 + z_2 + z_3 + z_4)))
        Cells.volume.append(1./6*np.fabs(det_vec_3d(b-a, c-a, d-a)))
        
        
    tmp = [[] for i in range(nbnodes)]
    longn = [0]*nbnodes
    for i in range(nbelements):
        for j in range(4):
            tmp[Cells.nodeid[i][j]].append(i)
            longn[Cells.nodeid[i][j]] = longn[Cells.nodeid[i][j]] + 1
    
    longc = [0]*nbelements
    tmp2 = [[] for i in range(nbelements)]
    for i in range(nbelements):
        for j in range(4):
            for k in range(len(tmp[Cells.nodeid[i][j]])):
                if (tmp[Cells.nodeid[i][j]][k] not in tmp2[i] and  tmp[Cells.nodeid[i][j]][k] != i):
                    tmp2[i].append(tmp[Cells.nodeid[i][j]][k])
                    longc[i] = longc[i] + 1

    maxlongn = max(longn)
    Nodes.cellid = [[-1]*maxlongn for i in range(nbnodes)]

    for i in range(len(tmp)):
        for j in range(len(tmp[i])):
            Nodes.cellid[i][j] = tmp[i][j]
            
    for i in range(nbnodes):
        Nodes.cellid[i].append(longn[i])
    
    maxlongc = max(longc)
    Cells.cellnid = [[-1]*maxlongc for i in range(nbelements)]
    
    for i in range(len(tmp2)):
        for j in range(len(tmp2[i])):
            Cells.cellnid[i][j] = tmp2[i][j]

    #Création des faces
    cellule = Cells.nodeid
    cellf = []
    faces = []
    
    k = 0
    for i in range(len(cellule)):
        faces.append([cellule[i][0], cellule[i][1], cellule[i][2]])
        faces.append([cellule[i][2], cellule[i][3], cellule[i][0]])
        faces.append([cellule[i][0], cellule[i][1], cellule[i][3]])
        faces.append([cellule[i][3], cellule[i][1], cellule[i][2]])
        cellf.append([faces[k], faces[k+1], faces[k+2], faces[k+3]])
        k = k+4
    for i in range(len(faces)):
        faces[i].sort()
    faces = set(tuple(x) for x in faces)
    faces = list(faces)
    
    nbfaces = len(faces)
    
    facesdict = OrderedDict()
    for i in range(nbfaces):
        facesdict[faces[i]] = i
        Faces.nodeid.append(faces[i])
        Faces.cellid.append((-1, -1))

    #Création des 4 faces de chaque cellule
    for i in range(len(cellf)):
        Cells.faceid.append([facesdict.get(tuple(cellf[i][0])), facesdict.get(tuple(cellf[i][1])),
                             facesdict.get(tuple(cellf[i][2])), facesdict.get(tuple(cellf[i][3]))])

    #Faces.cellid = [(-1, -1)]*len(faces)
    for i in range(nbelements):
        for j in range(4):
            if Faces.cellid[Cells.faceid[i][j]] == (-1, -1):
                Faces.cellid[Cells.faceid[i][j]] = (i, -1)

            if Faces.cellid[Cells.faceid[i][j]][0] != i:
                Faces.cellid[Cells.faceid[i][j]] = (Faces.cellid[Cells.faceid[i][j]][0], i)

    Cells.cellfid = [[] for i in range(nbelements)]
    #Création des 3/4 triangles voisins par face
    for i in range(nbelements):
        for j in range(4):
            f = Cells.faceid[i][j]
            if Faces.cellid[f][1] != -1:
                if i == Faces.cellid[f][0]:
                    Cells.cellfid[i].append(Faces.cellid[f][1])
                else:
                    Cells.cellfid[i].append(Faces.cellid[f][0])
            else:
                Cells.cellfid[i].append(-1)
    
    for i in range(nbfaces):
        Faces.name.append(0)
        if (Faces.cellid[i][1] == -1):
            Faces.name[i]=1
            
        normal, mesure = normalvector3d(np.asarray(Nodes.vertex[faces[i][0]]),
                                         np.asarray(Nodes.vertex[faces[i][1]]),
                                         np.asarray(Nodes.vertex[faces[i][2]]),
                                         Cells.center[Faces.cellid[i][0]])
        Faces.normal.append(normal)
        Faces.mesure.append(mesure)
        Faces.center.append([1./3 * (Nodes.vertex[faces[i][0]][0] + Nodes.vertex[faces[i][1]][0] + Nodes.vertex[faces[i][2]][0]),
                             1./3 * (Nodes.vertex[faces[i][0]][1] + Nodes.vertex[faces[i][1]][1] + Nodes.vertex[faces[i][2]][1]),
                             1./3 * (Nodes.vertex[faces[i][0]][2] + Nodes.vertex[faces[i][1]][2] + Nodes.vertex[faces[i][2]][2])
                             ])
    
    
    
    for i in range(nbfaces):
        if Faces.name[i] != 0:
            Faces.bound += 1
        
    #compute the outgoing normal faces for each cell
    for i in range(nbelements):
        ss = np.zeros((4, 3))

        G = np.asarray(Cells.center[i])

        for j in range(4):
            f = Cells.faceid[i][j]
            c = Faces.center[f]

            if np.dot(G-c, Faces.normal[f]) < 0.0:
                ss[j] = Faces.normal[f]
            else:
                ss[j] = -1.0*Faces.normal[f]
                
        Cells.nf.append(ss)
        
    #TODO improve the opp node creation
    Faces.oppnodeid = [[] for i in range(nbfaces)]
    for i in range(nbelements):
        f1 = Cells.faceid[i][0]; f2 = Cells.faceid[i][1]; f3 = Cells.faceid[i][2]; f4 = Cells.faceid[i][3]; 
        n1 = Cells.nodeid[i][0]; n2 = Cells.nodeid[i][1]; n3 = Cells.nodeid[i][2]; n4 = Cells.nodeid[i][3]; 
        
        if n1 not in Faces.nodeid[f1] :
            Faces.oppnodeid[f1].append(n1)
        if n1 not in Faces.nodeid[f2] :
            Faces.oppnodeid[f2].append(n1)
        if n1 not in Faces.nodeid[f3] :
            Faces.oppnodeid[f3].append(n1)
        if n1 not in Faces.nodeid[f4] :
            Faces.oppnodeid[f4].append(n1)
        
        if n2 not in Faces.nodeid[f1] :
            Faces.oppnodeid[f1].append(n2)
        if n2 not in Faces.nodeid[f2] :
            Faces.oppnodeid[f2].append(n2)
        if n2 not in Faces.nodeid[f3] :
            Faces.oppnodeid[f3].append(n2)
        if n2 not in Faces.nodeid[f4] :
            Faces.oppnodeid[f4].append(n2)
        
        if n3 not in Faces.nodeid[f1] :
            Faces.oppnodeid[f1].append(n3)
        if n3 not in Faces.nodeid[f2] :
            Faces.oppnodeid[f2].append(n3)
        if n3 not in Faces.nodeid[f3] :
            Faces.oppnodeid[f3].append(n3)
        if n3 not in Faces.nodeid[f4] :
            Faces.oppnodeid[f4].append(n3)
        
        if n4 not in Faces.nodeid[f1] :
            Faces.oppnodeid[f1].append(n4)
        if n4 not in Faces.nodeid[f2] :
            Faces.oppnodeid[f2].append(n4)
        if n4 not in Faces.nodeid[f3] :
            Faces.oppnodeid[f3].append(n4)
        if n4 not in Faces.nodeid[f4] :
            Faces.oppnodeid[f4].append(n4)
        
            
    for i in range(nbfaces):
        if len(Faces.oppnodeid[i]) < 2:
            Faces.oppnodeid[i].append(-1)
            
    cells = Cells(np.asarray(Cells.nodeid), np.asarray(Cells.faceid), np.asarray(Cells.center),
                  np.asarray(Cells.volume), np.asarray(Cells.cellfid), 
                  np.asarray(Cells.cellnid), np.asarray(Cells.nf))
    
    nodes = Nodes(np.asarray(Nodes.vertex), np.asarray(Nodes.name), np.asarray(Nodes.cellid),
                  np.asarray(Nodes.ghostcenter))

    faces = Faces(np.asarray(Faces.nodeid), np.asarray(Faces.cellid), np.asarray(Faces.name),
                  np.asarray(Faces.normal), np.asarray(Faces.mesure), np.asarray(Faces.center), 
                  Faces.bound, np.asarray(Faces.ghostcenter), np.asarray(Faces.oppnodeid) )
    
    return cells, nodes, faces

def generate_structure(dim):
    
    filename = 'mesh.txt'
    txt_file = open(filename)

    if dim == 2:
        cells, nodes, faces = Create2DStructure(txt_file)
    elif dim == 3:
        cells, nodes, faces = Create3DStructure(txt_file)

    grid = {}

    grid["cells"] = cells
    grid["nodes"] = nodes
    grid["faces"] = faces

    txt_file.close()

    return grid

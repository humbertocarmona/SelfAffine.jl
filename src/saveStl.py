#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import meshio as msh
import pandas as pd


# In[2]:


def csv2stl(csvfile, stlfile):
    z = pd.read_csv(csvfile).to_numpy()
    points = []
    nx = z.shape[0]
    ny = z.shape[1]

    for j in range(ny):
        for i in range(nx):
            points.append([i-1,j-1,z[j,i]])


    def ij2n(i,j):
        return i+(nx)*(j)

    tri=[]
    for j in range(ny-1):
        for i in range(nx-1):
            t = ij2n(i,j)
            u = ij2n(i+1,j+1)
            v = ij2n(i, j+1)
            tri.append([t,u,v])

            t = ij2n(i,j)
            u = ij2n(i+1,j)
            v = ij2n(i+1, j+1)
            tri.append([t,u,v])

    cells=dict()
    cells["triangle"]=np.array(tri)
    msh.write_points_cells(stlfile, np.array(points), cells, point_data={'pos':np.array(points)})
    print("ok")


# In[3]:


csv2stl("04.csv", "04.stl")


# In[5]:


csv2stl("05.csv", "n05.stl")


# In[ ]:





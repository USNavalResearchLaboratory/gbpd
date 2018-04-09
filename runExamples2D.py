'''
/****************************************************************************
 *****                                                                  *****
 *****                   Classification: UNCLASSIFIED                   *****
 *****                    Classified By:                                *****
 *****                    Declassify On:                                *****
 *****                                                                  *****
 ****************************************************************************
 *
 *
 * Developed by: Naval Research Laboratory,
 *               Multifunctional Materials Branch, Code 6353
 *               4555 Overlook Ave.
 *               Washington, D.C. 20375-5339
 *
 *
 * The U.S. Government retains all rights to use, duplicate, distribute,
 * disclose, or release this software.
 *
 ''' 

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 14:22:50 2018

@author: kteferra
"""

import os
#import importlib
import numpy as np
import matplotlib.pyplot as plt



'============================================================================'
'This part checks the fit with data and GBPD: beta- Ti'

os.chdir('/Users/kteferra/Documents/research/projects/SIMPL/codes/Tessellation/')

import Tessellation2D
os.chdir('/Users/kteferra/Documents/research/projects/SIMPL/data/Titanium/')

zdata = Tessellation2D.cellidData2D()

fil1='bTi_slice.npy'
zdata.voxelID=np.load(fil1)
zdata.dims = np.asarray(zdata.voxelID.shape)
zdata.spacing=np.array([.665,.665])
pbounds = np.array([[0,zdata.dims[0]],[0,zdata.dims[1]]])
nd = zdata.dims
dx = zdata.spacing
zgbpd = Tessellation2D.GBPD()
zgbpd.fitParameters(zdata)
zgbpd.generateMicrostructure(nd,dx)
tmp = np.gradient(zdata.voxelID)
z1 = np.zeros(zdata.dims[0:2])
z1[np.where(tmp[0]!=0)] = 1
z1[np.where(tmp[1]!=0)] = 1
tmp = np.gradient(zgbpd.voxelID)
zval =1
z2=np.ones(np.concatenate((zdata.dims[0:2],np.array([3]))))*zval
z2a = np.ones(zdata.dims[0:2])*zval
z2a[np.where(tmp[0]!=0)] = 0.0
z2a[np.where(tmp[1]!=0)] = 0.0
z2[:,:,1] = z2a
z2[:,:,2] = z2a
myfig=plt.figure(figsize=(12,10))
myimg0=plt.imshow(z1,cmap='Greys')#,extent=(2,50,2,50))
plt.savefig('bTi_slice.png',dpi=300)
myimg1=plt.imshow(z2,alpha=.5,cmap='Reds')#,extent=(2,50,2,50))
plt.tick_params(axis='both',which='both',labelbottom='off',labelleft='off',bottom='off',left='off')
plt.show()
plt.savefig('GBPDfit_bTi_slice70.png',dpi=300)
plt.close()

'compute percent voxel mismatch'
z1 = np.copy(zdata.voxelID)
z2 = np.copy(zgbpd.voxelID)
z1 = np.reshape(z1,zdata.dims[0]*zdata.dims[1])
z2 = np.reshape(z2,zdata.dims[0]*zdata.dims[1])
b = zgbpd.dataID[np.where(zgbpd.gID==1)]
numMatch= 0.0
for j in range(0,len(np.where(zgbpd.gID==1)[0])):  
    k1 = np.where(z1==b[j])[0]
    k2 = np.where(z2==b[j])[0]
    numMatch = numMatch + len(np.intersect1d(k1,k2))
percentVoxelMatch =numMatch/zdata.dims[0]/zdata.dims[1]    
    

zgbpd.computeNeighbors()
zdata.computeNeighbors()
neighborStats = np.zeros([zdata.Ngrain,3]) # col1: if all correct
                                             # col2: if only 1 wrong
                                             # col3: difference in # of grains
bbad = np.where(zgbpd.gID==0)[0]
for j in range(0,zdata.Ngrain):
    k1 = np.setdiff1d(zdata.neighbors[j],bbad)
    k2 = np.setdiff1d(zgbpd.neighbors[j],bbad)
    if (len(np.setdiff1d(k1,k2))==0): neighborStats[j,0]=1.0
    if (len(np.setdiff1d(k1,k2))<=1): neighborStats[j,1]=1.0
    neighborStats[j,2] = len(k2)-len(k1)

k1 = neighborStats[np.where(zgbpd.gID==1)[0],0]
percentAllNeighbor = np.sum(k1)/len(np.where(zgbpd.gID==1)[0])
k1 = neighborStats[np.where(zgbpd.gID==1)[0],1]
percentOneNeighbor = np.sum(k1)/len(np.where(zgbpd.gID==1)[0])

k1 = np.where( (neighborStats[:,2]>0) & (zgbpd.gID==1)      )[0]
percentNeighborExcess = np.sum(neighborStats[k1,2])/ \
    len(np.where(zgbpd.gID==1)[0])
k1 = np.where( (neighborStats[:,2]<0) & (zgbpd.gID==1)      )[0]
percentNeighborLess = np.abs(np.sum(neighborStats[k1,2]))/ \
    len(np.where(zgbpd.gID==1)[0])

print('beta-Ti results:')
print('Percent voxel match=',percentVoxelMatch)
print('Percent grains with all neighbors correct=',percentAllNeighbor)
print('Percent grains with at least 1 neighbor correct=',percentOneNeighbor)
print('#Erroneous extra neighbors per grain=',percentNeighborExcess)
print('#Erroneous missing neighbors  per grain=',percentNeighborLess)
      

'This part checks the fit with data and GBPD: beta- Ti'
'============================================================================'

'============================================================================'
'This part checks the fit with data and GBPD: IN100'


zdata = Tessellation2D.cellidData3D()

fil='Small_IN100_20140317.am'

zdata.readAMFormatIN100(fil)

nSlice = np.arange(10,100,10)
j = np.floor(np.random.rand()*len(nSlice)).astype(int)
pbounds = np.floor(np.array([[5, 42],[7, 43]])/zdata.spacing[0:2]).astype(int)
zgbpd = [None]*len(nSlice)
jstar = 2 # np.copy(j)
zdata2D = Tessellation2D.cellidData2D()
zdata2D.convertFrom3D(zdata,nSlice[jstar],pbounds)
nd = np.array(zdata2D.voxelID.shape)
dx =zdata2D.spacing

'''
for jstar in range(0,len(nSlice)):
    zgbpd[jstar] = Tessellation2D.GBPD()
    zgbpd[jstar].fitParameters(zdata2D)
    zgbpd[jstar].generateMicrostructure(nd,dx)
'''
zgbpd[jstar] = Tessellation2D.GBPD()
zgbpd[jstar].fitParameters(zdata2D)
zgbpd[jstar].generateMicrostructure(nd,dx)

tmp = np.gradient(zdata2D.voxelID)
z1 = np.zeros(zdata2D.dims[0:2])
z1[np.where(tmp[0]!=0)] = 1
z1[np.where(tmp[1]!=0)] = 1
tmp = np.gradient(zgbpd[jstar].voxelID)
zval =1
z2=np.ones(np.concatenate((zdata2D.dims[0:2],np.array([3]))))*zval
z2a = np.ones(zdata2D.dims[0:2])*zval
z2a[np.where(tmp[0]!=0)] = 0.0
z2a[np.where(tmp[1]!=0)] = 0.0
z2[:,:,1] = z2a
z2[:,:,2] = z2a


#z1 = z1[pbounds[0,0]:pbounds[0,1],pbounds[1,0]:pbounds[1,1]]
#z2a = np.copy(z2)
#z2= np.zeros(np.concatenate((z1.shape,np.array([3]))))
#z2[:,:,0]=z2a[pbounds[0,0]:pbounds[0,1],pbounds[1,0]:pbounds[1,1],0]
#z2[:,:,1]=z2a[pbounds[0,0]:pbounds[0,1],pbounds[1,0]:pbounds[1,1],1]
#z2[:,:,2]=z2a[pbounds[0,0]:pbounds[0,1],pbounds[1,0]:pbounds[1,1],2]

myfig=plt.figure
#plt.subplot(121)
myimg0=plt.imshow(z1,cmap='Greys',extent=(2,50,2,50))
#plt.subplot(122)
myimg1=plt.imshow(z2,alpha=.5,cmap='Reds',extent=(2,50,2,50))
plt.tick_params(axis='both',which='both',labelbottom='off',labelleft='off',bottom='off',left='off')

plt.show()
plt.savefig('GBPDfit_IN100_slice'+str(nSlice[j])+'.png',dpi=300)


'compute percent voxel mismatch'
z1 = np.copy(zdata2D.voxelID)
z2 = np.copy(zgbpd[jstar].voxelID)
z1 = np.reshape(z1,z1.shape[0]*z1.shape[1])
z2 = np.reshape(z2,z2.shape[0]*z2.shape[1])
b = zgbpd[jstar].dataID[np.where(zgbpd[jstar].gID==1)]
numMatch= 0.0
for j in range(0,len(np.where(zgbpd[jstar].gID==1)[0])):  
    k1 = np.where(z1==b[j])[0]
    k2 = np.where(z2==b[j])[0]
    numMatch = numMatch + len(np.intersect1d(k1,k2))
percentVoxelMatch =numMatch/z1.shape[0]    
    

zgbpd[jstar].computeNeighbors()
zdata2D.computeNeighbors()
neighborStats = np.zeros([zdata2D.Ngrain,3]) # col1: if all correct
                                             # col2: if only 1 wrong
                                             # col3: difference in # of grains
bbad = np.where(zgbpd[jstar].gID==0)[0]
for j in range(0,zdata2D.Ngrain):
    k1 = np.setdiff1d(zdata2D.neighbors[j],bbad)
    k2 = np.setdiff1d(zgbpd[jstar].neighbors[j],bbad)
    if (len(np.setdiff1d(k1,k2))==0): neighborStats[j,0]=1.0
    if (len(np.setdiff1d(k1,k2))<=1): neighborStats[j,1]=1.0
    neighborStats[j,2] = len(k2)-len(k1)

k1 = neighborStats[np.where(zgbpd[jstar].gID==1)[0],0]
percentAllNeighbor = np.sum(k1)/len(np.where(zgbpd[jstar].gID==1)[0])
k1 = neighborStats[np.where(zgbpd[jstar].gID==1)[0],1]
percentOneNeighbor = np.sum(k1)/len(np.where(zgbpd[jstar].gID==1)[0])

k1 = np.where( (neighborStats[:,2]>0) & (zgbpd[jstar].gID==1)      )[0]
percentNeighborExcess = np.sum(neighborStats[k1,2])/ \
    len(np.where(zgbpd[jstar].gID==1)[0])
k1 = np.where( (neighborStats[:,2]<0) & (zgbpd[jstar].gID==1)      )[0]
percentNeighborLess = np.abs(np.sum(neighborStats[k1,2]))/ \
    len(np.where(zgbpd[jstar].gID==1)[0])

print('IN100 results:')
print('Percent voxel match=',percentVoxelMatch)
print('Percent grains with all neighbors correct=',percentAllNeighbor)
print('Percent grains with at least 1 neighbor correct=',percentOneNeighbor)
print('#Erroneous extra neighbors per grain=',percentNeighborExcess)
print('#Erroneous missing neighbors  per grain=',percentNeighborLess)
      

'This part checks the fit with data and GBPD: IN100'
'============================================================================'


'============================================================================'
'This part checks the fit with data and GBPD: Cr-Al'

zdata = Tessellation2D.cellidData2D()

fil1='Elongated_2.png'
fil2='Elongated_celliddat.txt'


zdata.readImage(fil1)
pbounds = np.array([[20,500],[0,485]])
zdata.readCellid(fil2,pbounds)

nd = zdata.dims
dx = np.array([1.0,1.0])

zgbpd = Tessellation2D.GBPD()
zgbpd.fitParameters(zdata)
zgbpd.generateMicrostructure(nd,dx)

tmp = np.gradient(zdata.voxelID)
z1 = np.zeros(zdata.dims[0:2])
z1[np.where(tmp[0]!=0)] = 1
z1[np.where(tmp[1]!=0)] = 1
tmp = np.gradient(zgbpd.voxelID)
zval =1
z2=np.ones(np.concatenate((zdata.dims[0:2],np.array([3]))))*zval
z2a = np.ones(zdata.dims[0:2])*zval
z2a[np.where(tmp[0]!=0)] = 0.0
z2a[np.where(tmp[1]!=0)] = 0.0
z2[:,:,1] = z2a
z2[:,:,2] = z2a


myfig=plt.figure
#plt.subplot(121)
myimg0=plt.imshow(z1,cmap='Greys')#,extent=(2,50,2,50))
#plt.subplot(122)
myimg1=plt.imshow(z2,alpha=.5,cmap='Reds')#,extent=(2,50,2,50))
plt.tick_params(axis='both',which='both',labelbottom='off',labelleft='off',bottom='off',left='off')
plt.show()
plt.savefig('GBPDfit_CRAl.png',dpi=300)


'compute percent voxel mismatch'
z1 = np.copy(zdata.voxelID)
z2 = np.copy(zgbpd.voxelID)
z1 = np.reshape(z1,z1.shape[0]*z1.shape[1])
z2 = np.reshape(z2,z2.shape[0]*z2.shape[1])
b = zgbpd.dataID[np.where(zgbpd.gID==1)]
numMatch= 0.0
for j in range(0,len(np.where(zgbpd.gID==1)[0])):  
    k1 = np.where(z1==b[j])[0]
    k2 = np.where(z2==b[j])[0]
    numMatch = numMatch + len(np.intersect1d(k1,k2))
percentVoxelMatch =numMatch/z1.shape[0]    
    

zgbpd.computeNeighbors()
zdata.computeNeighbors()
neighborStats = np.zeros([zdata.Ngrain,3]) # col1: if all correct
                                             # col2: if only 1 wrong
                                             # col3: difference in # of grains
bbad = np.where(zgbpd.gID==0)[0]
for j in range(0,zdata.Ngrain):
    k1 = np.setdiff1d(zdata.neighbors[j],bbad)
    k2 = np.setdiff1d(zgbpd.neighbors[j],bbad)
    if (len(np.setdiff1d(k1,k2))==0): neighborStats[j,0]=1.0
    if (len(np.setdiff1d(k1,k2))<=1): neighborStats[j,1]=1.0
    neighborStats[j,2] = len(k2)-len(k1)

k1 = neighborStats[np.where(zgbpd.gID==1)[0],0]
percentAllNeighbor = np.sum(k1)/len(np.where(zgbpd.gID==1)[0])
k1 = neighborStats[np.where(zgbpd.gID==1)[0],1]
percentOneNeighbor = np.sum(k1)/len(np.where(zgbpd.gID==1)[0])

k1 = np.where( (neighborStats[:,2]>0) & (zgbpd.gID==1)      )[0]
percentNeighborExcess = np.sum(neighborStats[k1,2])/ \
    len(np.where(zgbpd.gID==1)[0])
k1 = np.where( (neighborStats[:,2]<0) & (zgbpd.gID==1)      )[0]
percentNeighborLess = np.abs(np.sum(neighborStats[k1,2]))/ \
    len(np.where(zgbpd.gID==1)[0])

print('Cr-Al results:')
print('Percent voxel match=',percentVoxelMatch)
print('Percent grains with all neighbors correct=',percentAllNeighbor)
print('Percent grains with at least 1 neighbor correct=',percentOneNeighbor)
print('#Erroneous extra neighbors per grain=',percentNeighborExcess)
print('#Erroneous missing neighbors  per grain=',percentNeighborLess)
      



'This part checks the fit with data and GBPD: Cr-Al'
'============================================================================'

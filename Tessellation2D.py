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
Created on Tue Jan  2 12:10:25 2018

@author: kteferra
"""
import numpy as np
import os as os


' ==========================================================================='
' Obtain the data'



' -------------------------------------------------------------------------'
' cellid of 3D data'    
class cellidData3D:

    def readAMFormatIN100(self,filname):
        # filname is string of path to file name
        self.dims=np.array([189,201,117],dtype=int)
        self.spacing = np.array([.25,.25,.25])
        self.nVoxel = self.dims[0]*self.dims[1]*self.dims[2]
        self.voxelID = np.zeros([self.nVoxel],dtype=int)
        b=np.loadtxt(filname,dtype=int,skiprows=23)
        self.voxelID =np.reshape(np.reshape(b,self.nVoxel),self.dims,order='F')
        self.Ngrain = len(np.unique(self.voxelID))-1
        # note subtract 1 above because 0 is background

    ' read in cell data from hdf5 file Dave gave me of Ti data'
    def readHDF5Ti(self,filename):
        import h5py
        f=h5py.File(filename,'r')
        z=f.get('VoxelData/SCALARS/GrainID')
        self.dims = np.array(f.get('VoxelData/DIMENSIONS'))
        self.spacing = np.array(f.get('VoxelData/SPACING'))
        self.voxelID = np.reshape(z,self.dims[[2,1,0]]).transpose()
        
class cellidData2D:

    def readImage(self,fil):
        from skimage import io
        self.img = io.imread(fil, as_grey=True)
        self.dims = np.array(self.img.shape)
        
    def readCellid(self,fil,pbounds):
        z=np.loadtxt(fil,dtype=float)
        z= np.reshape(z,np.array([self.dims[1],self.dims[0]])).transpose()
        self.voxelID = z[pbounds[0,0]:pbounds[0,1], \
                           pbounds[1,0]:pbounds[1,1]]
        self.Ngrain = len(np.unique(self.voxelID))
        self.spacing = np.array([1.0,1.0])
        self.dims = np.array(self.voxelID.shape)

        # note subtract 1 above because 0 is background
        
    def convertFrom3D(self,zdata3D,jslice,pbounds):
        self.spacing = zdata3D.spacing[0:2]
        self.voxelID = zdata3D.voxelID[pbounds[0,0]:pbounds[0,1], \
                           pbounds[1,0]:pbounds[1,1],jslice]
        self.dims = np.array(self.voxelID.shape)
        self.Ngrain = len(np.unique(self.voxelID))

    def computeNeighbors(self):
        b=np.unique(self.voxelID)
        self.neighbors = [np.array(-1)]*len(b)        
        for j1 in range(0,self.dims[0]-1):
            for j2 in range(0,self.dims[1]-1):        
                k = np.where(b==self.voxelID[j1,j2])[0][0]
                k1 = np.where(b==self.voxelID[j1+1,j2])[0][0]
                k2 = np.where(b==self.voxelID[j1,j2+1])[0][0]
                self.neighbors[k]=np.append(self.neighbors[k],k1)
                self.neighbors[k]=np.append(self.neighbors[k],k2)
                self.neighbors[k1]=np.append(self.neighbors[k1],k)
                self.neighbors[k2]=np.append(self.neighbors[k2],k)                
        for k in range(0,self.Ngrain):
            self.neighbors[k] = np.unique( \
                 self.neighbors[k][1:len(self.neighbors[k])] )
        
        
' cellid of 3D data'    
' -------------------------------------------------------------------------'




' Obtain the data'
' ==========================================================================='

' -------------------------------------------------------------------------'
' Define the 2D GBPD model (generalized balanced power diagrams) see Alpers'

class GBPD:

    def fitParameters(self,cellDat2D):    
        cellDat = cellDat2D.voxelID
        self.Ngrain = len(np.unique(cellDat))
        spacing =cellDat2D.spacing
        self.sites = np.zeros([self.Ngrain,2])
        self.A = np.zeros([self.Ngrain,2,2])
        self.sigma = np.zeros(self.Ngrain)
        self.EulerAngles = np.zeros([self.Ngrain,3])
        b= np.unique(cellDat)
        self.gID = np.ones(self.Ngrain,dtype=int)
        ' get centroids'
        '''as a result of taking 2D slices from 3D images, there are some 
        grains that do not have any volume. these will be removed from the
        fitting
        '''
        mincheck = 5e-3
        mincheck2 = 5
        for j in range(0,self.Ngrain):
            x = np.asarray(np.where(cellDat==b[j]))
            gvol = x.shape[1]*spacing[0]*spacing[1]
            coors = np.concatenate(([x[0,:]*spacing[0]], \
                   [x[1,:]*spacing[1]]))
            self.sites[j,:] = np.mean(coors,axis=1)
            ' get covariance'
            moi2nd = np.array([[1.0,0.0],[0.0,1.0]])
            if (coors.shape[1]>mincheck2): moi2nd = np.cov(coors) 
            [w,c] = np.linalg.eig(moi2nd)
            i1 = np.argsort(w) 
            w=4.0*w[i1]
            c=c[:,i1]
            la=np.zeros([2,2])
            if ((w[0]>mincheck) and (b[j]!=0)):
                la[np.diag_indices(2)]=w**(-1)
                self.A[j,:,:] = np.matmul(c,np.matmul(la,c.transpose()))
                self.sigma[j] = gvol/(np.pi*np.sqrt(w[0]*w[1]))
            if (coors.shape[1]<=mincheck2 or w[0]<mincheck or b[j]==0): self.gID[j]=0
        self.dataID = b
        
    def generateMicrostructure(self,Ndim,dx):           
        'Ndim - 2D np array of # of voxels in x & y'
        'dx - 2D np array of spacing in x & y'
        self.dims = Ndim
        self.spacing = dx
        self.voxelID = np.zeros(self.dims,dtype=int)
        for j1 in range(0,self.dims[0]):
            for j2 in range(0,self.dims[1]):
                node  = np.array([j1*self.spacing[0],j2*self.spacing[0]])
                x = self.sites - node
                d = x[:,0]*(self.A[:,0,0]*x[:,0]+self.A[:,0,1]*x[:,1]) + \
                    x[:,1]*(self.A[:,1,0]*x[:,0]+self.A[:,1,1]*x[:,1]) - \
                    self.sigma[:]
                d[np.where(self.gID==0)[0]] = np.max(d)
                self.voxelID[j1,j2] = self.dataID[np.where(d==np.min(d))[0][0]]
        for j in range(0,self.Ngrain):
            k = np.where(self.voxelID==self.dataID[j])[0]
            if (len(k)==0): self.gID[j] = 0
        
    def computeNeighbors(self):
        b = self.dataID
        self.neighbors = [np.array(-1)]*len(b)
        for j1 in range(0,self.dims[0]-1):
            for j2 in range(0,self.dims[1]-1):        
                k = np.where(b==self.voxelID[j1,j2])[0][0]
                k1 = np.where(b==self.voxelID[j1+1,j2])[0][0]
                k2 = np.where(b==self.voxelID[j1,j2+1])[0][0]
                self.neighbors[k]=np.append(self.neighbors[k],k1)
                self.neighbors[k]=np.append(self.neighbors[k],k2)
                self.neighbors[k1]=np.append(self.neighbors[k1],k)
                self.neighbors[k2]=np.append(self.neighbors[k2],k)                
        for k in range(0,self.Ngrain):
            if (self.gID[k] ==1):
                self.neighbors[k] = np.unique( \
                 self.neighbors[k][1:len(self.neighbors[k])] )
        
            

' class GBPD'
' Define the 2D GBPD model (generalized balanced power diagrams) see Alpers'
' -------------------------------------------------------------------------'

    
    

# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 15:19:17 2023

@author: lshen
"""
import numpy as np
import cmath

def croot(k, n):
    if n<=0:
        return None
    return cmath.exp((2 * cmath.pi * 1j * k) / n)


def boundary(complex,p):
    maxdim = len(max(complex, key=len))
    simplices = [sorted([spx for spx in complex if len(spx)==i]) for i in range(1,maxdim+1)]
    #print(simplices)
    
    bnd = [np.zeros((0,len(simplices[0])))]
    
    for spx_k, spx_kp1 in zip(simplices, simplices[1:]):
        #print(spx_k,spx_kp1)
        mtx = []
        for sigma in spx_kp1:
            faces = get_faces(sigma)
            mtx.append([get_coeff(spx, faces,p) for spx in spx_k])
            #print(mtx)
        bnd.append(np.array(mtx).T)

    return bnd




def get_faces(lst):
    return [lst[:i] + lst[i+1:] for i in range(len(lst))]


def get_coeff(simplex, faces,p):
    if simplex in faces:
        idxs = [i for i in range(len(faces)) if faces[i] == simplex]
        #print(idx,simplex)
        return sum([croot(idx,p) for idx in idxs])
    else:
        return 0
  
def matrix_rank(A):
    if A.size == 0:
        return 0
    else:
        return np.linalg.matrix_rank(A,tol = 1e-10)

def betti_number(f,g):#kerf/img
    rk_f = matrix_rank(f)
    
    dim_kerf = f.shape[1]-rk_f
    #print("f",f.size,f.shape,dim_kerf)
    rk_g = matrix_rank(g)
    #print("g",g.size,rk_g)
    return dim_kerf-rk_g



def mat_multiply(bnd_op,dim,q):
    f = bnd_op.get(dim-q+1)
    #print(dim-q+1)
    for i in range(1,q):
        #print(dim+i-q+1)
        f = f@bnd_op.get(dim+i-q+1)
        

    return f

def calculate_betti_laplacian(bnd_op,n,q,max_dim=1):#ker d^q/imd^(n-q) for q<n
    #print(q)
    #print(bnd)
    dim = bnd_op.max_dim
    #print(dim)
    betti=[]
    gamma_min=[]
    gamma_max = []
    gamma_mean = []
    gamma_std = []
    for i in range(max_dim+1):#only consider b0,b1,l0,l1 if max_dim =1
        #print("f")
        f = mat_multiply(bnd_op,i,q)#f = d^q
        #print("g")
        g = mat_multiply(bnd_op,i+n-q,n-q)# g = d^(n-q)
        #print(f.shape,g.shape)

        betti.append(betti_number(f, g))
    
        f = np.matrix(f)
        g = np.matrix(g)
        
        lap = f.H@f+g@g.H
        #print(i,lap.shape)
        eigvals = np.linalg.eigvalsh(lap)
        positive_eigvals = eigvals.real[eigvals.real>1e-6]
        gamma_min.append(np.min(positive_eigvals) if positive_eigvals.size>0 else 0)
        gamma_max.append(np.max(positive_eigvals) if positive_eigvals.size>0 else 0)
        gamma_mean.append(np.mean(positive_eigvals)if positive_eigvals.size>0 else 0)
        gamma_std.append(np.std(positive_eigvals) if positive_eigvals.size>0 else 0)

    return betti,gamma_min,gamma_max,gamma_mean,gamma_std

class boundary_operators():
    def __init__(self,complex,p):
        self.boundaries = {}
        bnd = boundary(complex,p)
        self.max_dim = len(bnd)
        for i in range(self.max_dim):
            self.boundaries[i] = bnd[i]
    def get(self,index):
        if index in self.boundaries:
            return self.boundaries[index]
        elif index == self.max_dim:
            return np.zeros((self.boundaries[self.max_dim-1].shape[1],0))
        else:
            return np.zeros((0,0))
    
    

def betti_laplacian(X,n,max_dim=1):
    Bettis=[]
    Gmin=[]
    Gmax=[]
    Gmean=[]
    Gstd=[]
    bnd_op = boundary_operators(X,n)
    for q in range(1,n):
        betti,gmin,gmax,gmean,gstd = calculate_betti_laplacian(bnd_op,n,q,max_dim=max_dim)
        Bettis.append(betti)
        Gmin.append(gmin)
        Gmax.append(gmax)
        Gmean.append(gmean)
        Gstd.append(gstd)
    return Bettis,Gmin,Gmax,Gmean,Gstd


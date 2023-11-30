
from Mayer_homology_Laplacian import betti_laplacian
from rips import rips_complex
import numpy as np
import time


def calculate(pointcloud,distances,n,name):
     B = {}
     Gmin = {}
     Gmax = {}
     Gmean = {}
     Gstd = {}
     pre_complex_num = 0
     start = time.time()
     rips = rips_complex(pointcloud,distances[-1],n)
     for i,distance in enumerate(distances):
         cur_complex = rips.get_complex(distance)
         if len(cur_complex) == pre_complex_num:
             B[distance] = B[distances[i-1]]
             Gmin[distance] = Gmin[distances[i-1]]
             Gmax[distance] = Gmax[distances[i-1]]
             Gmean[distance] = Gmean[distances[i-1]]
             Gstd[distance] = Gstd[distances[i-1]]
         else:
             #print(distance)
             B[distance],Gmin[distance],Gmax[distance],Gmean[distance],Gstd[distance] = betti_laplacian(cur_complex,n)
             #print("betti,lap,done",time.time()-start)
         pre_complex_num = len(cur_complex)
     for i in range(n-1):
         b0=[]
         b1=[]
         gmin0= []
         gmin1 = []
         gmax0 = []
         gmax1 = []
         gmean0=[]
         gmean1=[]
         gstd0=[]
         gstd1=[]
         for distance in distances:
             betti = B[distance][i]
             gmin = Gmin[distance][i]
             gmax = Gmax[distance][i]
             gmean = Gmean[distance][i]
             gstd = Gstd[distance][i]
             b0.append(betti[0] if len(betti)>0 else 0)
             b1.append(betti[1] if len(betti)>1 else 0)
             gmin0.append(gmin[0] if len(gmin)>0 else 0)
             gmin1.append(gmin[1] if len(gmin)>1 else 0)
             gmax0.append(gmax[0] if len(gmax)>0 else 0)
             gmax1.append(gmax[1] if len(gmax)>1 else 0)
             gmean0.append(gmean[0] if len(gmean)>0 else 0)
             gmean1.append(gmean[1] if len(gmean)>1 else 0)
             gstd0.append(gstd[0] if len(gstd)>0 else 0)
             gstd1.append(gstd[1] if len(gstd)>1 else 0)
         pname = './result/'+name
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_b0.npy',np.array(b0))
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_b1.npy',np.array(b1))
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_gmin0.npy',np.array(gmin0)) 
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_gmin1.npy',np.array(gmin1))
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_gmax0.npy',np.array(gmax0))
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_gmax1.npy',np.array(gmax1))
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_gmean0.npy',np.array(gmean0))
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_gmean1.npy',np.array(gmean1))
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_gstd0.npy',np.array(gstd0))
         np.save(pname+"_n="+str(n)+'_q='+str(i+1)+'_gstd1.npy',np.array(gstd1))
     return None


                



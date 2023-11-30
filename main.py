# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:52:03 2023

@author: lshen
"""
import numpy as np
from persistence import calculate
from functions import read_file,plot_betti_lap,plot_betti_lap_channels

##chose your parameters
file = 'X1.xy' # the 'xy' file, 'xyz' file or 'mol2' file you want to study,including 'X1.xy','X2.xyz','C20.xyz','C60.xyz','CB7.mol2'.
radius = np.linspace(0,2,201)# filtration radius.
n = 7 # N value for Mayer homology/Laplacian.
####
def process(file,radius,n):

    filename = file.split('.')[0]
    fileformat= file.split('.')[1]
    filepath = './data/'+file
    pointcloud = read_file(filepath)
    distances = radius*2

    
    calculate(pointcloud,distances,n,filename)
    fig1,axes1 = plot_betti_lap(filename,radius,n)
    fig2,axes2 =plot_betti_lap_channels(filename,radius,n)
    
    axes1[0,0].set_ylabel('$\\beta_0$',rotation=0)
    axes1[0,1].set_ylabel('$\\beta_1$',rotation=0)
    axes1[1,0].set_ylabel('$\\lambda_0$',rotation=0)
    axes1[1,1].set_ylabel('$\\lambda_1$',rotation=0)
    
    axes2[0,0].set_ylabel('$\\beta_0$',rotation=0)
    axes2[0,1].set_ylabel('$\\beta_1$',rotation=0)
    axes2[1,0].set_ylabel('$\\lambda_0$',rotation=0)
    axes2[1,1].set_ylabel('$\\lambda_1$',rotation=0)
    
    return fig1,axes1,fig2,axes2



fig1,axes1,fig2,axes2 = process(file,radius,n)


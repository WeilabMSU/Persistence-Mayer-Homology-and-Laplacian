# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:52:03 2023

@author: lshen
"""
import numpy as np
from persistence import calculate
from functions import read_file,plot_betti_lap,plot_betti_lap_channels,process

##chose your parameters
file = 'X1.xy' # the 'xy' file, 'xyz' file or 'mol2' file you want to study,including 'X1.xy','X2.xyz','C20.xyz','C60.xyz','CB7.mol2'.
radius = np.linspace(0,2,201)# filtration radius.
n = 7 # N value for Mayer homology/Laplacian.
##process the plot the results
fig1,axes1,fig2,axes2 = process(file,radius,n)


# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 15:30:32 2023

@author: lshen
"""

import gudhi
import numpy as np
from biopandas.mol2 import PandasMol2

      
class rips_complex():
    def __init__(self,pointcloud,max_distance,max_dim):
        self.simplex_list=[]
        self.filtration_list=[]
        comp = gudhi.RipsComplex(points=pointcloud,
                                     max_edge_length=max_distance)
        simplex_tree = comp.create_simplex_tree(max_dimension=max_dim)
        for filtered_value in simplex_tree.get_filtration():
            self.simplex_list.append(filtered_value[0])
            self.filtration_list.append(filtered_value[1])
        
    
    def get_complex(self,filtration_distance):
        fil_list = np.array(self.filtration_list)
        index = np.array(fil_list)<=filtration_distance
        simplexes = np.array(self.simplex_list,dtype=object)[index]
        return simplexes.tolist()

# Persistence-Mayer-Homology-and-Laplacian
The code for persistence Mayer homology and laplacian for N-chain complexes
---

## Data

The source data used in this work is provided.

https://github.com/WeilabMSU/PathHom/tree/main/data

### Requirments

Python Dependencies
  - python      (>=3.7)
  - numpy       (1.23.5)
  - biopandas   (0.4.1)
  - matplotlib  (3.6.2)
  - gudhi       (3.8.0)
    
### Tutorial

1. Basic examples, for calculating Mayer Bettis and Mayer Laplacians:
   
```
from Mayer_homology_Laplacian import betti_laplacian
##define your simplicial complex
X = [[0], [1], [2], [3],
[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3],
[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
##define your N value and the maximal dimensions of homology and laplacian
N = 3
max_dim = 2 #default to be 1
#####
Betti,Gmin,_,_,_ = betti_laplacian(X,N,max_dim=max_dim)

str_betti = "The Mayer Betti {} at stage {} is:"
str_lap = "The smallest positive eigenvalue of Mayer Laplacian {} at stage {} is:"

for q in range(1,N):
    for i,(b,l) in enumerate(zip(Betti[q-1],Gmin[q-1])):
        print(str_betti.format(i,q+1),b)
        print(str_betti.format(i,q+1),l)
```
2. Persistent Mayer homology and Laplacian:
   the following example is the same as the 'main.py' file
"""
import numpy as np
from persistence import calculate
from functions import read_file,plot_betti_lap,plot_betti_lap_channels,process

##chose your example
file = 'X1.xy' # the 'xy' file, 'xyz' file or 'mol2' file you want to study,including 'X1.xy','X2.xyz','C20.xyz','C60.xyz','CB7.mol2'.
radius = np.linspace(0,2,201)# filtration radius.
n = 7 # N value for Mayer homology/Laplacian.
##
fig1,axes1,fig2,axes2 = process(file,radius,n)

### Citing
 To be added...
 
### Contributors
The persistence Mayer homology and Laplacian code was created by [Li Shen](https://github.com/shenli0920)) and is maintained by [Li Shen](https://github.com/shenli0920)), Jian Liu and [Weilab](https://github.com/msuweilab) at MSU Math.

### License
All codes released in this study is under the MIT License.

---

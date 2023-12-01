def read_file(filename):
    file_format = filename.split('.')[-1]
    pointcloud=[]
    if file_format == 'mol2':
        from biopandas.mol2 import PandasMol2
        data =PandasMol2().read_mol2(filename).df
        for x,y,z in zip(data["x"],data["y"],data["z"]):
              pointcloud.append(np.array([x,y,z]))
    if file_format == 'xyz':
        with open(filename,'r') as f:
            for item in f.readlines():
                x,y,z = item.split()
                pointcloud.append(np.array([float(x),float(y),float(z)]))
    if file_format == 'xy':
        with open(filename,'r') as f:
            for item in f.readlines():
                x,y = item.split()
                pointcloud.append(np.array([float(x),float(y)]))
    
    return pointcloud

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, BoundaryNorm
def plot(filepath,radius,axe):
    data = np.load(filepath)
    points_x =[]
    points_y=[]
    for i,x in enumerate(radius):
        points_x.append(x)
        points_y.append(data[i])
        if i <len(radius)-1:
            points_x.append(radius[i+1])
            points_y.append(data[i])
            
    axe.plot(points_x,points_y)
    #axe.set_xticks(radius)

def plot_channels(name,radius,axe,n):
    datas=[]
    for q in range(1,n):
        filepath = name.format(q)
        data = np.load(filepath)
        datas.append(data)
    X = np.concatenate([data.reshape(-1,1) for data in datas],axis=1)
    #print(X)
    percentiles = np.percentile(np.unique(X), np.arange(10, 100, 10)).round(2)
    breakpoints = [0] + list(percentiles) + [np.max(X)]
    #print(breakpoints)
    # if name.split('.')[1][-2] == 'b':
    #     breakpoints = [int(value) for value in breakpoints]
    #     breakpoints = np.unique(np.array(breakpoints)).tolist()
    norm = BoundaryNorm(breakpoints, len(breakpoints), clip=True)
    cmap = plt.get_cmap('viridis', len(breakpoints) - 1)
    pos = axe.imshow(X.T,aspect='auto',interpolation='None',cmap=cmap,norm=norm)
    
    axe.set_yticks([])
    return pos


def plot_betti_lap(filename,radius,n):
    fig,axes = plt.subplots(nrows=2,ncols=2)
    plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)
    name = './result/'+filename+'_n={}_q={}_{}.npy'
    for q in range(1,n):
        plot(name.format(n,q,'b0'),radius,axes[0,0])
        plot(name.format(n,q,'b1'),radius,axes[0,1])
        plot(name.format(n,q,'gmin0'),radius,axes[1,0])
        plot(name.format(n,q,'gmin1'),radius,axes[1,1])
    return fig,axes
        
def plot_betti_lap_channels(filename,radius,n):
    fig,axes = plt.subplots(nrows=2,ncols=2)
    plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)
    name = './result/'+filename+'_n={}_q={}_{}.npy'
    pos1 = plot_channels(name.format(n,'{}','b0'),radius,axes[0,0],n)
    pos2 = plot_channels(name.format(n,'{}','b1'),radius,axes[0,1],n)
    pos3 = plot_channels(name.format(n,'{}','gmin0'),radius,axes[1,0],n)
    pos4 = plot_channels(name.format(n,'{}','gmin1'),radius,axes[1,1],n)   
    fig.colorbar(pos1,ax=axes[0,0])
    fig.colorbar(pos2,ax=axes[0,1])
    fig.colorbar(pos3,ax=axes[1,0])
    fig.colorbar(pos4,ax=axes[1,1])
    return fig,axes
from persistence import calculate
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

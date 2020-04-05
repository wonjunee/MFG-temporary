import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import csv
import dask.dataframe as dd
import pandas as pd

directory="./data"

#--------------------------------------------------
#   Getting n1 and nt
#--------------------------------------------------

with open("{}/parameters.csv".format(directory)) as F:
    csvReader = csv.reader(F)
    for i in csvReader:
        n1 = int(i[0])
        n2 = int(i[1])
        nt = int(i[2])

print(n1,n2,nt)
tau = 1.0/(nt)
#--------------------------------------------------
#   Getting Rho Data
#--------------------------------------------------

def open_csv(filename,nt,n1,n2):
    A=dd.read_csv(filename,header=None)
    return np.array(A).reshape((nt,n2,n1))


rho = open_csv("{}/rho.csv".format(directory),nt,n1,n2)
Phi = open_csv("{}/Phi.csv".format(directory),nt,n1,n2)



#--------------------------------------------------
#   Create animation
#--------------------------------------------------
px=[]
py=[]

num_particles=20

RHO_MAX = np.max(rho[0])
while len(px)<num_particles:
    # x = 0.25 + (np.random.random()-0.5)*0.2
    # y = 0.25 + (np.random.random()-0.5)*0.2

    x = np.random.random()
    y = np.random.random()

    argx=(int)(x*n1-0.5)
    argy=(int)(y*n2-0.5)

    # if(rho[0,argy,argx]/RHO_MAX*np.random.random()>0.5):
    if(rho[0,argy,argx]/RHO_MAX>0.6):
        add=True
        # for i in range(len(px)):
        #     dis = np.sqrt((px[i]-x)**2+(py[i]-y)**2)
        #     if dis<0.05:
        #         add=False
        #         break
        if add:
            px.append(x)
            py.append(y)



        
print("starting")
    # px.append(np.random.random())
    # py.append(np.random.random())

px=np.array(px)
py=np.array(py)
# First set up the figure, the axis, and the plot element we want to animate


# animation function.  This is called sequentially

gradxPhi=np.zeros((n2,n1))
gradyPhi=np.zeros((n2,n1))

import os
os.system("rm ./frames/*")

for n in range(nt):
    fig, ax = plt.subplots(1,1,figsize=(6,5))
    
    cax = ax.scatter(px,py,c='r',alpha=0.5,s=15)
    cax1= ax.imshow(np.flipud(rho[n]), zorder=0, extent=[0.0, 1.0, 0.00, 1.], aspect='auto')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    plt.suptitle("t={:.4}".format(n/(nt-1)))
    # ax.set_axis_off()

    # ax.set_facecolor((1.0, 0.47, 0.42))

    plt.xticks([])
    plt.yticks([])

    filename='frames/figure-{:0>3d}.png'.format(n)
    plt.savefig(filename, dpi=96)
    

    # plt.gca()
    plt.close()

    if(n<nt-1):
        argpx = (px*n1 - 0.5)
        argpy = (py*n2 - 0.5)

        argpx = argpx.astype(int)
        argpy = argpy.astype(int)

        gradxPhi[:,:-1] = (Phi[n+1,:,1:] - Phi[n+1,:,:-1])*(1.0*n1)
        gradyPhi[:-1,:] = (Phi[n+1,1:,:] - Phi[n+1,:-1,:])*(1.0*n2)

        px = px - tau * gradxPhi[argpy,argpx]
        py = py - tau * gradyPhi[argpy,argpx]


os.system("ffmpeg -r 30 -f image2 -s 1280x1280 -i ./frames/figure-%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test2.mp4")


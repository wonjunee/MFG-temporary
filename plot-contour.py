import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import animation
import numpy as np
import csv
import dask.dataframe as dd
import pandas as pd
from numpy.ma import masked_array


import matplotlib.patches as patches

def get_circle(center, radius):
    angle = np.linspace( 0 , 2 * np.pi , 150 )     
    x = radius * n1 * np.cos( angle ) + center[0] * n1 - 0.5
    y = radius * n2 * np.sin( angle ) + center[1] * n2 - 0.5
    return x,y


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

# var = 0.02
#--------------------------------------------------
#   Getting Rho Data
#--------------------------------------------------

def open_csv(filename,nt,n1,n2):
    A = np.fromfile(filename, dtype=np.float64)
    return A.reshape((nt,n2,n1))

rho = open_csv("{}/rho.csv".format(directory),nt,n1,n2)
print(rho.shape)

#--------------------------------------------------
#   Create animation
#--------------------------------------------------


import sys

type_video = sys.argv[1]

N = 6
tlist = [0,0.2,0.4,0.6,0.8,1.0]

def save_plot_contour(visual=True):
    fig, ax = plt.subplots(1,N,figsize=(12,3))

    fig.subplots_adjust(bottom=0, top=0.95, left=0, right=1, wspace=0, hspace=0.144)

    vmax = np.max(rho)
    
    for i in range(N):
        n = int((nt-1) * tlist[i]);
        ax[i].imshow(rho[n], cmap='inferno').set_clim(0, vmax)
        ax[i].set_axis_off()
        ax[i].set_title("t = {:.2f}\nsum = {:.3e}".format(1.0*n/(nt-1), np.sum(rho[n])/(n1*n2)))
        plt.savefig("figures/SIR-with-mass.png")
        
    if(visual==True):
        plt.show()
    plt.close();

save_plot_contour();

def save_animation():
    # First set up the figure, the axis, and the plot element we want to animate
    fig, ax = plt.subplots(1,1,figsize=(6,7))
    fig.subplots_adjust(bottom=0, top=0.8, right=1, left=0, hspace=0.1, wspace=0)

    cax = ax.imshow(rho[0], cmap='inferno', origin='lower')
    
    ax.set_axis_off()
    plt.tight_layout()

    vmin = np.min(rho)
    vmax = np.max(rho)
    
    # animation function.  This is called sequentially
    def animate(n):
        # fig.clear()
        cax.set_array(rho[n])
        
        # cax.set_clim(0, vmax)
        # cax.set_clim(vmin, vmax)
        cax.set_clim(np.min(rho[n]), np.max(rho[n]))
        
        

        ax.set_title("Video: {}\nMin: {:.4}\nMax: {:.4}".format(n,np.min(rho[n]),np.max(rho[n])))
        
        # cax0.set_clim(0, 10)
        # plt.suptitle("$\\beta$={:.2}, $\\gamma$={:.2}, $\\sigma$={:.3}\nt={:.2f}".format(beta,gamma,var,n/nt))
        return cax, 

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, 
                                   frames=nt, interval=100, blit=True)
    anim.save("video.mp4", fps=10)

save_animation()
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import csv
import dask.dataframe as dd
import pandas as pd
from numpy.ma import masked_array

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
        c0 = float(i[3])
        c1 = float(i[4])
        c2 = float(i[5])
        beta = float(i[6])
        gamma = float(i[7])
        var   = float(i[8])

# var = 0.02
#--------------------------------------------------
#   Getting Rho Data
#--------------------------------------------------

def open_csv(filename,nt,n1,n2):
    A = np.fromfile(filename, dtype=np.float64)

    # A=dd.read_csv(filename,header=None)
    return A.reshape((nt,n2,n1))

    



rho0 = open_csv("{}/rho0.csv".format(directory),nt,n1,n2)
rho1 = open_csv("{}/rho1.csv".format(directory),nt,n1,n2)
rho2 = open_csv("{}/rho2.csv".format(directory),nt,n1,n2)
rho3 = open_csv("{}/rho3.csv".format(directory),nt,n1,n2)
obstacle = open_csv("{}/obstacle.csv".format(directory),1,n1,n2)
mask_obstacle = masked_array(obstacle,obstacle==0)

#--------------------------------------------------
#   Create animation
#--------------------------------------------------


import sys

type_video = sys.argv[1]

N = 5
tlist = [0,0.25,0.5,0.75,1.0]

def save_plot_contour(visual=True, SIR=True, title_type=0):
    fig, ax = plt.subplots(4,N,figsize=(9,6))

    fig.subplots_adjust(bottom=0, top=0.9, right=1, left=0, wspace=0, hspace=0.2)

    vmax0 = np.max(rho0)
    print("vmax0 : ",vmax0)
    vmax0 = 0.8
    vmax1 = np.max(rho1)
    vmax2 = np.max(rho2)

    vmax = max(vmax0, vmax1, vmax2)

    for i in range(N):
        n = int((nt-1) * tlist[i]);
        ax[0,i].imshow(rho0[n], cmap='inferno').set_clim(0, vmax0)
        ax[1,i].imshow(rho1[n], cmap='inferno').set_clim(0, vmax0)
        ax[2,i].imshow(rho2[n], cmap='inferno').set_clim(0, vmax0)
        ax[3,i].imshow(rho3[n], cmap='inferno').set_clim(0, vmax0)


        ax[0,i].set_axis_off()
        ax[1,i].set_axis_off()
        ax[2,i].set_axis_off()
        ax[3,i].set_axis_off()

        if title_type == 0:
            ax[0,i].set_title("t = {:.2f}\nsum = {:.3f}".format(1.0*n/(nt-1), np.sum(rho0[n])/(n1*n2)))
            ax[1,i].set_title("sum = {:.3f}".format(np.sum(rho1[n])/(n1*n2)))
            ax[2,i].set_title("sum = {:.3f}".format(np.sum(rho2[n])/(n1*n2)))
            # ax[3,i].set_title("sum = {:.3f}".format(np.sum(rho3[n])/(n1*n2)))
            ax[3,i].set_title("max = {:.3f}".format(np.max(rho3[n])))
            plt.savefig("figures/SIR-with-mass.png")
        elif title_type == 1:
            ax[0,i].set_title("t = {:.2f}".format(1.0*n/(nt-1)))
            plt.savefig("figures/SIR-no-mass.png")


    
    if(visual==True):
        plt.show()
    plt.close();



def save_plot(visual=True, total=False):

    fig, ax = plt.subplots(1,2,figsize=(8,3))

    fig.subplots_adjust(bottom=0.1, top=0.95, right=1, left=0.1, wspace=0.2, hspace=0.2)

    xx = np.linspace(0,1,nt)
    max0 = np.sum(rho0[0])
    yy0 = np.sum(np.sum(rho0,axis=1),axis=1)/(n1*n2)
    yy1 = np.sum(np.sum(rho1,axis=1),axis=1)/(n1*n2)
    yy2 = np.sum(np.sum(rho2,axis=1),axis=1)/(n1*n2)
    yy3 = np.sum(np.sum(rho3,axis=1),axis=1)/(n1*n2)

    # ax.plot(xx,yy0,'.-',label="S")
    # ax.plot(xx,yy1,'.-',label="I")

    ax[0].plot(xx,yy0,label="S")
    ax[0].plot(xx,yy1,label="I")
    ax[0].plot(xx,yy2,label="R")
    ax[0].plot(xx,yy2+yy1+yy0,label="Total")
    ax[0].set_xlabel("Time")
    # ax.set_ylim(200,700)

    ax[1].plot(xx,yy3,label="V")

    V_0 = np.zeros_like(rho3)
    V_1 = np.zeros_like(rho3)

    X,Y = np.meshgrid(np.linspace(0.5/n1,1-0.5/n1,n2),np.linspace(0.5/n1,1-0.5/n1,n1))

    V_0[:,X-Y>0.5] = 1
    r = 0.075

    # ax[1].plot(xx,np.sum(np.sum(rho3[:,:,:n1//2],axis=1),axis=1)/(n1*n2),label="V left")
    # ax[1].plot(xx,np.sum(np.sum(rho3[:,:,n1//2:],axis=1),axis=1)/(n1*n2),label="V right")

    print(rho3[:,(X-0.3)**2+(Y-0.3)**2<r**2].shape)
    ax[1].plot(xx,np.sum(rho3[:,(X-0.5)**2+(Y-0.2)**2<r**2],axis=1)/(n1*n2),label="V bottom")
    ax[1].plot(xx,np.sum(rho3[:,(X-0.5)**2+(Y-0.5)**2<r**2],axis=1)/(n1*n2),label="V middle")
    ax[1].plot(xx,np.sum(rho3[:,(X-0.5)**2+(Y-0.8)**2<r**2],axis=1)/(n1*n2),label="V top")

    plt.legend(loc='upper right', framealpha=0.5)
    ax[0].grid()
    ax[1].grid()
    plt.savefig("figures/SIR-plot.eps")
    if(visual==True):
        plt.show()
    plt.close()

fig, ax = plt.subplots(1,N,figsize=(9,6))

fig.subplots_adjust(bottom=0, top=0.9, right=1, left=0, wspace=0, hspace=0.2)

vmax = np.max(rho0+rho1+rho2)

for i in range(N):
    n = int((nt-1) * tlist[i]);
    ax[i].imshow(rho0[n]+rho1[n]+rho2[n], cmap='inferno').set_clim(0, vmax)
    ax[i].set_axis_off()
    ax[i].set_title("{:.4f}".format(np.max(rho0[n]+rho1[n]+rho2[n])))
    
plt.show()
plt.close();

save_plot_contour(visual=True, SIR=True, title_type = 0)
save_plot(visual=True,  total=True)

def save_animation():
    # First set up the figure, the axis, and the plot element we want to animate
    num = 4; w = 16
    
    fig, ax = plt.subplots(1,num,figsize=(w,6))
    fig.subplots_adjust(bottom=0, top=0.8, right=1, left=0, hspace=0.1, wspace=0)

    cax0 = ax[0].imshow(rho0[0], cmap='inferno', origin='lower')
    cax1 = ax[1].imshow(rho1[0], cmap='inferno', origin='lower')
    cax2 = ax[2].imshow(rho2[0], cmap='inferno', origin='lower')
    cax3 = ax[3].imshow(rho3[0], cmap='inferno', origin='lower')

    ax[0].set_axis_off()
    ax[1].set_axis_off()
    ax[2].set_axis_off()
    ax[3].set_axis_off()
    plt.tight_layout()

    vmax0 = np.max(rho0)
    vmax1 = np.max(rho1)
    vmax2 = np.max(rho2)
    vmax3 = np.max(rho3)
    vmax  = max(vmax0, vmax1, vmax2)
    obstacle[obstacle > 0] = 1
    obs_sum = np.sum(obstacle) * 100 / (n1*n2)

    # animation function.  This is called sequentially
    def animate(n):
        # fig.clear()
        

        rho0[np.maximum(0,n-1)][obstacle[0] > 0] = 100
        rho1[np.maximum(0,n-1)][obstacle[0] > 0] = 100
        rho2[np.maximum(0,n-1)][obstacle[0] > 0] = 100
        rho3[np.maximum(0,n-1)][obstacle[0] > 0] = 100
        cax0.set_array(rho0[np.maximum(0,n-1)])
        cax1.set_array(rho1[np.maximum(0,n-1)])
        cax2.set_array(rho2[np.maximum(0,n-1)])
        cax3.set_array(rho3[np.maximum(0,n-1)])


        rho0sum = np.sum(rho0[np.maximum(0,n-1)])/(n1*n2) - obs_sum
        rho1sum = np.sum(rho1[np.maximum(0,n-1)])/(n1*n2) - obs_sum
        rho2sum = np.sum(rho2[np.maximum(0,n-1)])/(n1*n2) - obs_sum
        rho3sum = np.sum(rho3[np.maximum(0,n-1)])/(n1*n2) - obs_sum

        # cax0.set_clim(0, vmax)
        cax0.set_clim(0, vmax0)
        cax1.set_clim(0, vmax1)
        cax2.set_clim(0, vmax0)
        cax3.set_clim(0, vmax3*0.8)

        

        ax[0].set_title("Susceptible: {:.4f}".format(rho0sum))
        ax[1].set_title("Infected: {:.4f}".format(rho1sum))
        ax[2].set_title("Recovered: {:.4f}".format(rho2sum))
        ax[3].set_title("Vaccine: {:.4f}".format(rho3sum))

        # cax0.set_clim(0, 10)
        plt.suptitle("$\\beta$={:.2}, $\\gamma$={:.2}, $\\sigma$={:.3}\nt={:.2f}".format(beta,gamma,var,n/nt))
        return cax0, 

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, 
                                   frames=nt+1, interval=100, blit=True)
    anim.save("video.mp4", fps=10)

save_animation()



print(np.max(rho0[n]+rho1[n]+rho2[n]))
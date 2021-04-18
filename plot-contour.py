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

N = 6
tlist = [0,0.25,0.5,0.55,0.65,1.0]

def save_plot_contour(visual=True, SIR=True, title_type=0):
    fig, ax = plt.subplots(4,N,figsize=(8,6))

    fig.subplots_adjust(bottom=0, top=0.95, left=0, right=1, wspace=0, hspace=0.144)

    vmax0 = np.max(rho0)
    vmax1 = np.max(rho1)
    vmax2 = np.max(rho2*2)
    vmax3 = np.max(rho3*0.8)

    vmax = max(vmax0, vmax1, vmax2)

    for i in range(N):
        n = int((nt-1) * tlist[i]);
        ax[0,i].imshow(rho0[n]+obstacle[0]*100, cmap='inferno').set_clim(0, vmax0)
        ax[1,i].imshow(rho1[n]+obstacle[0]*100, cmap='inferno').set_clim(0, vmax1)
        ax[2,i].imshow(rho2[n]+obstacle[0]*100, cmap='inferno').set_clim(0, vmax2)
        ax[3,i].imshow(rho3[n]+obstacle[0]*100, cmap='inferno').set_clim(0, vmax3)

        # factory area
        x,y = get_circle((0.5,0.5), 0.075)
        ax[3,i].plot( x, y , 'tab:green') 
        x,y = get_circle((0.5,0.2), 0.075)
        ax[3,i].plot( x, y , 'tab:green') 
        x,y = get_circle((0.5,0.8), 0.075)
        ax[3,i].plot( x, y , 'tab:green') 


        ax[0,i].set_axis_off()
        ax[1,i].set_axis_off()
        ax[2,i].set_axis_off()
        ax[3,i].set_axis_off()

        if title_type == 0:
            ax[0,i].set_title("t = {:.2f}\nsum = {:.3e}".format(1.0*n/(nt-1), np.sum(rho0[n])/(n1*n2)))
            ax[1,i].set_title("sum = {:.2e}".format(np.sum(rho1[n])/(n1*n2)))
            ax[2,i].set_title("sum = {:.2e}".format(np.sum(rho2[n])/(n1*n2)))
            ax[3,i].set_title("sum = {:.2e}".format(np.sum(rho3[n])/(n1*n2)))
            ax[3,i].set_title("max = {:.3f}".format(np.max(rho3[n])))
            plt.savefig("figures/SIR-with-mass.png")
        elif title_type == 1:
            ax[0,i].set_title("t = {:.2f}".format(1.0*n/(nt-1)))
            plt.savefig("figures/SIR-no-mass.png")


    
    if(visual==True):
        plt.show()
    plt.close();



def save_plot_initial_densities():
    fig, ax = plt.subplots(1,2,figsize=(5,2.5))

    fig.subplots_adjust(bottom=0, top=0.9, right=1, left=0, wspace=0, hspace=0.2)

    vmax0 = np.max(rho0)
    vmax1 = np.max(rho1)
    vmax2 = np.max(rho2*10)
    vmax3 = np.max(rho3)

    vmax = max(vmax0, vmax1, vmax2)

    n = 0
    ax[0].imshow(rho0[n], cmap='inferno',origin='lower').set_clim(0, vmax0)
    ax[1].imshow(rho1[n], cmap='inferno',origin='lower').set_clim(0, vmax1)
    # ax[2].imshow(rho2[n], cmap='inferno').set_clim(0, vmax2)
    # ax[3].imshow(rho3[n], cmap='inferno').set_clim(0, vmax3)

    # factory area
    x,y = get_circle((0.5,0.5), 0.075)
    ax[0].plot( x, y , 'tab:green') 
    ax[1].plot( x, y , 'tab:green') 
    x,y = get_circle((0.5,0.2), 0.075)
    ax[0].plot( x, y , 'tab:green') 
    ax[1].plot( x, y , 'tab:green') 
    x,y = get_circle((0.5,0.8), 0.075)
    ax[0].plot( x, y , 'tab:green') 
    ax[1].plot( x, y , 'tab:green') 


    ax[0].set_axis_off()
    ax[1].set_axis_off()
    # ax[2].set_axis_off()
    # ax[3].set_axis_off()

    ax[0].set_title("S")
    ax[1].set_title("I")
    # ax[2].set_title("R")
    # ax[3].set_title("V")
    plt.savefig("figures/SIR-initial-densities.png")
    
    plt.show()
    plt.close();

save_plot_initial_densities()

def save_plot(visual=True, total=False):

    fig, ax = plt.subplots(1,2,figsize=(9,3))

    fig.subplots_adjust(bottom=0.1, top=0.9, right=0.95, left=0.1, wspace=0.2, hspace=0.2)

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

    # ax[1].plot(xx,yy3,label="V")

    V_0 = np.zeros_like(rho3)
    V_1 = np.zeros_like(rho3)

    X,Y = np.meshgrid(np.linspace(0.5/n1,1-0.5/n1,n2),np.linspace(0.5/n1,1-0.5/n1,n1))

    V_0[:,X-Y>0.5] = 1
    r = 0.1

    # ax[1].plot(xx,np.sum(np.sum(rho3[:,:,:n1//2],axis=1),axis=1)/(n1*n2),label="V left")
    # ax[1].plot(xx,np.sum(np.sum(rho3[:,:,n1//2:],axis=1),axis=1)/(n1*n2),label="V right")

    # ax[1].plot(xx,np.sum(rho3[:,(X-0.5)**2+(Y-0.2)**2<r**2],axis=1)/(n1*n2),label="V bottom")
    # ax[1].plot(xx,np.sum(rho3[:,(X-0.5)**2+(Y-0.5)**2<r**2],axis=1)/(n1*n2),label="V middle")
    # ax[1].plot(xx,np.sum(rho3[:,(X-0.5)**2+(Y-0.8)**2<r**2],axis=1)/(n1*n2),label="V top")

    ax[1].plot(xx[:nt//2],np.sum(rho3[:nt//2,(X-0.5)**2+(Y-0.5)**2<r**2],axis=1)/(n1*n2),label="Factory")

    plt.legend(loc='upper right', framealpha=0.5)
    ax[0].grid()
    ax[1].grid()
    plt.savefig("figures/SIR-plot.eps")
    if(visual==True):
        plt.show()
    plt.close()



def save_plot_plots_for_rho_V_production_and_delivery(visual=True):

    fig, ax = plt.subplots(1,1,figsize=(5,3))

    fig.subplots_adjust(bottom=0.15, top=0.98, right=0.97, left=0.2, wspace=0.3, hspace=0.3)
    xx = np.linspace(0,1,nt)
    directory1 = "./exp2-data/data-no-obs"
    rho00 = open_csv("{}/rho0.csv".format(directory1),nt,n1,n2)
    rho10 = open_csv("{}/rho1.csv".format(directory1),nt,n1,n2)
    rho20 = open_csv("{}/rho2.csv".format(directory1),nt,n1,n2)
    rho30 = open_csv("{}/rho3.csv".format(directory1),nt,n1,n2)


    X,Y = np.meshgrid(np.linspace(0.5/n1,1-0.5/n1,n2),np.linspace(0.5/n1,1-0.5/n1,n1))
    r = 0.1

    ax.plot(xx[:nt//2],np.sum(rho30[:nt//2,(X-0.5)**2+(Y-0.5)**2<r**2],axis=1)/(n1*n2),label="Without obstacle")
    ax.plot(xx[:nt//2],np.sum(rho3[:nt//2,(X-0.5)**2+(Y-0.5)**2<r**2],axis=1)/(n1*n2),label="With obstacle")
    
    ax.legend(loc='upper left', framealpha=0.5)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.grid()
    ax.set_xlabel("Time $t$")

    ax.set_ylabel("Total mass of produced vaccines")

    plt.savefig("figures/exp2-comparison-plots1.eps")
    plt.show()

    fig, ax = plt.subplots(1,1,figsize=(5,3))

    fig.subplots_adjust(bottom=0.15, top=0.98, right=0.97, left=0.2, wspace=0.3, hspace=0.3)
    xx = np.linspace(0,1,nt)
    directory1 = "./exp2-data/data-no-obs"
    rho00 = open_csv("{}/rho0.csv".format(directory1),nt,n1,n2)
    rho10 = open_csv("{}/rho1.csv".format(directory1),nt,n1,n2)
    rho20 = open_csv("{}/rho2.csv".format(directory1),nt,n1,n2)
    rho30 = open_csv("{}/rho3.csv".format(directory1),nt,n1,n2)


    X,Y = np.meshgrid(np.linspace(0.5/n1,1-0.5/n1,n2),np.linspace(0.5/n1,1-0.5/n1,n1))
    r = 0.1

    ax.plot(xx[nt//2:],np.sum(np.sum(rho30[nt//2:,:,:n1//2-20],axis=1),axis=1)/(n1*n2),'--',label="Left without obstacle")
    ax.plot(xx[nt//2:],np.sum(np.sum(rho3[nt//2:,:,:n1//2-20],axis=1),axis=1)/(n1*n2),label="Left with obstacle")
    ax.plot(xx[nt//2:],np.sum(np.sum(rho30[nt//2:,:,n1//2+20:],axis=1),axis=1)/(n1*n2),'--',label="Right without obstacle")
    ax.plot(xx[nt//2:],np.sum(np.sum(rho3[nt//2:,:,n1//2+20:],axis=1),axis=1)/(n1*n2),label="Right with obstacle")
    
    ax.legend(loc='lower right', framealpha=0.5)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.grid()
    ax.set_xlabel("Time $t$")

    ax.set_ylabel("Total mass of delivered vaccines")

    plt.savefig("figures/exp2-comparison-plots2.eps")
    plt.show()
    plt.close()

save_plot_plots_for_rho_V_production_and_delivery()

def save_plot2(visual=True, total=False):

    directory0 = "data"
    rho01 = open_csv("{}/rho0.csv".format(directory0),nt,n1,n2)
    rho11 = open_csv("{}/rho1.csv".format(directory0),nt,n1,n2)
    rho21 = open_csv("{}/rho2.csv".format(directory0),nt,n1,n2)
    rho31 = open_csv("{}/rho3.csv".format(directory0),nt,n1,n2)

    directory1 = "data-large"
    rho00 = open_csv("{}/rho0.csv".format(directory1),nt,n1,n2)
    rho10 = open_csv("{}/rho1.csv".format(directory1),nt,n1,n2)
    rho20 = open_csv("{}/rho2.csv".format(directory1),nt,n1,n2)
    rho30 = open_csv("{}/rho3.csv".format(directory1),nt,n1,n2)

    fig, ax = plt.subplots(1,4,figsize=(15,3))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.08, right=0.99, wspace=0.3, hspace=0.2)

    xx = np.linspace(0,1,nt)
    max0 = np.sum(rho0[0])
    yy01 = np.sum(np.sum(rho01,axis=1),axis=1)/(n1*n2)
    yy11 = np.sum(np.sum(rho11,axis=1),axis=1)/(n1*n2)
    yy21 = np.sum(np.sum(rho21,axis=1),axis=1)/(n1*n2)
    yy31 = np.sum(np.sum(rho31,axis=1),axis=1)/(n1*n2)

    yy00 = np.sum(np.sum(rho00,axis=1),axis=1)/(n1*n2)
    yy10 = np.sum(np.sum(rho10,axis=1),axis=1)/(n1*n2)
    yy20 = np.sum(np.sum(rho20,axis=1),axis=1)/(n1*n2)
    yy30 = np.sum(np.sum(rho30,axis=1),axis=1)/(n1*n2)

    # ax.plot(xx,yy0,'.-',label="S")
    # ax.plot(xx,yy1,'.-',label="I")

    ax[0].plot(xx,yy01,label="Simulation 1")
    ax[0].plot(xx,yy00,label="Simulation 2")

    ax[1].plot(xx,yy11,label="Simulation 1")
    ax[1].plot(xx,yy10,label="Simulation 2")

    ax[2].plot(xx,yy21,label="Simulation 1")
    ax[2].plot(xx,yy20,label="Simulation 2")

    ax[0].set_ylabel("The Total Mass of Densities")
    for k in range(4):
        ax[k].set_xlabel("Time")
        
    ax[0].set_title("S")
    ax[1].set_title("I")
    ax[2].set_title("R")
    ax[3].set_title("V")
        

    V_0 = np.zeros_like(rho3)
    V_1 = np.zeros_like(rho3)

    X,Y = np.meshgrid(np.linspace(0.5/n1,1-0.5/n1,n2),np.linspace(0.5/n1,1-0.5/n1,n1))

    V_0[:,X-Y>0.5] = 1
    r = 0.1
    for k in range(4):
        ax[k].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

    # ax[1].plot(xx,np.sum(np.sum(rho3[:,:,:n1//2],axis=1),axis=1)/(n1*n2),label="V left")
    # ax[1].plot(xx,np.sum(np.sum(rho3[:,:,n1//2:],axis=1),axis=1)/(n1*n2),label="V right")

    # ax[1].plot(xx,np.sum(rho3[:,(X-0.5)**2+(Y-0.2)**2<r**2],axis=1)/(n1*n2),label="V bottom")
    # ax[1].plot(xx,np.sum(rho3[:,(X-0.5)**2+(Y-0.5)**2<r**2],axis=1)/(n1*n2),label="V middle")
    # ax[1].plot(xx,np.sum(rho3[:,(X-0.5)**2+(Y-0.8)**2<r**2],axis=1)/(n1*n2),label="V top")

    ax[3].plot(xx[:nt//2],np.sum(rho31[:nt//2,(X-0.5)**2+(Y-0.5)**2<r**2],axis=1)/(n1*n2),label="Simulation 1")
    ax[3].plot(xx[:nt//2],np.sum(rho30[:nt//2,(X-0.5)**2+(Y-0.5)**2<r**2],axis=1)/(n1*n2),label="Simulation 2")

    # plotting all
    # ax[3].plot(xx,yy31,label="Sim 1")
    # ax[3].plot(xx,yy30,label="Sim 2")

    ax[0].legend(loc='lower left', framealpha=0.5)
    ax[1].legend(loc='lower right', framealpha=0.5)
    ax[2].legend(loc='upper left', framealpha=0.5)
    ax[3].legend(loc='upper left', framealpha=0.5)
    for k in range(4):
        ax[k].grid()
    plt.savefig("figures/SIR-plot2.eps")
    plt.show()
    plt.close()

save_plot_contour(visual=True, SIR=True, title_type = 1)
save_plot(visual=True,  total=True)
# save_plot2(visual=True,  total=True)

def get_rect(x,y,w,h):
    return x*n1-0.5-w*n1, y*n2-0.5-h*n2, w*n1*2, h*n2*2
def save_animation():
    # First set up the figure, the axis, and the plot element we want to animate
    num = 4; w = 16
    
    fig, ax = plt.subplots(1,num,figsize=(w,6))
    fig.subplots_adjust(bottom=0, top=0.8, right=1, left=0, hspace=0.1, wspace=0)

    cax0 = ax[0].imshow(rho0[0], cmap='inferno', origin='lower')
    cax1 = ax[1].imshow(rho1[0], cmap='inferno', origin='lower')
    cax2 = ax[2].imshow(rho2[0], cmap='inferno', origin='lower')
    cax3 = ax[3].imshow(rho3[0], cmap='inferno', origin='lower')

    angle = np.linspace( 0 , 2 * np.pi , 150 ) 
 
    # factory area
    x,y = get_circle((0.5,0.5), 0.075)
    ax[3].plot( x, y , 'tab:green') 
    x,y = get_circle((0.5,0.2), 0.075)
    ax[3].plot( x, y , 'tab:green') 
    x,y = get_circle((0.5,0.8), 0.075)
    ax[3].plot( x, y , 'tab:green')

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
        
        rho3max = np.max(rho3[np.maximum(0,n-1)])

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
        cax2.set_clim(0, vmax2)
        cax3.set_clim(0, vmax3*0.99)

        

        ax[0].set_title("Susceptible: {:.4f}".format(rho0sum))
        ax[1].set_title("Infected: {:.4f}".format(rho1sum))
        ax[2].set_title("Recovered: {:.4f}".format(rho2sum))
        # ax[3].set_title("Vaccine: {:.4f}".format(rho3sum))
        ax[3].set_title("Vaccine: {:.4f}".format(rho3max))

        # cax0.set_clim(0, 10)
        plt.suptitle("$\\beta$={:.2}, $\\gamma$={:.2}, $\\sigma$={:.3}\nt={:.2f}".format(beta,gamma,var,n/nt))
        return cax0, 

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, 
                                   frames=nt+1, interval=100, blit=True)
    anim.save("video.mp4", fps=10)

save_animation()
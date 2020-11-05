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
        c0 = float(i[3])
        c1 = float(i[4])
        c2 = float(i[5])
        beta = float(i[6])
        gamma = float(i[7])
        var   = float(i[8])

var = 0.02
#--------------------------------------------------
#   Getting Rho Data
#--------------------------------------------------

def open_csv(filename,nt,n1,n2):
    A=dd.read_csv(filename,header=None)
    return np.array(A).reshape((nt,n2,n1))


rho0 = open_csv("{}/rho0.csv".format(directory),nt,n1,n2)
rho1 = open_csv("{}/rho1.csv".format(directory),nt,n1,n2)
rho2 = open_csv("{}/rho2.csv".format(directory),nt,n1,n2)

directory_2 = "data-2"
rho0_2 = open_csv("{}/rho0.csv".format(directory_2),nt,n1,n2)
rho1_2 = open_csv("{}/rho1.csv".format(directory_2),nt,n1,n2)
rho2_2 = open_csv("{}/rho2.csv".format(directory_2),nt,n1,n2)



#--------------------------------------------------
#   Create animation
#--------------------------------------------------


import sys

type_video = sys.argv[1]

N = 5
tlist = [0,0.25,0.5,0.75,1.0]

def save_plot_contour(visual=True, SIR=True, title_type=0):
    fig, ax = plt.subplots(3,N,figsize=(9,6))

    fig.subplots_adjust(bottom=0, top=0.9, right=1, left=0, wspace=0, hspace=0.2)

    vmax0 = np.max(rho0)
    vmax1 = np.max(rho1)
    vmax2 = np.max(rho2)

    vmax = max(vmax0, vmax1, vmax2)

    for i in range(N):
        n = int((nt-1) * tlist[i]);
        ax[0,i].imshow(rho0[n], cmap='inferno').set_clim(0, vmax)
        ax[1,i].imshow(rho1[n], cmap='inferno').set_clim(0, vmax)
        ax[2,i].imshow(rho2[n], cmap='inferno').set_clim(0, vmax)


        ax[0,i].set_axis_off()
        ax[1,i].set_axis_off()
        ax[2,i].set_axis_off()

        if title_type == 0:
            ax[0,i].set_title("t = {:.2f}\nsum = {:.3f}".format(1.0*n/(nt-1), np.sum(rho0[n])/(n1*n2)))
            ax[1,i].set_title("sum = {:.3f}".format(np.sum(rho1[n])/(n1*n2)))
            ax[2,i].set_title("sum = {:.3f}".format(np.sum(rho2[n])/(n1*n2)))
            plt.savefig("figures/SIR-with-mass.png")
        elif title_type == 1:
            ax[0,i].set_title("t = {:.2f}".format(1.0*n/(nt-1)))
            plt.savefig("figures/SIR-no-mass.png")


    
    if(visual==True):
        plt.show()
    plt.close();



def save_plot(visual=True, total=False):

    fig, ax = plt.subplots(1,1,figsize=(5,4))

    fig.subplots_adjust(bottom=0.1, top=0.95, right=1, left=0.1, wspace=0, hspace=0.2)

    xx = np.linspace(0,1,nt)
    max0 = np.sum(rho0[0])
    yy0 = np.sum(np.sum(rho0,axis=1),axis=1)/(n1*n2)
    yy1 = np.sum(np.sum(rho1,axis=1),axis=1)/(n1*n2)
    yy2 = np.sum(np.sum(rho2,axis=1),axis=1)/(n1*n2)

    # ax.plot(xx,yy0,'.-',label="S")
    # ax.plot(xx,yy1,'.-',label="I")

    ax.plot(xx,yy0,label="S")
    ax.plot(xx,yy1,label="I")
    ax.plot(xx,yy2,label="R")

    if(total==True):
        ax.plot(xx,yy2+yy1+yy0,label="Total")

    ax.set_xlabel("t")
    # ax.set_ylim(200,700)

    plt.legend(loc='lower left', framealpha=0.5)
    plt.savefig("figures/SIR-plot.eps")
    if(visual==True):
        plt.show()
    plt.close()

def save_plot_2(visual=True, total=False):

    fig, ax = plt.subplots(1,1,figsize=(5,4))

    fig.subplots_adjust(bottom=0.1, top=0.95, right=1, left=0.1, wspace=0, hspace=0.2)

    xx = np.linspace(0,1,nt)
    max0 = np.sum(rho0[0])
    yy0 = np.sum(np.sum(rho0,axis=1),axis=1)/(n1*n2)
    yy1 = np.sum(np.sum(rho1,axis=1),axis=1)/(n1*n2)
    yy2 = np.sum(np.sum(rho2,axis=1),axis=1)/(n1*n2)

    yy0_2 = np.sum(np.sum(rho0_2,axis=1),axis=1)/(n1*n2)
    yy1_2 = np.sum(np.sum(rho1_2,axis=1),axis=1)/(n1*n2)
    yy2_2 = np.sum(np.sum(rho2_2,axis=1),axis=1)/(n1*n2)

    # ax.plot(xx,yy0,'.-',label="S")
    # ax.plot(xx,yy1,'.-',label="I")

    ax.plot(xx,yy0,'.',label="S with velocity")
    ax.plot(xx,yy1,'.',label="I with velocity")
    ax.plot(xx,yy2,'.',label="R with velocity")

    ax.plot(xx,yy0_2,label="S without velocity")
    ax.plot(xx,yy1_2,label="I without velocity")
    ax.plot(xx,yy2_2,label="R without velocity")

    if(total==True):
        ax.plot(xx,yy2+yy1+yy0,label="Total")

    ax.set_xlabel("t")
    # ax.set_ylim(200,700)

    plt.legend(loc='best', bbox_to_anchor=(0.5, 0.5))
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
save_plot_contour(visual=False, SIR=True, title_type = 1)
save_plot(visual=True,  total=True)
save_plot(visual=False,  total=False)
save_plot_2(visual=False,  total=False)


def save_animation(SIR=True):
    # First set up the figure, the axis, and the plot element we want to animate
    if SIR==True:
        num = 3; w = 12
    else:
        num = 2; w = 8

    fig, ax = plt.subplots(1,num,figsize=(w,4))
    fig.subplots_adjust(bottom=0, top=0.8, right=1, left=0, hspace=0.1, wspace=0)

    cax1 = ax[0].imshow(rho0[0], cmap='inferno')
    cax2 = ax[1].imshow(rho1[0], cmap='inferno')

    ax[0].set_axis_off()
    ax[1].set_axis_off()

    if(SIR==True):
        cax3 = ax[2].imshow(rho2[0], cmap='inferno')
        ax[2].set_axis_off()
    
    vmax = np.max(rho1)

    # animation function.  This is called sequentially
    def animate(n):
        # fig.clear()
        cax1.set_array(np.flipud(rho0[n]))
        cax2.set_array(np.flipud(rho1[n]))
        
                
        vmax0 = np.max(rho0)
        vmax1 = np.max(rho1)
        vmax2 = np.max(rho2)

        vmax = max(vmax0, vmax1, vmax2)

        # cax1.set_clim(0, vmax)
        cax1.set_clim(0, vmax)
        cax2.set_clim(0, vmax)
    
        max0 = n1*n2

        ax[0].set_title("{:.4f}".format(np.sum(rho0[n])/max0))
        ax[1].set_title("{:.4f}".format(np.sum(rho1[n])/max0))
        

        if(SIR==True):
            max1 = 2*np.max(rho2)
            cax3.set_array(np.flipud(rho2[n]))
            cax3.set_clim(0, vmax)
            ax[2].set_title("{:.4f}".format(np.sum(rho2[n])/max0))


        ax[0].set_title("{:.4f}\n{:.4f}".format(np.sum(rho0[n])/max0,np.max(rho0[n])))
        ax[1].set_title("{:.4f}\n{:.4f}".format(np.sum(rho1[n])/max0,np.max(rho1[n])))
        ax[2].set_title("{:.4f}\n{:.4f}".format(np.sum(rho2[n])/max0,np.max(rho2[n])))

        # cax1.set_clim(0, 10)
        plt.suptitle("c0={:.2}, c1={:.2}, c2={:.2}, beta={:.2}, gamma={:.2}, var={:.3}\nt={:.2f},{:.2f}".format(c0,c1,c2,beta,gamma,var,n/(nt-1),np.max(rho0[n]+rho1[n]+rho2[n])))
        return cax1, 

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, 
                                   frames=nt, interval=10, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    anim.save("video.mp4", fps=30)


def save_animation2():
    # First set up the figure, the axis, and the plot element we want to animate
    num = 3; w = 12
    
    fig, ax = plt.subplots(1,num,figsize=(w,4))
    fig.subplots_adjust(bottom=0, top=0.8, right=1, left=0, hspace=0.1, wspace=0)

    cax1 = ax[0].imshow(rho0[0], cmap='inferno')
    cax2 = ax[1].imshow(rho1[0], cmap='inferno')

    ax[0].set_axis_off()
    ax[1].set_axis_off()

    cax3 = ax[2].imshow(rho2[0], cmap='inferno')
    ax[2].set_axis_off()

    vmax = np.max(rho1)

    # animation function.  This is called sequentially
    def animate(n):
        # fig.clear()
        cax1.set_array(np.flipud(rho0[n]))
        cax2.set_array(np.flipud(rho1[n]))
        
        vmax0 = np.max(rho0)
        vmax1 = np.max(rho1)
        vmax2 = np.max(rho2)

        vmax  = max(vmax0, vmax1, vmax2)

        # cax1.set_clim(0, vmax)
        cax1.set_clim(0, vmax)
        cax2.set_clim(0, vmax)

        cax3.set_array(np.flipud(rho2[n]))
        cax3.set_clim(0, vmax)

        ax[0].set_title("S: {:.4f}".format(np.sum(rho0[n])/(n1*n2)))
        ax[1].set_title("I: {:.4f}".format(np.sum(rho1[n])/(n1*n2)))
        ax[2].set_title("R: {:.4f}".format(np.sum(rho2[n])/(n1*n2)))

        # cax1.set_clim(0, 10)
        plt.suptitle("$\\beta$={:.2}, $\\gamma$={:.2}, $\\sigma$={:.3}\nt={:.2f}".format(beta,gamma,var,n/(nt-1)))
        return cax1, 

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, 
                                   frames=nt, interval=10, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    anim.save("video.mp4", fps=10)



if type_video=="0":
    save_animation(SIR=True)
elif type_video=="1":
    save_animation2()
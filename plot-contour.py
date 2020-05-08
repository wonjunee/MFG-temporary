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

#--------------------------------------------------
#   Getting Rho Data
#--------------------------------------------------

def open_csv(filename,nt,n1,n2):
    A=dd.read_csv(filename,header=None)
    return np.array(A).reshape((nt,n2,n1))


rho0 = open_csv("{}/rho0.csv".format(directory),nt,n1,n2)
rho1 = open_csv("{}/rho1.csv".format(directory),nt,n1,n2)
rho2 = open_csv("{}/rho2.csv".format(directory),nt,n1,n2)



#--------------------------------------------------
#   Create animation
#--------------------------------------------------


import sys

type_video = sys.argv[1]

N = 5
tlist = [0,0.25,0.5,0.75,1.0]



# fig, ax = plt.subplots(2,N,figsize=(9,4))

# fig.subplots_adjust(bottom=0, top=0.9, right=1, left=0, wspace=0, hspace=0.2)

# for i in range(N):
#     n = int((nt-1) * tlist[i]);
#     ax[0,i].imshow(rho0[n], cmap='inferno').set_clim(0, np.max(rho0))
#     ax[1,i].imshow(rho1[n], cmap='inferno').set_clim(0, np.max(rho1))


#     ax[0,i].set_axis_off()
#     ax[1,i].set_axis_off()

#     ax[0,i].set_title("t = {:.2f}\nsum = {:.2f}".format(1.0*n/(nt-1), np.sum(rho0[n])/(n1*n2)))
#     ax[1,i].set_title("sum = {:.2f}".format(np.sum(rho1[n])/(n1*n2)))

# plt.savefig("figures/SI.png")
# plt.show()


# fig, ax = plt.subplots(1,1,figsize=(4,3))

# fig.subplots_adjust(bottom=0.1, top=1.0, right=1, left=0.1, wspace=0, hspace=0.2)

# xx = np.linspace(0,1,nt)
# yy0 = np.sum(np.sum(rho0,axis=1),axis=1)/(n1*n2)
# yy1 = np.sum(np.sum(rho1,axis=1),axis=1)/(n1*n2)

# # ax.plot(xx,yy0,'.-',label="S")
# # ax.plot(xx,yy1,'.-',label="I")

# ax.plot(xx,yy0,label="S")
# ax.plot(xx,yy1,label="I")

# ax.set_xlabel("t")
# # ax.set_ylim(200,700)

# plt.legend()
# plt.savefig("figures/SI-plot.png")
# plt.show()





fig, ax = plt.subplots(3,N,figsize=(9,6))

fig.subplots_adjust(bottom=0, top=0.9, right=1, left=0, wspace=0, hspace=0.2)

for i in range(N):
    n = int((nt-1) * tlist[i]);
    ax[0,i].imshow(rho0[n], cmap='inferno').set_clim(0, np.max(rho0))
    ax[1,i].imshow(rho1[n], cmap='inferno').set_clim(0, np.max(rho1))
    ax[2,i].imshow(rho2[n], cmap='inferno').set_clim(0, np.max(rho2)*2)


    ax[0,i].set_axis_off()
    ax[1,i].set_axis_off()
    ax[2,i].set_axis_off()

    ax[0,i].set_title("t = {:.2f}\nsum = {:.2f}".format(1.0*n/(nt-1), np.sum(rho0[n])/(n1*n2)))
    ax[1,i].set_title("sum = {:.2f}".format(np.sum(rho1[n])/(n1*n2)))
    ax[2,i].set_title("sum = {:.2f}".format(np.sum(rho2[n])/(n1*n2)))

plt.savefig("figures/SIR.png")
plt.show()


fig, ax = plt.subplots(1,1,figsize=(5,3))

fig.subplots_adjust(bottom=0.1, top=1.0, right=1, left=0.1, wspace=0, hspace=0.2)

xx = np.linspace(0,1,nt)
yy0 = np.sum(np.sum(rho0,axis=1),axis=1)/(n1*n2)
yy1 = np.sum(np.sum(rho1,axis=1),axis=1)/(n1*n2)
yy2 = np.sum(np.sum(rho2,axis=1),axis=1)/(n1*n2)

# ax.plot(xx,yy0,'.-',label="S")
# ax.plot(xx,yy1,'.-',label="I")

ax.plot(xx,yy0,label="S")
ax.plot(xx,yy1,label="I")
ax.plot(xx,yy2,label="R")

ax.set_xlabel("t")
# ax.set_ylim(200,700)

plt.legend()
plt.savefig("figures/SIR-plot.png")
plt.show()




if type_video=="0":

    # First set up the figure, the axis, and the plot element we want to animate
    fig, ax = plt.subplots(1,3,figsize=(12,4.5))
    fig.subplots_adjust(bottom=0, top=0.85, right=1, left=0)

    cax1 = ax[0].imshow(rho0[0], cmap='inferno')
    cax2 = ax[1].imshow(rho1[0], cmap='inferno')
    cax3 = ax[2].imshow(rho2[0], cmap='inferno')
    ax[0].set_axis_off()
    ax[1].set_axis_off()
    ax[2].set_axis_off()

    # animation function.  This is called sequentially
    def animate(n):
        # fig.clear()
        cax1.set_array(np.flipud(rho0[n]))
        cax2.set_array(np.flipud(rho1[n]))
        cax3.set_array(np.flipud(rho2[n]))
        cax1.set_clim(0, np.max(rho0))
        cax2.set_clim(0, np.max(rho1))
        cax3.set_clim(0, 2*np.max(rho2))

        ax[0].set_title("{:.4f}".format(np.sum(rho0[n])))
        ax[1].set_title("{:.4f}".format(np.sum(rho1[n])))
        ax[2].set_title("{:.4f}".format(np.sum(rho2[n])))
        # cax1.set_clim(0, 10)
        plt.suptitle("t={:.4}\n{}".format(n/(nt-1),np.sum(rho0[n]+rho1[n]+rho2[n])))
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

elif type_video=="1":

    # First set up the figure, the axis, and the plot element we want to animate
    fig, ax = plt.subplots(1,1,figsize=(8,6))

    cax = ax.imshow(rho[0], cmap='inferno')
    fig.colorbar(cax)
    plt.axis('off')

    # animation function.  This is called sequentially
    def animate(n):
        # fig.clear()
        cax.set_array(np.flipud(rho[n]))
        cax.set_clim(0, np.max(rho[n]))
        # cax1.set_clim(0, 10)
        plt.suptitle("t={:.4}\n{}".format(n/(nt-1),np.sum(rho[n])))
        return cax, 

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, 
                                   frames=nt, interval=10, blit=True)

    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    anim.save("video.mp4", fps=30)
import matplotlib.pyplot as plt
import numpy as np
import dask.dataframe as dd
import pandas as pd
import sys
import csv

directory = "data"

def open_csv(filename,n1,n2):
    A=dd.read_csv(filename,header=None)
    return np.array(A).reshape((n2,n1))

with open("{}/parameters.csv".format(directory)) as F:
    csvReader = csv.reader(F)
    for i in csvReader:
        n1 = int(i[0])
        n2 = int(i[1])
        nt = int(i[2])

if len(sys.argv)==2:
    nt = int(sys.argv[1])

# #--------------------------------------------------
# #   Getting Rho Data
# #--------------------------------------------------

def open_and_reshape(filename, n1, n2):
    X = np.fromfile(filename, dtype=np.float64)
    X.shape = (n1,n2)
    return X

print(n1,n2,nt)

#--------------------------------------------------
#   Create animation
#--------------------------------------------------
from matplotlib import animation

x = np.linspace(0, 1, n1)
y = np.linspace(0, 1, n2)
x, y = np.meshgrid(x, y)

# First set up the figure, the axis, and the plot element we want to animate
rho1=open_and_reshape("./data/rho1-{}.csv".format(0),n1,n2)
rho2=open_and_reshape("./data/rho2-{}.csv".format(0),n1,n2)
rho3=open_and_reshape("./data/rho3-{}.csv".format(0),n1,n2)
fig, ax = plt.subplots(1,3,figsize=(12,4))
cax1 = ax[0].imshow(rho1, cmap='inferno')
cax2 = ax[1].imshow(rho2, cmap='inferno')
cax3 = ax[2].imshow(rho3, cmap='inferno')

# fig.colorbar(cax2)
ax[0].set_axis_off()
ax[1].set_axis_off()
ax[2].set_axis_off()

# max_rho = np.max(rho)
max_rho1 = np.max(rho1)
max_rho2 = np.max(rho2)
max_rho3 = np.max(rho3)
# animation function.  This is called sequentially
def animate(n):

    print("\rProcessing frame %d/%d..." % (n, nt), end='')

    # fig.clear()
    rho1=open_and_reshape("./data/rho1-{}.csv".format(n),n1,n2)
    rho2=open_and_reshape("./data/rho2-{}.csv".format(n),n1,n2)
    rho3=open_and_reshape("./data/rho3-{}.csv".format(n),n1,n2)
    # cax.set_array(np.flipud(rho))
    cax1.set_array(rho1)
    cax1.set_clim(np.min(0), np.max(rho1))
    cax2.set_array(rho2)
    cax2.set_clim(np.min(0), np.max(rho2))
    cax3.set_array(rho3)
    cax3.set_clim(np.min(0), np.max(rho3))
    
    # cax1.set_clim(np.min(0), np.max(rho))
    plt.suptitle("n={}, sum={:.4f}".format(n,np.sum(rho2)/(n1*n2)),fontsize=6)
    return cax1, cax2, cax3,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, 
                               frames=nt+1, interval=10, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save("video.mp4", fps=30)
# sorry Julia

import matplotlib.pyplot as plt
import numpy as np

def save(fig, name):
    print("Saving %s.png." % (name))
    fig.savefig(name+".png", dpi=200)

characteristics = True

if characteristics:
    def phi (XX, YY, c):
        return (YY**2)/2.0 + c * XX**4

    [XX,YY] = np.meshgrid (*[np.linspace(-1.0,1.0,100) for _ in [0,1]])
    ZZ1 = phi (XX, YY,  1.0/2.0)
    ZZ2 = phi (XX, YY, -1.0/2.0)

    levels=np.linspace(-1.0,1.0,21)
    levels = np.sign(levels) * (abs(levels))**(3/2)+0.00006
    fig, (ax1,ax2) = plt.subplots (1,2,figsize=(8,3))
    for (ax,ZZ,mul,tag) in zip([ax1,ax2],[ZZ1,ZZ2],[100,10],[r"$f_e$",r"$f_i$"]):
        # ax.set_xticks([])
        # ax.set_yticks([])
        # ax.set_xlabel("x")
        # ax.set_ylabel("v")
        ax.contourf(XX,mul*YY,ZZ, levels=levels, cmap="jet")
        ax.contour (XX,mul*YY,ZZ, linewidths=2, colors="k", levels=levels)
        ax.set_title("characteristics for %s" % (tag))
    # plt.show ()
    save(fig, "characteristics")

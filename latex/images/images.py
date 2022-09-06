# sorry Julia

import matplotlib.pyplot as plt
import numpy as np

def save(fig, name):
    print("Saving %s.png." % (name))
    fig.savefig(name+".png", dpi=200)

characteristics = False
malkov_solution = True

if characteristics:
    def phi (XX, YY, c):
        # return (YY**2)/2.0 + c * XX**4
        # return (YY**2)/2.0 + XX**2 / 2.0
        return (YY**2)/2.0 + (XX**2 / 2.0) * (np.abs(XX) <= 1.0) + (np.abs(XX) - 1/2) * (np.abs(XX) > 1.0)
        # return (YY**2)/2.0 + np.min([XX**2 / 2.0, np.abs(XX)/2.0])

    [XX,YY] = np.meshgrid (*[np.linspace(-2.0,2.0,100) for _ in [0,1]])
    ZZ1 = phi (XX, YY,  1.0/2.0)
    ZZ2 = phi (XX, YY, -1.0/2.0)

    levels=np.linspace(-1.0,3.0,21)
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
    plt.show ()
    # save(fig, "characteristics")

if malkov_solution:
    def ff (XX, VV):
        res = 0.0 * XX
        mask = (XX**2 + VV**2 < 1.0)
        res[mask] = 1 / np.sqrt(1.0 - XX[mask]**2 - VV[mask]**2)
        return res

    classics = False
    extended = True

    Nx = 400; Nv = 400
    vv = np.linspace(-2.0, 2.0, Nv)
    fsizeE = (10,5)
    fsizeF = (10,5)
    
    if (classics):
        xxc = np.linspace(-1.0,1.0,Nx)
        EEc = - xxc
        figEc, axec = plt.subplots(figsize=fsizeE)
        axec.plot(xxc, EEc)
        plt.show ()

        [XX, VV] = np.meshgrid(xxc, vv)
        FFc = ff(XX, VV)
        figFc, axfc = plt.subplots(figsize=fsizeF)
        axfc.imshow(FFc, origin="lower", extent=(xxc[0], xxc[Nx-1], vv[Nv-1], vv[0]))
        plt.show()


    if (extended):
        eps = 0.5
        xxe = np.linspace(-1.0-eps,1.0+eps,Nx)
        EEe = np.minimum(1.0, np.maximum(-1.0, - xxe))
        figEe, axee = plt.subplots(figsize=fsizeE)
        axee.plot(xxe, EEe)
        save(figEe, "malkov_solution_Ee")
        # plt.show ()

        # [XX, VV] = np.meshgrid(xxe, vv)
        # FFe = ff(XX, VV)
        # FFe = np.minimum(10, FFe)
        # # figFe, axfe = plt.subplots(figsize=fsizeF)
        # # axfe.imshow(FFe, origin="lower", extent=(xxe[0], xxe[Nx-1], vv[Nv-1], vv[0]))
        # plt.figure()
        # plt.imshow(FFe, origin="lower", extent=(xxe[0], xxe[Nx-1], vv[Nv-1], vv[0]))
        # # figFe.colorbar(plt.cm.ScalarMappable()) # norm=norm, cmap=cmap
        # plt.colorbar()
        # # plt.show()
        # plt.savefig ("malkov_solution_fe.png")





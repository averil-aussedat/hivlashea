# sorry Julia
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

plt.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'

def save(fig, name):
    print("Saving %s.png." % (name))
    fig.savefig(name+".png", dpi=200)

characteristics = False
malkov_solution = False
fixed_point_char_maps = True

if fixed_point_char_maps:
    def phi(x):
        # return - 5 * x**2 / 2.0 - 2.0 * x**4
        return - 5 * x**2 

    def phim1(y):
        return np.sqrt(- y / 5)

    fig, ax = plt.subplots()
    # fig.tight_layout()

    lspineoffset = 0.1
    bspineoffset = 0.9
    ax.spines['left'].set_position(('axes',lspineoffset))
    ax.spines['bottom'].set_position(('axes',bspineoffset))
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')

    for spine in [ax.spines['left'],ax.spines['bottom']]:
        spine.set_linewidth(2)

    ax.plot(1, bspineoffset, ">k", transform=ax.transAxes, clip_on=False)
    ax.plot(lspineoffset, 1, "^k", transform=ax.transAxes, clip_on=False)

    fs=20
    ax.text(lspineoffset+0.03, 1.0, "velocity", fontsize=fs, va='top', ha='left', transform=ax.transAxes)
    ax.text(1.0+0.02, bspineoffset-0.03, "space", fontsize=fs, va='top', ha='center', transform=ax.transAxes)

    lowerboundni_1      = False # oldie
    domain_map_1        = False # oldie
    lowerboundni_2      = False # oldie
    lowerboundni_3      = False
    domain_map_2        = False
    global_bound_fi     = True

    if global_bound_fi:
        xx = np.linspace(0.0,1.0,400)
        vstar = -2.0 
        minv = 1.5*vstar
        mu = 1.0
        charestar = -np.sqrt(np.maximum(0,vstar**2 + 2/mu * phi(xx)))
        xstar = phim1(-mu/2 * vstar**2)
        thev = 0.3 * vstar
        chariv = -np.sqrt(thev**2 - 2*phi(xx))
        inter_x = phim1((thev**2/2 - vstar**2/2)/(1+1/mu))
        inter_v = -np.sqrt(thev**2 - 2*phi(inter_x))

        domcolor = "palegreen"
        linecolor = "k"
        charicolor = "saddlebrown"
        charecolor = "seagreen"

        ax.set_xticks([inter_x, 1.0])
        ax.set_xticklabels([r"$a(x,v)$", "1"])
        ax.tick_params('x', direction='in', width=2, pad=-25.0, labelsize=fs)
        maxx = 1.1
        ax.set_xlim([- lspineoffset * maxx / (1.0 - lspineoffset),maxx])

        ax.set_yticks([thev, vstar])
        ax.set_yticklabels([r"$(x,v)$", r"$-v_*$"])
        ax.tick_params('y', width=2, pad=4.0, labelsize=fs)
        ax.set_ylim([minv,minv * (1 - 1.0 / bspineoffset)])

        ax.fill_between(xx, -10, charestar, color=domcolor)
        ax.plot(xx[xx <= xstar+1e-2], charestar[xx <= xstar+1e-2], linewidth=2, linestyle="--", color=charecolor)
        ax.plot(xx, chariv, linewidth=2, linestyle="--", color=charicolor)
        ax.plot([inter_x, inter_x],[0.0,inter_v],linewidth=2,linestyle="--",color="lightgray")

        save(fig, "global_bound_fi")

    if lowerboundni_3:
        xx = np.linspace(0.0,1.0,400)
        mu = 1.0
        aa = 0.4
        va = -np.sqrt(-2*phi(aa))
        vlow = -np.sqrt(va**2 - 2/mu * phi(aa))
        vupp = vlow-1.0
        minv = vupp-1.5
        mu = 1.0
        charelow = -np.sqrt(np.maximum(0,vlow**2 + 2/mu * phi(xx)))
        chareupp = -np.sqrt(np.maximum(0,vupp**2 + 2/mu * phi(xx)))
        lower_chari = -np.sqrt(np.maximum(0.0,vupp**2 - 2 * phi(xx)))
        critichar = -np.sqrt(- 2*phi(xx))
        thex = 0.75
        gx0 = -np.sqrt(-2*phi(thex))
        gxvupp = -np.sqrt(vupp**2-2*phi(thex))

        domcolor = "palegreen"
        linecolor = "k"
        charicolor = "saddlebrown"
        charecolor = "seagreen"
        hatchcolor = "royalblue"

        ax.set_xticks([aa, thex, 1.0])
        ax.set_xticklabels([r"$a$", r"$x$", "1"])
        ax.tick_params('x', direction='in', width=2, pad=-25.0, labelsize=fs)
        maxx = 1.1
        ax.set_xlim([- lspineoffset * maxx / (1.0 - lspineoffset),maxx])

        ax.set_yticks([vlow, vupp])
        ax.set_yticklabels([r"$-\underline{v}$", r"$-\overline{v}$"])
        ax.tick_params('y', width=2, pad=4.0, labelsize=fs)
        ax.set_ylim([minv,minv * (1 - 1.0 / bspineoffset)])

        # domains
        plt.rcParams['hatch.linewidth'] = 3.5
        ax.fill_between(xx, chareupp, charelow, color=domcolor)
        ax.fill_between(xx[xx<=aa], chareupp[xx<=aa], charelow[xx<=aa], hatch="///", color=domcolor, fc=hatchcolor)

        # lines other than char
        ax.plot([thex,thex],[minv,0.0],linestyle="--",linewidth=2,color=linecolor) # vertical line x
        ax.plot([aa, aa],[0.0,va],linewidth=2,linestyle="--",color="lightgray")

        # chars
        ax.plot(xx, charelow, linewidth=2, linestyle="--", color=charecolor)
        ax.plot(xx, chareupp, linewidth=2, linestyle="--", color=charecolor)
        ax.plot(xx, lower_chari, linewidth=2, linestyle="--", color=charicolor)
        ax.plot(xx, critichar, linewidth=2, linestyle="-", color="k")

        # annotations
        ax.scatter(thex, gx0,    s=12, color=linecolor)
        ax.scatter(thex, gxvupp, s=12, color=linecolor)
        ax.text(thex+0.03, gx0, r"$-g_x(0)$", fontsize=fs, ha="left", va="bottom")
        ax.text(thex-0.03, gxvupp, r"$-g_x(-\overline{v})$", fontsize=fs, ha="right", va="top")

        save(fig, "fpcharmaps_lowerboundni")

    if lowerboundni_1:
        thex = 0.8
        domfel = 0.5
        domfeu = 1.0
        mu = 0.4
        xx = np.linspace(0.0,1.0,200)
        upper_chare = -np.sqrt(np.maximum(0.0,domfel**2 + 2/mu * phi(xx)))
        lower_chare = -np.sqrt(np.maximum(0.0,domfeu**2 + 2/mu * phi(xx)))
        upper_chari = -np.sqrt(np.maximum(0.0,domfel**2 - 2 * phi(xx)))
        lower_chari = -np.sqrt(np.maximum(0.0,domfeu**2 - 2 * phi(xx)))
        wlow = -np.sqrt(np.maximum(0.0,domfel**2 - 2 * phi(thex)))
        wupp = -np.sqrt(np.maximum(0.0,domfeu**2 - 2 * phi(thex)))
        minv = 2.0 * np.min(lower_chare)

        domcolor = "palegreen"
        linecolor = "k"
        charicolor = "saddlebrown"
        charecolor = "seagreen"

        ax.set_xticks([thex, 1.0])
        ax.set_xticklabels([r"$x$", "1"])
        ax.tick_params('x', direction='in', width=2, pad=-22.0, labelsize=fs)
        maxx = 1.1
        ax.set_xlim([- lspineoffset * maxx / (1.0 - lspineoffset),maxx])

        ax.set_yticks([-domfeu, -domfel])
        ax.set_yticklabels([r"$-\overline{v}$", r"$-\underline{v}$"])
        ax.tick_params('y', width=2, pad=4.0, labelsize=fs)
        ax.set_ylim([minv,minv * (1 - 1.0 / bspineoffset)])

        ax.plot([thex,thex],[minv,0.0],linestyle="--",linewidth=2,color=linecolor) # vertical line x
        ax.plot(xx, lower_chare, linewidth=2, linestyle="--", color=charecolor)
        ax.plot(xx, upper_chari, linewidth=2, linestyle="--", color=charicolor)
        ax.plot(xx, lower_chari, linewidth=2, linestyle="--", color=charicolor)
        ax.fill_between(xx, lower_chare, upper_chare, color=domcolor)

        ax.scatter(thex,wlow,s=50,color=charicolor)
        ax.scatter(thex,wupp,s=50,color=charicolor)
        ax.text(thex+0.03, wlow, r"$-\underline{w}$", fontsize=fs, ha="left", va="bottom")
        ax.text(thex-0.03, wupp, r"$-\overline{w}$", fontsize=fs, ha="right", va="top")
        # thetext = r"$\begin{pmatrix}\overline{y}(w) \\ w \end{pmatrix} = \begin{pmatrix} y \\ h(y) \end{pmatrix}$"
        thetext = r"$\begin{pmatrix}\overline{y}(w) \\ w \end{pmatrix}$"
        # thetext = r"$(\overline{y}(w), w ) = (y, h(y))$"
        ax.text(0.54, -0.6*domfel, thetext, fontsize=fs, ha="left", va="center", color=charecolor)

        save(fig, "fpcharmaps_lowerboundni")

    if lowerboundni_2:
        thex = 0.8
        domfel = 0.5
        domfeu = 1.0
        mu = 0.4
        xx = np.linspace(0.0,1.0,200)
        upper_chare = -np.sqrt(np.maximum(0.0,domfel**2 + 2/mu * phi(xx)))
        lower_chare = -np.sqrt(np.maximum(0.0,domfeu**2 + 2/mu * phi(xx)))
        upper_chari = -np.sqrt(np.maximum(0.0,domfel**2 - 2 * phi(xx)))
        lower_chari = -np.sqrt(np.maximum(0.0,domfeu**2 - 2 * phi(xx)))
        criti_char  = -np.sqrt(- 2 * phi(xx))
        wlow = -np.sqrt(np.maximum(0.0,domfel**2 - 2 * phi(thex)))
        criw = -np.sqrt(-2*phi(thex))
        wupp = -np.sqrt(np.maximum(0.0,domfeu**2 - 2 * phi(thex)))
        minv = 2.0 * np.min(lower_chare)

        domcolor = "palegreen"
        linecolor = "k"
        charicolor = "saddlebrown"
        charecolor = "seagreen"

        ax.set_xticks([thex, 1.0])
        ax.set_xticklabels([r"$x$", "1"])
        ax.tick_params('x', direction='in', width=2, pad=-22.0, labelsize=fs)
        maxx = 1.1
        ax.set_xlim([- lspineoffset * maxx / (1.0 - lspineoffset),maxx])

        ax.set_yticks([-domfeu, -domfel])
        ax.set_yticklabels([r"$-\overline{v}$", r"$-\underline{v}$"])
        ax.tick_params('y', width=2, pad=4.0, labelsize=fs)
        ax.set_ylim([minv,minv * (1 - 1.0 / bspineoffset)])

        ax.plot([thex,thex],[minv,0.0],linestyle="--",linewidth=2,color=linecolor) # vertical line x
        ax.plot(xx, upper_chare, linewidth=2, linestyle="--", color=charecolor)
        ax.plot(xx, lower_chare, linewidth=2, linestyle="--", color=charecolor)
        ax.plot(xx, lower_chari, linewidth=2, linestyle="--", color=charicolor)
        ax.plot(xx, criti_char, linewidth=2, linestyle="-", color=linecolor)
        ax.fill_between(xx, lower_chare, upper_chare, color=domcolor)

        ax.scatter(thex,criw,s=30,color=linecolor)
        # ax.scatter(thex,wlow,s=30,color=charicolor)
        ax.scatter(thex,wupp,s=30,color=charicolor)
        ax.text(thex+0.03, criw, r"$-g_x(0)$", fontsize=fs, ha="left", va="bottom")
        ax.text(thex-0.03, wupp, r"$-g_x(-\overline{v})$", fontsize=fs, ha="right", va="top")
        # thetext = r"$\begin{pmatrix}\overline{y}(w) \\ w \end{pmatrix} = \begin{pmatrix} y \\ h(y) \end{pmatrix}$"
        # thetext = r"$\begin{pmatrix}\overline{y}(w) \\ w \end{pmatrix}$"
        # thetext = r"$(\overline{y}(w), w ) = (y, h(y))$"
        # ax.text(0.54, -0.6*domfel, thetext, fontsize=fs, ha="left", va="center", color=charecolor)

        save(fig, "fpcharmaps_lowerboundni")

    if domain_map_1:

        linecolor = "k"
        domUpL = r"$\mathcal{D}_1$"; colordomUpL = "gold"
        domUpR = r"$\mathcal{D}_2$"; colordomUpR = "violet"
        domLow = r"$\mathcal{D}_3$"; colordomLow = "palegreen"

        thex = 0.7
        xx = np.linspace(0.0,1.0,200)
        xxchar = - np.sqrt(np.maximum(0.0, 2*(phi(thex)-phi(xx))))
        critical_char = -np.sqrt(-2*phi(xx))
        minv = 1.5 * np.min(critical_char)
        v0 = -np.sqrt(-2*phi(thex))

        ax.set_xticks([thex, 1.0])
        ax.set_xticklabels([r"$x$", "1"])
        ax.tick_params('x', direction='in', width=2, pad=-22.0, labelsize=fs)
        maxx = 1.1
        ax.set_xlim([- lspineoffset * maxx / (1.0 - lspineoffset),maxx])

        ax.set_yticks([v0])
        ax.set_yticklabels([r"$v_0(x)$"])
        ax.tick_params('y', width=2, pad=4.0, labelsize=fs)
        ax.set_ylim([minv,minv * (1 - 1.0 / bspineoffset)])
        
        ax.fill_between(xx, minv * np.ones_like(xx), critical_char, color=colordomLow)
        ax.fill_between(xx, critical_char * (xx <= thex), xxchar * (xx <= thex), color=colordomUpL)
        ax.fill_between(xx, critical_char * (xx >= thex), xxchar * (xx >= thex), color=colordomUpR)

        ax.text(0.58, minv*0.13, domUpL, fontsize=fs, ha='center', va='center')
        ax.text(0.78, minv*0.33, domUpR, fontsize=fs, ha='center', va='center')
        ax.text(0.30, minv*0.50, domLow, fontsize=fs, ha='center', va='center')

        ax.plot(xx,critical_char,color=linecolor,linewidth=2) # critical char
        ax.plot([0.0,thex],[v0,v0],linestyle="--",linewidth=2,color="mediumseagreen") # horiz line v0
        ax.plot([thex,thex],[minv,0.0],linestyle="--",linewidth=2,color=linecolor) # vertical line x
    
        save(fig, "fpcharmaps_domainmap")

    if domain_map_2:

        linecolor = "k"
        domUpL = r"$\mathcal{I}_1$"; colordomUpL = "gold"
        domUpR = r"$\mathcal{I}_2$"; colordomUpR = "violet"
        domLow = r"$\mathcal{I}_3$"; colordomLow = "palegreen"
        charcolorupl="peru"
        charcolorlow="mediumseagreen"

        thex = 0.7
        they = 0.6
        xx = np.linspace(0.0,1.0,200)
        xxchar = - np.sqrt(np.maximum(0.0, 2*(phi(thex)-phi(xx))))
        critical_char = -np.sqrt(-2*phi(xx))
        minv = 1.2 * np.min(critical_char)
        v0 = -np.sqrt(-2*phi(thex))
        vy = -np.sqrt(2*phi(they)-2*phi(thex))
        xx_y = np.linspace(they, thex, 100)
        chary = - np.sqrt(np.maximum(0.0, 2*(phi(they)-phi(xx_y))))
        charlows = [-np.sqrt((multiple*minv)**2 / 2.0 - 2 * phi(xx)) for multiple in [0.5,0.85,1.15]]

        ax.set_xticks([they, thex, 1.0])
        ax.set_xticklabels([r"$y$", r"$x$", "1"])
        ax.tick_params('x', direction='in', width=2, pad=-22.0, labelsize=fs)
        maxx = 1.1
        ax.set_xlim([- lspineoffset * maxx / (1.0 - lspineoffset),maxx])

        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.tick_params('y', width=2, pad=4.0, labelsize=fs)
        ax.set_ylim([minv,minv * (1 - 1.0 / bspineoffset)])
        
        ax.fill_between(xx, minv * np.ones_like(xx), critical_char, color=colordomLow)
        ax.fill_between(xx, critical_char * (xx <= thex), xxchar * (xx <= thex), color=colordomUpL)
        ax.fill_between(xx, critical_char * (xx >= thex), xxchar * (xx >= thex), color=colordomUpR)

        ax.text(0.40, minv*0.13, domUpL, fontsize=1.2*fs, ha='center', va='center')
        ax.text(0.89, minv*0.60, domUpR, fontsize=1.2*fs, ha='center', va='center')
        ax.text(0.12, minv*0.23, domLow, fontsize=1.2*fs, ha='center', va='center')

        ax.plot(xx,critical_char,color=linecolor,linewidth=2) # critical char
        # ax.plot([0.0,thex],[v0,v0],linestyle="--",linewidth=2,color=charcolor) # horiz line v0
        ax.plot([thex,thex],[minv,0.0],linestyle="--",linewidth=2,color=linecolor) # vertical line x

        ax.plot(xx_y,chary,linestyle="--",linewidth=2,color=charcolorupl)
        for charlow in charlows:
            ax.plot(xx,charlow,linestyle="--",linewidth=2,color=charcolorlow)

        ax.scatter(thex,v0,s=15,color="k",zorder=10) # -g_x(0)
        ax.scatter(thex,vy,s=15,color="k",zorder=10) # -g_x(0)
        ax.text(thex+0.02,v0,r"$-g_x(0)$", fontsize=fs, ha="left", va="bottom") # -g_x(0)
        ax.text(thex+0.02,vy,r"$-g_x(y)$", fontsize=fs, ha="left", va="center") # -g_x(y)
        
        save(fig, "fpcharmaps_domainmap")

    plt.show()

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






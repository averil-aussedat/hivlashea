import h5py
import numpy as np
import matplotlib.pyplot as plt

from os import listdir # browsing current directory 
import re # regex for filename matching
from datetime import datetime # date.
from mpl_toolkits.axes_grid1.inset_locator import inset_axes # colorbar

import matplotlib.animation as animation # mp4 videos

####################################################################
#          Bank of plotting functions ---- CEMRACS 2022            #
####################################################################

print("Welcome in diag_sheath.")
plt.rcParams.update({'font.size': 10})
plt.rcParams['font.family'] = 'monospace' # fixed-case font *-*

#################################
# Parameters                    
#################################

loadrep = "../data/"
saverep_img = "../images/"
saverep_mp4 = "../videos/"
yamlfilename = "test_2sp.yaml"

informations_on_graphs = True # add parameters on graphs. Keep true please...

plot_fi = True # density of ions 
plot_fe = True # density of electrons
plot_EE = True # electric field

video_fi = True # density of ions 
video_fe = True # density of electrons
video_EE = True # electric field

thefiles = listdir (loadrep)
times = [] # strong assumption : E, fe an fi are all plotted on same times

show_graphs = False # only for graphs, not videos
save_graphs = True # only for graphs, not videos

#################################
# Loading data
#################################

### yaml
yamlfile = open(yamlfilename, "r")
yamldata = yamlfile.readlines ()
yamlfile.close ()

def readyamlnode (indent, name, type, node):
    nocom = node[len(indent*'  '+name+':'):].split("#")[0]
    nocom = nocom.replace(' ','')
    nocom = nocom.replace('\n','')
    return type(nocom)

# print("yamldata : ", yamldata)

def readmesh (meshname, yamldata):
    res = dict()
    for iy, y in enumerate(yamldata):
        if meshname+":" in y:
            res["name"] = meshname
            res["min"] = readyamlnode (2,"min",float,yamldata[iy+2])
            res["max"] = readyamlnode (2,"max",float,yamldata[iy+3])
            res["N"] = readyamlnode (2,"N",int,yamldata[iy+4])
    return res

meshx  = readmesh("meshx", yamldata)
meshve = readmesh("meshve", yamldata)
meshvi = readmesh("meshvi", yamldata)
# print("meshx : ", meshx)
# print("meshve : ", meshve)
# print("meshvi : ", meshvi)

def readadv (advname, yamldata): # only non-periodic ones
    res = dict()
    for iy, y in enumerate(yamldata):
        if advname+":" in y:
            res["name"] = advname
            res["d"] = readyamlnode (2,"d",int,yamldata[iy+2])
            res["kb"] = readyamlnode (2,"kb",int,yamldata[iy+3])
            res["v"] = readyamlnode (1,"v",float,yamldata[iy+4])
    return res

adv_x  = readadv("adv_x", yamldata)
adv_ve = readadv("adv_ve", yamldata)
adv_vi = readadv("adv_vi", yamldata)
# print("adv_x : ", adv_x)
# print("adv_ve : ", adv_ve)
# print("adv_vi : ", adv_vi)

for y in yamldata:
    if "dt:" in y:
        dt = readyamlnode (1,"dt",float,y)
    if "num_iteration:" in y:
        num_iteration = readyamlnode (1,"num_iteration",int,y)
    if "lambda:" in y:
        llambda = readyamlnode (1,"lambda",float,y)
    if "nu:" in y:
        nu = readyamlnode (1,"nu",float,y)
# print("dt : ", dt, ", num iterations : ", num_iteration, ", lambda : ", llambda, ", nu : ", nu)

if (informations_on_graphs):
    Infos = ""
    # unstructered information
    Infos += "SL for Vlasov-Poisson, " + datetime.now().strftime("%d-%m-%Y %Hh%M") + "\n"
    Infos += "T=%6.3f, %d time iterations (dt=%6.3f)\n" % (dt*num_iteration, num_iteration, dt)
    Infos += "lambda = %6.3f, nu = %6.3f\n" % (llambda, nu)
    # mesh infos
    for mesh in [meshx, meshve, meshvi]:
        Infos += "%8s : %4d points in [%7.3f, %7.3f]\n" % (mesh["name"],mesh["N"]+1,mesh["min"],mesh["max"])
    # advection infos
    for adv in [adv_x, adv_ve, adv_vi]:
        Infos += "%8s : v = %7.2f, d = %2d, kb = %2d\n" % (adv["name"], adv["v"], adv["d"], adv["kb"])

### Electric fields E
if (plot_EE or video_EE or video_fe or video_fi): # get times.
    # first line : simulation time
    # all the others : x_i E_i
    Efields = []; Emin=1e6; Emax=-1e6; E_its = []
    Efieldsnames = [filename for filename in thefiles if re.match("E[0-9]*\.dat", filename)] 
    Efieldsnames = np.sort(Efieldsnames)
    if (plot_EE and not video_EE):
        Efieldsnames = [Efieldsnames[0],Efieldsnames[-1]]
    print("Reading electric field files%s\n" % ("" if (plot_EE or video_EE) else" (just for times)"), Efieldsnames)
    for Efieldname in Efieldsnames:
        Efile = open(loadrep+Efieldname, "r")
        if (plot_EE or video_EE):
            Edata = Efile.readlines ()
            times.append (float(Edata.pop(0).replace('\n','')))
            Edata = [[float(x.replace('\n','')) for x in line.split(' ')] for line in Edata]
            Efields.append ([ee[1] for ee in Edata]) 
            xxE = [ee[0] for ee in Edata]
            E_its.append(int(Efieldname[-len("0000.dat"):-len(".dat")])) 
            Emin = np.min([Emin,np.min(Efields[-1])])
            Emax = np.max([Emax,np.max(Efields[-1])])
        else:
            times.append (float(Efile.readline().replace('\n',''))) # just read first line
        Efile.close ()

fis=[]; fes=[] # values 
fi_its = []; fe_its = [] # number of iterations

# Densities
for (plot,vid,tag,ffs,its) in zip([plot_fe,plot_fi],[video_fe,video_fi],["fe","fi"],[fes,fis],[fe_its,fi_its]):
    if (plot or vid):
        fsnames = [filename for filename in thefiles if re.match("%s[0-9]*-values\.h5" % tag, filename)] 
        fsnames = np.sort(fsnames)
        if (plot and not vid):
            fsnames = [fsnames[0], fsnames[-1]]
        print("Reading %s density files\n" % tag, fsnames)
        for fname in fsnames:
            fdict = h5py.File(loadrep+fname,'r')
            ffs.append(np.transpose(fdict["values"]))
            its.append(int(fname[-len("0000-values.h5"):-len("-values.h5")])) 

maxe = max([np.max(ff) for ff in fes])
maxi = max([np.max(ff) for ff in fis])

if (maxe<=0.0 or maxi<=0):
    print("WARNING : one of the two species is absent (max f_{i/e} <=0.0).")
    maxe=max([maxe,0.5]); maxi=max([maxi,0.5])

#################################
# Plotting figures
#################################

feextent = [meshx["min"],meshx["max"],meshve["min"],meshve["max"]] # left, right, bottom, top
fiextent = [meshx["min"],meshx["max"],meshvi["min"],meshvi["max"]] # left, right, bottom, top

def writeInfos (axInfos, Infos):
    axInfos.text(-0.1,0.5,Infos,ha='left',va='center',fontsize=14)
    for pos in ['left', 'right', 'top', 'bottom']:
        axInfos.spines[pos].set_color('none')
    axInfos.set_xticks([])
    axInfos.set_yticks([])

def createAxesInfos (fig, Infos):
    ax0 = fig.add_subplot (1,3,1)
    axT = fig.add_subplot (1,3,2,sharex=ax0,sharey=ax0)
    axInfos = fig.add_subplot (1,3,3)
    writeInfos (axInfos, Infos)
    return [ax0, axT, axInfos]

def show_andor_save (fig, name):
    if (show_graphs):
        plt.show ()
    if (save_graphs):
        title = saverep_img+name+".png"
        print("Saving %s..." % title)
        fig.savefig(title,dpi=200)
        if not show_graphs:
            plt.close(fig)
        print("... Done saving %s." % title)

if (plot_EE):
    if (informations_on_graphs):
        figE = plt.figure (figsize=(15,4))
        [axE0,axET,axEInfos] = createAxesInfos (figE, Infos)
    else:
        figE, (axE0, axET) = plt.subplots (1,2,sharex=True,sharey=True,figsize=(11,4))

    for iax, ax, tt, Efield in zip([0,1],[axE0,axET], [times[0],times[-1]], [Efields[0],Efields[-1]]):
        ax.plot(xxE, Efield, marker=".", color="steelblue")
        ax.set_title("Time t=%6.3f" % tt)
        ax.set_xlabel("x")
        if (iax==0):
            ax.set_ylabel("E")
        ax.set_ylim([Emin-0.05*(Emax-Emin),Emax+0.01*(Emax-Emin)])

    figE.suptitle("Electric field at initial and final time")
    show_andor_save (figE, "E")

for (plot,ffs,extent,species,maxx) in zip([plot_fi, plot_fe],[fis,fes],[fiextent,feextent],["Ion", "Electron"],[maxi,maxe]):
    if (plot):
        if (informations_on_graphs):
            fig = plt.figure (figsize=(15,4))
            [ax0,axT,axInfos] = createAxesInfos (fig, Infos)
        else:
            fig, (ax0, axT) = plt.subplots (1,2,sharex=True,sharey=True,figsize=(11,4))

        for iax, ax, ff in zip([0,1], [ax0, axT], [ffs[0],ffs[-1]]):
            im = ax.imshow (ff,vmin=0.0,vmax=maxx,cmap='jet',origin='lower',extent=extent,aspect="auto",interpolation='spline36')
            ax.set_xlabel("x")
            if (iax==0):
                ax.set_ylabel("v")
                # source : https://stackoverflow.com/questions/13310594/positioning-the-colorbar
                axins = inset_axes(ax, width="5%", height="100%", loc='center left', borderpad=-10)
                fig.colorbar(im, cax=axins, orientation="vertical")

        fig.suptitle("%s density at initial and final time" % species)
        show_andor_save (fig, species)

#################################
# Creating videos
#################################

# Source : https://stackoverflow.com/questions/34975972/how-can-i-make-a-video-from-array-of-images-in-matplotlib

factor = 1.0 # 1.0 = real time, less = slow motion, more = accelerated
video_length = factor * dt * num_iteration * 1000 # in seconds

if (video_EE):
    name = saverep_mp4+"E.mp4"
    print("Creating video %s... (please wait, it can take a few seconds)" % name)
    if (informations_on_graphs):
        figvE, (axvE, axvEInfos) = plt.subplots (1,2,figsize=(11,4))
        writeInfos (axvEInfos, Infos)
    else:
        figvE, axvE = plt.subplots (figsize=(8,4))
    axvE.set_xlabel("x")
    axvE.set_ylabel("E")

    frames = []
    for tt, it, E in zip(times, E_its, Efields):
        label = "t = %6.3f / %6.3f, it %d / %d" % (tt, dt*num_iteration, it, num_iteration)
        txt = axvE.text(0.0,-0.1*Emin+1.1*Emax,label,fontsize=14,ha="center",va="bottom")
        im, = axvE.plot (xxE, E, marker=".", color="steelblue", animated=True)
        frames.append([im, txt])

    videoE = animation.ArtistAnimation (figvE, frames, interval=video_length / len(frames), repeat=False)
    videoE.save (name)
    print("...Done creating video %s." % name)
    # plt.show ()

for (vid,its,ffs,meshv,extent,tag,maxx) in zip([video_fe,video_fi],[fe_its,fi_its],[fes,fis],[meshve,meshvi],[feextent,fiextent],["fe","fi"],[maxe,maxi]):
    if (vid):
        name = saverep_mp4+tag+".mp4"
        print("Creating video %s... (please wait, it can take a few seconds)" % name)
        if (informations_on_graphs):
            figv, (axv, axvInfos) = plt.subplots (1,2,figsize=(11,4))
            writeInfos (axvInfos, Infos)
        else:
            figv, axv = plt.subplots (figsize=(8,4))
        axv.set_xlabel("x")
        axv.set_ylabel("v")

        frames = []
        for tt, it, ff in zip(times, its, ffs):
            label = "t = %6.3f / %6.3f, it %d / %d" % (tt, dt*num_iteration, it, num_iteration)
            txt = axv.text(0.0,-0.01*meshv["min"]+1.01*meshv["max"],label,fontsize=14,ha="center",va="bottom")
            im = axv.imshow (ff,vmin=0.0,vmax=maxx,cmap='jet',origin='lower',extent=extent,aspect="auto",interpolation='spline36', animated=True)
            if (it==its[0]):
                # source : https://stackoverflow.com/questions/13310594/positioning-the-colorbar
                axins = inset_axes(axv, width="5%", height="100%", loc='center left', borderpad=-10)
                figv.colorbar(im, cax=axins, orientation="vertical")
            frames.append([im, txt])            

        video = animation.ArtistAnimation (figv, frames, interval=video_length / len(frames), repeat=False)
        video.save (name)
        print("...Done creating video %s." % (name))
        # plt.show ()

####################################################################
#                                                       the end :) #
####################################################################

print("Bye")

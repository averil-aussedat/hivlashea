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

# loadrep = "data_output/"
# root = "run_Malkov_1sp_d2_Nx800_Nt1000/"
# root = "run_comp_short_time_2sp_Nx1000_Nvi2001_Nve2001_Nt6250/"
root = "" # current directory (default)

loadrep = root+"data_output/"
saverep_img = root+"python_diags/"
#saverep_img = "images/"
# rep = "nu1_T25_Yann"
# loadrep = "../%s/data/" % rep
# saverep_img = "../%s/"  % rep
#saverep_mp4 = "videos/"
#saverep_diag = "diags/" # where to save stationary tests, for instance
saverep_mp4 = root+"python_diags/"
saverep_diag = root+"python_diags/" # where to save stationary tests, for instance

# yamlfilename = "test_2sp.yaml"
# yamlfilename = "Badsi_2022.yaml"
# yamlfilename = "Badsi_2022_apr_test.yaml"
# yamlfilename = "fpcase2_256.yaml"
# yamlfilename = "fpcase2_512.yaml"
# yamlfilename = "fpcase2_1024.yaml"
# yamlfilename = "compMehdi.yml"
# yamlfilename = "badsiBerthonCrestetto2014.yml"

# yamlfilename = "Malkov_1sp.yaml"
# yamlfilename = "comp_short_time.yaml"
# yamlfilename = "comp_long_time.yaml"
yamlfilename = "fun.yaml"

if (root==""):
    yamlfilename = "yaml/" + yamlfilename # add path to yaml
else:
    yamlfilename = root + yamlfilename # add path to yaml

informations_on_graphs = False # add parameters on graphs. Keep true please...
fullImages = False # if True, initial + final + infos on images. Else, only final time

plot_fi  = True # density of ions 
plot_fe  = True # density of electrons
plot_EE  = True # electric field
plot_rho = True # \int_v f_i - f_e dv

video_fi  = True # density of ions 
video_fe  = True # density of electrons
video_EE  = True # electric field
video_rho = True # \int_v f_i - f_e dv

# compute and display stationary test. Saved in file as lines + plotted
# norm of (E_k - E_{k-1} / dt)
var_and_drifts = True
# variation_EE  = True # electric field
# variation_rho = True # \int_v f_i - f_e dv
# # norm of (E_k - E_0)
# drift_EE = True
# drift_rho = True

thefiles = listdir (loadrep)
times = [] # strong assumption : all data are all plotted on same times

show_graphs = False # only for graphs, not videos
save_graphs = True # only for graphs, not videos

factor = 500.0 # 1.0 = real time, more = slow motion, less = accelerated
cmap="gnuplot"

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
            if ("_non_" in yamldata[iy+1]):
                res["kb"] = readyamlnode (2,"kb",int,yamldata[iy+3])
                res["v"] = readyamlnode (1,"v",float,yamldata[iy+4])
            else:
                res["kb"] = 0
                res["v"] = readyamlnode (1,"v",float,yamldata[iy+3])
                
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
    Infos += "T=%6.3f, %d time iterations (dt=%2.2e)\n" % (dt*num_iteration, num_iteration, dt)
    Infos += "lambda = %6.3f, nu = %6.3f\n" % (llambda, nu)
    # mesh infos
    for mesh in [meshx, meshve, meshvi]:
        Infos += "%8s : %4d points in [%7.2f, %7.2f]\n" % (mesh["name"],mesh["N"]+1,mesh["min"],mesh["max"])
    # advection infos
    for adv in [adv_x, adv_ve, adv_vi]:
        Infos += "%8s : v = %7.2f, d = %2d, kb = %2d\n" % (adv["name"], adv["v"], adv["d"], adv["kb"])

Efields = []; Emin=1e6; Emax=-1e-6; E_its=[]; 
rhos = []; rhomin=1e6; rhomax=-1e-6; rho_its=[]; 
fes=[]; fe_its = []; minfe=1.0; maxfe=-1.0
fis=[]; fi_its = []; minfi=1.0; maxfi=-1.0
var_EE = []; var_rho = []; dri_EE = []; dri_rho = [] # measures of stationarity

### Electric fields E + we do it to get times and xx
# first line : simulation time
# all the others : x_i E_i
Efieldsnames = [filename for filename in thefiles if re.match("E[0-9]*\.dat", filename)] 
Efieldsnames = np.sort(Efieldsnames)
if (plot_EE and not (video_EE and var_and_drifts)):
    Efieldsnames = [Efieldsnames[0],Efieldsnames[-1]]
print("Reading electric field files\n", Efieldsnames)
for Efieldname in Efieldsnames:
    Efile = open(loadrep+Efieldname, "r")
    if (plot_EE or video_EE or var_and_drifts):
        Edata = Efile.readlines ()
        times.append (float(Edata.pop(0).replace('\n','')))
        Edata = [[float(x.replace('\n','')) for x in line.split(' ')] for line in Edata]
        Efields.append ([ee[1] for ee in Edata]) 
        E_its.append(int(Efieldname[-len("000000.dat"):-len(".dat")])) 
        Emin = np.min([Emin,np.min(Efields[-1])])
        Emax = np.max([Emax,np.max(Efields[-1])])
    else:
        times.append (float(Efile.readline().replace('\n',''))) # just read first line
    Efile.close ()

### rho
if (plot_rho or video_rho or var_and_drifts):
    # first line : simulation time
    # all the others : x_i rho_i
    rhosnames = [filename for filename in thefiles if re.match("rho[0-9]*\.dat", filename)] 
    rhosnames = np.sort(rhosnames)
    if (plot_rho and not (video_rho or var_and_drifts)):
        rhosnames = [rhosnames[0],rhosnames[-1]]
    print("Reading rho files\n", rhosnames)
    for rhoname in rhosnames:
        rhofile = open(loadrep+rhoname, "r")
        rhodata = rhofile.readlines ()
        rhodata.pop(0) # simu time
        rhodata = [[float(x.replace('\n','')) for x in line.split(' ')] for line in rhodata]
        rhos.append ([ee[1] for ee in rhodata]) 
        rho_its.append(int(rhoname[-len("000000.dat"):-len(".dat")])) 
        rhomin = np.min([rhomin,np.min(rhos[-1])])
        rhomax = np.max([rhomax,np.max(rhos[-1])])
        rhofile.close ()

# Densities
for (plot,vid,tag,ffs,its) in zip([plot_fe,plot_fi],[video_fe,video_fi],["fe","fi"],[fes,fis],[fe_its,fi_its]):
    if (plot or vid):
        fsnames = [filename for filename in thefiles if re.match("%s[0-9]*-values\.h5" % tag, filename)] 
        fsnames = np.sort(fsnames)
        # if (plot and not (vid or statio)):
        if (plot and not vid):
            fsnames = [fsnames[0], fsnames[-1]]
        print("Reading %s density files\n" % tag, fsnames)
        for fname in fsnames:
            fdict = h5py.File(loadrep+fname,'r')
            ffs.append(np.transpose(fdict["values"]))
            its.append(int(fname[-len("000000-values.h5"):-len("-values.h5")])) 

# Computing the stationarity
if (var_and_drifts):
    dx = (meshx["max"]-meshx["min"])/meshx["N"]
    var_EE_i = 0.0; var_EE_1 = 0.0; dri_EE_i = 0.0; dri_EE_1 = 0.0
    for i in range(len(Efields)-1):
        var_EE_i = np.max([var_EE_i,             np.max ([(Efields[i+1][j] - Efields[i][j]) / (times[i+1] - times[i]) for j in range(len(Efields[i]))])])
        var_EE_1 = np.max([var_EE_1, dx * np.linalg.norm([(Efields[i+1][j] - Efields[i][j]) / (times[i+1] - times[i]) for j in range(len(Efields[i]))])])
        dri_EE_i = np.max([dri_EE_i,             np.max ([ Efields[i+1][j] - Efields[0][j] for j in range(len(Efields[i]))])])
        dri_EE_1 = np.max([dri_EE_1, dx * np.linalg.norm([ Efields[i+1][j] - Efields[0][j] for j in range(len(Efields[i]))])])
    var_rho_i = 0.0; var_rho_1 = 0.0; dri_rho_i = 0.0; dri_rho_1 = 0.0
    for i in range(len(rhos)-1):
        var_rho_i = np.max([var_rho_i,             np.max ([(rhos[i+1][j] - rhos[i][j]) / (times[i+1] - times[i]) for j in range(len(rhos[i]))])])
        var_rho_1 = np.max([var_rho_1, dx * np.linalg.norm([(rhos[i+1][j] - rhos[i][j]) / (times[i+1] - times[i]) for j in range(len(rhos[i]))])])
        dri_rho_i = np.max([dri_rho_i,             np.max ([ rhos[i+1][j] - rhos[0][j] for j in range(len(rhos[i]))])])
        dri_rho_1 = np.max([dri_rho_1, dx * np.linalg.norm([ rhos[i+1][j] - rhos[0][j] for j in range(len(rhos[i]))])])

# if variation_EE:
#     [var_EE.append(dx * np.linalg.norm([(Efields[i+1][j] - Efields[i][j]) / (times[i+1] - times[i]) for j in range(len(Efields[i]))])) for i in range(len(Efields)-1)]
# if variation_rho:
#     [var_rho.append(dx * np.linalg.norm([(rhos[i+1][j] - rhos[i][j]) / (times[i+1] - times[i]) for j in range(len(rhos[i]))])) for i in range(len(rhos)-1)]
# if drift_EE:
#     [dri_EE.append(dx * np.linalg.norm([Efields[i+1][j] - Efields[0][j] for j in range(len(Efields[i]))])) for i in range(len(Efields)-1)]
# if drift_rho:
#     [dri_rho.append(dx * np.linalg.norm([rhos[i+1][j] - rhos[0][j] for j in range(len(rhos[i]))])) for i in range(len(rhos)-1)]

# print("var_EE : ", var_EE)
# print("dri_EE : ", dri_EE)
# print("var_rho : ", var_rho)
# print("dri_rho : ", dri_rho)

def get_minmax(dowe, datas, tag):
    themin = 0.0; themax = 0.5
    if dowe:
        themax = max([np.max(dd) for dd in datas[-min([20,len(datas)]):]])
        themin = min([np.min(dd) for dd in datas[-min([20,len(datas)]):]])
        if (themax <= 0):
            print("WARNING : max %s <= 0.0)." % tag)
            themax=0.5
    return themin, themax

minfi, maxfi = get_minmax(plot_fi or video_fi, fis, "fi")
minfe, maxfe = get_minmax(plot_fe or video_fe, fes, "fe")

# print("WARNING setting manual values for min/max fi/e")
# maxfe = 10.0 # for Malkov

# maxfe = 0.04 # short time comparison
# maxfi = 0.5  # short time comparison

#################################
# Verifications
#################################

# EE = Efields[-1]
# rr = rhos[-1]
# dx = (meshx["max"] - meshx["min"])/meshx["N"]
# # for i in range(1,meshx["N"]):
# #     print("err : ", llambda**2 * (EE[i] - EE[i-1])/dx - rr[i])

# err = max([llambda**2 * (EE[i] - EE[i-1])/dx - rr[i] for i in range(1,meshx["N"])])
# print("err : ", err)

# print("EE : ", EE)
# print("rho : ", rr)

#################################
# Plotting figures
#################################

xx = np.linspace (meshx["min"],meshx["max"],meshx["N"]+1)
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
    axT = fig.add_subplot (1,3,2)
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

if (fullImages):
    # functions from (t,x) -> IR (E, rho)
    for (plot, ffs, col, name, tag, minn, maxx) in zip ([plot_EE, plot_rho], [Efields,rhos],["steelblue","coral"],["Electric field","Rho"],["E","rho"],[Emin,rhomin],[Emax,rhomax]):
        if (plot):
            if (informations_on_graphs):
                fig = plt.figure (figsize=(15,4))
                [ax0,axT,axInfos] = createAxesInfos (fig, Infos)
            else:
                fig, (ax0, axT) = plt.subplots (1,2,figsize=(11,4))

            for iax, ax, tt, ff in zip([0,1],[ax0,axT], [times[0],times[-1]], [ffs[0],ffs[-1]]):
                ax.plot(xx, ff, marker=".", linestyle="none", markersize=4, color=col)
                ax.set_title("Time t=%6.3f" % tt)
                ax.set_xlabel("x")
                if (iax==0):
                    ax.set_ylabel(tag)
                ax.set_ylim([minn-0.05*(maxx-minn),maxx+0.01*(maxx-minn)])

            fig.suptitle("%s at initial and final time" % (name))
            show_andor_save (fig, tag)

    # functions from (t,x,v) -> IR (fe, fi)
    for (plot,ffs,extent,species,minn,maxx) in zip([plot_fi, plot_fe],[fis,fes],[fiextent,feextent],["Ion", "Electron"],[minfi,minfe],[maxfi,maxfe]):
        if (plot):
            if (informations_on_graphs):
                fig = plt.figure (figsize=(15,4))
                [ax0,axT,axInfos] = createAxesInfos (fig, Infos)
            else:
                fig, (ax0, axT) = plt.subplots (1,2,figsize=(11,4))

            for iax, ax, ff in zip([0,1], [ax0, axT], [ffs[0],ffs[-1]]):
                im = ax.imshow (ff,cmap=cmap,vmin=minn,vmax=maxx,origin='lower',extent=extent,aspect="auto",interpolation='spline36')
                ax.set_xlabel("x")
                if (iax==1):
                    ax.set_ylabel("v")
                if (iax==0):
                    # source : https://stackoverflow.com/questions/13310594/positioning-the-colorbar
                    axins = inset_axes(ax, width="4%", height="100%", loc='center left', borderpad=-9)
                    fig.colorbar(im, cax=axins, orientation="vertical")

            fig.suptitle("%s density at initial and final time" % species)
            show_andor_save (fig, species)
else:
    # functions from (t,x) -> IR (E, rho)
    for (plot, ffs, col, name, tag, minn, maxx) in zip ([plot_EE, plot_rho], [Efields,rhos],["steelblue","coral"],["Electric field","Rho"],["E","rho"],[Emin,rhomin],[Emax,rhomax]):
        if (plot):
            fig, ax = plt.subplots (figsize=(6,5))
            ax.plot(xx, ffs[-1], marker=".", linestyle="none", markersize=4, color=col)
            ax.set_xlabel("x")
            ax.set_ylabel(tag)
            ax.set_ylim([minn-0.05*(maxx-minn),maxx+0.01*(maxx-minn)])

            fig.suptitle("%s at final time T=%6.3f" % (name, times[-1]))
            show_andor_save (fig, tag)

    # functions from (t,x,v) -> IR (fe, fi)
    for (plot,ffs,extent,species,tag,minn,maxx) in zip([plot_fi, plot_fe],[fis,fes],[fiextent,feextent],["Ion", "Electron"],["fi", "fe"],[minfi,minfe],[maxfi,maxfe]):
        if (plot):
            fig, ax = plt.subplots (figsize=(7,5))
            im = ax.imshow (ffs[-1],cmap=cmap,vmin=minn,vmax=maxx,origin='lower',extent=extent,aspect="auto",interpolation='spline36')
            # im = ax.imshow (ffs[0],cmap=cmap,origin='lower',extent=extent,aspect="auto",interpolation='spline36')
            # ax.set_xlabel("x"); ax.set_ylabel("v")
            fig.colorbar(im)
            fig.tight_layout()
            # fig.suptitle("%s density at final time T=%6.3f" % (species, times[-1]))
            show_andor_save (fig, tag)
            # show_andor_save (fig, tag+"_init")

if (var_and_drifts):
    thefilename = saverep_diag+"var_and_drifts.txt"
    print("Saving %s..." % (thefilename))
    thefile = open(thefilename, "w")
    towrite = []; labels = []
    # physics
    towrite.append("%f" % llambda); labels.append("lambda")
    towrite.append("%f" % nu); labels.append("nu")
    # time mesh
    towrite.append("%f" % dt); labels.append("dt")
    towrite.append("%i" % num_iteration); labels.append("Nt")
    towrite.append("%f" % (num_iteration*dt)); labels.append("T")
    # spatial mesh
    towrite.append("%f" % meshx["min"]); labels.append("xmin")
    towrite.append("%f" % meshx["max"]); labels.append("xmax")
    towrite.append("%d" % meshx["N"]); labels.append("Nx")
    # speed mesh fi
    towrite.append("%f" % meshvi["min"]); labels.append("vimin")
    towrite.append("%f" % meshvi["max"]); labels.append("vimax")
    towrite.append("%d" % meshvi["N"]); labels.append("Nvi")
    # speed mesh fe
    towrite.append("%f" % meshve["min"]); labels.append("vemin")
    towrite.append("%f" % meshve["max"]); labels.append("vemax")
    towrite.append("%d" % meshve["N"]); labels.append("Nve")
    # advection fi
    towrite.append("%f" % adv_vi["d"]); labels.append("di")
    towrite.append("%f" % adv_vi["kb"]); labels.append("kbi")
    towrite.append("%f" % adv_vi["v"]); labels.append("vi")
    # advection fe
    towrite.append("%f" % adv_ve["d"]); labels.append("de")
    towrite.append("%f" % adv_ve["kb"]); labels.append("kbe")
    towrite.append("%f" % adv_ve["v"]); labels.append("ve")
    # variations (max over n of |g_{i+1} - g_i|/dt)
    towrite.append("%f" % var_EE_1); labels.append("var_E_1")
    towrite.append("%f" % var_EE_i); labels.append("var_E_i")
    towrite.append("%f" % var_rho_1); labels.append("var_rho_1")
    towrite.append("%f" % var_rho_i); labels.append("var_rho_i")
    # drifts (max over n of |g_i - g_0|)
    towrite.append("%f" % dri_EE_1); labels.append("dri_E_1")
    towrite.append("%f" % dri_EE_i); labels.append("dri_E_i")
    towrite.append("%f" % dri_rho_1); labels.append("dri_rho_1")
    towrite.append("%f" % dri_rho_i); labels.append("dri_rho_i")
    # writing everything
    labels = [ll+"(%d)"%il for (il, ll) in enumerate(labels)]
    thefile.write("# " + "  ".join(labels) + "\n") 
    thefile.write("\t".join(towrite) + "\n")
    thefile.close()
    print("Done saving %s." % (thefilename))

# for (thebool, thedata, thetag, thetitle) in zip(
#     [variation_EE,variation_rho,drift_EE,drift_rho],
#     [var_EE,var_rho,dri_EE,dri_rho],
#     ["varEE","varrho","driEE","drirho"],
#     [r"$|E_{k+1}-E_{k}|$",r"$|\rho_{k+1}-\rho_{k}|$",r"$|E_{k+1}-E_{0}|$",r"$|\rho_{k+1}-\rho_{0}|$"]
#     ):
#     if thebool:
#         thefile = open(saverep_diag+thetag+".txt", "a")
#         thefile.write("\t".join(["%f" % dd for dd in thedata]) + "\n")
#         thefile.close()

#         fig,ax = plt.subplots()
#         ax.plot(thedata,marker=".")
#         show_andor_save(fig, thetag)

#################################
# Creating videos
#################################

# Source : https://stackoverflow.com/questions/34975972/how-can-i-make-a-video-from-array-of-images-in-matplotlib

# video_length = factor * dt * num_iteration * 1000 # in seconds
video_length = 10000.0 # in seconds

# functions from (t,x) -> IR (E, rho)
for (vid,tag,its,ffs,minn,maxx,col) in zip([video_EE,video_rho],["E","rho"],[E_its,rho_its],[Efields,rhos],[Emin,rhomin],[Emax,rhomax],["steelblue","coral"]):
    if (vid):
        name = saverep_mp4+tag+".mp4"
        print("Creating video %s... (please wait, it can take a few seconds)" % name)
        if (informations_on_graphs):
            figv, (axv, axvInfos) = plt.subplots (1,2,figsize=(11,4))
            writeInfos (axvInfos, Infos)
        else:
            figv, axv = plt.subplots (figsize=(8,4))
        axv.set_xlabel("x")
        axv.set_ylabel(tag)

        frames = []
        for tt, it, ff in zip(times, its, ffs):
            label = "t = %6.3f / %6.3f, it %d / %d" % (tt, dt*num_iteration, it, num_iteration)
            txt = axv.text(0.0,-0.1*minn+1.1*maxx,label,fontsize=14,ha="center",va="bottom")
            im, = axv.plot (xx, ff, marker=".", color=col,linestyle="none", markersize=4, animated=True)
            frames.append([im, txt])

        video = animation.ArtistAnimation (figv, frames, interval=video_length / len(frames), repeat=False)
        video.save (name)
        print("...Done creating video %s." % name)
        # plt.show ()

# functions from (t,x,v) -> IR (fe, fi)
for (vid,its,ffs,meshv,extent,tag,minn,maxx) in zip([video_fe,video_fi],[fe_its,fi_its],[fes,fis],[meshve,meshvi],[feextent,fiextent],["fe","fi"],[minfe,minfi],[maxfe,maxfi]):
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
            im = axv.imshow (ff,cmap=cmap,vmin=minn,vmax=maxx,origin='lower',extent=extent,aspect="auto",interpolation='spline36', animated=True)
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

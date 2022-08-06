##############################################################
# Bank of visual representations for unitary tests
##############################################################

import numpy as np
import matplotlib.pyplot as plt

test_adv = False
test_E = True

def readdata (filename, sep=" "):
    f = open(filename, "r")    
    data = [[float(zz.replace('\n','')) for zz in line.split(sep)] for line in f.readlines()]
    data.pop(0) # skip first line (simu time)
    f.close ()
    return data

if (test_adv):
    def u0 (x):
        # return np.cos(x)-np.cos(-1.0)
        return np.cos(x)-np.cos(3.0)
        # return np.exp(-x**2/0.1)
        # return np.ones_like(x)
        # return x+1.0

    T = 0.5
    speed = -2.0

    fig, (ax0,ax1) = plt.subplots(2,1,figsize=(12,10),sharey=True)
    data0 = readdata ("test_adv0000.dat")
    data1 = readdata ("test_adv0001.dat")

    xx0 = [ll[0] for ll in data0]
    ff0 = [ll[1] for ll in data0]

    xx1 = [ll[0] for ll in data1]
    ff1 = [ll[1] for ll in data1]
    if (speed >= 0):
        sol1 = [0.0 if x - speed * T < -1.0 else u0(x - speed * T) for x in xx1]
    else:
        sol1 = [0.0 if x - speed * T > 3.0 else u0(x - speed * T) for x in xx1]
    
    print("Erreur max : ", np.max([np.abs(ff1[i] - sol1[i]) for i in range(len(sol1))]))
    
    ax0.plot(xx0, [u0(x) for x in xx0], label="condition initiale")
    ax0.plot(xx0, ff0, marker="o", label="approximation")
    ax0.legend()
    
    ax1.plot(xx1, sol1, label="solution exacte")
    ax1.plot(xx1, ff1, marker="o", label="approximation")
    ax1.legend()

    plt.show ()

if (test_E):
    thecase = 0
    plotvalue = True
    plotorder = True

    if (plotvalue):
        if (thecase==0):
            def rho (x):
                # return np.sin(x)
                return np.cos(x)
            dt = 1.0
            speed = -2.0
            def realE (x):
                # return 1.0 - np.cos(x)
                return np.sin(x)

        fig, (ax0,ax1) = plt.subplots(2,1,figsize=(12,10),sharey=True)
        therhodata = readdata("test_E_rho0000.dat")
        xxrho = [zz[0] for zz in therhodata]
        therho = [zz[1] for zz in therhodata]

        E0data = readdata("test_E_E0000.dat")
        E1data = readdata("test_E_E0001.dat")
        xxE = [zz[0] for zz in E0data]
        E0 = [zz[1] for zz in E0data]
        E1 = [zz[1] for zz in E1data]
        E1sol = [realE(x) for x in xxE]
        
        print("Erreur max : ", np.max([np.abs(E1[i] - E1sol[i]) for i in range(len(E1))]))
        
        ax0.plot(xxrho, therho, label="rho")
        ax0.plot(xxE, E0, marker="o", label="E0")
        ax0.legend()
        
        ax1.plot(xxE, E1sol, label="E1 exact")
        ax1.plot(xxE, E1, marker="o", label="E1 approx")
        ax1.legend()

        plt.show ()
	    
    if (plotorder):
        errdata = readdata ("test_errors.dat", sep="\t")
        hh = [zz[0] for zz in errdata]
        ee = [zz[1] for zz in errdata]

        pp = np.polyfit(np.log(hh),np.log(ee),1)
        print("pp : ", pp)
        
        fig, ax = plt.subplots ()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("h")
        ax.set_ylabel("erreur")
        ax.plot (hh, ee, marker="o", linestyle="none", label="erreur", zorder=5)
        ax.plot (hh, [np.exp(pp[0]*np.log(h) + pp[1]) for h in hh], linestyle="--", label="regression")
        ax.legend ()
        ax.set_title ("Order %f" % (pp[0]))
        plt.show ()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

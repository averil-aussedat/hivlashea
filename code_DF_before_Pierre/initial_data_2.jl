""" Initial datum and boundary """

# uncomment the block corresponding to the testcase of interest
# block comment begins with "#=" and ends with "=#"
# have one uncommented, even if you initialize from a data file (for allocation)

#############
# testcase=0: full 2-species maxwellian case
#=
const testcase=0
const mu = 0.01
fi_0(x, v) = ones(size(x)) * exp.(-0.5 .* (v').^2)./sqrt(2pi);
fe_0(x, v) = ones(size(x)) * sqrt(mu).*exp.(-0.5 * mu .* (v').^2)./sqrt(2*pi);
E0(x) = 0.0*x; # to be checked before using Ampere solver
# end testcase=0
=#

#############
# testcase=1: full 2-species masked maxwellian case
#=
const testcase=1
const mu = 0.01#0.5
mask(x;d=0.1,xl=-0.1,xr=0.1) = 0.5*(tanh((x-xl)/d)-tanh((x-xr)/d))
fe_0(x, v) = mask.(x) .* sqrt(abs(mu)).*exp.(-0.5 * abs(mu) .* (v').^2)./sqrt(2*pi);
fi_0(x, v) = mask.(x) .* exp.(-0.5  .* (v').^2)./sqrt(2*pi);
E0(x) = 0.0*x; # to be checked before using Ampere solver
# end testcase=1
=#

#############
# testcase=2: constant 1-species validation (paper Malkov-Kudryavtsev)

const testcase=2
const mu = -1.
fi_0(x,v)=zeros((size(x,1),size(v,1)))
rho0=1.
z0=1.
alpha=1.
a=1.
ap=0.
function fe_0(xx,vv)
    Nx = size(xx,1)
    Nv = size(vv,1)
    fe = zeros((Nx,Nv))
    for j in 1:Nv
        for i in 1:Nx
            coef=alpha*alpha*(z0*z0-xx[i]*xx[i]/a/a)-a*a*(vv[j]-ap*xx[i]/a)*(vv[j]-ap*xx[i]/a)
            mask=0.5*(tanh((coef)/1.e-10)-tanh((coef-10.)/1.e-10))
            fe[i,j]=rho0/pi/sqrt(abs(coef))*mask;
        end
    end
    return fe
end
function E_0(xx)
    E=zeros(size(xx))
    for i in 1:Nx
        if xx[i]<-1.
            E[i]=1.
        elseif xx[i]>1.
            E[i]=-1.
        else
            E[i]=-xx[i]
        end
    end
    return E
end  
E_ex(x)=E_0(x)
fe_ex(x,v)=fe_0(x,v)
# end testcase=2


#############
# testcase=3: time variable 1-species validation (paper Malkov-Kudryavtsev)
#=
include("exact_1species_RK4solver.jl")
const testcase=3
const mu = -1.
fi_0(x,v)=zeros((size(x,1),size(v,1)))
rho0=1.
z0=1.
alpha=sqrt(0.98)
a=1.08
ap=0.
function fe_0(xx,vv)    
    Nx = size(xx,1)
    Nv = size(vv,1)
    fe = zeros((Nx,Nv))
    for j in 1:Nv
        for i in 1:Nx
            coef=alpha*alpha*(z0*z0-xx[i]*xx[i]/a/a)-a*a*(vv[j]-ap*xx[i]/a)*(vv[j]-ap*xx[i]/a)
            mask=0.5*(tanh((coef)/1.e-10)-tanh((coef-10.)/1.e-10))
            fe[i,j]=rho0/pi/sqrt(abs(coef))*mask;
        end
    end
    return fe
end
function E_0(xx)
     E=zeros(size(xx))
    for i in 1:Nx
        if xx[i]<-1.
            E[i]=1.
        elseif xx[i]>1.
            E[i]=-1.
        else
            E[i]=-rho0*xx[i]/a
        end
    end
    return E
end
# values of a and ap=dt_a at final time given by a python RK4 solver
# WARNING!!! values depend on final time
#aT=0.9889868233720464#T=10#0.9160923237643877#T=20#0.9405182218848853#T=5#0.9889868233720464#T=1
#apT=-0.14237509909915097#T=10#0.016234206996915204#T=20#-0.10619953933896482#T=5#-0.14237509909915097#T=1
function E_ex(xx)
    # values of a and ap=dt_a at final time given by a python RK4 solver
    # WARNING!!! values depend on final time
    aT,apT=rk4solver(T,alpha,rho0,z0,a,ap)
    E=zeros(size(xx))
    for i in 1:Nx
        if xx[i]<-1.
            E[i]=1.
        elseif xx[i]>1.
            E[i]=-1.
        else
            E[i]=-rho0*xx[i]/aT
        end
    end
    return E
end
function fe_ex(xx,vv)
    aT,apT=rk4solver(T,alpha,rho0,z0,a,ap)
    Nx = size(xx,1)
    Nv = size(vv,1)
    feex=zeros((Nx,Nv))
    for j in (1:Nv)
        for i in (1:Nx)    
            x=xx[i]
            v=vv[j]
            coef=alpha*alpha*(z0*z0-x*x/aT/aT)-a*a*(v-apT*x/aT)*(v-apT*x/aT)
            mask=0.5*(tanh((coef)/1.e-10)-tanh((coef-10.)/1.e-10))
            feex[i,j]=rho0/pi/sqrt(abs(coef))*mask
        end
    end
    return feex
end
# end testcase=3
=#


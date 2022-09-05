"""

Two-species Vlasov-Poisson solver on [-1,1] \\times \\RR using
 
 with :
 - non-emitting boundary conditions on the distribution functions
 - Source term for the ions
 - Dynamic electric field at the boundary
 
 
 Simulation of the symmetric sheath problem with ionization term.
 
 Discretization :
 - Upwind for the transport
 - Integration for the Poisson equation -\\lambda^2 dx E = rho with trapeze formula for rho.  The integration  preserves the symmetry of rho if it is even.
 
 
 Initial data :
 
 - Initial data are in the file initial_data.hpp
 
 Author : Mehdi BASDI. 
 Date : 27/07/2022.

 Translation in Julia : Pierre Navaro & Averil Prost.
 Date : 05/08/2022.
 
"""

####################### Memo Julia
# Ouvrir une session : Alt+j Alt+o
# Marche à suivre : aller dans packages, taper instantiate
# aller dans package : ]
# aller dans shell : ;
# taper activate (! location)
# taper instantiate
# revenir : backspace
# lancer un main : include("main_file.jl")
# quit : exit()

#= 
    this is a block comment (not doc)
=#

using JLD
using Plots
#using ProgressMeter
using Printf
using DelimitedFiles

include("initial_data_2.jl")  # constant mu and analytic initializations
# functions fi_0, fe_0, E0 

# Physical parameters
const nu = 2.#0.0; # Collision frequency
const lambda = 0.1#1;#0.5; # Debye length

# Discretization parameters
const Nx = 512 #512;#24;#1024;#301; # Number of points of the space mesh
const Nv = 513 #513;#24;#1024#500; # Number of points of the speed mesh
const xmin = -1.# -1.5#-1.0; # lower bound space mesh
const xmax = 1.#1.5#1.0; # upper bound space mesh
# The grid in velocity corresponds to the electronic one which contains the support of boths ions and electron density function
const vmin = -5.0 / sqrt(mu)# -2.#-10.#-5.0 / sqrt(mu); # lower bound speed mesh
const vmax = 5.0 / sqrt(mu)#2.#10.#5.0 / sqrt(mu); # upper bound speed mesh
const dx = (xmax - xmin) / Nx; # space step
const dv = (vmax - vmin) / Nv; # speed step
const CFL_x = 0.5 * dx / vmax;
const CFL_v = abs(mu) * dv / 10;

# Solver and other functionality choices
const init_file=0 # =0: initialization thanks to analytic expression, =1: initialization from file
const field_solver=1 # =0: Ampere solver, =1: Poisson solver E(Nx/2=0), =2: Poisson solver with current
# testcases are implemented in file initial_data.jl :
# 0: full 2-species, 1: cst 1-species validation, 2: time variable 1-species validation
const animate=0 # =0: no gif animation, =1: gif animation
const bar_prog=0 # =1: see the progression, else: nothing

# plot parameters
const savedir = "jlimages/"
const thecmap = "gnuplot"
const fislice = floor(Int, Nv):ceil(Int, Nv) # ions will be plotted on these colums


# Main function with final time in argument
function main(T)

    # Iteration number and time step respecting CFL
    Nt = (floor(Int, T / min(CFL_x, CFL_v)) + 1) # 'Int' argument for int output
    dt = T / Nt
    
    println("Welcome in main.")

    println("Parameters : ")
    println("nu = $nu, lambda = $lambda, mu = $mu")
    println("T = $T, Nt = $Nt, dt = $dt (CFL_x = $CFL_x, CFL_v = $CFL_v)")
    println("[xmin,xmax] = [$xmin,$xmax], Nx = $Nx, dx = $dx")
    println("[vmin,vmax] = [$vmin,$vmax], Nv = $Nv, dv = $dv")

    # The discretized quantities 
    xx = collect(LinRange(xmin, xmax, Nx + 1)) # vector (vertical) of spatial grid points (shared by e and i)
    vv = collect(LinRange(vmin, vmax, Nv + 1)) # vector (vertical) of speeds (shared by e and i)
    vv_plus = vv .* (vv .> 0.0)
    vv_minus = vv .* (vv .< 0.0) # signed positive and negative parts

    #############################################################
    # Initialisation with analytic expression or datas from files
    # Initialisation from analytical expression (and array declaration)
    EE = E0(xx)         # electric field (has to be initialized for Ampere solver)
    EEtemp = E0(xx)
    fi = fi_0(xx, vv)    # ion density 
    fe = fe_0(xx, vv)    # electron density
    rho = 0.0 * xx        # charge density (will be computed later)
  
    # Initialization of fi and fe from a data file 
    # WARNING!!! Files have to contain one value per line, we could probably
    # read datas presented in a different way, but it is not implemented
    # one value per line means : for x=0 values of f_{i/e} for all v, then for x=dx, etc.
    if init_file == 1
        # opening data files
        # be sure there are the right ones
        file_e = open("fe_phase_space0.dat", "r")
        file_i = open("fi_phase_space0.dat", "r")
        file_EE = open("E_space0.dat", "r")

        # reading datas and initializing
        for i in 1:Nx+1
            for j in 1:Nv+1
                
                res=readline(file_e)
                fe[i,j]=parse(Float64,res)
                res=readline(file_i)
                fi[i,j]=parse(Float64,res)
                
            end
            res=readline(file_EE)
            EE[i]=parse(Float64,res)

            # comments: if the data files have void lines separating x=0, x=dx, etc.
            #res=readline(file_e)
            #res=readline(file_i)

        end
        close(file_e)
        close(file_i)
        close(file_EE)
    end

    # store initial distribution functions (for diagnostics purpose)
    fi_init=fi*1.
    fe_init=fe*1.         
    
    # WARNING! not the right way to initialize E for the Ampere solver...
    # but works if the analytic E(t=0,x) is such that E(0,xmin)=-E(0,xmax)
    # better way: implement E(t=0,x) in the initial_data.jl file (done for some testcases only)
    #=
    if field_solver==0
        # (electronic charge density and) total mass
        rho .= vec(sum(fi .- fe, dims=2)) * dv
        # symmetrization of rho
        for i in 1:floor(Int, (Nx-1)/2) 
            rho[i] = (rho[i] + rho[Nx-i])*0.5;
            rho[Nx-i] = rho[i];
        end
        # compute current at boundaries : rectangle integration of \int_{v} v * (fi(t,+-1,v) - fe(t,+-1,v)) dv
        @views J_l = dv * sum(vv .* (fi[begin, :] - fe[begin, :]))
        @views J_r = dv * sum(vv .* (fi[end, :] - fe[end, :]))
        
        # update electric field : use formula E(x) = 0.5*(E(1) + E(-1)) + nu/2 * ( \int_{-1}^{x} rho(y)dy - \int_{x}^{1} rho(y)dy )
        # with E(t^{n+1},1) + E(t^{n+1},-1) approximated by E(t^n,1) + E(t^n,-1) + dt/lambda^2 * (J_l + J_r) (Ampère boundary conditions)
        EE .= 0.5 * (EE[begin] + EE[end] + dt / (lambda^2) * (J_l + J_r)) # the +- Mass_e cancels out (?????)
        @views EE[begin+1:end] .+= 0.5 / (lambda^2) * 0.5 * dx .* cumsum(rho[1:end-1] .+ rho[2:end])                        # add forward  trapeze integral
        @views EE[begin:end-1] .-= 0.5 / (lambda^2) * 0.5 * dx .* cumsum(rho[end-1:-1:1] .+ rho[end:-1:2])[end:-1:begin]    # add backward trapeze integral
        
    end
    =#
    #############################################################

    println("init done")
    
    # preparing time loop
    t = 0.0
    iplot = 0

    # Boundary conditions (all 0, never updated afterwards)
    fi[begin, :] .*= 0.0
    fi[end, :] .*= 0.0 # speed distribution is almost 0 
    fi[:, begin] .*= 0.0
    fi[:, end] .*= 0.0 # non-emmiting wall
    fe[begin, :] .*= 0.0
    fe[end, :] .*= 0.0 # speed distribution is almost 0 
    fe[:, begin] .*= 0.0
    fe[:, end] .*= 0.0 # non-emmiting wall

    # plots at time 0
    title = @sprintf(" at t=%7.3f", t)
    contourf(xx[2:end-1], vv[2:end-1], fi[begin+1:end-1, begin+1:end-1]', colormap=thecmap, title="Ions"*title, xlabel="x", ylabel="v", linewidth=0)
    png("jlimages/fi$iplot.png")
    contourf(xx[2:end-1], vv[2:end-1], fe[begin+1:end-1, begin+1:end-1]', colormap=thecmap, title="Electrons"*title, xlabel="x", ylabel="v", linewidth=0)
    png("jlimages/fe$iplot.png")
    plot(xx, EE, legend=false, title="E"*title, xlabel="x", ylabel="E")
    png("jlimages/E$iplot.png")
    plot(xx, rho, legend=false, title="rho"*title, xlabel="x", ylabel="rho")
    png("jlimages/rho$iplot.png")

    # nice functionalities (thanks to Pierre Navaro)
    # to see a progression bar in the execution terminal
    if bar_prog==1
        bar = Progress(Nt)
    end

    # animate=1 will construct a gif from plots
    if (animate==1)
        anim = Animation()
    end
    rho_results = Vector{Float64}[]

    ######################
    # Main temporal loop #
    ######################
    for n in 1:Nt
        if bar_prog==1
            next!(bar)
        end
        
        ######################################################
        # (electronic charge density and) total mass at time n
        #rho .= vec(sum(fi .- fe, dims=2)) * dv
        for i in 1:Nx+1
            rho[i]=0.
            for j in 1:Nv
                rho[i]+=0.5*dv*(fi[i,j]+fi[i,j+1]-fe[i,j]-fe[i,j+1])
            end
        end
        # symmetrization of rho
        #for i in 1:floor(Int, (Nx-1)/2) 
        #    rho[i] = (rho[i] + rho[Nx-i])*0.5;
        #    rho[Nx-i] = rho[i];
        #end
        ######################################################

        ####################################
        # Electric field evolution
        if field_solver==0
        # first way of computing EE: Ampere
        # --> we get E at time n+1  

            # compute current at boundaries : rectangle integration of \int_{v} v * (fi(t,+-1,v) - fe(t,+-1,v)) dv
            @views J_l = dv * sum(vv .* (fi[begin, :] - fe[begin, :]))
            @views J_r = dv * sum(vv .* (fi[end, :] - fe[end, :]))
    
            # update electric field : use formula E(x) = 0.5*(E(1) + E(-1)) + nu/2 * ( \int_{-1}^{x} rho(y)dy - \int_{x}^{1} rho(y)dy )
            # with E(t^{n+1},1) + E(t^{n+1},-1) approximated by E(t^n,1) + E(t^n,-1) + dt/lambda^2 * (J_l + J_r) (Ampère boundary conditions)
            EE .= 0.5 * (EE[begin] + EE[end] + dt / (lambda^2) * (J_l + J_r)) # the +- Mass_e cancels out (?????)
            @views EE[begin+1:end] .+= 0.5 / (lambda^2) * 0.5 * dx .* cumsum(rho[1:end-1] .+ rho[2:end])                        # add forward  trapeze integral
            @views EE[begin:end-1] .-= 0.5 / (lambda^2) * 0.5 * dx .* cumsum(rho[end-1:-1:1] .+ rho[end:-1:2])[end:-1:begin]    # add backward trapeze integral
            
        elseif field_solver==1
        # second way of computing EE update electric field with Poisson equation and E(Nx/2) fixed to 0
        # --> we get E at time n
        	
            i0 = floor(Int, Nx/2)+1
            if 2*i0==Nx+2
                # Boundary condition
                #EEtemp[i0]=0.
                EE[i0]=0.
                # Poisson solved on two subdomains
                #EEtemp[begin:end]=EE[begin:end]
                for i in i0+1:Nx+1
                    EE[i]=EE[i-1]+ 0.5 * dx / (lambda*lambda) * (rho[i] + rho[i-1])
                end
                for i in i0-1:-1:1
                    EE[i]=EE[i+1]- 0.5 * dx / (lambda*lambda) * (rho[i+1] + rho[i])
                end
            else 
                println("Warning!!! 0 is not a mesh point...")
            end

        else
        # third way of computing EE update electric field with Poisson equation and E(Nx/2) fixed to 0
        # --> we get E at time n
            @views J_l = dv * sum(vv .* (fi[begin, :] - fe[begin, :]))
            @views J_r = dv * sum(vv .* (fi[end, :] - fe[end, :]))
            Mass_e = 0.
            for i in 1:Nx
                for j in 1:Nv
                    Mass_e+=fe[i,j]*dv*dx
                end
            end
            E_xmax = EE[end] + 0.5 * (nu * dt/(lambda*lambda)) * Mass_e  - (dt/(lambda*lambda)) * J_r;
            E_xmin = EE[begin]  - 0.5 * (nu * dt/(lambda*lambda)) * Mass_e  - (dt/(lambda*lambda)) * J_l;
            
            for i in 1:Nx+1
                mass_rho_plus= 0; 
                mass_rho_minus= 0;
                # Quadrature for int_{-1}^{x_i} rho(x)dx.
                for k in 1:i-1
                    mass_rho_minus += 0.5*dx * (rho[k+1] + rho[k])
                end
                # Quadrature for int_{x_i}^{1} rho(x)dx
                for  k in i:Nx
                    mass_rho_plus += 0.5 *dx * (rho[k+1]+rho[k])
                end
                EE[i] = (1/(2*lambda*lambda)) * (mass_rho_minus - mass_rho_plus) + 0.5 * (E_xmax + E_xmin)
            end

        end
        ####################################
 

        # plots every Nplots iterations (and at initial and final times)
        Nplots=Nt/100 #10
        if ((mod(n+1,floor(Int,Nt/Nplots))==0) || (n==Nt) || (n==1))
            println("Iteration $n / $Nt, time t = ", dt * (n-1))
            title = @sprintf(" at t=%7.3f", t)

			open("jlimages/rho_$iplot.txt", "w") do io
    			writedlm(io, [xx rho EE])
			end


            #contourf(xx[2:end-1],vv[fislice],fi[begin+1:end-1,fislice]',colormap=thecmap,title="Ions"*title,xlabel="x",ylabel="v",linewidth=0)
            contourf(xx[2:end-1],vv[2:end-1],fi[begin+1:end-1,begin+1:end-1]',colormap=thecmap,title="Ions"*title,xlabel="x",ylabel="v",linewidth=0)
            png("jlimages/fi$iplot.png")
            save("jlfiles/fi$iplot.jld", "fi", fi)
            contourf(xx[2:end-1],vv[2:end-1],fe[begin+1:end-1,begin+1:end-1]',colormap=thecmap,title="Electrons"*title,xlabel="x",ylabel="v",linewidth=0);
            png("jlimages/fe$iplot.png")
            save("jlfiles/fe$iplot.jld", "fe", fe)
            plot(xx,EE,legend=false,title="E"*title,xlabel="x",ylabel="E");
            png("jlimages/E$iplot.png")
            save("jlfiles/E$iplot.jld", "E", EE)
            plot(xx,rho,legend=false,title="rho"*title,xlabel="x",ylabel="rho");
            png("jlimages/rho$iplot.png")
            save("jlfiles/rho$iplot.jld", "rho", rho)
            iplot += 1;
        end

        # for a gif animation of rho
        if animate ==1 && (n == 1 || mod(n,Nplots) == 0)
            title = @sprintf(" at t=%7.3f", t)
            p = plot(xx, rho, legend=false, 
            #title="ρ"*title, xlabel="x", ylabel="rho", xlim=[-0.8, 0.8],ylim=[-0.05, 0.5])
            title="ρ"*title, xlabel="x", ylabel="rho")
            frame(anim, p)
            push!(rho_results, copy(rho))
        end
        
        t = n * dt

        ####################################
        # Update ion and electron densities
        # dt fi + v dx fi + E dv fi = nu * fe
        # upwind coefficient, depending on field sign (multiplied by mu)
        if mu>0.
            EE_plus = max.(EE, 0.0)
            EE_minus = min.(EE, 0.0)
        else
            EE_plus = min.(EE, 0.0)
            EE_minus = max.(EE, 0.0)
        end

        # time evolution
        @views fi[2:Nx, 2:Nv] .+= (
            -(dt / dx) .* (vv_plus[2:Nv]' .* (fi[2:Nx, 2:Nv] .- fi[1:Nx-1, 2:Nv]) .+ vv_minus[2:Nv]' .* (fi[3:Nx+1, 2:Nv] .- fi[2:Nx, 2:Nv]))
            .-
            (dt / dv) .* (EE_plus[2:Nx] .* (fi[2:Nx, 2:Nv] .- fi[2:Nx, 1:Nv-1]) .+ EE_minus[2:Nx] .* (fi[2:Nx, 3:Nv+1] .- fi[2:Nx, 2:Nv]))
            .+
            dt * nu .* fe[2:Nx, 2:Nv])
        # dt fe + v dx fe - 1/mu E dv fe = 0
        @views fe[2:Nx, 2:Nv] .+=  (
            -(dt / dx) .* (vv_plus[2:Nv]' .* (fe[2:Nx, 2:Nv] .- fe[1:Nx-1, 2:Nv]) .+ vv_minus[2:Nv]' .* (fe[3:Nx+1, 2:Nv] .- fe[2:Nx, 2:Nv]))
            .+
            (dt / (dv * mu)) .* (EE_minus[2:Nx] .* (fe[2:Nx, 2:Nv] .- fe[2:Nx, 1:Nv-1]) .+ EE_plus[2:Nx] .* (fe[2:Nx, 3:Nv+1] .- fe[2:Nx, 2:Nv])))
        ####################################
       
   
    end
    #catch e
    #println(e.msg) # just so that ctrl+c works

    ###############################
    # diagnostics, data files, etc.
    # constructing gif animation
    if animate ==1
        gif(anim; loop=1)
    end

    # saving x, v, fe at final time (for gnuplot vizualisation) 
    file_e = open("fe_phase_spaceFinal.dat", "w")
    for i in 1:Nx+1
        for j in 1:Nv+1
            println(file_e,string(xx[i]), " ",string(vv[j])," ",string(fe[i,j])," ")
        end
        println(file_e," ")
    end
    close(file_e)

    # diagnostics, 1 species case with exact solution
    if testcase==2 || testcase==3
        # exact E and fe
        Eex=E_ex(xx)
        feex=fe_ex(xx,vv)

        # saving error on fe at final time (for gnuplot vizualisation) 
        file_err = open("fex-fe_phase_spaceFinal.dat", "w")
        for i in 1:Nx+1
            for j in 1:Nv+1
                println(file_err,string(xx[i]), " ",string(vv[j])," ",string(abs(feex[i,j]-fe[i,j]))," ")
            end
            println(file_err," ")
        end
        close(file_err)

        # L1 norm on electric field
        L1sum=0.
        for i in 1:Nx-1
            L1sum=L1sum+dx*abs(Eex[i]-EE[i])
        end
        println("dt, dx, dv, L1 error on E ",dt, " ", dx, " ", dv, " ",L1sum)
    end

    # Datas saved in files, ready to be read again
    # it let go further in time with the "init_file=1" parameter
    # opening data files
    file_fe_end = open("fe_for_continuation.dat", "w")
    file_fi_end = open("fi_for_continuation.dat", "w")
    file_E_end = open("E_for_continuation.dat", "w")
    for i in 1:Nx+1
        for j in 1:Nv+1
            println(file_fe_end,string(fe[i,j]))        
            println(file_fi_end,string(fi[i,j]))        
        end
        println(file_E_end,string(EE[i]))        
    end
    close(file_fe_end)
    close(file_fi_end)
    close(file_E_end)
    ###############################


    #return rho_results, fe_init, fi_init, fe, fi
end # try-catch

####################################
# Final time choice and call to main
T = 0.2 #1.#1.#0.1; # Time horizon

#@time rho_results, fe_init, fi_init, fe, fi = main(T)
@time main(T)
####################################


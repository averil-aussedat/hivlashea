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

using Plots

include("initial_data.jl")  # constant mu
                            # functions fi_0, fe_0, E0 

# Physical parameters
const nu     =  0.0; # Collision frequency
const lambda =  0.1; # Debye length
const T      =  0.5; # Time horizon

# Discretization parameters
const Nx     =  20; # Number of points of the space mesh
const Nv     =  40; # Number of points of the speed mesh
const xmin   = -1.0; # lower bound space mesh
const xmax   =  1.0; # upper bound space mesh
# The grid in velocity corresponds to the electronic one which contains the support of boths ions and electron density function
const vmin   = -10.0/sqrt(mu); # lower bound speed mesh
const vmax   =  10.0/sqrt(mu); # upper bound speed mesh
const dx     =  (xmax - xmin)/Nx; # space step
const dv     =  (vmax - vmin)/Nv; # speed step
const CFL_x  =  0.5*dx/vmax;
const CFL_v  =  mu * dv/10;
const Nt     =  floor(Int,T/min(CFL_x,CFL_v))+1; # 'Int' argument for int output
const dt     =  T/Nt;

println("Welcome in main.")

# """
#  Compute the electric field with the symmetric formula:
#  ```math
#  E(x) = 1/(2 Lambda^2) * (int_{-1}^{x} rho ds - int_{x}^{1} rho ds) + 0.5 (E(1) + E(-1))
#  ```
 
#  where
 
# - E^n+1(1) = E^n(1)   + (\nu dt/2 Lambda^2) * M_e^n   -(dt/Lambda^2) J^n(1),
# - E^n+1(-1) = E^n(-1) - (\nu dt/2 Lambda^2) * M_e^n   -(dt/Lambda^2) J^n(-1)
# """
# function compute_electric_field!(E,rho,J_l,J_r,Mass_e)
#     E[begin] += dt/(lambda^2) * (- nu * 0.5 * Mass_e + J_l);
#     E[end]   += dt/(lambda^2) * (  nu * 0.5 * Mass_e + J_r);
#     forward_int = 0.5 *  dx .* cumsum(rho[1: 1:end-1] .+ rho[2: 1:end]);
#     backward_int = 0.5 * dx .* cumsum(rho[end-1:-1:1] .+ rho[end:-1:2])[end:-1:begin];
#     E[begin+1:end-1] =  0.5/(lambda^2) * (forward_int - backward_int) .+ 0.5*(E[begin]+E[end]);
# end 

# function main()   
    println("Parameters : ")
    println("nu = $nu, lambda = $lambda, mu = $mu")
    println("T = $T, Nt = $Nt, dt = $dt (CFL_x = $CFL_x, CFL_v = $CFL_v)")
    println("[xmin,xmax] = [$xmin,$xmax], Nx = $Nx, dx = $dx")
    println("[vmin,vmax] = [$vmin,$vmax], Nv = $Nv, dv = $dv")

    # The discretized quantities 
    xx = collect(LinRange(xmin,xmax,Nx+1)) # vector (vertical) of spatial grid points (shared by e and i)
    vv = collect(LinRange(vmin,vmax,Nv+1)) # vector (vertical) of speeds (shared by e and i)
    vv_plus = vv .* (vv .> 0.0); vv_minus = vv .* (vv .< 0.0) # signed positive and negative parts

    # Initialisation
    EE = E0(xx)         # electric field
    fi = fi_0(xx,vv)    # electron density 
    fe = fe_0(xx,vv)    # ion density
    
    try 
        # Main temporal loop
        for n in 1:Nt
            if (mod(n,floor(Int,Nt/20))==0)
                println("Iteration $n / $Nt, time t = ", dt * n)
                println("Difference fi : ", maximum(abs.(fi - fi_0(xx,vv))));
                println("Difference fe : ", maximum(abs.(fe - fe_0(xx,vv))));
            end

            # electronic charge density and total mass
            rho = vec(sum(fi.-fe,dims=2)) * dv  
            # Mass_e = sum(fe) * dx * dv              

            # compute current at boundaries : rectangle integration
            J_l = dv * sum(vv .* (fi[begin,:] - fe[begin,:]))
            J_r = dv * sum(vv .* (fi[end  ,:] - fe[end  ,:]))

            # update electric field
            # EEmax = EE[begin] + dt/(lambda^2) * (- nu * 0.5 * Mass_e + J_l);
            # EEmin = EE[end]   + dt/(lambda^2) * (  nu * 0.5 * Mass_e + J_r);
            # EE .= 0.5 * (EEmax + EEmin)

            EE .= 0.5 * (EE[begin]+EE[end] + dt/(lambda^2) * (J_l + J_r)); # the +- Mass_e cancels out (?????)
            EE[begin+1:end] .+= 0.5/(lambda^2) * 0.5*dx.*cumsum(rho[1:end-1] .+ rho[2:end]);                        # add forward integral
            EE[begin:end-1] .-= 0.5/(lambda^2) * 0.5*dx.*cumsum(rho[end-1:-1:1] .+ rho[end:-1:2])[end:-1:begin];    # add backward integral

            # EE .*= 0.0; EE[begin] = 0.5 * EEmin; EE[end] = 0.5 * EE[end];
            # EE[begin+1:end] .+= 0.5/(lambda^2) .* 

            # EE[begin] += dt/(lambda^2) * (- nu * 0.5 * Mass_e + J_l);
            # EE[end]   += dt/(lambda^2) * (  nu * 0.5 * Mass_e + J_r);
            # forward_int = 0.5 *  dx .* cumsum(rho[1:end-1] .+ rho[2:end]);
            # backward_int = 0.5 * dx .* cumsum(rho[end-1:-1:1] .+ rho[end:-1:2])[end:-1:begin];
            # EE[begin+1:end-1] = 0.5/(lambda^2) * (forward_int - backward_int) .+ 0.5*(EE[begin]+EE[end]);
            # compute_electric_field!(EE,rho,J_l,J_r,Mass_e); 
            EE_plus = EE .* (EE .> 0.0); EE_minus = EE .* (EE .< 0.0)

            # Update ion and electron densities
            # dt fi + v dx fi + E dv fi = nu * fe
            fi[2:Nx,2:Nv] += (0.0
                .- (dt/dx) .* (vv_plus[2:Nv]' .* (fi[2:Nx,2:Nv].-fi[1:Nx-1,2:Nv]) .+ vv_minus[2:Nv]' .* (fi[3:Nx+1,2:Nv].-fi[2:Nx,2:Nv]))
                .- (dt/dv) .* (EE_plus[2:Nx]  .* (fi[2:Nx,2:Nv].-fi[2:Nx,1:Nv-1]) .+ EE_minus[2:Nx]  .* (fi[2:Nx,3:Nv+1].-fi[2:Nx,2:Nv]))
                .+ dt * nu .* fe[2:Nx,2:Nv]);

            # dt fe + v dx fe - 1/mu E dv fe = 0
            fe[2:Nx,2:Nv] += (0.0
                .- (dt/dx)      .* (vv_plus[2:Nv]' .* (fe[2:Nx,2:Nv].-fe[1:Nx-1,2:Nv]) .+ vv_minus[2:Nv]' .* (fe[3:Nx+1,2:Nv].-fe[2:Nx,2:Nv]))
                .- (dt/(dv*mu)) .* (EE_minus[2:Nx] .* (fe[2:Nx,2:Nv].-fe[2:Nx,1:Nv-1]) .+ EE_plus[2:Nx]   .* (fe[2:Nx,3:Nv+1].-fe[2:Nx,2:Nv])));
        end # while n
    catch e
        println(e.msg) # just so that ctrl+c works
    end # try-catch

#     return rho,fi,fe;
# end


# @time rho,fi,fe = main()
# plot(rho)
# plot(fe)
p = plot(layout=(3,1))
contour!(p[1,1],fe') # transpose
contour!(p[2,1],fi') # transpose
plot!(p[3,1],EE)
println("Bye")
#=

""" 
    Compute the charge density the total eletronic mass :
    ```math
    rho(x) = \\int_{\\RR} (fi - fe)(x,v)dv;
    Me     = \\int_{-1}^{1} \\int_{\\RR} fe(x,v)dvdx
    ```
"""
function compute_rho(fi, fe)
    rho = vec(sum(fi.-fe,dims=2)) * dv
    Me = sum(fe) * dx * dv
    return rho, Me
end 



""" Compute the current at x = xmax and x = xmin 
     J(t,x) = \\int_{\\RR} (f_i - f_e)(t,x,v)vdv
"""
 
function compute_current_at_boundaries(fi, fe)
    J_r = J_l = 0
    for j=1:Nv+1 # change to sum
        J_r += dv * (fi[end,j] - fe[end,j]) * (vmin_e + j * dv);
        J_l += dv* (fi[begin,j] * (vmin_e + j * dv) - fe[begin,j] * (vmin_e + j *dv));
    end # for
    return J_l, J_r
end

# """
#  Compute the electric field with the symmetric formula:
#  ```math
#  E(x) = 1/(2 Lambda^2) * (int_{-1}^{x} rho ds - int_{x}^{1} rho ds) + 0.5 (E(1) + E(-1))
#  ```
 
#  where
 
# - E^n+1(1) = E^n(1)   + (\nu dt/2 Lambda^2) * M_e^n   -(dt/Lambda^2) J^n(1),
# - E^n+1(-1) = E^n(-1) - (\nu dt/2 Lambda^2) * M_e^n   -(dt/Lambda^2) J^n(-1)
# """
# function compute_electric_field! (E,rho,J_l,J_r,Mass_e)
#     E_xmin = E[begin]  + (dt/(lambda*lambda)) * (- nu * 0.5 * Mass_e + J_l);
#     E_xmax = E[end]    + (dt/(lambda*lambda)) * (  nu * 0.5 * Mass_e + J_r);
#     E .*= 0.0 # reset electric field without changing memory

#     mass_rho_minus = 0.0; # Quadrature for int_{-1}^{x_i} rho(x)dx.
#     mass_rho_plus  = 0.0; # Quadrature for int_{x_i}^{1} rho(x)dx
#     E[begin] = 0.5*E_xmin; E[end] = 0.5*E_xmax # will be changed later
#     for i=2:Nx+1
#         mass_rho_minus += 0.5 * dx * (rho[i-1]+rho[i])          #  forward integration
#         mass_rho_plus  += 0.5 * dx * (rho[Nx+3-i]+rho[Nx+2-i])  # backward integration 

#         E[i] += 1/(2*lambda^2) * mass_rho_minus + 0.5*E_xmin;
#         E[Nx+2-i] += -1/(2*lambda^2) * mass_rho_plus + 0.5*E_xmax;
#     end # for i
# end 

"""
    Performs one step of the upwind scheme for
    ```math
        \\partial_t f + v_x \\partial_x f + v_v \\partial_v f = \text{source}
    ```
    Parameters:
        * ff [in/out] : 2D array of values at time tn, augmented with (here null) boundary conditions
        * source [in] : 2D array of the source at time tn (same shape at ff)
        * v{x/v}_{minus/plus} [in] : 2D arrays of negative/positive part of speed along x/v (same shape as ff)
        
"""
function update_f! (ff, source, vx_plus, vx_minus, vv_plus, vv_minus)
    ff[begin+1:end-1,begin+1:end-1] += 
        - (dt/dx) * (vx_plus .* () + vx_minus .* ())
        - (dt/dv) * (vv_plus .* () + vv_minus .* ())
        + dt * source;
end # update_f!
 
""" Update the distribution function """
function update_fi!(fi_0, fe_0, E)
    new_fi = zeros(Nx+1,Nv+1);
    v_plus = 0.0; v_minus = 0.0;
    E_plus  =  E[begin] >= 0.0 ? E[begin] : 0.0
    E_minus =  E[end] <= 0.0 ? E[end] : 0.0

    E_plus = 0.0; E_minus = 0.0;
    #  Evole the discrete density in the interior domain
    for i in 2:Nx # was [1,Nx-1] in c++, so skip the first and last
        E_plus  = E[i] >= 0 ? E[i] : 0.0
        E_minus = E[i] <= 0 ? E[i] : 0.0
        for j=2:Nx # skip first and last
            v_plus  = vmin_e + (j-1) * dv >= 0.0 ? vmin_e + (j-1) * dv : 0.0
            v_minus = vmin_e + (j-1) * dv <= 0.0 ? vmin_e + (j-1) * dv : 0.0

            new_fi[i,j] = fi_0[i,j]
            - (dt/dx) * ( v_plus * (fi_0[i,j] - fi_0[i-1,j]) + v_minus * (fi_0[i+1,j] - fi_0[i,j]) )
            - (dt/dv) * ( E_plus * (fi_0[i,j] - fi_0[i,j-1]) + E_minus * (fi_0[i,j+1] - fi_0[i,j]) )
            + nu * dt * fe_0[i,j];
        end # for j
    end # for i

    #  Apply the boundary condition at x = xmin
    for j=2:Nv # skip the first and last
        v_plus  = min_e + (j-1) * dv >= 0.0 ? min_e + (j-1) * dv : 0.0
        v_minus = min_e + (j-1) * dv <= 0.0 ? min_e + (j-1) * dv : 0.0
        if min_e + (j-1) * dv > 0.0
            new_fi[begin,j] = 0.;
        else 
            new_fi[begin,j] = fi_0[begin,j]
            - (dt/dx) * ( v_minus *(fi_0[begin+1,j] - fi_0[begin,j]) )
            - (dt/dv) * ( E_plus * (fi_0[begin,j] - fi_0[begin,j-1]) + E_minus * (fi_0[begin,j+1] - fi_0[begin,j]) )
            + nu * dt * fe_0[1,j];
        end # if upper boundary (known)
    end # for j

    #  Apply the boundary condition at x = xmax
    for j=2:Nv # skip the first and last
        v_plus  = min_e + (j-1) * dv >= 0.0 ? min_e + (j-1) * dv : 0.0
        v_minus = min_e + (j-1) * dv <= 0.0 ? min_e + (j-1) * dv : 0.0
        if min_e + (j-1) * dv < 0.0
            new_fi[end,j] = 0.;
        else 
            new_fi[end,j] = fi_0[end,j]
            - (dt/dx) * ( v_plus * (fi_0[end,j] - fi_0[end-1,j]) )
            - (dt/dv) * ( E_plus * (fi_0[end,j] - fi_0[end,j-1]) + E_minus * (fi_0[end,j+1] - fi_0[end,j]) )
            + nu * dt * fe_0[1,j];
        end # if lower boundary (known)
    end # for j

    # the upper and lower rows stays, already initialized at 0
    return new_fi;
end
#=

# Udpate the distribution function

Matrix update_fe(Matrix & fe_0, Vector & E)
{
    Matrix new_fe(Nx+1,Nv+1);
    double v_plus, v_minus;
    double E_plus, E_minus;
    #  Evole the discrete density in the interior domain
    for(int i = 1 ; i < Nx ; i++)
    {
        E_plus =  0.5  * (E[i] + abs(E[i]));
        E_minus = 0.5  * (E[i] - abs(E[i]));
        for(int j =1 ; j < Nv; j++)
        {
            v_plus = 0.5 * (vmin_e + j * dv + abs(vmin_e + j *dv));
            v_minus = 0.5 * (vmin_e + j * dv - abs(vmin_e + j *dv));
            new_fe[i][j] = fe_0[i][j]
            - (dt/dx) * ( v_plus * (fe_0[i][j] - fe_0[i-1][j])  +v_minus *(fe_0[i+1][j] - fe_0[i][j]) )
            - (dt/dv) * ( E_plus * (fe_0[i][j] - fe_0[i][j-1]) + E_minus * (fe_0[i][j+1] - fe_0[i][j]) );
            
                                                
        }
    }
    
    #  Apply the boundary condition at x = 0
    for(int j = 0; j < Nv  ; j++)
    {
        v_plus =  0.5 * (vmin_e + j * dv + abs(vmin_e + j *dv));
        v_minus = 0.5 * (vmin_e + j * dv - abs(vmin_e + j *dv));
        if( vmin_e + j*dv > 0)
        {
            new_fe[0][j] = 0.;
        }
        else
        {
            E_plus =  0.5  * (E[0] + abs(E[0]));
            E_minus = 0.5  * (E[0] - abs(E[0]));
            new_fe[0][j] = fe_0[0][j]
            - (dt/dx) * ( v_minus *(fe_0[1][j] - fe_0[0][j]) )
            - (dt/dv) * ( E_plus * (fe_0[0][j] - fe_0[0][j-1]) + E_minus * (fe_0[0][j+1] - fe_0[0][j]) );
            
        }
    }
    #  Apply the boundary condition at  x= 1
    for(int j = 0; j < Nv ; j++)
    {
        v_plus = 0.5 * (vmin_e + j * dv + abs(vmin_e + j *dv));
        v_minus = 0.5 * (vmin_e + j * dv - abs(vmin_e + j *dv));
        if( vmin_e + j*dv < 0)
        {
            new_fe[Nx][j] = 0.;
        }
        else
        {
            E_plus =  0.5  * (E[Nx] + abs(E[Nx]));
            E_minus = 0.5  * (E[Nx] - abs(E[Nx]));
            new_fe[Nx][j] = fe_0[Nx][j]
            - (dt/dx) * ( v_plus * (fe_0[Nx][j] - fe_0[Nx-1][j])  )
            - (dt/dv) * ( E_plus * (fe_0[Nx][j] - fe_0[Nx][j-1]) + E_minus * (fe_0[Nx][j+1] - fe_0[Nx][j]) );
            
        }
        
        #  Apply the boundary condition at v = v_max and v = vmin
        for(int i = 1 ; i < Nx ; i++)
        {
            
            
            new_fe[i][Nv] = 0;
            new_fe[i][0]  = 0;
         }
    }
    

    return new_fe;
}

=#

=#
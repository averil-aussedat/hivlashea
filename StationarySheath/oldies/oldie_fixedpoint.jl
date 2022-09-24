###########################################################################
# Test code for fixed-point algorithm
# Vlasov-Poisson with nonperidic boundary conditions
# CEMRACS 2022 
###########################################################################

using LinearAlgebra 
using Plots
using DelimitedFiles # writedlm, readdlm

################################
# Inputs (in conf file)
################################

########## Parameters ##########
# physics
const lambda    = 0.1; # Debye length
# const mu        = 1.0/3672; # ratio mass_electron / mass_ion
const mu        = 0.01; # ratio mass_electron / mass_ion
const nu        = 42.0; # collision frequency

# const xmax = 1.0;
fact = 1.0; # the largest fact, the tighest the distributions
# electronic distribution at x=0 (expression and latex-printable description)
f_eb(v) = sqrt(mu*fact) .* exp.(-mu * fact * (v.^2)/2) ./ sqrt(2pi); f_eb_str = "\$ $fact \\mu \\frac{- $fact \\mu \\frac{v^2}{2}}{2 \\pi}\$";
# f_eb(v) = 1.0.*(abs.(v).<=1.0); f_eb_str = "\$1_{\left\{|v|\leqslant 1\right\}}\$";
# f_eb(v) = max.(1.0.+v,0.0); f_eb_str = "\$(1+v)_{+}\$";
# f_eb(v) = 0.0005 .* max.(0.0, v.^2 .* (100.0 .- v.^2)); f_eb_str = "\$5e-4\\times (v^2 (100-v^2))_{+}\$"; 
# f_eb(v) = exp.(-mu .* (abs.(v) .- 2.0).^2 ./(2*fact)) ./ sqrt(mu * 2pi*fact); f_eb_str = "\$\\frac{- \\mu \\frac{(|v|-2)^2}{2\\times $fact}}{$fact \\times 2 \\pi \\mu}\$";
# f_eb(v) = 0.01 .* v.^2 .* exp.(-mu * fact .* (abs.(v) .- 10.0).^2 ./2) ./ sqrt(mu * 2pi*fact); f_eb_str = "\$v^2\\frac{- $fact \\mu \\frac{(|v|-10)^2}{2}}{$fact \\times 2 \\pi \\mu}\$";
# f_eb(v) = 1.0 .* (7 .<= abs.(v) .<= 12); f_eb_str = "\$1_{\\left\\{7<|v|<12\\right\\}}\$";

# numerics
const cardinal  = 5 # number of basis functions for phi
# const basis     = "canonical"; # choice of representation
const basis     = "pair"; # choice of representation
# Nx              = 100; # number of mesh points in [0,1]
Nit_max         = 300; # maximal number of fixed-point iterations

# expectations on the solution
const ve_bound = 50.0; # electron speed ∈ [-ve_bound, ve_bound]
const vi_bound = 10.0; #      ion speed ∈ [-vi_bound, vi_bound]

coeff = [(-30.0)*(i==1) for i=1:cardinal]; # initial coefficients in the basis

# lambda = 0.1, mu = 0.01, nu = 25, cardinal 3, base paire
# coeff = [-29.275193716990835, 9.970222191354972, -2.9013822611511237, 0.0, 0.0]
# lambda = 0.1, mu = 0.01, nu = 25, cardinal 5, base paire
# coeff = [-23.598907425704247, -32.53313797857284, 57.34577553580175, -43.941298025794595, 12.723603170563207] # 200 it
# coeff = [-23.506721799618767, -32.56394245948412, 57.395524709951175, -43.984835123067164, 12.73761089330094] # 400 it
# coeff = [-23.502145358576897, -32.565589032948814, 57.39814412865666, -43.98709531643629, 12.738332411753726] # 600 it
# coeff = [-23.502136645441738, -32.565586352306724, 57.398140716014325, -43.98709365154615, 12.738332128080874] # 1600 it

# coeff = [-23.502150074610846, -32.56558316308049, 57.39813556967014, -43.98708889237952, 12.738330543533573] # 600 it, test case 256
# coeff = [-23.224480416200326, -33.38923429223604, 58.73999579734637, -45.097067780949935, 13.09492160234383] # 600 it, test case 512
# coeff = [-23.006121534107088, -34.06829975649596, 59.81915146207058, -45.95973509844458, 13.363114940523264] # 400 it, test case 1024... unreachable.

# coeff = [-22.60208623599732, -35.15117074584164, 61.10948126984762, -46.661726647373, 13.490994216825257] # avec mehdi's integral, 256 

################ lambda = 0.1, mu = 0.01, nu = 42, cardinal 5, base paire, maxwellienne
coeff = [-42.74868771747942, -10.323119689802478, 31.225554184355076, -27.10895958278373, 8.33076483840167] # 256
# coeff = [-42.66063631441532, -10.711261741127007, 31.9741690031144, -27.788216769478087, 8.561672720273881] # 512
# coeff = [-42.647067207473654, -10.814841942994336, 32.21592060277268, -28.03279554275593, 8.650916996617742] # 1024

################ lambda = 0.1, mu = 0.01, nu = 42, cardinal 5, base paire, nulle en 0


struct Polybase
    "Number of elements in the basis"
    cardinal    :: Int64 
    "Chosen basis"
    basis       :: String 
    "Array of Horner coefficients for x -> p(x)"
    evals       :: Vector{Vector{Float64}} 
    "Array of Horner coefficients for x -> p'(x)"
    evals_dx    :: Vector{Vector{Float64}} 
    # "Array of Horner coefficients for x -> p''(x)"
    # evals_dxx   :: Vector{Vector{Float64}} 
    "Matrix of the linear system for the Poisson solver"
    Vandermonde :: Matrix{Float64} 
end

mutable struct Phi_data
    "Basis for polynomial representation"
    polybase :: Polybase;      
    "Coefficients in the representation system" 
    coeff :: Vector{Float64};  
end

################################
# Small functions (to mooooove)
################################

function trapezes(values, dx)
    return dx * (0.5 * values[1] + sum(values[2:end-1]) + 0.5 .* values[end])
end

""" Variable step """
function trapezes_inh(values, mesh)
    return 0.5 * sum((values[1:end-1] .+ values[2:end]) .* (mesh[2:end] .- mesh[1:end-1]))
end

"""
    Returns the negative v such that
     v²/2 - 1/μ φ(x) = v₀²/2 - 1/μ φ(x₀)     
"""
function get_v_char_e(v0, phix0, phix)
    val = v0.^2 .- 2.0/mu .* (phix0 .- phix)
    if (any(val .< -1e-10))
        # throw(DomainError("[get_v_char_e] char issued from v=$v0 produces val=$val."))
        throw(DomainError("negative value in get_v_char_e"))
    else
        res = - sqrt.(max.(0.0,val));
    end
    return res
end

"""
    Returns the negative v such that
     v²/2 + φ(x) = v₀²/2 + φ(x₀)
"""
function get_v_char_i(v0, phix0, phix)
    val = v0.^2 .+ 2.0.*(phix0 .- phix)
    if (any(val .< -1e-10))
        # throw(DomainError("[get_v_char_i] char issued from v=$v0 produces val=$val."))
        throw(DomainError("negative value in get_v_char_i"))
    else
        res = - sqrt.(max.(0.0,val));
    end
    return res
end

"""
    Newton descent to find the x such that
     v²/2 - 1/μ φ(x) = v₀²/2 - 1/μ φ(x₀) 
"""
function get_x_char_e(v0, x0, phix0, vbar, phi_data)
    # initialization with the intersection of v=vbar and the tangent to the characteristic crossing (x,v)
    phixprime = eval_phi_dx(phi_data, x0)
    x = abs(phixprime)>1e-8 ? x0 + mu * v0*(vbar - v0)/phixprime : x0;
    Jx = vbar.^2/2 -1/mu * eval_phi(phi_data, x) - (v0^2/2 -1/mu * phix0) # we minimize Jx
    k=0; 
    while ((k<=10) && (abs(Jx)>1e-14))
        Jxprime = - 1/mu * eval_phi_dx(phi_data, x);
        x = abs(Jxprime)>1e-8 ? x - Jx / Jxprime : x;
        # x -= Jx / Jxprime;
        Jx = vbar.^2 /2.0 -1.0/mu * eval_phi(phi_data, x) - (v0^2 /2.0 -1.0/mu * phix0);
        # println("iteration $k, Jx = $Jx, x = $x")
        k += 1;
    end
    return x
end

"""
    Newton descent to find the x such that
    vbar²/2 + φ(x) = v₀²/2 + φ(x₀) 
"""
function get_x_char_i(v0, x0, phix0, vbar, phi_data, x_init)
    # initialization with the intersection of v=vbar and the tangent to the characteristic crossing (x,v)
    if (x_init < 0.0) # if no particular init
        phixprime = eval_phi_dx(phi_data, x0)
        x = abs(phixprime)>1e-8 ? x0 - v0*(vbar - v0)/phixprime : x0 + 1e-2;
    else
        x = x_init;
    end
    Jx = (v0^2/2 + phix0) - (vbar.^2/2 + eval_phi(phi_data, x)) # we minimize Jx
    k=0; 
    while ((k<=10) && (abs(Jx)>1e-14))
        Jxprime = - eval_phi_dx(phi_data, x);
        x = abs(Jxprime)>1e-8 ? x - Jx / Jxprime : x;
        Jx = (v0^2 /2.0 + phix0) - (vbar.^2 /2.0 + eval_phi(phi_data, x));
        # println("iteration $k, Jx = $Jx, x = $x")
        k += 1;
    end
    return x
end

# @doc raw"""
#     Returns ``x`` such that
#     ```math
#         \frac{\overline{v}^2}{2} - \frac{1}{\mu} \phi(x) = \frac{v_0^2}{2} - \frac{1}{\mu} \phi(x_0)
#     ```
#     using Newton descent.
# """
# function get_x_char_e(phi::Phidata, params::Params, v0, x0, phix0, vbar, xinit=[])
#     if (length(xinit) == 0)
#         # initialization with the intersection of v=vbar and the tangent to the characteristic crossing (x,v)
#         phixprime = eval_phi_dx(phi, params, x0)
#         x = abs(phixprime)>1e-8 ? x0 + mu * v0*(vbar - v0)/phixprime : x0
#     else
#         x = xinit
#     end
#     Jx = vbar.^2/2 - 1/mu * eval_phi(phi, params, x) - (v0^2/2 - 1/mu * phix0) # we minimize Jx
#     k=0
#     while ((k<=10) && (abs(Jx)>1e-14))
#         Jxprime = - 1/mu * eval_phi_dx(phi, params, x)
#         x = abs(Jxprime)>1e-8 ? x - Jx / Jxprime : x
#         Jx = vbar.^2/2 - 1/mu * eval_phi(phi, params, x) - (v0^2/2 - 1/mu * phix0)
#         k += 1
#     end
#     return x
# end

# @doc raw"""
#     Returns ``x`` such that
#     ```math
#         \frac{\overline{v}^2}{2} + \phi(x) = \frac{v_0^2}{2} + \phi(x_0)
#     ```
#     using Newton descent.
# """
# function get_x_char_i(phi::Phidata, params::Params, v0, x0, phix0, vbar, xinit=[])
#     if (length(xinit) == 0)
#         # initialization with the intersection of v=vbar and the tangent to the characteristic crossing (x,v)
#         phixprime = eval_phi_dx(phi, params, x0)
#         x = abs(phixprime)>1e-8 ? x0 - v0*(vbar - v0)/phixprime : x0 + 1e-2
#     else
#         x = xinit
#     end
#     Jx = (v0^2/2 + phix0) - (vbar.^2/2 + eval_phi(phi, params, x)) # we minimize Jx
#     k=0
#     while ((k<=10) && (abs(Jx)>1e-14))
#         Jxprime = eval_phi_dx(phi, params, x)
#         x = abs(Jxprime)>1e-8 ? x - Jx / Jxprime : x
#         Jx = (v0^2/2 + phix0) - (vbar.^2/2 + eval_phi(phi, params, x))
#         k += 1
#     end
#     return x
# end

################################
# Polynomial tools
################################

# wrapper
function eval_phi(phi_data, x)
    return sum((phi_data.coeff .* [evaluate_poly(eval, x) for eval = phi_data.polybase.evals])')'
end

# wrapper
function eval_phi(polyb, coeff, x)
    return sum((coeff .* [evaluate_poly(eval, x) for eval = polyb.evals])')'
end

function eval_phi_dx(phi_data, x)
    return sum((phi_data.coeff .* [evaluate_poly(eval, x) for eval = phi_data.polybase.evals_dx])')'
end

function eval_phi_dxx(phi_data, x)
    return sum((phi_data.coeff .* [evaluate_poly(eval, x) for eval = phi_data.polybase.evals_dxx])')'
end

function evaluate_poly(eval, x)
    res = 0.0*x; # higher degree coefficient
    for ee=eval
        res = x .* res .+ ee;
    end
    return res;
end

# constructor
function Polybase(cardinal, basis, meshx)
    if (basis=="canonical")
        # p(x) = a_n + x (a_{n-1} + x (a_{n-2} + ... x (a_2 + a_1 x)))
        # starts from degree 2
        evals = [[1.0*(i==j-2) for j=cardinal+2:-1:1] for i=1:cardinal]
        evals_dx = [[1.0*j*(i==j-1) for j=cardinal+2:-1:1] for i=1:cardinal]
        evals_dxx = [[1.0*j*(j+1)*(i==j) for j=cardinal+2:-1:1] for i=1:cardinal]
    elseif (basis=="pair")
        evals = [[1.0*(i==j-2) for j=2*cardinal+2:-1:1] for i=1:2:2*cardinal]
        evals_dx = [[1.0*j*(i==j-1) for j=2*cardinal+2:-1:1] for i=1:2:2*cardinal]
        evals_dxx = [[1.0*j*(j+1)*(i==j) for j=2*cardinal+2:-1:1] for i=1:2:2*cardinal]
    else
        throw(ArgumentError("Unknown phi basis '$basis'."))
    end
    Vandermonde = zeros(length(meshx), cardinal)
    [Vandermonde[:,j] = evaluate_poly(eval_dxx,meshx) for (j,eval_dxx) = enumerate(evals_dxx)]
    return Polybase(cardinal, basis, evals, evals_dx, Vandermonde)
end 


################################
# Solvers for ni and ne
################################

"""
    Computes nₑ = ∫fₑ(x,v)dv for x in meshx.
    Uses a v mesh with roughly the same step as meshx.
"""
function update_ne!(ne, f_eb, phix)
    meshv = collect(LinRange(-ve_bound,0.0,length(phix)))
    for (i,phixi) in enumerate(phix)
        # use the characteristics to determine the values of fₑ
        fe_values = f_eb(get_v_char_e(meshv,phixi,0.0))
        # integral (using the assumption of symmetry over fₑ)
        ne[i] = trapezes(vcat(fe_values,fe_values[end-1:-1:begin]), meshv[2]-meshv[1])
    end
end

"""
    Integrates nu * f_e over the characteristic issued from (meshx[end], v_1) stopping at x=meshx[stop_index].
    The characteristic is meshed by (x_i,v_i), where x_i ∈ meshx and v_i is computed.
    Returns the integral (scalar) and the vector of v_i.
"""
function integrate_char(f_eb, meshx, phix, v_1, stop_index)
    integral = 0.0; 
    v_is = get_v_char_i(v_1, phix[end], phix[stop_index:end])
    fe_previous = f_eb(get_v_char_e(v_1,phix[end],0.0))

    for ix in length(meshx)-1:-1:stop_index
        dt = (meshx[ix]-meshx[ix+1]) * 2.0 / (v_is[ix-stop_index+1] + v_is[ix-stop_index+2]) # screw numbering from 1 -.-
        fe_current = f_eb(get_v_char_e(v_is[ix-stop_index+1],phix[ix],0.0))
        integral += dt * 0.5 * (fe_current + fe_previous)
        fe_previous = fe_current
    end # for ix
    return nu * integral, v_is
end

"""
    Integrates nu * f_e over the characteristic issued from (meshx[1], v_0).
    Uses the formula (when the characteristic intersects (x=0,v<0))
      fᵢ(x,v) = ∫[1,x] ν fₑ(a,v(a)) / v(a) da,   v(a) := - √(v² + 2(φ(x)-φ(a)))  
    or (when the characteristic intersects (x>=0,v=0))
      fᵢ(x,v) = ∫[vₘᵢₙ,0] ν fₑ(x(b),b) / (-∂ₓφ(b)) db,   x(b) := φ⁻¹(v²/2 + φ(x) - b²/2)  
"""
function integrate_char_mehdi(f_eb, meshx, phix, phidx, ix0, v_0)
    myzero = 1e-10
    if ((v_0^2)/2.0 + phix[ix0] > 0.0) # x is the reference variable
        v_is = min.(-myzero, get_v_char_i(v_0, phix[ix0], phix[ix0:end]))
        fes = f_eb(get_v_char_e(v_is,phix[ix0:end],0.0)) ./ v_is
        integral = nu * trapezes(fes, meshx[1] - meshx[2])
    else # v is the reference variable
        v_is = get_v_char_i(v_0, phix[ix0], phix[ix0:end])
        fes = f_eb(get_v_char_e(v_is,phix[ix0:end],0.0)) ./ max.(myzero,-phidx[ix0:end])
        integral = nu * trapezes_inh(fes, -v_is)
    end
    return integral, v_is # integration over [1,x] with x<1
end

"""
    Integrates nu * f_e over the characteristic issued from (meshx[ix0], v_0).
    Uses the formula
    fᵢ(x,v) = ∫[vₘᵢₙ,0] ν fₑ(x(b),b) / (-∂ₓφ(b)) db,   x(b) := φ⁻¹(v²/2 + φ(x) - b²/2)  
    THE v IS NEGATIVE
"""
function integrate_char_outofmeshx_mehdi(f_eb, meshx, phi_data, phix, phidx, ix0, v_0, ind)
    if ((v_0^2)/2.0 + phix[ix0] >= 0.0)
        # println("liap : ", (v_0^2)/2.0 + phix[ix0], ", v0 $v_0")
        char, = integrate_char_mehdi(f_eb, meshx, phix, phidx, 1, get_v_char_i(v_0,phix[ix0],phix[1])) 
    else
        alpha = (v_0                                   - get_v_char_i(0.0,phix[ind+1],phix[ix0])) / 
                (get_v_char_i(0.0,phix[ind],phix[ix0]) - get_v_char_i(0.0,phix[ind+1],phix[ix0]))
        x_init = alpha * meshx[ind] + (1.0-alpha) * meshx[ind+1]
        # println("aa 3")
        x_targ = get_x_char_i(v_0,meshx[ix0], phix[ix0], 0.0, phi_data, x_init)

        newmeshx = collect(LinRange(x_targ, 1.0, ceil(Int,2+(1.0-x_targ)*length(meshx))))
        newphix  = eval_phi(phi_data, newmeshx)
        newphidx = eval_phi_dx(phi_data, newmeshx)
        char, = integrate_char_mehdi(f_eb, newmeshx, newphix, newphidx, 1, 0.0)

        # if (ind<length(phix))
        #     # we know that x_target ∈ [meshx[ind], meshx[ind+1][
        #     # println("aa 1")
        #     char1, = integrate_char_mehdi(f_eb, meshx, phix, phidx, ind+1, get_v_char_i(v_0,phix[ix0],phix[ind+1])) # first part
        #     # determine x_target
        #     # println("aa 2")
        #     # println("ind  $ind, ix0 : $ix0")
        #     alpha = (v_0                                   - get_v_char_i(0.0,phix[ind+1],phix[ix0])) / 
        #             (get_v_char_i(0.0,phix[ind],phix[ix0]) - get_v_char_i(0.0,phix[ind+1],phix[ix0]))
        #     x_init = alpha * meshx[ind] + (1.0-alpha) * meshx[ind+1]
        #     # println("aa 3")
        #     x_targ = get_x_char_i(v_0,meshx[ix0], phix[ix0], 0.0, phi_data, x_init)
        #     phixtarg = eval_phi(phi_data,x_targ)
        #     # println("aa 4")
        #     v_prec = get_v_char_i(v_0,phix[ix0],phix[ind+1])
        #     fe_targ = f_eb(get_v_char_e(0.0,   phixtarg,0.0)) / (-eval_phi_dx(phi_data, x_targ))
        #     fe_prec = f_eb(get_v_char_e(v_prec,phix[ind+1],0.0)) / (-phidx[ind+1])
        #     # println("alpha : $alpha, x_targ : $x_targ, norm to init : $(abs(x_targ-x_init)), char1 : $char1")
        #     char = char1 + nu * (0 - v_prec) * 0.5 * (fe_targ + fe_prec)
        # else
        #     char = 0.0
        # end
    end
    return char
end

"""
    Integrates nu * f_e over the characteristic issued from (meshx[end], v_1) stopping at the boundary.
    Difficulty: the endpoint (located on v=0 or x=0) may not belong to meshx.
"""
function integrate_char_outofmeshx(f_eb, meshx, phi_data, phix, v_1, ind)
    if (v_1 <= get_v_char_i(0.0,0.0,phix[end])) # easy case : endpoint on x=0
        char, = integrate_char(f_eb, meshx, phix, v_1, 1)
    else
        # we know that x_target ∈ [meshx[ind], meshx[ind+1][
        char1, = integrate_char(f_eb, meshx, phix, v_1, ind+1) # first part
        if (ind==length(phix))
            char = char1
        else
            # determine x_target
            alpha = (v_1                                   - get_v_char_i(0.0,phix[ind+1],phix[end])) / 
                    (get_v_char_i(0.0,phix[ind],phix[end]) - get_v_char_i(0.0,phix[ind+1],phix[end]))
            x_init = alpha * meshx[ind] + (1.0-alpha) * meshx[ind+1]
            x_targ = get_x_char_i(v_1, meshx[end], phix[end], 0.0, phi_data, x_init)
            phixtarg = eval_phi(phi_data,x_targ)
            v_prec = get_v_char_i(v_1,phix[end],phix[ind+1])
            # print("xl ", meshx[ind], ", xr ", meshx[ind+1], ", x_init $x_init, x_targ $x_targ")
            dt = (x_targ - meshx[ind+1]) * 2.0 / (0.0 + v_prec)
            fe_targ = f_eb(get_v_char_e(0.0,   phixtarg,0.0))
            fe_prec = f_eb(get_v_char_e(v_prec,phix[ind+1],0.0))
            char = char1 + nu * dt * 0.5 * (fe_targ + fe_prec)
        end
    end
    return char
end

"""
    Given a space mesh, returns an adapted mesh of (x=1,v<=0).
    Outputs:
        * upper_mesh : Vector{Float64}, intersections of x=1 and the characteristics issued from meshx
        * lower_mesh : Vector{Float64}, starting of characteristics below upper_mesh
"""
function get_fi_char_meshes(phix)
    # we suppose that fi = 0 for v < lowerbound
    lowerbound = min(get_v_char_i(get_v_char_e(0.0,phix[end],0.0), 0.0, phix[end]),-vi_bound) 
    upper_mesh = get_v_char_i(0.0,phix,phix[end]) 
    N_lm = ceil(Int,(upper_mesh[1] - lowerbound)*Nx) # proportional to Nx and the length of the interval
    lower_mesh = collect(LinRange(lowerbound,upper_mesh[1],N_lm)) # ends at the charac. crossing (0.0,0.0)
    return upper_mesh[2:end], lower_mesh
end

function update_ni!(ni, f_eb, phix, meshx)
    ni .*= 0.0
    # get the meshes
    upper_mesh, lower_mesh = get_fi_char_meshes(phix)
    # lower mesh part : all coordinates of nᵢ are updated
    integral_char_previous, v_is_previous = integrate_char(f_eb, meshx, phix, lower_mesh[1], 1) # first characteristic's integral
    for vlm in lower_mesh[2:end]
        integral_char, v_is = integrate_char(f_eb, meshx, phix, vlm, 1) # stop_index still at 1
        ni .+= (v_is .- v_is_previous) .* (integral_char + integral_char_previous) # times 2 is embedded
        integral_char_previous = integral_char; v_is_previous .= v_is # update
    end # for v lower mesh

    # upper mesh part : some characteristics stops 
    for (stop_index, vum) in enumerate(upper_mesh)
        integral_char, v_is = integrate_char(f_eb, meshx, phix, vum, stop_index+1) # v_is is of size length(meshx)-stop_index now
        ni[stop_index+1:end] .+= (v_is .- v_is_previous[2:end]) .* (integral_char + integral_char_previous) # times 2 is embedded
        integral_char_previous = integral_char; v_is_previous = 1*v_is # update
    end # for v upper mesh
end

function update_ni_mehdi!(ni, f_eb, phi_data, phix, phidx, meshx)
    ni .*= 0.0

    # limit of the mesh of (x=0,v<0) that we do
    lowerbound = min(get_v_char_e(0.0,phix[end],0.0),-vi_bound) 
    N_lm = ceil(Int,(-lowerbound)*Nx) # proportional to the discretization of meshx
    lower_mesh = collect(LinRange(lowerbound,0.0,N_lm)) # CAREFUL this time, mesh of x=0,v<0 (unlike previously)

    # mesh of (x=0, v<=0) : all coordinates of nᵢ are updated
    integral_char_previous, v_is_previous = integrate_char_mehdi(f_eb, meshx, phix, phidx, 1, lowerbound) # first characteristic's integral
    for vlm in lower_mesh[2:end-1]
        integral_char, v_is = integrate_char_mehdi(f_eb, meshx, phix, phidx, 1, vlm) # stop_index still at 1
        ni .+= (v_is .- v_is_previous) .* (integral_char + integral_char_previous) # times 2 is embedded
        integral_char_previous = integral_char; v_is_previous = 1*v_is # update
    end # for v lower mesh

    # special case of (0,0)
    integral_char, v_is = integrate_char_mehdi(f_eb, meshx, phix, phidx, 2, get_v_char_i(0.0,0.0,phix[2])) 
    factor = sqrt(2*abs(phi_data.coeff[1]))
    dx = meshx[2]-meshx[1]
    dv = abs(v_is[1])
    mKsur2 = - atan(-dv/(factor * dx)) / 2.0
    localint = factor * dx * (log(sin(mKsur2)+cos(mKsur2)) - log(cos(mKsur2)-sin(mKsur2))) + dv * (log(sin(mKsur2)) - log(cos(mKsur2)))
    v_is = vcat([0.0],v_is)
    integral_char += localint * nu * f_eb(0.0) 

    # println("factor : $factor, x : $dx, dv : $dv, mKsur2 : $mKsur2, localint : $localint")

    # upper mesh part : some characteristics stops. We start from (x ∈ meshx[2:end], v = 0) 
    for ix0 in 2:length(meshx)
        integral_char, v_is = integrate_char_mehdi(f_eb, meshx, phix, phidx, ix0, 0.0) # v_is is of size length(meshx)-ix0 now
        ni[ix0:end] .+= (v_is .- v_is_previous[2:end]) .* (integral_char + integral_char_previous) # times 2 is embedded
        integral_char_previous = integral_char; v_is_previous = 1*v_is # update
    end # for v upper mesh
end

################################
# Poisson solver
################################

"""
    Solves the problem 
        - λ² Δφ(x) = nᵢ(x) - nₑ(x)
              φ(0) = 0
            ∂ₓφ(0) = 0
    with φ decomposed in a polynomial basis.
"""
function update_phi!(phi_data, ni, ne)
    phi_data.coeff .= - 1.0/(lambda^2) .* (phi_data.polybase.Vandermonde \ (ni - ne)) # least squares if Vandermonde is not square
end

"""
    Returns the solution of the problem 
        - λ² Δφ(x) = nᵢ(x) - nₑ(x)
              φ(0) = 0
            ∂ₓφ(0) = 0
    with phi decomposed in a polynomial basis.
"""
function update_phi(phi_data, ni, ne)
    return - 1.0/(lambda^2) .* (phi_data.polybase.Vandermonde \ (ni - ne)) # least squares if Vandermonde is not square
end

################################
# Building solutions
################################

"""
    Evaluates fe on a given cartesian mesh.
    Outputs:
        * eval_fe : Matrix (Nx,Nv)
"""
function eval_cartmesh_fe(f_eb, phi_data, meshx, meshv)
    eval_fe = zeros(length(meshv), length(meshx))
    for (iv, v) in enumerate(meshv)
        for (ix, x) in enumerate(meshx)
            eval_fe[iv,ix] = f_eb(get_v_char_e(v, eval_phi(phi_data,x), 0.0))
        end
    end
    return eval_fe
end

"""
    Evaluates fi on a given cartesian mesh. 
    Outputs:
        * eval_fi : Matrix (Nx,Nv)
"""
function eval_cartmesh_fi(f_eb, phi_data, meshx, meshv)
    eval_fi = zeros(length(meshv), length(meshx))
    phix  = eval_phi(phi_data, meshx)

    # mesh of x=0, v<0
    meshv0 = collect(LinRange(min(get_v_char_e(0.0,phix[end],0.0),meshv[1],-vi_bound),0.0,length(meshx)))
    # for each x, we build the vi's 
    for (ix, x) in enumerate(meshx)
        lowervis = get_v_char_i(meshv0, 0.0, phix[ix])      # mesh of (x=x_target, v<0)
        uppervis = get_v_char_i(0.0, phix[2:ix], phix[ix])  # mesh of (0<x<=x_target, v=0)
        vis = vcat(lowervis, uppervis)

        for (iv, v) in enumerate(meshv)
            if (v<=0)
                char, vi = integrate_char(f_eb, meshx, phix, get_v_char_i(v,phix[ix],phix[end]), ix)
                eval_fi[iv,ix] = char
            else
                ind = findlast(vis .<= -v) # v \in [vis[ind],vis[ind+1][
                stop_index = max(1, ind - length(lowervis) + 1)
                char1 = integrate_char_outofmeshx(f_eb, meshx, phi_data, phix, get_v_char_i(-v,phix[ix],phix[end]), stop_index)
                char2, = integrate_char(f_eb, meshx, phix, get_v_char_i(-v,phix[ix],phix[end]), ix)
                eval_fi[iv,ix] = 2.0*char1 - char2
            end
        end
    end
    return eval_fi
end

function eval_cartmesh_fi_mehdi(f_eb, phi_data, meshx, meshv)
    eval_fi = zeros(length(meshv), length(meshx))
    phix  = eval_phi(phi_data, meshx)
    phidx = eval_phi_dx(phi_data, meshx)

    # mesh of x=0, v<0
    meshv0 = collect(LinRange(min(get_v_char_e(0.0,phix[end],0.0),meshv[1],-vi_bound),0.0,length(meshx)))
    # for each x, we build the vi's 
    for (ix, x) in enumerate(meshx)
        lowervis = get_v_char_i(meshv0, 0.0, phix[ix])      # mesh of (x=x_target, v<0)
        uppervis = get_v_char_i(0.0, phix[2:ix], phix[ix])  # mesh of (0<x<=x_target, v=0)
        vis = vcat(lowervis, uppervis)

        for (iv, v) in enumerate(meshv)
            if (v<=0)
                char, = integrate_char_mehdi(f_eb, meshx, phix, phidx, ix, v)
                eval_fi[iv,ix] = char
            else
                ind = findlast(vis .<= -v) # v \in [vis[ind],vis[ind+1][
                stop_index = max(1, ind - length(lowervis) + 1)
                # println("bruh, ix $ix, iv $iv, ind $ind, stop_index $stop_index")
                char1 = integrate_char_outofmeshx_mehdi(f_eb, meshx, phi_data, phix, phidx, ix, -v, stop_index)
                char2, = integrate_char_mehdi(f_eb, meshx, phix, phidx, ix, -v)
                eval_fi[iv,ix] = 2.0*char1 - char2
            end
        end
    end
    return eval_fi
end

################################
# Saving solutions
################################

"""
    Extension by radial symmetry to [-1,1] (instead of [0,1])
    Transpose to have the right way.
"""
function format_for_yann(meshx, meshv, ff)
    # extension
    ext_eval = zeros(length(meshv), 2*length(meshx)-1)
    ext_eval[:,length(meshx):end] .= ff;
    ext_eval[:,begin:length(meshx)-1] .= ff[end:-1:begin,end:-1:begin+1];

    return ext_eval' # from variables (v,x) to (x,v)
end

function save_cart(f_eb, phi_data, meshx, Nve, Nvi, folder, filetag)
    meshve = collect(LinRange(-ve_bound,ve_bound,Nve))
    meshvi = collect(LinRange(-vi_bound,vi_bound,Nvi))

    eval_fe = eval_cartmesh_fe(f_eb, phi_data, meshx, meshve)
    # eval_fi = eval_cartmesh_fi(f_eb, phi_data, meshx, meshvi)
    eval_fi = eval_cartmesh_fi_mehdi(f_eb, phi_data, meshx, meshvi)

    ext_eval_fe = format_for_yann(meshx, meshve, eval_fe)
    ext_eval_fi = format_for_yann(meshx, meshvi, eval_fi)

    display(plot(heatmap(ext_eval_fe',size=(400,400)), heatmap(ext_eval_fi',size=(400,400)), layout=[1,1]))

    # for (name, data) in zip(["$folder/fe_$(filetag)_$Nve-$Nvi.dat", "$folder/fi_$(filetag)_$Nve-$Nvi.dat"], [ext_eval_fe, ext_eval_fi])
    #     println("Saving $name, size ", length(data[:,1]), "×", length(data[1,:]))
    #     writedlm(name, data) 
    # end

    return meshve, meshvi, eval_fe, eval_fi
end

function evaluate_fi(f_eb, phi_data, Ncharmesh, x, v)
    # in all cases, we need the value of fᵢ(x,-|v|)
    meshx = collect(LinRange(x,1.0,2+ceil(Int,Ncharmesh*(1.0-x))))
    phix  = eval_phi(phi_data, meshx)
    # v_1 = get_v_char_i(-abs(v),phix[begin],phix[end])
    v_1 = get_v_char_i(-abs(v),eval_phi(phi_data, x),phix[end])
    fi, = integrate_char(f_eb, meshx, phix, v_1, 1) # first index relatively to the local mesh

    if (v>0) # bad case : we also need to integrate the whole characteristic
        meshx = collect(LinRange(0.0,1.0,Ncharmesh))
        phix  = eval_phi(phi_data, meshx)
        ind = findlast(v^2/2.0 .+ eval_phi(phi_data,x) .- phix .<= 0.0)
        ind = (ind == nothing ? 1 : ind) # si toutes les coords de liap sont positives
        fitot, = integrate_char_outofmeshx(f_eb, meshx, phi_data, phix, v_1, ind)

        fi = 2.0 * fitot - fi; # uses fᵢ(x,v) + fᵢ(x,-v) = 2.0 * fᵢ(x₀,0)
    end
    return fi
end

function evaluate_fi_mehdi(f_eb, meshx, phi_data, phix, phidx, Ncharmesh, x, v)
    # in all cases, we need the value of fᵢ(x,-|v|)
    newmeshx = collect(LinRange(x,1.0,2+ceil(Int,Ncharmesh*(1.0-x))))
    newphix  = eval_phi(phi_data, meshx)
    newphidx = eval_phi_dx(phi_data, meshx);
    fi, = integrate_char_mehdi(f_eb, newmeshx, newphix, newphidx, 1, v) # first index relatively to the local mesh

    if (v>0)
        liap = v^2/2.0 .+ eval_phi(phi_data,x)
        if (liap >= 0) # integration touches (x=0,v<0)
            fitot, = integrate_char_mehdi(f_eb, meshx, phix, phidx, 1, 0.0)
        else # need to find the intersection with (x>0, v=0)
            ind = findlast(liap .- phix .<= 0.0)
            fitot, = integrate_char_outofmeshx_mehdi(f_eb, meshx, phi_data, phix, phidx, 1, 0.0, ind)
        end
        fi = 2.0 * fitot - fi; # uses fᵢ(x,v) + fᵢ(x,-v) = 2.0 * fᵢ(x₀,0)
    end
    return fi
end

"""
    Draws points in (x,v) ∈ [0,1]×[-ve_bound,0], and checks
    whether the advection equation on fₑ is satisfied.
"""
function test_df_fe(f_eb, phi_data, Nsamples, eps)
    errs  = zeros(Nsamples,Nsamples);

    max_err = 0.0;   mean_err = 0.0;
    max_dx_fe = 0.0; mean_dx_fe = 0.0;
    max_dv_fe = 0.0; mean_dv_fe = 0.0;
    Neffective = 0
    for (ix,x) in enumerate(LinRange(0.0,1.0,Nsamples))
        phix   = eval_phi(phi_data, x)
        phixpe = eval_phi(phi_data, x+eps)
        phixme = eval_phi(phi_data, x-eps)
        for (iv,v) in enumerate(LinRange(-vi_bound, vi_bound, Nsamples))
            # if we are in the right bounds
            if ((eps < x < 1.0 - eps) && (- ve_bound+eps < v < - eps))
                dx_fe = (f_eb(get_v_char_e(v    ,phixpe,0.0)) - f_eb(get_v_char_e(v    ,phixme,0.0))) / (2.0 * eps)
                dv_fe = (f_eb(get_v_char_e(v+eps,phix  ,0.0)) - f_eb(get_v_char_e(v-eps,phix  ,0.0))) / (2.0 * eps)
                err = abs(v * dx_fe + 1/mu * eval_phi_dx(phi_data, x) * dv_fe)

                mean_err += err; max_err = max(max_err, err)
                mean_dx_fe += abs(dx_fe); max_dx_fe = max(max_dx_fe, abs(dx_fe))
                mean_dv_fe += abs(dv_fe); max_dv_fe = max(max_dv_fe, abs(dv_fe))

                Neffective += 1
                errs[iv,ix] = err;
            end
        end
    end
    if (Neffective>0)
        mean_err = mean_err / Neffective
        mean_dx_fe = mean_dx_fe / Neffective
        mean_dv_fe = mean_dv_fe / Neffective
    end

    println("$Neffective / $(Nsamples^2) effective evaluations")
    println("err : mean $mean_err, max $max_err")
    println("dx  : mean $mean_dx_fe, max $max_dx_fe")
    println("dv  : mean $mean_dv_fe, max $max_dv_fe")

    file = open("test_df_fe.txt", "a")
    print(file, "$f_eb_str\t$lambda\t$mu\t$nu\t$cardinal\t$Nx\t$eps\t$Nsamples\t$mean_err\t$max_err\n")
    close(file)

    return errs
end

"""
    Draws uniform points in (x,v) ∈ [0,1]×[-vi_bound,vi_bound], and checks
    whether the advection equation on fᵢ is satisfied.
"""
function test_df_fi(f_eb, phi_data, Nsamples, eps, Ncharmesh, thresh_liap)
    errs  = zeros(Nsamples,Nsamples);

    max_err = 0.0;   mean_err = 0.0;
    max_dx_fi = 0.0; mean_dx_fi = 0.0;
    max_dv_fi = 0.0; mean_dv_fi = 0.0;
    Neffective = 0
    meshx = collect(LinRange(0.0,1.0,Nsamples))
    meshv = collect(LinRange(-vi_bound, vi_bound, Nsamples))
    phix = eval_phi(phi_data, meshx)
    phidx = eval_phi_dx(phi_data, meshx)
    for (ix,x) in enumerate(meshx)
        for (iv,v) in enumerate(meshv)
            liap = v^2 / 2.0 + phix[ix]
            # if we are in the right bounds
            if ((eps < x < 1.0 - eps) && (- vi_bound + eps < v < vi_bound - eps) && (abs(liap) >= thresh_liap))
                # dx_fi = (evaluate_fi(f_eb,phi_data,Ncharmesh,x+eps,v) - evaluate_fi(f_eb,phi_data,Ncharmesh,x-eps,v)) / (2.0 * eps)
                # dv_fi = (evaluate_fi(f_eb,phi_data,Ncharmesh,x,v+eps) - evaluate_fi(f_eb,phi_data,Ncharmesh,x,v-eps)) / (2.0 * eps)
                dx_fi = (evaluate_fi_mehdi(f_eb,meshx,phi_data,phix,phidx,Ncharmesh,x+eps,v) - evaluate_fi_mehdi(f_eb,meshx,phi_data,phix,phidx,Ncharmesh,x-eps,v)) / (2.0 * eps)
                dv_fi = (evaluate_fi_mehdi(f_eb,meshx,phi_data,phix,phidx,Ncharmesh,x,v+eps) - evaluate_fi_mehdi(f_eb,meshx,phi_data,phix,phidx,Ncharmesh,x,v-eps)) / (2.0 * eps)
                source = nu * f_eb(get_v_char_e(v,phix[ix],0.0))
                err = abs(v * dx_fi - phidx[ix] * dv_fi - source)

                # lhs = v * dx_fi - eval_phi_dx(phi_data, x) * dv_fi
                # rhs = source
                # println("lhs : $lhs, rhs : $rhs, liap : $liap")

                mean_err += err; max_err = max(max_err, err)
                mean_dx_fi += abs(dx_fi); max_dx_fi = max(max_dx_fi, abs(dx_fi))
                mean_dv_fi += abs(dv_fi); max_dv_fi = max(max_dv_fi, abs(dv_fi))

                Neffective += 1
                errs[iv,ix] = err;
            end
        end
    end
    if (Neffective>0)
        mean_err = mean_err / Neffective
        mean_dx_fi = mean_dx_fi / Neffective
        mean_dv_fi = mean_dv_fi / Neffective
    end

    println("$Neffective / $(Nsamples^2) effective evaluations")
    println("err : mean $mean_err, max $max_err")
    println("dx  : mean $mean_dx_fi, max $max_dx_fi")
    println("dv  : mean $mean_dv_fi, max $max_dv_fi")
    # println("$eps\t$Nsamples\t$Ncharmesh\t$thresh_liap\t$mean_err\t$max_err")

    file = open("test_df_fi.txt", "a")
    print(file, "$f_eb_str\t$lambda\t$mu\t$nu\t$cardinal\t$Nx\t$eps\t$Nsamples\t$Ncharmesh\t$thresh_liap\t$mean_err\t$max_err\n")
    close(file)

    return errs
end

################################
# Diagnostics
################################

"""
    Computes Jᵢ(x) := ∫ᵥ v fᵢ(x,v) dv for x ∈ meshx.
"""
function get_Ji(f_eb, meshx, phix)
    Ji = 0.0 * phix; 

    filoc_old = 0.0*Ji; viloc_old = 0.0*phix;
    filoc = 0.0*Ji; viloc = 0.0*phix;

    # get the meshes
    upper_mesh, lower_mesh = get_fi_char_meshes(phix)

    for (ivm, vm) in enumerate(vcat(lower_mesh,upper_mesh))
        stop_index = max(1,ivm-length(lower_mesh)+1)

        filoc .*= 0.0;
        viloc[end] = vm;
        # compute fi and vi along the characteristic 
        for (ix, x) in enumerate(meshx[end-1:-1:stop_index])
            jx = length(phix)-ix # real index in arrays
            viloc[jx] = get_v_char_i(vm, phix[end], phix[jx])
            dt = (x-meshx[jx+1]) * 2.0 / (viloc[jx] + viloc[jx+1])
            filoc[jx] = filoc[jx+1] + nu * dt * 0.5 .* (f_eb(get_v_char_e(viloc[jx],phix[jx],0.0)) + f_eb(get_v_char_e(viloc[jx+1],phix[jx+1],0.0)))
        end

        # update the current
        if (ivm>2)
            # negative speeds : Jᵢ += Δvₖ * (vₖ fᵢₖ + vₖ₋₁ fᵢₖ₋₁) / 2
            Ji[stop_index:end] .+= (viloc[stop_index:end] .- viloc_old[stop_index:end]) .* 0.5 .* 
                                   (viloc[stop_index:end] .* filoc[stop_index:end] .+ viloc_old[stop_index:end] .* filoc_old[stop_index:end]) 
            # positive speeds : fᵢ(x,v) + fᵢ(x,-v) = cte = 2 * fᵢ(x₀,0)
            Ji[stop_index:end] .+= (viloc_old[stop_index:end] .- viloc[stop_index:end]) .* 0.5 .* 
                                   (viloc[stop_index:end] .* (2.0*filoc[stop_index] .- filoc[stop_index:end]) .+ 
                                   viloc_old[stop_index:end] .* (2.0*filoc_old[stop_index] .- filoc_old[stop_index:end])) 
        end

        filoc_old .= filoc; viloc_old .= viloc;
    end
    return Ji
end

function check_Ji(f_eb, phi_data, Nxs)
    errors = zeros(length(Nxs), 2)
    for (iNx, Nx) in enumerate(Nxs)
        println("Computing error for Nx=$Nx")
        meshx = collect(LinRange(0.0,1.0,Nx))
        phix = eval_phi(phi_data, meshx)
        Ji = get_Ji(f_eb, meshx, phix)
        ne = zeros(length(phix))
        update_ne!(ne, f_eb, phix)

        # Again trapezes integral 
        int_nu_ne = nu .* vcat([0.0], cumsum(0.5 .* (ne[begin+1:end] + ne[begin:end-1]) .* (meshx[begin+1:end] - meshx[begin:end-1])))
        # p = plot(meshx, Ji, label="Jᵢ")
        # plot!(p, meshx, int_nu_ne, label="∫νnₑ")
        # display(p)

        errors[iNx, 1] = maximum(abs.(Ji - int_nu_ne)) / maximum(max.(abs.(Ji), abs.(int_nu_ne)))
        errors[iNx, 2] = sum(abs.(Ji - int_nu_ne)) / sum(max.(abs.(Ji), abs.(int_nu_ne)))
    end
    println("Erreur relative L∞ entre Jᵢ et ∫[0,x] ν nₑ(y)dy : ", errors[:,1])
    println("Erreur relative L1 entre Jᵢ et ∫[0,x] ν nₑ(y)dy : ", errors[:,2])

    # p = plot(meshx, Ji, label="Jᵢ", legend=:bottomright, xlabel="x")
    # plot!(p, meshx, int_nu_ne, label="∫νnₑ")
    # display(p)
    errors = hcat(Nxs, errors)

    return errors
end

function plot_variations(folder, variations, starts)
    len = length(variations)
    for start in starts
        p = plot(start:len, variations[start:end],marker=true, label="|ϕₖ - ϕₖ₋₁|", xlabel="k", ylabel="variation")
        png(p, "$folder/variations$start.png")
    end
end

################################
# Main
################################

function fixed_point(Nx::Int)
    println("\nBegin of Nx=$Nx")
    meshx = LinRange(0.0,1.0,Nx) |> collect; # space mesh.
    polybase = Polybase(cardinal, basis, meshx);
    phi_data = Phi_data(polybase, coeff);
    variations = []
    
    ni = zeros(size(meshx)); #      ion density (depends on space only)
    ne = zeros(size(meshx)); # electron density (depends on space only)

    it=1;
    nolds = 4;
    old_coeffs = [1.0*phi_data.coeff for _ in 1:nolds] # old coefficients
    phix = eval_phi(phi_data, meshx);
    phidx = eval_phi_dx(phi_data, meshx);
    variation = 1.0; # Linf-norm of (phix - phix_old) in the space of coefficients

    # var_file = open("variations.txt", "w")

    while ((it <= Nit_max) && (variation > 1e-9))
        println("Iteration $it / $Nit_max")

        update_ne!(ne, f_eb, phix);      
        # p=plot(meshx, ne, label="nₑ")
        # update_ni!(ni, f_eb, phix, meshx); 
        # p=plot(meshx, ni, label="mine")
        # println("coord en 0 : ", ni[1])
        update_ni_mehdi!(ni, f_eb, phi_data, phix, phidx, meshx); 
        # plot!(p,meshx, ni, label="nᵢ")
        # display(p)
        # # println("coord en 0 : ", ni[1])
        # readline()
        # update_phi!(phi_data, ni, ne); 
        newcoeff = update_phi(phi_data, ni, ne);
        # newcoeff = update_phi(phi_data, min.(1e3,ni), ne);

        newphix = eval_phi(phi_data.polybase, newcoeff, meshx)
        if (all(newphix .<= 0.0))
            phi_data.coeff[:] = newcoeff[:];
        else
            println("\tNew coeffs make phi positive -> incomplete iteration")
            prop = 0.5*phix[end] / newphix[end];
            phi_data.coeff[:] .+= prop .* newcoeff[:];
        end
        phix = eval_phi(phi_data, meshx);
        phidx = eval_phi_dx(phi_data, meshx);

        variation = max(abs.(phi_data.coeff - old_coeffs[begin])...)
        push!(variations, variation)
        # print(var_file, "$Nx\t$it\t$variation\n")
        for iold in nolds:-1:2 # not good, i know
            println("\tNorm (∞,coeff) of phi_k - phi_{k-$(iold)} : ", max(abs.(phi_data.coeff - old_coeffs[iold])...))
            old_coeffs[iold] .= old_coeffs[iold-1]
        end
        println("Norm (∞,coeff) of phi_k - phi_{k-1} : ", variation, ", coeffs : ", phi_data.coeff)
        old_coeffs[begin] .= phi_data.coeff;
        it += 1;
    end

    # close(var_file)

    println("End of Nx=$Nx\n")
    return meshx, phi_data, variations
end 

function my_heatmap(meshx, meshv, ext_eval, thresh)
    ext_meshx = vcat(-meshx[end:-1:begin+1],meshx)
    return heatmap(ext_meshx, meshv, sqrt.(min.(thresh,ext_eval)))
end

#=
Reminder of what to do : 
    - clean the folder
    - run the code haha
    - save variations
    - export errors for Ji
    - export errors for df_fi and df_fe
=#

println("Welcome in fixedpoint.")

Nracine = 128
# thefolder = "../latex/test_cases/case_2/julia$(2*Nracine)"

Nx  = Nracine    # !!! This is the mesh on [0,1]. Will be extended to a mesh of Nx' := 2 Nx - 1 points.
Nve = Nracine*2   # Number of points of the velocity mesh for electrons
Nvi = Nracine*2   # Number of points of the velocity mesh for ions

# rm("test_df_fe.txt", force=true)
# rm("test_df_fi.txt", force=true)
println("Run for $(2*Nx - 1) points in x, $Nve in ve, $Nvi in vi (so Nx=$(2*Nx - 2), Nve=$(Nve-1), Nvi=$(Nvi-1)).")
meshx, phi_data, variations = fixed_point(Nx)
print("coeff = ["); [print("$a, ") for a in phi_data.coeff[begin:end-1]]; a=phi_data.coeff[end]; println("$a]")
# save_cart(f_eb, phi_data, meshx, Nve, Nvi, "../hivlashea/build/input_data", "case2")
# writedlm("$thefolder/variations.txt", variations)
# plot_variations(thefolder, variations, 1)
# writedlm("$thefolder/err_Ji.txt", check_Ji(f_eb, phi_data, 100:50:300))
# for eps in [0.001,0.0001,0.00001]
#     errfe = test_df_fe(f_eb, phi_data, 100, eps)
#     errfi = test_df_fi(f_eb, phi_data, 100, eps, 1000, 0.0)
#     writedlm("$thefolder/errfe_eps$eps.dat", errfe)
#     writedlm("$thefolder/errfi_eps$eps.dat", errfi)
# end
# mv("test_df_fe.txt", "$thefolder/test_df_fe.txt")
# mv("test_df_fi.txt", "$thefolder/test_df_fi.txt")

# plot(heatmap(ext_meshx, meshvi, ext_eval_fi,size=(400,500)), heatmap(ext_meshx, meshve, ext_eval_fe,size=(400,500)), layout=[1,1])

println("Bye")

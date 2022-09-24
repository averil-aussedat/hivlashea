###########################################################################
# Test code for fixed-point algorithm
# Vlasov-Poisson with nonperidic boundary conditions
# CEMRACS 2022 
###########################################################################

using LinearAlgebra 
using Plots

################################
# Inputs (in conf file)
################################

########## Parameters ##########
# physics
const lambda    = 0.1; # Debye length
# const mu        = 1.0/3672; # ratio mass_electron / mass_ion
const mu        = 1.0; # ratio mass_electron / mass_ion
const nu        = 25.0; # collision frequency

# 1sp parameters 
const G         = 1.0 # gravitational force
const rho_0     = 1.0 
const z0        = 1.0 # keep != 0 for non-null solution

# const xmax = 1.0;
const xmax = 0.2;
fact = 10.0;
# electronic distribution at x=0 (expression and latex-printable description)
f_eb(v) = sqrt(mu*fact) .* exp.(-mu * fact * (v.^2)/2) ./ sqrt(2pi); f_eb_str = "\$\frac{- $fact \\mu \\frac{v^2}{2}}{$fact \\times 2 \\pi \\mu}\$";
# f_eb(v) = 1.0.*(abs.(v).<=1.0); f_eb_str = "\$1_{\left\{|v|\leqslant 1\right\}}\$";
# f_eb(v) = max.(1.0.+v,0.0); f_eb_str = "\$(1+v)_{+}\$";
# f_eb(v) = 0.0005 .* max.(0.0, v.^2 .* (100.0 .- v.^2)); f_eb_str = "\$5e-4\\times (v^2 (100-v^2))_{+}\$"; 
# f_eb(v) = exp.(-mu .* (abs.(v) .- 2.0).^2 ./(2*fact)) ./ sqrt(mu * 2pi*fact); f_eb_str = "\$\\frac{- \\mu \\frac{(|v|-2)^2}{2\\times $fact}}{$fact \\times 2 \\pi \\mu}\$";
# f_eb(v) = 0.01 .* v.^2 .* exp.(-mu * fact .* (abs.(v) .- 10.0).^2 ./2) ./ sqrt(mu * 2pi*fact); f_eb_str = "\$v^2\\frac{- $fact \\mu \\frac{(|v|-10)^2}{2}}{$fact \\times 2 \\pi \\mu}\$";

# sol_1sp(x,v) = min.(500.0,(v.^2 .< z0^2 .- x.^2) .* rho_0 ./ (pi .* sqrt.(max.(1e-14,z0^2 .- x.^2 .- v.^2))));
sol_1sp(x,v) = (v.^2 .< z0^2 .- x.^2) .* rho_0 ./ (pi .* sqrt.(max.(1e-14,z0^2 .- x.^2 .- v.^2)));
f_b(v) = sol_1sp(0.0,v);
phi_sol_coeff = 2*pi*G*rho_0;

# numerics
const cardinal  = 3 # number of basis functions for phi
# const basis     = "canonical"; # choice of representation
const basis     = "pair"; # choice of representation
Nx              = 100; # number of mesh points in [0,1]
Nit_max         = 200; # maximal number of fixed-point iterations

# expectations on the solution
const ve_bound = 50.0; # electron speed ∈ [-ve_bound, ve_bound]
const vi_bound = 10.0; #      ion speed ∈ [-vi_bound, vi_bound]

coeff = [(-30.0)*(i==1) for i=1:cardinal]; # initial coefficients in the basis
# coeff = [phi_sol_coeff*(i==1) for i=1:cardinal]; # initial coefficients in the basis
# coeff = [(i==1) for i=1:cardinal]; # initial coefficients in the basis
# coeff = [-0.702674, -0.0986662, 2.29555, -12.3569, 37.3819, -70.7272, 85.2078, -63.5163, 26.7158, -4.84969] # l=1, mu=1, nu=4
# coeff = [-0.531322, -0.093826, 1.96247, -10.5839, 32.2984, -61.8136, 75.366, -56.8386, 24.1725, -4.43374] # l=1, mu=1, nu=3
# coeff = [-0.348437, -0.084665, 1.57194, -8.49328, 26.2306, -50.9633, 63.0867, -48.2639, 20.7982, -3.86115] # l=1, mu=1, nu=2
# coeff = [-0.155825, -0.0641009, 1.0371, -5.61937, 17.6735, -35.0784, 44.3161, -34.5319, 15.1248, -2.84867] # l=1, mu=1, nu=1

# l=0.5, mu=0.5, nu=1.0
# coeff = [-0.0822392744643385, -0.02597504749200164, 0.4008622169690375, -2.1153265411723816, 6.601625660930767, -13.084139950035043, 16.544191860500963, -12.913062134436489, 5.6664499683479495, -1.0692160805200983]
# l=0.1, mu=0.1, nu=1.0
# coeff = [-0.014381754590819622, -0.0027868196094572974, 0.03855291983935587, -0.1891983516966756, 0.5753986954106687, -1.130137452211309, 1.4248756776249012, -1.1114751636064732, 0.4878498099847744, -0.09210151771684331]
# l=0.01, mu=0.01, nu=1.0
# coeff = [-0.0008930558884573836, -0.00010193942488896136, 0.0011906227713942333, -0.005266630614402255, 0.015453758851181726, -0.029835124109966353, 0.03719374716574583, -0.02876007885100652, 0.012531960762369478, -0.0023511728582132053]
# l=0.05, mu=0.01, nu=1.0
# coeff = [-0.00925028981173814, -0.004008546647795378, 0.05275881822154073, -0.2930938048642883, 0.9723294439140623, -2.029872306882787, 2.673623348493392, -2.154843667570554, 0.9701290285376363, -0.1869366011207519]
# l=0.06, mu=0.01, nu = 2.0
# coeff = [-0.028473754748786738, -0.009759246689599594, 0.13464188443881345, -0.7508651377889196, 2.475818017354825, -5.134697208592925, 6.7272114622655925, -5.399603106596326, 2.4231607594353206, -0.46573202669516944]
# l=0.1, mu=0.01, nu = 25
# coeff = [-1.0393848102983225, -0.06554124852998945, 1.502418082294132, -8.425095981888637, 26.33389164725233, -51.26782622343658, 63.33414210911024, -48.26671303441629, 20.704159475122253, -3.825376369312274]
# l=0.1, mu=0.01, nu=20
# coeff = [-0.8388996865373262, -0.06777453615344116, 1.3934545222020371, -7.87279589659249, 24.897164464566128, -49.10049547884866, 61.40228162932173, -47.31795438084298, 20.501599750667378, -3.822293224067575]
# l=0.1, mu=0.01, nu=15
# coeff = [-0.6262692012638873, -0.06864288066340456, 1.2670807209704122, -7.217991832736354, 23.152740548821292, -46.345813801301624, 58.75056393313815, -45.82205948817675, 20.06427320628274, -3.7757606388802154]
# l=0.1, mu=0.01, nu=10
# coeff = [-0.40140227853143445, -0.06600391146589527, 1.0955794032826784, -6.293872351683455, 20.52605741210693, -41.77810635269354, 53.742599041416284, -42.44745140780811, 18.788290951505658, -3.5687470368201133] 
# l=0.1, mu=0.01, nu=5
# coeff = [-0.17115092517535477, -0.05341070662788475, 0.7964956016151127, -4.610508387729514, 15.32610410986369, -31.778781792209955, 41.53100195309726, -33.23763542531345, 14.874625301271875, -2.8517936354789537]
# l=0.1, mu=0.01, nu=4
# coeff = [-0.12698669450854022, -0.04828993595939699, 0.7036475048168308, -4.077924350835785, 13.614595785627833, -28.348054974849564, 37.17691265163232, -29.837996649785303, 13.38449319037601, -2.5711167278408027]
# l=0.1, mu=0.01, nu=3
# coeff = [-0.08501900556979783, -0.04161405076633612, 0.5914276371381305, -3.4316543410211207, 11.512659565946235, -24.08225245912275, 31.70290605999349, -25.522915691244883, 11.477500474892226, -2.2093539093756354]
# l=0.1, mu=0.01, nu=3, fact=10
# coeff = [-0.04444563540657112, -0.011074976413198387, 0.15735878877669554, -0.8689113452677147, 2.8367980232355703, -5.83567103586794, 7.598019376479059, -6.069304120270537, 2.7134232575028125, -0.5199309589314749]
# l=0.1, mu=0.1, nu=20.0
# coeff = [-0.3488412026940758, -0.01689698787582626, 0.38559745987391586, -1.875181277288754, 5.334673470772697, -9.593524267103227, 11.061731030230055, -7.9308966736852575, 3.220240405855763, -0.5658390602527037]

# lambda = 0.1, mu = 0.1, nu = 20, double max fact = 0.1
# coeff = [-0.5953343710367655, -0.010128445371518463, 0.06733707897607889, -0.5104885207588676, 1.3198118373651397, -1.9374951549928499, 1.7699159196292717, -0.9401760841539892, 0.2423268721572459, -0.016950850384284812]
# lambda = 0.1, mu = 0.01, nu = 20, double max fact = 0.1
# coeff = [-1.1596945276801447, -0.05642332746058975, 2.0678278022919843, -11.82283497475025, 31.08691162271426, -48.57023319308005, 49.00395771466013, -31.677642152701136, 11.962922547353696, -2.003003721196997]
# lambda = 0.01, mu = 0.01, nu = 20, double max fact = 0.1
# coeff = [-1.3370441664862496, 0.04105143880635959, -0.5330543809644839, 4.42195611307296, -17.449622938983467, 41.7487142275851, -58.71673716320553, 46.72632159102738, -19.56350822226556, 3.3703307271558645]
# lambda = 0.01, mu = 0.008, nu = 20, double max fact = 0.1
# coeff = [-1.45257946563025, -0.006402158690387889, 0.14581504439581083, 0.9977735121764416, -6.956119550021013, 26.69837584323909, -51.87555542704614, 51.138793366683394, -24.92062908906574, 4.812082508380441]
# lambda = 0.01, mu = 0.005, nu = 25, double max fact = 0.1
# coeff = [-2.011700471778535, -0.4026982491676736, 6.385923404446249, -39.269743333674064, 158.163281828977, -361.7784146171964, 463.13144213414165, -332.69115471540067, 125.7414591443631, -19.45926688012755]

# 0.01, mu = 0.005, nu = 25, test 2.0
# coeff = [-3.891692146609259, 3.6001427631509273, -48.012861693570905, 308.05811203563667, -1023.3700805078103, 1919.9473961415167, -2142.2266110075166, 1416.9201427244395, -514.4806406212646, 79.13866669650508]

# après correction du bug de parenthèses
# lambda = 0.1, mu = 0.01, nu = 25, cardinal 3, base paire
coeff = [-29.275193716990835, 9.970222191354972, -2.9013822611511237]

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
        println("iteration $k, Jx = $Jx, x = $x")
        k += 1;
    end
    return x
end

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
    Computes n = ∫f(x,v)dv for x in meshx.
    Uses a v mesh with roughly the same step as meshx.
"""
function update_n!(n, f_b, phix)
    meshv = collect(LinRange(-2.0, 2.0, 500000))
    for (i,phixi) in enumerate(phix)
        f_values = f_b(get_v_char_i(meshv,phixi,0.0))
        n[i] = trapezes(f_values, meshv[2]-meshv[1])
        # n[i] = 1.0
    end

    println("erreur nᵢ à 1 : Linf ", maximum(abs.(1.0 .- n)), ", L1 ", sum(abs.(1.0 .- n))/length(n))
end

"""
    Integrates nu * f_e over the characteristic issued from (meshx[end], v_1) stopping at x=meshx[stop_index].
    The characteristic is meshed by (x_i,v_i), where x_i ∈ meshx and v_i is computed.
    Returns the integral (scalar) and the vector of v_i.
"""
function integrate_char(f_eb, meshx, phix, v_1, stop_index)
    # println("coucou, stop_index : ", stop_index)
    integral = 0.0; 
    v_is = get_v_char_i(v_1, phix[end], phix[stop_index:end])
    # println("length v_is : ", length(v_is))
    # println("length phix : ", length(phix))
    # print("mesh of char : ", meshx[stop_index:end], v_is)
    fe_previous = f_eb(get_v_char_e(v_1,phix[end],0.0))

    for ix in length(meshx)-1:-1:stop_index

        # x1 = meshx[ix]; x0 = meshx[ix+1]; v1 = v_is[ix-stop_index+1]; v0 = v_is[ix-stop_index+2]
        # dttest = acosh((x1*x0 - v1*v0) / (x0^2 - v0^2))

        dt = (meshx[ix]-meshx[ix+1]) * 2.0 / (v_is[ix-stop_index+1] + v_is[ix-stop_index+2]) # screw numbering from 1 -.-
        fe_current = f_eb(get_v_char_e(v_is[ix-stop_index+1],phix[ix],0.0))
        # fe_current = 1.0 * (-0.8 <= get_v_char_e(v_is[ix-stop_index+1],phix[ix],0.0) <= -0.4)

        integral += dt * 0.5 * (fe_current + fe_previous)
        # integral += dttest * 0.5 * (fe_current + fe_previous)

        # println("ix = $ix, x = ", meshx[ix], ", v = ", v_is[ix-stop_index+1], ", fe = $fe_current, dt = $dt, dttest = $dttest, err ", abs(dt - dttest))
        # println("ix = $ix, x = ", meshx[ix], ", v = ", v_is[ix-stop_index+1], ", x^2+v^2 = ", x1^2 + v1^2, ", fe = ", fe_current)
        # if (v_is[ix-stop_index+1] > -1.5)
        #     println("\tvar = ", sqrt(meshx[ix]^2 + v_is[ix-stop_index+1]^2), ", fe = ", fe_current, ", diff ", meshx[ix]^2 + v_is[ix-stop_index+1]^2 - get_v_char_e(v_is[ix-stop_index+1],phix[ix],0.0)^2)
        # end

        fe_previous = 1*fe_current
    end # for ix
    return nu * integral, v_is
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
            # print("xl ", meshx[ind], ", xr ", meshx[ind+1], ", x_init $x_init, x_targ $x_targ")
            dt = (x_targ - meshx[ind+1]) * 2.0 / (0.0 + get_v_char_i(v_1,phix[end],phix[ind+1]))
            fe_targ = f_eb(get_v_char_e(0.0,eval_phi(phi_data,x_targ),0.0))
            fe_prec = f_eb(get_v_char_e(meshx[ind+1],phix[ind+1],0.0))
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

function update_ni!(ni, f_eb, phix, meshx, beta1, beta2)

    # meshh = collect(LinRange(-1.0,1.0,100))
    # display(plot(meshh, f_eb(meshh)))

    ni .*= 0.0

    abserr = 0.0; sumerr = 0.0;

    # get the meshes
    upper_mesh, lower_mesh = get_fi_char_meshes(phix)

    # lower mesh part : all coordinates of nᵢ are updated
    integral_char_previous, v_is_previous = integrate_char(f_eb, meshx, phix, lower_mesh[1], 1) # first characteristic's integral
    for vlm in lower_mesh[2:end]
        integral_char, v_is = integrate_char(f_eb, meshx, phix, vlm, 1) # stop_index still at 1
        ni .+= (v_is .- v_is_previous) .* (integral_char + integral_char_previous) # times 2 is embedded
        # ni .+= 0.01 .* (integral_char + integral_char_previous) # times 2 is embedded

        # println("vis - vis previous : ", v_is .- v_is_previous)

        v00 = get_v_char_i(vlm, phix[end], 0.0)
        sol = nu * (asinh(sqrt(max((beta1^2 - v00^2)/ (2.0 * v00^2), 0.0))) - asinh(sqrt(max((beta2^2 - v00^2)/ (2.0 * v00^2), 0.0))))
        if (vlm < lower_mesh[end])
            abserr = max(abserr, abs(integral_char - sol))
            sumerr += abs(integral_char - sol)
        end
        # println("integral char issued from $vlm : $integral_char, sol : $sol, error : ", abs(integral_char - sol))

        integral_char_previous = integral_char; v_is_previous .= v_is # update
    end # for v lower mesh

    # upper mesh part : some characteristics stops 
    for (stop_index, vum) in enumerate(upper_mesh)
        integral_char, v_is = integrate_char(f_eb, meshx, phix, vum, stop_index+1) # v_is is of size length(meshx)-stop_index now
        ni[stop_index+1:end] .+= (v_is .- v_is_previous[2:end]) .* (integral_char + integral_char_previous) # times 2 is embedded
        # ni[stop_index+1:end] .+= 0.01 .* (integral_char + integral_char_previous) # times 2 is embedded

        x00 = sqrt(1.0 - vum^2)
        # println("x00 : $x00, supposedly ", meshx[stop_index+1])
        sol = nu * (asinh(sqrt(max((beta1^2 - x00^2)/ (2.0 * x00^2), 0.0))) - asinh(sqrt(max((beta2^2 - x00^2)/ (2.0 * x00^2), 0.0))))
        sumerr += abs(integral_char - sol)
        abserr = max(abserr, abs(integral_char - sol))
        # println("integral char issued from $vum : $integral_char, sol : $sol, error : ", abs(integral_char - sol))

        integral_char_previous = integral_char; v_is_previous = 1*v_is # update
    end # for v upper mesh

    println("maximum error : $abserr, sum : ", sumerr / (length(upper_mesh) + length(lower_mesh) - 2.0))
end

function update_ni_graphical!(ni, f_eb, phix, meshx)
    ni .*= 0.0
    v_iss = []; integrals = []

    # get the meshes
    upper_mesh, lower_mesh = get_fi_char_meshes(phix)

    # lower mesh part : all coordinates of nᵢ are updated
    integral_char_previous, v_is_previous = integrate_char(f_eb, meshx, phix, lower_mesh[1], 1) # first characteristic's integral
    push!(integrals, integral_char_previous); push!(v_iss, v_is_previous)
    for vlm in lower_mesh[2:end]
        integral_char, v_is = integrate_char(f_eb, meshx, phix, vlm, 1) # stop_index still at 1
        # ni .+= (v_is .- v_is_previous) .* 0.5 .* (integral_char + integral_char_previous)
        ni .+= (v_is .- v_is_previous) .* (integral_char + integral_char_previous) # times 2 embedded

        integral_char_previous = 1*integral_char; v_is_previous = 1*v_is # update
        push!(integrals, integral_char); push!(v_iss, v_is)
        # test = v_is.^2 / 2 + phix
        # println("test min : ", min(test...), ", max : ", max(test...))
    end # for v lower mesh

    # upper mesh part : some characteristics stops 
    for (stop_index, vum) in enumerate(upper_mesh)
        integral_char, v_is = integrate_char(f_eb, meshx, phix, vum, stop_index+1) # v_is is of size length(meshx)-stop_index now
        # ni[stop_index+1:end] .+= (v_is .- v_is_previous[2:end]) .* 0.5 .* (integral_char + integral_char_previous)
        ni[stop_index+1:end] .+= (v_is .- v_is_previous[2:end]) .* (integral_char + integral_char_previous) # times 2 embedded

        integral_char_previous = 1*integral_char; v_is_previous = 1*v_is # update
        push!(integrals, integral_char); push!(v_iss, v_is)
    end # for v upper mesh

    # ni .*= 2.0;
    return integrals, v_iss, upper_mesh, lower_mesh
end

# fuck
function update_ni_bruteforce!(ni, f_eb, phix, meshx, beta1, beta2)
    ni .*= 0.0

    # get the meshes
    upper_mesh, lower_mesh = get_fi_char_meshes(phix)

    xx = []; vv = []; fi = []

    for (ix, x) in enumerate(meshx)
        # vmax = get_v_char_i(0.0,phix[ix],phix[end])
        # integrate fᵢ(x,v) for v <= 0
        for (iv,v) in enumerate(lower_mesh[begin:end-1]) # rectangle integrals
            # println("low, v = $v")
            char, = integrate_char(f_eb, meshx, phix, v, ix)
            dv = get_v_char_i(lower_mesh[iv+1],phix[end],phix[ix]) - get_v_char_i(v,phix[end],phix[ix])
            ni[ix] += dv * char

            push!(xx, x); push!(vv, get_v_char_i(v,phix[end],phix[ix])); push!(fi, char)
        end
        if (ix>1)
            char, = integrate_char(f_eb, meshx, phix, lower_mesh[end], ix)
            dv = get_v_char_i(upper_mesh[1],phix[end],phix[ix]) - get_v_char_i(lower_mesh[end],phix[end],phix[ix])
            ni[ix] += dv * char

            push!(xx, x); push!(vv, get_v_char_i(lower_mesh[end],phix[end],phix[ix])); push!(fi, char)
            if (ix>2)
                for (iv,v) in enumerate(upper_mesh[begin:ix-2]) # rectangle integrals
                    # println("up, v = $v")
                    char, = integrate_char(f_eb, meshx, phix, v, ix)
                    dv = get_v_char_i(upper_mesh[iv+1],phix[end],phix[ix]) - get_v_char_i(v,phix[end],phix[ix])
                    ni[ix] += dv * char
                    push!(xx, x); push!(vv, get_v_char_i(v,phix[end],phix[ix])); push!(fi, char)
                end
            end
        end
        # integrate fᵢ(x,v) for v > 0
        if (ix>1)
            if (ix > 2)
                for (iv,v) in enumerate(upper_mesh[ix-1:-1:begin+1])
                    stop_index = ix - iv + 1 
                    # println("attempt to go from (x=$x,v=$v) to stop_index=", stop_index, ", i.e. x=", meshx[stop_index])
                    char1, = integrate_char(f_eb, meshx, phix, v, stop_index)
                    char2, = integrate_char(f_eb, meshx, phix, v, ix)
                    dv = -get_v_char_i(upper_mesh[ix-1-iv],phix[end],phix[ix]) + get_v_char_i(v,phix[end],phix[ix])
                    ni[ix] += dv * (2.0 * char1 - char2)
                    push!(xx, x); push!(vv, -get_v_char_i(v,phix[end],phix[ix])); push!(fi, 2.0 * char1 - char2)
                end
            end
            # println("transition")
            char1, = integrate_char(f_eb, meshx, phix, upper_mesh[begin], 2)
            char2, = integrate_char(f_eb, meshx, phix, upper_mesh[begin], ix)
            dv = -get_v_char_i(lower_mesh[end],phix[end],phix[ix]) + get_v_char_i(upper_mesh[begin],phix[end],phix[ix])
            ni[ix] += dv * (2.0 * char1 - char2)
            push!(xx, x); push!(vv, -get_v_char_i(upper_mesh[begin],phix[end],phix[ix])); push!(fi, 2.0 * char1 - char2)
        end
        for (iv,v) in enumerate(lower_mesh[end:-1:begin+1]) # rectangle integrals
            # println("second lower wesh, iv=$iv")
            char1, = integrate_char(f_eb, meshx, phix, v, 1)
            char2, = integrate_char(f_eb, meshx, phix, v, ix)
            dv = -get_v_char_i(lower_mesh[end-iv],phix[end],phix[ix]) + get_v_char_i(v,phix[end],phix[ix])
            # if (iv==1)
            #     ni[ix] += dv * (2.0 * char1 - char2) # the characteristic never reaches 0
            #     push!(xx, x); push!(vv, -get_v_char_i(v,phix[end],phix[ix])); push!(fi, char1 - char2)
            # else
                ni[ix] += dv * (2.0 * char1 - char2)
                push!(xx, x); push!(vv, -get_v_char_i(v,phix[end],phix[ix])); push!(fi, 2.0 * char1 - char2)
            # end
        end
    end
    return xx, vv, fi
end


################################
# Poisson solver
################################

"""
    Solves the problem 
        - λ² Δφ(x) = nᵢ(x) - nₑ(x)
              φ(0) = 0
            ∂ₓφ(0) = 0
    with phi decomposed in a polynomial basis.
"""
function update_phi!(phi_data, ni, ne)
    phi_data.coeff = - 1.0/(lambda^2) .* (phi_data.polybase.Vandermonde \ (ni - ne)) # least squares if Vandermonde is not square
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

"""
    Solves the problem 
            Δφ(x) = 4πG n(x)
              φ(0) = 0
            ∂ₓφ(0) = 0
    with phi decomposed in a polynomial basis.
"""
function update_phi_1sp!(phi_data, n)
    phi_data.coeff = 4*pi*G .* (phi_data.polybase.Vandermonde \ n) # least squares if Vandermonde is not square
end

################################
# Building solutions
################################

"""
    Evaluates f_e on a mesh given by the characteristics (only for negative v).
    Outputs:
        * xx_fe the vector of x-coordinates
        * vv_fe the vector of v-coordinates
        * fe the vector of fe(x,v) for (x,v) in zip(xx,vv)
"""
function build_fe(f_eb, phix, meshx)
    xx_fe = []; vv_fe = []; fe = [];
    for (ix,x) in enumerate(meshx)
        char_x = meshx[1:ix];
        char_v = get_v_char_e(0.0, phix[ix], phix[1:ix]);
        char_f = f_eb(char_v[1] * ones(length(char_v)))

        xx_fe = vcat(xx_fe,char_x);
        vv_fe = vcat(vv_fe,char_v);
        fe    = vcat(   fe,char_f);
    end

    return xx_fe, vv_fe, fe
end

"""
    Evaluates f_i on a mesh given by the characteristics.
    Outputs:
        * xx_fi the vector of x-coordinates
        * vv_fi the vector of v-coordinates
        * fi the vector of fi(x,v) for (x,v) in zip(xx,vv)
"""
function build_fi(f_eb, phix, meshx)
    xx_fi = []; vv_fi = []; fi = [];

    # get the meshes
    upper_mesh, lower_mesh = get_fi_char_meshes(phix)

    for (ivm,vm) in enumerate(vcat(lower_mesh,upper_mesh)) # for each characteristic start points (x=1,v=vm)
        stop_index = max(1,ivm-length(lower_mesh)+1)
        # println("Computing char for vm = $vm, ivm = $ivm, stop_index = $stop_index, length lm = ", length(lower_mesh))
        for (ix,x) in enumerate(meshx[stop_index:end])
            v = get_v_char_i(vm, phix[end], phix[ix+stop_index-1]); # current speed
            floc, vi_s = integrate_char(f_eb, meshx, phix, vm, ix+stop_index-1); # fᵢ(x,v)
            ftot, vi_s = integrate_char(f_eb, meshx, phix, vm, stop_index); # 

            # push!(xx_fi, x);
            # push!(vv_fi, v);
            # push!(fi, floc);
            xx_fi = vcat(xx_fi, [x, x]);
            vv_fi = vcat(vv_fi, [v,-v]);
            fi = vcat(fi, [floc, 2*ftot-floc]);
        end 
    end

    return xx_fi, vv_fi, fi
end 

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
        lowervis = get_v_char_i(meshv0, 0.0, phix[ix])      # mesh of (x=0, v<0)
        uppervis = get_v_char_i(0.0, phix[2:ix], phix[ix])  # mesh of (0<x<=x_target, v=0)
        vis = vcat(lowervis, uppervis)

        for (iv, v) in enumerate(meshv)
            if (v<=0)
                # println("x : $x, neg v : $v")
                char, vi = integrate_char(f_eb, meshx, phix, get_v_char_i(v,0.0,phix[end]), ix)
                eval_fi[iv,ix] = char
            else
                # println("x : $x, pos v : $v")
                ind = findlast(vis .<= -v) # v \in [vis[ind],vis[ind+1][
                # println("ind : ", ind, ", vis : ", vis, ", -v : ", -v)
                stop_index = max(1, ind - length(lowervis) + 1)
                # char1, vi = integrate_char(f_eb, meshx, phix, get_v_char_i(-v,0.0,phix[end]), stop_index)
                char1 = integrate_char_outofmeshx(f_eb, meshx, phi_data, phix, get_v_char_i(-v,0.0,phix[end]), stop_index)
                char2, vi = integrate_char(f_eb, meshx, phix, get_v_char_i(-v,0.0,phix[end]), ix)
                # println("ix $ix, x $x, stop_index $stop_index, char2 $char2, char1 $char1")
                eval_fi[iv,ix] = 2.0*char1 - char2
                # eval_fi[iv,ix] = char1 - char2
            end
        end
    end
    return eval_fi
end

"""
    Extension by radial symmetry to [-1,1] (instead of [0,1])
"""
function extend(meshx, meshv, ff)
    ext_meshx = vcat(-meshx[end:-1:begin+1], meshx)
    ext_eval = zeros(length(meshv), 2*length(meshx)-1)
    ext_eval[:,length(meshx):end] = ff;
    ext_eval[:,begin:length(meshx)-1] = ff[end:-1:begin,end:-1:begin+1];
    return ext_meshx, ext_eval
end

################################
# Plotting
################################

function plot_fe(f_eb, phi_data, meshx)
    meshve = collect(LinRange(-ve_bound,ve_bound,100));
    # eval_fe, ext_meshxe, ext_eval_fe = eval_cartmesh_fe(f_eb, phi_data, meshx, meshve);
    eval_fe = eval_cartmesh_fe(f_eb, phi_data, meshx, meshve);
    # display(surface(ext_meshxe, meshve, ext_eval_fe));
    display(surface(meshx, meshve, eval_fe));
    # return ext_meshxe, meshve, ext_eval_fe
    return meshve, eval_fe
end

function plot_fi(f_eb, phi_data, meshx)
    meshvi = collect(LinRange(-vi_bound,vi_bound,100));
    # eval_fi, ext_meshxi, ext_eval_fi = eval_cartmesh_fi(f_eb, phi_data, meshx, meshvi);
    eval_fi = eval_cartmesh_fi(f_eb, phi_data, meshx, meshvi);
    # display(surface(ext_meshxi, meshvi, ext_eval_fi));
    display(surface(meshx, meshvi, eval_fi));
    # return ext_meshxi, meshvi, ext_eval_fi
    return meshvi, eval_fi
end

function plot_char(phix)
    meshve = collect(LinRange(-ve_bound,ve_bound,100))
    liape = meshve.^2 ./ 2.0 .- 1/mu .* phix'
    meshvi = collect(LinRange(-vi_bound,vi_bound,100))
    liapi = meshvi.^2 ./ 2.0 .+ phix'

    plot(
        contourf(liape, size=(500,600), title="electrons"),
        contourf(liapi, size=(500,600), title="ions"),
        layout=[1,1]
    )
end

function plot_char(phix, meshve, meshvi)
    liape = meshve.^2 ./ 2.0 .- 1/mu .* phix'
    liapi = meshvi.^2 ./ 2.0 .+ phix'

    plot(
        contourf(liape, size=(500,600), title="electrons"),
        contourf(liapi, size=(500,600), title="ions"),
        layout=[1,1]
    )
end

################################
# Diagnostics
################################

"""
    Returns an approximation of the L1-norm of 
        |v ∂ₓf + speed(v) ∂ᵥf - source(x,v)|.
    Discretization by cartesian centered differences (only inside the domain)
"""
function norm_equation_advection(meshx, meshv, f, speed, source, mask)
    dx = meshx[begin+2:end] - meshx[begin:end-2] # to be horizontal
    dv = meshv[begin+2:end] - meshv[begin:end-2] # to be vertical

    Dxf = (f[begin+1:end-1, begin+2:end] - f[begin+1:end-1, begin:end-2]) ./ dx'
    Dvf = (f[begin+2:end, begin+1:end-1] - f[begin:end-2, begin+1:end-1]) ./ dv
    
    to_integrate = meshv[begin+1:end-1] .* Dxf + speed[begin+1:end-1]' .* Dvf - source[begin+1:end-1,begin+1:end-1]
    return sum(dx' .* abs.(to_integrate) .* mask[begin+1:end-1,begin+1:end-1] .* dv)
end

"""
    Computes Jₑ(x) := ∫ᵥ v fₑ(x,v) dv for x ∈ meshx.
"""
function get_Je(f_eb, phix)
    Je = 0.0 * phix;
    for (i,phixi) in enumerate(phix)
        # mesh of the segment such that fe_value is not null
        lowerbound = min(get_v_char_e(0.0,phix[end],phixi), -ve_bound);
        Ne = ceil(Int, -lowerbound*length(phix))+2; # proportional to Nx (with |[0,1]|=1)
        vv = collect(LinRange(lowerbound,0.0,Ne))
        # use the characteristics to determine the values of fₑ
        to_integrate = f_eb(get_v_char_e(vv,phixi,0.0))
        to_integrate = vcat(vv .* to_integrate, (-vv[end-1:-1:begin]) .* to_integrate[end-1:-1:begin])
        # integral (using the assumption of symmetry over fₑ)
        Je[i] = trapezes(to_integrate, vv[2]-vv[1])
    end
    return Je
end

"""
    Computes Jᵢ(x) := ∫ᵥ v fᵢ(x,v) dv for x ∈ meshx.
"""
function get_Ji(f_eb, phix)
    Ji = 0.0 * phix; 

    Jiloc_old = 0.0*Ji; viloc_old = 0.0*phix;
    Jiloc = 0.0*Ji; viloc = 0.0*phix;

    # get the meshes
    upper_mesh, lower_mesh = get_fi_char_meshes(phix)

    # lower mesh part : all coordinates of nᵢ are updated
    for (ivm, vm) in enumerate(vcat(lower_mesh,upper_mesh))
        stop_index = max(1,ivm-length(lower_mesh)+1)

        Jiloc .*= 0.0;
        viloc[end] = vm;
        for (ix, x) in enumerate(meshx[end-1:-1:stop_index])
            jx = length(phix)-ix # real index in arrays
            viloc[jx] = get_v_char_i(vm, phix[end], phix[jx])
            dt = (x-meshx[jx+1]) * 2.0 / (viloc[jx] + viloc[jx+1])
            Jiloc[jx] = Jiloc[jx+1] + nu * dt * 0.5 .* (f_eb(get_v_char_e(viloc[jx],phix[jx],0.0)) + f_eb(get_v_char_e(viloc[jx+1],phix[jx+1],0.0)))
        end

        if (ivm>2)
            # negative speeds : Jᵢ += Δvₖ * (vₖ fᵢₖ + vₖ₋₁ fᵢₖ₋₁) / 2
            Ji[stop_index:end] .+= (viloc[stop_index:end] .- viloc_old[stop_index:end]) .* 0.5 .* 
                                   (viloc[stop_index:end] .* Jiloc[stop_index:end] .+ viloc_old[stop_index:end] .* Jiloc_old[stop_index:end]) 
            # positive speeds : fᵢ(x,v) + fᵢ(x,-v) = cte = 2 * fᵢ(x₀,0)
            Ji[stop_index:end] .+= (viloc_old[stop_index:end] .- viloc[stop_index:end]) .* 0.5 .* 
                                   (viloc[stop_index:end] .* (2.0*Jiloc[stop_index] .- Jiloc[stop_index:end]) .+ 
                                   viloc_old[stop_index:end] .* (2.0*Jiloc_old[stop_index] .- Jiloc_old[stop_index:end])) 
        end

        Jiloc_old = 1*Jiloc; viloc_old = 1*viloc;
    end
    return Ji
end

"""
    By integration of the advection equation on fᵢ, we obtain
        ∂ₓJᵢ(x) = ν nₑ(x)   ⟹   Jᵢ(x) - Jᵢ(0) = ∫_[0,x] ν nₑ(y) dy.
    Since Jᵢ(0) = ∫ᵥ v fᵢ(0,v) dv = 0 (as vfᵢ(0,v) is an impair function),
    this function compares Jᵢ and ∫_[0,x] ν nₑ(y) dy. 
"""
function check_Ji(f_eb, meshx, phix)
    Ji = get_Ji(f_eb, phix)
    ne = zeros(length(phix))
    update_ne!(ne, f_eb, phix)

    # Again trapezes integral 
    int_nu_ne = nu .* vcat([0.0], cumsum(0.5 .* (ne[begin+1:end] + ne[begin:end-1]) .* (meshx[begin+1:end] - meshx[begin:end-1])))
    p = plot(meshx, Ji, label="Jᵢ")
    plot!(p, meshx, int_nu_ne, label="∫νnₑ")
    display(p)

    println("Erreur relative entre Jᵢ et ∫_[0,x] ν nₑ(y)dy : ", maximum(abs.(Ji - int_nu_ne)) / maximum(max.(abs.(Ji), abs.(int_nu_ne))))
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

    println("$Neffective / $Nsamples effective evaluations")
    println("err : mean $mean_err, max $max_err")
    println("dx  : mean $mean_dx_fe, max $max_dx_fe")
    println("dv  : mean $mean_dv_fe, max $max_dv_fe")

    file = open("test_df_fe.txt", "a")
    print(file, "$f_eb_str\t$lambda\t$mu\t$nu\t$cardinal\t$Nx\t$eps\t$Nsamples\t$mean_err\t$max_err\n")
    close(file)

    return errs
end

function evaluate_fi(f_eb, phi_data, Ncharmesh, x, v)
    # in all cases, we need the value of fᵢ(x,-|v|)
    meshx = collect(LinRange(x,1.0,2+ceil(Int,Ncharmesh*(1.0-x))))
    phix  = eval_phi(phi_data, meshx)
    phidx = eval_phi_dx(phi_data, meshx)
    v_1 = get_v_char_i(-abs(v),phix[begin],phix[end])
    fi, = integrate_char(f_eb, meshx, phix, v_1, 1)

    if (v>0) # bad case : we also need to integrate the whole characteristic
        # println("bruh")
        meshx = collect(LinRange(0.0,1.0,Ncharmesh))
        phix  = eval_phi(phi_data, meshx)
        phidx = eval_phi_dx(phi_data, meshx)
        ind = findlast(v^2/2.0 .+ eval_phi(phi_data,x) .- phix .<= 0.0)
        ind = (ind == nothing ? 1 : ind) # si toutes les coords de liap sont positives
        # println("ind : $ind")
        fitot, = integrate_char_outofmeshx(f_eb, meshx, phi_data, phix, v_1, ind)
        # println("fitot : $fitot")

        fi = 2.0 * fitot - fi; # uses fᵢ(x,v) + fᵢ(x,-v) = 2.0 * fᵢ(x₀,0)
    end
    return fi
end

"""
    Draws uniform points in (x,v) ∈ [0,1]×[-vi_bound,vi_bound], and checks
    whether the advection equation on fᵢ is satisfied.
"""
function test_df_fi(f_eb, phi_data, Nsamples, eps, Ncharmesh, thresh_liap)
    xrand = zeros(Nsamples,Nsamples); 
    vrand = zeros(Nsamples,Nsamples); 
    errs  = zeros(Nsamples,Nsamples);

    max_err = 0.0;   mean_err = 0.0;
    max_dx_fi = 0.0; mean_dx_fi = 0.0;
    max_dv_fi = 0.0; mean_dv_fi = 0.0;
    Neffective = 0
    for (ix,x) in enumerate(LinRange(0.0,1.0,Nsamples))
        for (iv,v) in enumerate(LinRange(-vi_bound, vi_bound, Nsamples))
            liap = v^2 / 2.0 + eval_phi(phi_data, x)
            # if we are in the right bounds
            if ((eps < x < 1.0 - eps) && (- vi_bound + eps < v < vi_bound - eps) && (abs(liap) >= thresh_liap))
                dx_fi = (evaluate_fi(f_eb,phi_data,Ncharmesh,x+eps,v) - evaluate_fi(f_eb,phi_data,Ncharmesh,x-eps,v)) / (2.0 * eps)
                dv_fi = (evaluate_fi(f_eb,phi_data,Ncharmesh,x,v+eps) - evaluate_fi(f_eb,phi_data,Ncharmesh,x,v-eps)) / (2.0 * eps)
                source = nu * f_eb(get_v_char_e(v,eval_phi(phi_data,x),0.0))
                err = abs(v * dx_fi - eval_phi_dx(phi_data, x) * dv_fi - source)

                # lhs = v * dx_fi - eval_phi_dx(phi_data, x) * dv_fi
                # rhs = source
                # println("lhs : $lhs, rhs : $rhs, liap : $liap")

                mean_err += err; max_err = max(max_err, err)
                mean_dx_fi += abs(dx_fi); max_dx_fi = max(max_dx_fi, abs(dx_fi))
                mean_dv_fi += abs(dv_fi); max_dv_fi = max(max_dv_fi, abs(dv_fi))

                Neffective += 1
                xrand[iv,ix] = x; vrand[iv,ix] = v; errs[iv,ix] = err;
            end
        end
    end
    if (Neffective>0)
        mean_err = mean_err / Neffective
        mean_dx_fi = mean_dx_fi / Neffective
        mean_dv_fi = mean_dv_fi / Neffective
    end

    println("$Neffective / $Nsamples^2 effective evaluations")
    println("err : mean $mean_err, max $max_err")
    println("dx  : mean $mean_dx_fi, max $max_dx_fi")
    println("dv  : mean $mean_dv_fi, max $max_dv_fi")
    # println("$eps\t$Nsamples\t$Ncharmesh\t$thresh_liap\t$mean_err\t$max_err")

    file = open("test_df_fi.txt", "a")
    print(file, "$f_eb_str\t$lambda\t$mu\t$nu\t$cardinal\t$Nx\t$eps\t$Nsamples\t$Ncharmesh\t$thresh_liap\t$mean_err\t$max_err\t")
    close(file)

    return xrand, vrand, errs
end

################################
# Test case data
################################

function get_test_ni(meshx, beta1, beta2, Npoints)
    ni = 0.0 * meshx;
    for (ix, x) in enumerate(meshx)
        meshv = collect(LinRange(-sqrt(beta1^2+x^2), 0.0, Npoints))
        sum_value = 0.0 * meshv; # f_i(x,v) + f_i(x,-v)
        for (iv, v) in enumerate(meshv)
            var = abs(x^2-v^2)
            sum_value[iv] = 2.0 * (asinh.(sqrt.(0.5 .* max.(0.,beta1^2 .- var) ./ var)) .- asinh.(sqrt.(0.5 .* max.(0.,beta2^2 .- var) ./ var))) # integral along the whole char, also f_i(x,v) + f_i(x,-v)
        end
        ni[ix] = trapezes(sum_value, meshv[2] - meshv[1])
    end
    return nu .* ni
end

function unitary_test_fi(NN, NN2)
    meshx = collect(LinRange(0.0,1.0,NN));
    polybase = Polybase(1,"canonical",meshx);
    phi_data = Phi_data(polybase,[-1/2]);

    phix = eval_phi(phi_data, meshx);
    beta1 = 0.8; beta2 = 0.4;
    f_ebb = (v) -> 1.0 * (beta2^2 .<= v.^2 .<= beta1^2);

    # computed one
    ni = zeros(size(meshx));
    update_ni!(ni, f_ebb, phix, meshx, beta1, beta2);
    # xx, vv, fi = update_ni_bruteforce!(ni, f_ebb, phix, meshx, beta1, beta2);
    # display(scatter(xx,vv,marker_z=fi))

    # # analytical one
    ni_test = get_test_ni(meshx, beta1, beta2, NN2);

    # display(plot( meshx[2:end-20], ni[2:end-20], markers=true, ms=0.5, label="numerical"))
    # display(plot!(meshx[2:end-20], ni_test[2:end-20], markers=true, ms=0.5, label="analytic"))
    # # return meshx, ni, ni_test

    println("erreur : ", max(abs.(ni[2:end] .- ni_test[2:end])...))

    # x_1 = 1.0;
    # v_1 = -(sqrt(2.0) + 1.0)/2;
    # int_code, = integrate_char(f_eb, meshx, phix, v_1, 1);
    # int_test = nu * asinh(sqrt(0.5*(beta^2 - (v_1^2 - x_1^2)) / (v_1^2 - x_1^2)));

    # error = abs(int_code - int_test)
    # println("int code : $int_code, int test : $int_test, error : ", error)
    # return error
    # return xx, vv, fi
end

################################
# Main
################################

function fixed_point_2sp(Nx::Int)
    println("\nBegin of Nx=$Nx")
    meshx = LinRange(0.0,1.0,Nx) |> collect; # space mesh.
    polybase = Polybase(cardinal, basis, meshx);
    phi_data = Phi_data(polybase, coeff);
    
    ni = zeros(size(meshx)); #      ion density (depends on space only)
    ne = zeros(size(meshx)); # electron density (depends on space only)

    it=1;
    phix_old = eval_phi(phi_data, meshx);
    phix = 1*phix_old;
    variation = 1.0; # Linf-norm of (phix - phix_old)

    while ((it <= Nit_max) && (variation > 1e-6))
        println("Iteration $it / $Nit_max")
        # phidx = eval_phi_dx(phi_data, meshx)

        update_ne!(ne, f_eb, phix);      
        # display(plot(meshx, ne))
        # readline()
        update_ni!(ni, f_eb, phix, meshx); 
        # display(plot(meshx, ni))
        # readline()
        # update_phi!(phi_data, ni, ne); 
        newcoeff = update_phi(phi_data, ni, ne);
        newphix = eval_phi(phi_data.polybase, newcoeff, meshx)
        if (newphix[end] <= 0.0)
            phi_data.coeff[:] = newcoeff[:];
        else
            prop = 0.5*phix[end] / newphix[end];
            phi_data.coeff[:] .+= prop .* newcoeff[:];
        end
        phix = eval_phi(phi_data, meshx);
        # display(plot(meshx, phix))
        # readline()

        variation = max(abs.(phix - phix_old)...)
        println("Norme ∞ de phi - phi_prec : ", variation)
        phix_old = 1*phix;
        it += 1;
    end

    # println("Error evaluation...")
    # meshv = collect(LinRange(-ve_bound,ve_bound,2*Nx))
    # EE = - eval_phi_dx(phi_data, meshx)
    # eval_fe = eval_cartmesh_fe(f_eb, phi_data, meshx, meshv)
    # eval_fi = eval_cartmesh_fi(f_eb, phi_data, meshx, meshv)
    # errfe = norm_equation_advection(meshx, meshv, eval_fe, - 1/mu * EE, 0. *eval_fe, eval_fe .== eval_fe)
    # liapi = meshv.^2 ./ 2.0 .+ eval_phi(phi_data, meshx)'
    # errfi = norm_equation_advection(meshx, meshv, eval_fi, EE, nu*eval_fe, abs.(liapi) .> 0.2)
    println("End of Nx=$Nx\n")

    return meshx, phi_data#, errfe, errfi
end 

function fixed_point_1sp(Nx::Int)
    println("\nBegin of Nx=$Nx")
    meshx = LinRange(0.0,xmax,Nx) |> collect; # space mesh.
    polybase = Polybase(cardinal, basis, meshx); 
    phi_data = Phi_data(polybase, coeff);
    
    n = zeros(size(meshx)); # density (depends on space only)

    it=1;
    phix_old = eval_phi(phi_data, meshx);
    phix = 1*phix_old;
    variation = 1.0; # Linf-norm of (phix - phix_old)

    while ((it <= Nit_max) && (variation > 1e-9))
        println("Iteration $it / $Nit_max")

        update_n!(n, f_b, phix); 
        update_phi_1sp!(phi_data, n);
        phix = eval_phi(phi_data, meshx);

        variation = max(abs.(phix - phix_old)...)
        println("Norme ∞ de phi - phi_prec : ", variation, ", erreur à la solution analytique : ", phi_data.coeff[1] - phi_sol_coeff)
        phix_old = 1*phix;
        it += 1;
    end

    println("End of Nx=$Nx\n")
    return meshx, phi_data
end 


# function main()
    println("Welcome in fixedpoint.")
    # Nxs = 20:30:200;
    Nxs = 100;
    meshxs = []
    phi_datas = []

    one_sp = false
    two_sp = true

    if (two_sp)
        errors = zeros(length(Nxs),2)
        for (iNx, Nx) in enumerate(Nxs)
            # meshx, phi_data, errors[iNx,1], errors[iNx,2] = fixed_point_2sp(Nx)
            meshx, phi_data = fixed_point_2sp(Nx)
            push!(meshxs, meshx)
            push!(phi_datas, phi_data)
        end
        # println("Errors : ")
        # println(errors)
        if (length(Nxs)>1)
            for iphi in 2:length(Nxs)
                println("norme phi $iphi - prec : ", maximum(abs.(phi_datas[iphi].coeff-phi_datas[iphi-1].coeff)))
            end
        end
    end

    if (one_sp)
        for (iNx, Nx) in enumerate(Nxs)
            meshx, phi_data, = fixed_point_1sp(Nx)
            push!(meshxs, meshx)
            push!(phi_datas, phi_data)
        end
    end

    println("Bye")
# end

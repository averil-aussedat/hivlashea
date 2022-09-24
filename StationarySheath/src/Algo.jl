
####################################
# Integrals along ion chars
####################################

"""
    $(SIGNATURES)

    Compute the integral of f_e along the ion characteristic reaching ``(x_b=0,v_b<0)`` at ``\\tau=0``.
    Uses 
    ```math
        \\int_{\\tau=-\\infty}^{0} f_e(x(\\tau),v(\\tau))
        = \\int_{z=0}^{1} \\frac{f_e(z,-\\sqrt{v_b^2 - 2*\\phi(z)})}{-\\sqrt{v_b^2 - 2*\\phi(z)}} dz
    ```
    and a trapeze integral on the points ``(z,-g_z(v_b))`` s.t. ``z\\in``meshx.
    Preallocated vectors vi and toint will be modified by the function.
"""
function integral_char_xb0!(vi::Vector{Float64}, toint::Vector{Float64}, params::Params, vb::Float64, phix::Vector{Float64}, f_eb::Feb)
    get_v_char_i!(vi, vb, 0.0, phix)
    toint .= eval_fe(f_eb, params, phix, vi) ./ (-vi)
    integral = params.meshx.step * (0.5 * (toint[begin]+toint[end]) + sum(toint[begin+1:end-1]))
    return integral
end

"""
    $(SIGNATURES)

    Special case of the critical characteristic:
    use of a rectangle integration scheme to avoid the equilibrium point ``(x=0,v=0)``.
"""
function integral_char_00!(vi::Vector{Float64}, toint::Vector{Float64}, params::Params, phix::Vector{Float64}, f_eb::Feb)
    get_v_char_i!(vi, 0.0, 0.0, phix)
    toint[2:end] .= eval_fe(f_eb, params, phix[2:end], vi[2:end]) ./ (-vi[2:end])
    integral = params.meshx.step * sum(toint[2:end])
    return integral
end

"""
    $(SIGNATURES)
    
    Compute the integral of f_e along the ion characteristic reaching ``(x_b>0,v_b=0)`` at ``\\tau=0``.
    Here ``x_b`` is known through ``\\phi(x_b)`` only. Uses 
    ```math
        \\int_{\\tau=-\\infty}^{0} f_e(x(\\tau),v(\\tau))
        = \\int_{v=vmin}^{0} \\frac{f_e(x(v),v)}{-\\phi'(x(v))} dv
    ```
    and a trapeze integral on the points ``(x(v),v)`` s.t. ``x(v)\\in``meshx.
    Preallocated vectors vi and toint will be modified by the function.
"""
function integral_char_vb0!(vi::Vector{Float64}, toint::Vector{Float64}, params::Params, ixb::Int, phix::Vector{Float64}, phidx::Vector{Float64}, f_eb::Feb)
    get_v_char_i!(vi[ixb:end], 0.0, phix[ixb], phix[ixb:end])
    toint[ixb:end] .= eval_fe(f_eb, params, phix[ixb:end], vi[ixb:end]) ./ (-phidx[ixb:end])
    integral = params.meshx.step * (0.5 * (toint[ixb]+toint[end]) + sum(toint[ixb+1:end-1]))
    return integral
end

####################################
# Densities
####################################

export compute_ne!, compute_ni!

"""
    $(SIGNATURES)
    
    Computation of 
    ```math
        n_e(x) = \\int_{v\\in\\mathbb{R}} f_e(x,v) dv = 2 \\int_{v \\leqslant 0} f_e(x,v) dv.
    ```
"""
function compute_ne!(ne::Vector, chare::Vector{Float64}, tointe::Vector{Float64}, params::Params, f_eb::Feb, phix::Vector{Float64})
    for (i,phixi) in enumerate(phix)
        # chare ranges from [-some bound, 0], so that chare[end]=0
        chare .= LinRange(get_v_char_e(params, params.meshx.mesh[i], phix[end], 0.0), 0.0, params.Nve+1)
        eval_fe!(tointe, f_eb::Feb, params::Params, phixi, chare)
        # trapeze integral using the assumption f_eb(v)=f_eb(-v)
        ne[i] = (chare[2]-chare[1]) * (tointe[begin] + 2.0 * sum(tointe[2:end-1]) + tointe[end])
    end
end

"""
    $(SIGNATURES)


    Computation of 
    ```math
        n_i(x) = \\int_{v\\in\\mathbb{R}} f_i(x,v) dv = 2 \\int_{v \\leqslant 0} f_i(x_b(x,v),v_b(x,v)) dv.
    ```
"""
function compute_ni!(ni::Vector, chari::Vector{Float64}, vi::Vector{Float64}, vip::Vector{Float64}, toint::Vector{Float64}, 
                     params::Params, f_eb::Feb, phix::Vector{Float64}, phidx::Vector{Float64})
    ni .= 0.0
    chari .= LinRange(get_v_char_e(params,0.0,phix[end],0.0), 0.0, params.Nvi+1)

    # first part: the end of the characteristic is on (x=0,v<0). p is for previous
    integralp = integral_char_xb0!(vip, toint, params, chari[1], phix, f_eb)
    for vb in chari[2:end-1]
        integral = integral_char_xb0!(vi, toint, params, vb, phix, f_eb)
        ni .+= (vi .- vip) .* (integral + integralp) # times 2 is embedded
        integralp = integral; vip .= vi # might do better by letting vi[p] in place and changing index
    end

    # special case of (x=0,v=0)
    integral = integral_char_00!(vi, toint, params, phix, f_eb)
    ni .+= (vi .- vip) .* (integral + integralp) # times 2 is embedded
    integralp = integral; vip .= vi 

    # last part: the end of the characteristic is on (x>0,v=0)
    for ixb=2:params.meshx.NN+1
        integral = integral_char_vb0!(vi, toint, params, ixb, phix, phidx, f_eb)
        ni[ixb:end] .+= (vi[ixb:end] .- vip[ixb:end]) .* (integral + integralp) # times 2 is embedded
        integralp = integral; vip[ixb:end] .= vi[ixb:end]
    end

    ni .*= params.ν
end

####################################
# Poisson problem
####################################

export solve_Poisson!

"""
    $(SIGNATURES)

    Solve of the Poisson problem 
    ```math
        \\begin{cases}
            - \\lambda^2 \\Delta \\phi(x) = n_i(x) - n_e(x), \\\\
            \\phi(0) = \\phi'(0) = 0.
        \\end{cases}
    ```
"""

"""
    $(SIGNATURES)

    Resolution by least squares when ``\\phi`` is decomposed in a polynomial basis.
"""
function solve_Poisson!(phi::Phidata_poly, params::Params, ni::Vector{Float64}, ne::Vector{Float64})
    phi.coeffs .= - 1.0/(params.λ^2) .* (phi.Vandermonde \ (ni - ne))
end

"""
    $(SIGNATURES)

    Resolution by quadrature when ``\\phi`` is approximated pointwise.
"""
function solve_Poisson!(phi::Phidata_grid, params::Params, ni::Vector{Float64}, ne::Vector{Float64})
    # trapeze integration
    toint = params.meshx.step * 0.5 * ((ni[begin:end-1] + ni[begin+1:end]) - (ne[begin:end-1] + ne[begin+1:end]))
    phi.values[begin] = 0.0
    phi.values[begin+1:end] .= - cumsum(cumsum(toint)) ./ (params.λ^2)
end

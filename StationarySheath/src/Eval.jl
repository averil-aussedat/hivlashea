####################################
# Infinitesimal energies
####################################

export 𝓛ₑ, 𝓛ᵢ

"""
$(SIGNATURES)

Electric infinitesimal energy.
"""
function 𝓛ₑ(params::Params, phix, v)
    return v.^2 / 2 .- 1/params.μ .* phix
end

"""
$(SIGNATURES)

Ion infinitesimal energy.
"""
function 𝓛ᵢ(phix, v) # bscrL
    return v.^2 / 2 .+ phix
end

####################################
# Polynomial phi
####################################

export evaluate_poly, eval_phi, eval_phi_dx

"""
$(SIGNATURES)

Return the polynomial 
```math
    a_1 + x (a_2 + x \\dots (a_{n-1} + x a_n))
```
with `eval```=[a_1,\\dots,a_n]``.
"""
function evaluate_poly(eval, x)
    res = 0.0*x; # higher degree coefficient
    for ee=eval
        res = x .* res .+ ee;
    end
    return res;
end

"""
$(SIGNATURES)

Value of ϕ at `x`.
"""
function eval_phi(phi::Phidata_poly, params::Params, x)
    return sum((phi.coeffs .* [evaluate_poly(horn, x) for horn = phi.horner])')'
end

"""
$(SIGNATURES)

Value of ϕ' at `x`.
"""
function eval_phi_dx(phi::Phidata_poly, params::Params, x)
    return sum((phi.coeffs .* [evaluate_poly(horn, x) for horn = phi.horner_dx])')'
end

####################################
# FD phi
####################################

"""
$(SIGNATURES)

Value of ϕ at `x`.
"""
function eval_phi(phi::Phidata_grid, params::Params, x)
    # linear interpolation 
    ix = max.(1,min.(params.meshx.NN,floor.(Int,x./params.meshx.step).+1))
    α = (x .- params.meshx.mesh[ix]) ./ params.meshx.step
    return phi.values[ix] .* (1 .- α) .+ phi.values[ix.+1] .* α
end

"""
$(SIGNATURES)

Value of ϕ' at `x`.
"""
function eval_phi_dx(phi::Phidata_grid, params::Params, x)
    # archibasic finite differences
    ix = max.(1,min.(params.meshx.NN,floor.(Int,x./params.meshx.step).+1))
    return (phi.values[ix.+1] .- phi.values[ix]) ./ params.meshx.step
end

####################################
# Characteristics
####################################

export get_v_char_e, get_v_char_e!

"""
$(SIGNATURES)

Return the negative ``v`` such that

```math
    \\frac{v^2}{2} - \\frac{1}{\\mu} \\phi(x) = \\frac{v_0^2}{2} - \\frac{1}{\\mu} \\phi(x_0)
```    
"""
function get_v_char_e(params::Params, v0, phix0, phix)
    res = v0.^2 .- 2.0./params.μ .* (phix0 .- phix)
    if (any(res .< -1e-10))
        throw(DomainError("negative value in get_v_char_e"))
    else
        res = - sqrt.(max.(0.0,res))
    end
    return res
end

function get_v_char_e!(res, params::Params, v0, phix0, phix)
    res .= v0.^2 .- 2.0./params.μ .* (phix0 .- phix)
    if (any(res .< -1e-10))
        throw(DomainError("negative value in get_v_char_e"))
    else
        res .= - sqrt.(max.(0.0,res))
    end
end

export get_v_char_i, get_v_char_i!

"""
$(SIGNATURES)

Return the negative ``v`` such that

```math
    \\frac{v^2}{2} + \\phi(x) = \\frac{v_0^2}{2} + \\phi(x_0)
```
"""
function get_v_char_i(v0, phix0, phix)
    res = v0.^2 .+ 2.0 .* (phix0 .- phix)
    if (any(res .< -1e-10))
        throw(DomainError("negative value in get_v_char_i"))
    else
        res = - sqrt.(max.(0.0,res))
    end
    return res
end

function get_v_char_i!(res, v0, phix0, phix)
    res .= v0.^2 .+ 2.0 .* (phix0 .- phix)
    if (any(res .< -1e-10))
        throw(DomainError("negative value in get_v_char_i"))
    else
        res .= - sqrt.(max.(0.0,res))
    end
end

####################################
# electron density
####################################

export eval_fe!, eval_fe

"""
$(SIGNATURES)

Return the value ``f_e(x,v)``. 

Note that the non-emitting boundary condition has precedence over
the prescribed value ``f_{e,b}``, and the result will be 0
whenever 𝓛ₑ(x,v) ⩾ 𝓛ₑ(1,0).
"""
function eval_fe!(res, f_eb::Feb, params::Params, phix, v)
    res .= 𝓛ₑ(params, phix, v) # temp: electrical energy
    res .= f_eb.value.(-sqrt.(2 .* res)) .* (res .<= f_eb.Le_threshold)
end

function eval_fe(f_eb::Feb, params::Params, phix, v)
    Le = 𝓛ₑ(params, phix, v) 
    return f_eb.value.(-sqrt.(2*Le)) .* (Le .<= f_eb.Le_threshold)
end

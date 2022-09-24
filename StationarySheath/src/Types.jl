export Mesh

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct Mesh
    "Minimal value"
    min         :: Float64
    "Maximal value"
    max         :: Float64
    "Step"
    step        :: Float64
    "Number of intervals"
    NN          :: Int
    "Mesh itself"
    mesh        :: Vector{Float64}

    function Mesh(min::Float64, max::Float64, NN::Int)
        new(min, max, (max-min)/NN, NN, collect(LinRange(min,max,NN+1)))
    end
end

export Params

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct Params 

    # Physical parameters
    "Debye length"
    λ           :: Float64 
    "Mass ratio"
    μ           :: Float64 
    "Ionization factor"
    ν           :: Float64

    # Numerical parameters
    "Space mesh"
    meshx       :: Mesh
    "Number of points in ion speed mesh"
    Nvi         :: Int
    "Number of points in electron speed mesh"
    Nve         :: Int
    "Maximal number of fixed-point iterations"
    Nit_max     :: Int

    function Params(;λ=1.0, μ=1.0, ν=42.0, Nx=42, Nvi=42, Nve=42, Nit_max=100)
        meshx  = Mesh(0.0,1.0,Nx)
        new(λ, μ, ν, meshx, Nvi, Nve, Nit_max)
    end

end

export Phidata, Phidata_poly, Phidata_grid

abstract type Phidata end

"""
$(TYPEDEF)

Polynomial representation of `ϕ`. Basis are
    * canonical: starting from degree 2. 
    * pair: only pair expononents of the canonical basis.

$(TYPEDFIELDS)
"""
struct Phidata_poly <: Phidata
    "Cardinal of the basis"
    cardinal    :: Int 
    "Basis choice"
    basis       :: String
    "Array of Horner coefficients for x -> p(x)"
    horner      :: Vector{Vector{Float64}} 
    "Array of Horner coefficients for x -> p'(x)"
    horner_dx   :: Vector{Vector{Float64}} 
    "Matrix of the linear system for the Poisson solver"
    Vandermonde :: Matrix{Float64} 
    "Coefficients in the basis"
    coeffs      :: Vector

    function Phidata_poly(cardinal::Int, basis::String, meshx::Mesh, coeffs::Vector)
        if (length(coeffs)!=cardinal)
            throw(ErrorException("initial coefficients of size $(length(coeffs)) for cardinal $cardinal"))
        end
        if (basis=="canonical")
            # p(x) = a_n + x (a_{n-1} + x (a_{n-2} + ... x (a_2 + a_1 x)))
            # starts from degree 2
            horner = [[1.0*(i==j-2) for j=cardinal+2:-1:1] for i=1:cardinal]
            horner_dx = [[1.0*j*(i==j-1) for j=cardinal+2:-1:1] for i=1:cardinal]
            horner_dxx = [[1.0*j*(j+1)*(i==j) for j=cardinal+2:-1:1] for i=1:cardinal]
        elseif (basis=="pair")
            horner = [[1.0*(i==j-2) for j=2*cardinal+2:-1:1] for i=1:2:2*cardinal]
            horner_dx = [[1.0*j*(i==j-1) for j=2*cardinal+2:-1:1] for i=1:2:2*cardinal]
            horner_dxx = [[1.0*j*(j+1)*(i==j) for j=2*cardinal+2:-1:1] for i=1:2:2*cardinal]
        else
            throw(ArgumentError("Unknown phi basis '$basis'"))
        end
        Vandermonde = zeros(meshx.NN+1, cardinal)
        [Vandermonde[:,j] = evaluate_poly(eval_dxx,meshx.mesh) for (j,eval_dxx) = enumerate(horner_dxx)]
        new(cardinal, basis, horner, horner_dx, Vandermonde, coeffs)
    end

    function Phidata_poly(cardinal::Int, basis::String, meshx::Mesh)
        Phidata_poly(cardinal, basis, meshx, -20.0 .* (1:cardinal .== 1))
    end
end

"""
$(TYPEDEF)

Representation of `ϕ` on a 1D regular mesh. 

$(TYPEDFIELDS)
"""
struct Phidata_grid <: Phidata
    "Pointwise values on the space mesh"
    values      :: Vector{Float64}

    function Phidata_grid(meshx::Mesh)
        Phidata_grid(meshx, - 5 .* meshx.mesh.^2)
    end

    function Phidata_grid(meshx::Mesh, values::Vector{Float64})
        if (length(values)!=meshx.NN+1)
            throw(ErrorException("initial coefficients of size $(length(values)) for $(meshx.NN+1) points"))
        end
        new(values)
    end
end

export Feb

"""
$(TYPEDEF)

Boundary condition ``f_{e,b}`` for the electron density at ``(x=0,v⩽0)``.

$(TYPEDFIELDS)
"""
mutable struct Feb
    "Maximal electric energy ``\\mathcal{L}_e(x,v)`` s.t. ``(x,v)\\in``supp``(f_{e,b})``"
    Le_threshold::Float64
    "Prescribed value at ``(x=0,v\\leqslant 0)``"
    value::Function

    function Feb(phi::Phidata, params::Params)
        # maxwellian distribution
        # new(𝓛ₑ(params,eval_phi(phi,params,1.0),0.0) , x -> sqrt(params.μ/(2π)) .* exp.(- params.μ .* x.^2 ./ 2))
        
        # constant distribution (not compactly supported ! Just for tests)
        # new(𝓛ₑ(params,eval_phi(phi,params,1.0),0.0) , x -> 0*x .+ 1.0)

        # # compacty supported constant distribution
        # new(𝓛ₑ(params,eval_phi(phi,params,1.0),0.0) , x -> 1.0.*(abs.(x) .<= 5.0))

        # # compacty supported (initially) continuous distribution
        # new(𝓛ₑ(params,eval_phi(phi,params,1.0),0.0) , x -> max.(0.0,(5.0-abs.(x))/5.0))

        # compacty supported null-around-0 distribution
        new(𝓛ₑ(params,eval_phi(phi,params,1.0),0.0) , x -> 1.0.*(3.0 .<= abs.(x) .<= 5.0))
    end
end
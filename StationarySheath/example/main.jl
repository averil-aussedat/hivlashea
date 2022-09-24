"""
    Implementation of a fixed-point algorithm for 
    the stationary Vlasov-Poisson system.

    ```math
        \\begin{cases}
            v \\partial_x f_e(x,v) + \\frac{1}{\\mu} \\phi''(x) \\partial_v f_e(x,v) = 0 \\\\
            v \\partial_x f_i(x,v) - \\phi'(x) \\partial_v f_i(x,v) = \\nu f_e(x,v) \\\\
            - \\lambda^2 \\phi''(x) = n_i(x) - n_e(x)
        \\end{cases}
    ```
    with the densities ``n_s``, ``s\\in\\{i,e\\}`` defined as
    ```math
        n_s(x) = \\int_{v\\in\\mathbb{R}} f_s(x,v) dv.
    ```

    Code initiated during CEMRACS 2022 by Michel Mehrenberger, Mehdi Badsi, 
    Anaïs Crestetto and Averil Prost.
"""

using StationarySheath

function run(; eps=1e-10)

    params = Params(Nx=100,Nve=100,Nvi=100) # only keywords arguments for this one
    # phi = Phidata_grid(params.meshx)
    phi = Phidata_poly(5, "pair", params.meshx)
    f_eb = Feb(phi,params)

    # physical quantities
    ni = zeros(params.meshx.NN+1)
    ne = zeros(params.meshx.NN+1)
    phix  = eval_phi(phi, params, params.meshx.mesh) # eval of phi on meshx
    phidx = eval_phi_dx(phi, params, params.meshx.mesh) # eval of phi' on meshx

    # allocations of temporary variables
    vi     = zeros(params.meshx.NN+1) # temp velocities along ion char
    vip    = zeros(params.meshx.NN+1) # temp velocities along ion char (previous)
    chari  = zeros(params.Nvi+1) # temp velocities at x=0 for ions
    chare  = zeros(params.Nve+1) # temp velocities at x=z for electrons
    tointe = zeros(params.Nve+1) # electron velocities to integrate
    toint  = zeros(params.meshx.NN+1) # to integrate on space
    phixp  = copy(phix) # previous phix, for stopping criterion

    goOn = true; n=0
    while (goOn && n<params.Nit_max)

        compute_ni!(ni, chari, vi, vip, toint, params, f_eb, phix, phidx)
        compute_ne!(ne, chare, tointe, params, f_eb, phix)
        solve_Poisson!(phi, params, ni, ne)
        
        phix  .=    eval_phi(phi, params, params.meshx.mesh)
        phidx .= eval_phi_dx(phi, params, params.meshx.mesh)

        goOn = diags_loop(params, phi, ni, ne, phix, phixp, n, eps)
        phixp .= phix
        f_eb.Le_threshold = 𝓛ₑ(params, phix[end], 0.0)

        n += 1

    end # main loop

    diags_end(params, phi, f_eb, ni, ne, phix, phidx, n)

    ploot(params.meshx, phix, "phi")
    ploot(params.meshx, ni, "ni")
    ploot(params.meshx, ne, "ne")
end

@time run()
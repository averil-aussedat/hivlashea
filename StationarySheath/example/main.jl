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
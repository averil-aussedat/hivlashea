####################################
# Reconstruction of fᵢₑ
####################################

export build_fe, build_fi

"""
    $(SIGNATURES)

    Evaluation of ``f_e`` on the regular grid meshx × meshve. 

    The speed mesh is computed to fit 0 ⩽ 𝓛ₑ ⩽ 𝓛ₑ(1,0).
"""
function build_fe!(fe::Matrix, params::Params, f_eb::Feb, phix::Vector{Float64})
    meshve = collect(LinRange(get_v_char_e(params, 0.0, phix[end], 0.0), 0.0, params.Nve+1))
    for (i, phixi) in enumerate(phix)
        @views eval_fe!(fe[:,i], f_eb, params, phixi, meshve)
    end
end

#   zᵢ   zᵢ₋₁
#   __↑_____._____.__  
#     |\
#   l | \_  upper 
#   o |   \_.
#   w x      \___
#   e |\__       `--- critical characteristic
#   r |   \_._       
#     |       \___.
#     |            \_ 

"""
    $(SIGNATURES)
    
    Approximation of ``f_i`` on a regular grid meshx × meshvi.

"""
function build_fi!(fi::Matrix, params::Params, f_eb::Feb, phix::Vector{Float64}, phidx::Vector{Float64})
    meshvi = collect(LinRange(get_v_char_e(params, 0.0, phix[end], 0.0), 0.0, params.Nvi+1))

    # lower part of the mesh (0 included): (Nvi+1)*(Nx+1) points
    # upper part of the mesh: Nx*(Nx+1)/2
    number_of_point = round(Int,(params.Nvi+1)*(params.meshx.NN+1)+(params.meshx.NN+1)*params.meshx.NN/2)
    println("number of points : $number_of_point")
    scatteredpoints = zeros(2, number_of_point) # (x,v)
    scattered_fi    = zeros(number_of_point)
    iscat = 1 # running index along data

    # lower part of the mesh (critical included)
    for vb in meshvi
        twoLi = vb^2 # + 2 * phi(0), which is equal to 0 
        wold = get_v_char_i(vb, 0.0, phix[end])
        tointold = eval_fe(f_eb, params, phix[end], wold) / (twoLi - 2*phix[end])
        scatteredpoints[:,iscat] = [params.meshx.mesh[end],wold]
        iscat += 1
        for iz = params.meshx.NN:-1:1
            w = get_v_char_i(vb, 0.0, phix[iz])
            scatteredpoints[:,iscat] = [params.meshx.mesh[iz],w]
            if (abs(vb)>1e-6 && iz>1)
                toint = eval_fe(f_eb, params, phix[iz], w) / (twoLi - 2*phix[iz])
                scattered_fi[iscat] = scattered_fi[iscat-1] + params.meshx.step * 0.5 * (toint + tointold) 
                tointold = toint
            else
                scattered_fi[iscat] = scattered_fi[iscat-1] + params.meshx.step * tointold 
            end
            iscat += 1
            wold = w
        end
    end

    # upper part of the mesh
    for ixb in 2:params.meshx.NN+1
        wold = get_v_char_i(0.0, phix[ixb], phix[end])
        tointold = eval_fe(f_eb, params, phix[end], wold) / (-phidx[ixb])
        scatteredpoints[:,iscat] = [params.meshx.mesh[end],wold]
        iscat += 1
        for iz = params.meshx.NN:-1:ixb
            w = get_v_char_i(0.0, phix[ixb], phix[iz])
            scatteredpoints[:,iscat] = [params.meshx.mesh[iz],w]
            toint = eval_fe(f_eb, params, phix[iz], w) / (-phidx[ixb])
            scattered_fi[iscat] = scattered_fi[iscat-1] + (w - wold) * 0.5 * (toint + tointold) 
            iscat += 1
            wold = w; tointold = toint
        end
    end

    scattered_fi .*= params.ν

    # interpolation 
    XX = repeat(params.meshx.mesh, params.Nvi+1)[:]
    VV = repeat(meshvi', params.meshx.NN+1)[:]
    gridpoints = [XX VV]'
    itp = interpolate(Shepard(), scatteredpoints, scattered_fi)
    interpolated = evaluate(itp, gridpoints)
    fi .= reshape(interpolated, params.meshx.NN+1, params.Nvi+1)'

    # return scatteredpoints, scattered_fi
end

####################################
# Plots
####################################

export ploot

"""
    $(SIGNATURES)
    
    Plot and wait for user input.
"""
function ploot(mesh::Mesh, value::Vector{Float64}, tag::String; erase::Bool=true)
    println("Plotting $tag (press enter to go on)")
    if (erase)
        display(plot(mesh.mesh,value,legend=false))
    else
        display(plot!(mesh.mesh,value,legend=false))
    end
    readline()
end

####################################
# Aggregation of diagnostics
####################################

export diags_loop, diags_end

"""
    $(SIGNATURES)
    
    Aggregation of diagnostics performed during the fixed-point loop.
"""
function diags_loop(params::Params, phi::Phidata, ni::Vector{Float64}, ne::Vector{Float64}, 
    phix::Vector{Float64}, phixp::Vector{Float64}, n::Int, eps::Float64)

    variation_∞ = maximum(abs.(phix - phixp)) / maximum(abs.(phixp)) 
    variation_1 = sum(abs.(phix - phixp)) / sum(abs.(phixp)) 
    goOn = ((variation_∞ > eps) && (variation_1 > eps))

    println("Iteration $n / $(params.Nit_max) : 
        variation ∞ $(variation_∞), 
        variation_1 $(variation_1)")

    return goOn
end

"""
    $(SIGNATURES)
    
    Aggregation of diagnostics performed at the end of the program.
"""
function diags_end(params::Params, phi::Phidata, f_eb::Feb, ni::Vector{Float64}, ne::Vector{Float64}, 
    phix::Vector{Float64}, phidx::Vector{Float64}, n::Int)

    fe = zeros(params.Nve+1, params.meshx.NN+1) # regular grid evaluation
    fi = zeros(params.Nvi+1, params.meshx.NN+1) # regular grid approximation

    build_fe!(fe, params, f_eb, phix)
    png(heatmap(fe), "data/fe.png")
    build_fi!(fi, params, f_eb, phix, phidx)
    png(heatmap(fi), "data/fi.png")

end
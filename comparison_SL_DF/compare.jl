using Plots 
using DelimitedFiles
using Printf

root_DF = "../../DynamicElectricSheath.jl/data/two_species"
root_SL = "../code_SL"

folders_DF = [
    # "run_comp_short_time_2sp_Nx1000_Nv2000_Nt6250"
    # "run_comp_short_time_2sp_Nx1000_Nv2000_Nt1563"
    "run_comp_long_time_2sp_Nx200_Nv400_Nt2500000"
]
folders_SL = [
    # "run_comp_short_time_2sp_Nx1000_Nvi2001_Nve2001_Nt6250"
    # "run_comp_short_time_2sp_Nx1000_Nvi2001_Nve2001_Nt1563"
    "run_comp_long_time_2sp_Nx1000_Nvi2001_Nve2001_Nt100000"
]

Elims = (-8.0,8.0)
ρlims = (-1.0,2.5)

for (folder_DF, folder_SL) in zip(folders_DF, folders_SL)

    Nt  = parse(Int,folder_SL[match(r"Nt",folder_SL).offset+2:end])

    # Electric field

    EDF = readdlm("$root_DF/$folder_DF/EE.dat")
    SLEData = readdlm("$root_SL/$folder_SL/data_output/E$(@sprintf("%06d",Nt)).dat")
    ESL = SLEData[2:end,2]

    println("min / max E DF : $(minimum(EDF)) / $(maximum(EDF))")
    println("min / max E SL : $(minimum(ESL)) / $(maximum(ESL))")

    xDF = LinRange(-1.0, 1.0, length(EDF))
    xSL = SLEData[2:end,1]
    p = plot(xDF, EDF, ylims=Elims, label="E (DF)", lw=4, legend=:bottomright, legendfontsize=14)
    plot!(p, xSL, ESL, label="E (SL)", lw=4)
    png(p, "$root_DF/$folder_DF/comp_E_$folder_SL.png")
    png(p, "$root_SL/$folder_SL/comp_E_$folder_DF.png")
    # display(p)

    # Density ρ

    ρDF = readdlm("$root_DF/$folder_DF/rho.dat")
    SLρData = readdlm("$root_SL/$folder_SL/data_output/rho$(@sprintf("%06d",Nt)).dat")
    ρSL = SLρData[2:end,2] 

    println("min / max ρ DF : $(minimum(ρDF)) / $(maximum(ρDF))")
    println("min / max ρ SL : $(minimum(ρSL)) / $(maximum(ρSL))")

    xDF = LinRange(-1.0, 1.0, length(ρDF))
    xSL = SLρData[2:end,1]
    p = plot(xDF, ρDF, ylims=ρlims, label="ρ (DF)", lw=4, legend=:bottomright, legendfontsize=14)
    plot!(p, xSL, ρSL, label="ρ (SL)", lw=4)
    png(p, "$root_DF/$folder_DF/comp_rho_$folder_SL.png")
    png(p, "$root_SL/$folder_SL/comp_rho_$folder_DF.png")
    # display(p)
end
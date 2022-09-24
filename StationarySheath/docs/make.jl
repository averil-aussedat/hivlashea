using Documenter
using StationarySheath

ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(StationarySheath, :DocTestSetup, :(using StationarySheath); recursive=true)

makedocs(;
    sitename = "Stationary Sheath",
    # repo="https://github.com/averil-prost/hivlashea/blob/{commit}{path}#{line}",
    # format = Documenter.HTML(;
    #     prettyurls=get(ENV, "CI", "false") == "true",
    #     canonical="https://averil-prost.github.io/StationarySheath.jl",
    #     edit_link="main",
    #     assets=String[]),
    # modules = [StationarySheath],
    # pages = ["index.md"]
)

# deploydocs(;
#     repo="github.com/averil-prost/StationarySheath.jl",
#     devbranch="main",
# )

mkdir run
cp -r jlfiles run/
cp -r jlimages run/
cp initial_data.jl main_julia.jl run/
rm jlfiles/*.jld
rm jlimages/*.png jlimages/*.mp4

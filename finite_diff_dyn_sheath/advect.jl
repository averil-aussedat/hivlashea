"""
    Advection tests
"""

using Plots
# using VideoIO
# using PerceptualColourMaps
# using Images

# params
mu  = 0.1  
N   = 50 
M   = 60
vx  = 1.0
vy  = 0.0
T   = 1.0

# init
xx = collect(LinRange(-1.0,1.0,N+2))
yy = collect(LinRange(-10.0,10.0,M+2))

ff = sqrt(mu).*exp.(-0.5 * mu .* (xx.^2 .+ yy'.^2))./sqrt(2*pi)
ff[begin,:] .*= 0.0; ff[end,:] .*= 0.0
ff[:,begin] .*= 0.0; ff[:,end] .*= 0.0

p = contourf(1.0*ff)
png("jlimages/advect0.png")

# advect : \partial_t f + vx \partial_x f + vy \partial_y f = 0
dx = xx[2]-xx[1]; dy = yy[2] - yy[1]
dt = 0.8 * min(dx / max(abs(vx),1e-6), dy / max(abs(vy),1e-6))
Nt = ceil(Int,T / dt)

println("T=$T, dt=$dt avec Nt=$Nt itérations")

for n in 1:Nt
    ff[begin+1:end-1,begin+1:end-1] .-= (
        (dt/dx) .* (max(vx,0.0) .* (ff[begin+1:end-1,begin+1:end-1] .- ff[begin:end-2,  begin+1:end-1]) 
                .+  min(vx,0.0) .* (ff[begin+2:end,  begin+1:end-1] .- ff[begin+1:end-1,begin+1:end-1]))
      + (dt/dy) .* (max(vy,0.0) .* (ff[begin+1:end-1,begin+1:end-1] .- ff[begin+1:end-1,begin:end-2]) 
                .+  min(vy,0.0) .* (ff[begin+1:end-1,begin+2:end]   .- ff[begin+1:end-1,begin+1:end-1]))
    )

    contourf!(p,1.0*ff)
    png("jlimages/advect$n.png")
end # for n

# then go to shell (type ;) and call ffmpeg -framerate 24 -i jlimages/advect%d.png jlimages/advect.mp4

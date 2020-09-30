using NBInclude
using Plots;plotly()
using LaTeXStrings
using LsqFit
using Printf
using HDF5, JLD

@nbinclude("fftFilter1d.ipynb");
@nbinclude("interpolatedFBM.ipynb")
@nbinclude("dfa.ipynb");

function series1d(hrange=0.2:0.1:0.8, prange=5:11)
    data = Dict()
    hursts = collect(hrange)
    powers = collect(prange)
    for hurst in hursts
        scales = Dict()
        # for each hurst <- scales 
        for scale in powers
            #for each hurst and scale <- samples
            samples = Dict()
            for sample=1:5
                # for each sample a zh(x) 
#                 hm, z = fbmfft(scale, hurst)
                hm, z = fbmmid(hurst, scale)
                samples["z$sample"]=z
            end
            L= convert(Int64,2.0^scale)
            scales["L$L"]=samples
        end
        println(hurst)
        data["h$hurst"]=scales
    end
    return data
end

function series1d(filename::String)
    # saved using save(filename, "data", dataFrame)
    data = load(filename)["data"]
end
data = series1d("1dseries.jld")
# data = series1d(0.2:0.1:0.8, 7:12)

y = data["h0.6"]["L1024"]["z5"]
pyplot()
plot(y, yaxis=[-50,50], grid=false)
savefig("t.pdf")

# save("1dseries.jld", "data", data)

# p1 = plot(xaxis=("L", :log10), yaxis=("Δz",:log10))
# for h in hursts
#     Δz = []
#     L = []
#     for scale in powers
#         Ls = convert(Int64,2.0^scale)
#         dz = 0
#         for s=1:5
#             dz += std(data["h$h"]["L$Ls"]["z$s"])
#         end
#         dz /= 5
#         push!(Δz, dz)
#         push!(L, Ls)
#     end
#     plot!(L, Δz, markershape=:circle, label=@sprintf "%0.2f" h)
#     plot!(L, (0.25/h).*L.^h, color=:black, label="")
# end
# plot(p1)


function perimenter(z, step=1)
    z2 =  z[1+step:step:end]
    z1 = z[1:step:end-step]
    δx = step.*ones(length(z1))
    Δz = z2 .- z1
    h = sqrt.(δx.^2 .+ Δz.^2)
    Lp = sum(h)/sum(δx)
end

p1=plot()
for p = 10:11
    L = convert(Int64, 2^p)
    Lp = []
    hurst = []
    for h=0.2:0.1:0.8
        per = 0.0
        for s=1:5
            z = data["h$h"]["L$L"]["z$s"]
            z = z./std(z)
            per += perimenter(z)
        end
        push!(Lp, per/5)
        push!(hurst, h)
    end
    plot!(hurst, (Lp.^2), markershape=:circle, label="$L")
end

plot(p1)



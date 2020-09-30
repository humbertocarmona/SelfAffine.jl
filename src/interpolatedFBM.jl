using Random
using Distributions


function midpointRecursion(z::Array{Float64,1}, i0::Int64, i1::Int64, level::Int64, maxlevel::Int64, Δn)
    i = (i0+i1)÷2
    z[i] = 0.5*(z[i0]+z[i1]) + rand(Normal(0.0, Δn[level]))
    if level < maxlevel
        midpointRecursion(z,i0,i, level+1, maxlevel,Δn)
        midpointRecursion(z,i,i1, level+1, maxlevel,Δn)
    end
end

function interpolatedFBM(H::Float64, maxlevel::Int64, σ::Float64=1.0)
    """
        generate a fractional brownian motion, L = 2^maxlevel
        using mid-point interpolation
    """
    L = 2^maxlevel
    f = σ*(L^H)*(1-2^(2*H-2))^0.5    # var(z(L)-z(0)) ≡ σ^2 L^2H
    Δn = [f*2^(-n*H) for n=1:maxlevel]
    z = zeros(L)
    z[L] = rand(Normal(0.0,σ*(L^H)))
    midpointRecursion(z,1,L,1,maxlevel, Δn)
    return z
end

include("dfa.jl")
function fbmmid(h, maxlevel)
    hurst = -1
    n = 1
    z = []
    while !isapprox(hurst, h; atol=0.01) && n<100
        z = interpolatedFBM(h, maxlevel)
        hurst, x, y = dfa(z; scal=4:0.25:maxlevel-2, offset=50)
        n+=1
    end
    return hurst, z
end

# using Plots; plotly()
# @nbinclude("dfa.ipynb")
# function interpolateFBMTest_(hurst, nsamples=20)
#     """
#     for a given hurst, plot chech if var(z) ≈ L^2h with L the size of the series
#     """
#     s = []
#     L = []
#     h = []
#     for n=8:14
#         sav = 0.0
#         hav = 0.0
#         for i=1:nsamples
#             z = interpolateFBM(hurst, n)
#             x, y = dfa(z; scal=3:0.25:n-2, offset=50)
#             p = curve_fit((x,p)-> p[1] .+ p[2].*x, x,y,[0.0,0.0])
#             sav += std(z)
#             hav += p.param[2] -1
#         end
#         push!(L, 2^n)
#         push!(h, hav/nsamples)
#         push!(s, sav/nsamples)
#     end
#     p1 = plot(L,s, markershape=:circle, axis=:log10)
#     plot!(L, 0.3.*L.^hurst)
# end





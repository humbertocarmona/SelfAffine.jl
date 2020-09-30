using FFTW
using Distributions
using Random

function fftFilter2D(maxlevel, h)
    #Random.seed!(2342)
    L = 2^maxlevel
    L2 = L ÷ 2
    β = (h+1)/2
    f = zeros(Complex{Float64}, L,L) 
    
    
    # upper left  and lower right 
    rad = [(i^2 + j^2)^(-β) for i=1:L2, j=1:L2]
    rad = rad .* rand(Normal(), L2, L2)
    ϕ = 2π .* rand(L2, L2)
    exp_iϕ = cos.(ϕ) .+ sin.(ϕ) .* im 
    Z = rad.*exp_iϕ
    Z[L2,L2] = real(Z[L2, L2])
    Zt = conj(Z[L2:-1:1, L2:-1:1])
    
    f[2:L2+1,2:L2+1] = Z     # line 1 is dealt with separetely
    f[(L2+1):L,(L2+1):L] = Zt  
    
    
    # lower left and  upper right
    rad = [(i^2 + j^2)^(-β) for i=1:(L2-1), j=1:(L2-1)]
    rad = rad .* rand(Normal(), (L2-1), (L2-1))
    ϕ = 2π .* rand((L2-1), (L2-1))
    exp_iϕ = cos.(ϕ) .+ sin.(ϕ) .* im 
    W = rad.*exp_iϕ
    Wt = conj(W[(L2-1):-1:1, (L2-1):-1:1])

    f[L:-1:(L2+2),2:L2] = W
    f[L2:-1:2,(L2+2):L] = Wt 
    
    
    # first line
    rad = [j^(-2*β) for j=1:L2]
    rad = rad .* rand(Normal(), L2)
    ϕ = (2*π) .* rand(L2)
    exp_iϕ = cos.(ϕ) .+ sin.(ϕ) .* im 
    line = rad.*exp_iϕ
    line[L2] = real(line[L2])
    linet = conj(line[L2:-1:1])
    
    f[1,2:(L2+1)] = line
    f[1,(L2+1):L] = linet

    # first column
    rad = [i^(-2*β) for i=1:L2]
    rad = rad .* rand(Normal(), L2)
    ϕ = (2*π) .* rand(L2)
    exp_iϕ = cos.(ϕ) .+ sin.(ϕ) .* im 
    column = rad.*exp_iϕ
    column[L2] = real(column[L2])
    columnt = conj(column[L2:-1:1])
    
    f[2:(L2+1),1] = column
    f[(L2+1):L,1] = columnt

    f[1,1] = 1

    g = ifft(f)
    # println(sum(abs.(imag.(g))))
    # h = fftshift(f)
    # display( heatmap(log.(abs.(h.*h)).+ 1.0,size=(600,500), color=ColorGradient([:black,:yellow,:red])))
    z = real(g)
    z = (z .- mean(z)).*L^2
    return z
end

include("dfa.jl")

function fbmSurf(maxlevel, h; atol=0.01)
    """
    - uses fftFilter2D to generate surfaces
    - checks Hurst using DFA
    - returns surface heights z with Hurst within atol
    """
    hurst = 50.0
    n=0
    z = []
    while !isapprox(hurst, h; atol=atol) && n<1e3
        z = fftFilter2D(maxlevel, h)
        hurst = dfa2d(z; scales=3:0.5:maxlevel-2)
        n+=1
    end
    println(n)
    return z
end

function make2DSurfs(;maxlevel=9, hursts=0.3:0.1:0.8, samples=1:1, fbm=nothing, atol=0.01)
    """
    returns a dictionary with the set of selfAffine surfaces
    """
    
    fbm = fbm == nothing ? Dict() : fbm
    fbm["hursts"]=collect(hursts)
    fbm["samples"]=collect(samples)
    fbm["maxlevel"] = maxlevel
    for h=hursts
        sample = Dict()
        for s=samples
            z = fbmSurf(maxlevel,h, atol=atol)
            sample["$s"]=z
        end
        fbm["$h"]=sample
        println("finished $h")
    end
    return fbm
end

function make2DSurfs(filename::String)
    fbm = load(filename)["fbm"]
end



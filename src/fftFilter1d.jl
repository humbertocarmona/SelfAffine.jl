using FFTW
using Random
using Distributions
using LsqFit

function fourierFiltering(pw, h)
    N= 2^pw
    N2 = N ÷ 2
    β = 2*h+1
    f = zeros(Complex{Float64}, N2)

    rad = [(i^(-β/2)) for i in 1:N2]
    rad = rad .* rand(Normal(0.0,1.0), N2)
    ϕ = 2π .* rand(N2)
    exp_iϕ = cos.(ϕ) .+ sin.(ϕ) .* im
    f[1:N2] = rad.*exp_iϕ
#     f[N2] = real(f[N2])
    g = irfft(f, 2*N2-1)  # irfft is nomalized by 1/N2
    zh = (g .- mean(g)).*N2^(1+h)
end

function fbmfft(pw, h)
    hurst=-1
    n=1
    z = fourierFiltering(pw,h)
    while !isapprox(hurst,h; atol=0.01) && n<10000
        z = fourierFiltering(pw,h)
        x, y = dfa(z; scal=3:0.25:(pw-2), offset=50)
        p = curve_fit((x,p)-> p[1] .+ p[2].*x, x,y,[0.0,0.0])
        
        hurst = p.param[2]-1
        n+=1
    end
    return hurst, z
end

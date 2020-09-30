using LsqFit

include("utils.jl")
function dfa(z; scales=3:0.5:8, offset=-1)
    # cummulative sum of (z - <z>)
    zcum = cumsum(z .- mean(z))
    
    # define the rulers for each x-window
    scales = collect(scales)
    rulers = convert.(Int64,round.(2 .^ scales))

    
    # each rulers will have rms averaged for each window in the data
    fluct = zeros(length(rulers))
    
    for (k, ℓ) in enumerate(rulers)        
        # x data is the same for all windows in this scale
        xdata = collect(1:ℓ)
        
        #each element contail the data in one window
        yvec = strided(zcum,ℓ,offset)

        #for each window one rms
        rms = zeros(length(yvec))

        for (i,ydata) in enumerate(yvec)
            fit = curve_fit((x, p) -> p[1] .+ p[2].*x, xdata, ydata, [1.0, 1.0])
            σ = stderror(fit)
            rms[i] = sqrt(mean(fit.resid.^2))
            # yfit = model(xdata, fit.param)
            # rms2 = sqrt(mean((ydata .- yfit).^2))
        end
        fluct[k] = mean(rms)
    end
    
    # now fit (fluct vs. scale) 
    x = log10.(rulers)
    y = log10.(fluct)
    fit = curve_fit((x, p) -> p[1] .+ p[2].*x ,x, y, [.0, .0])
    h = fit.param[2]-1.0
    return  h,x,y
end

function dfa2d(z; scales=3:0.5:8, offset=-1)
    L = size(z)[1]
    hurst = 0
    n = 0
    for i=10:10:L-10
        h, x, y = dfa(z[:,i]; scales=scales, offset=offset)
        hurst += h
        n+=1
    end
    for i=10:10:L-10
        h, x, y = dfa(z[i,:]; scales=scales, offset=offset)
        hurst += h
        n+=1
    end
    hurst = hurst/n
end



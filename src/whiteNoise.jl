using Random
using Distributions
function whiteNoise(pw)
    """
    inegrated white noise: H=0.5
    """
    N = 2^pw
    whitenoise = rand(Normal(0.0,1.0), N)
    z = cumsum(whitenoise)
end

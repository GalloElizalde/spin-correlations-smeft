using EFTfitter
using Distributions
using IntervalSets
using Plots

# ----------------------------------------------------
# Parameters
# ----------------------------------------------------
parameters = BAT.distprod(
    Qd8 = Uniform(-2.5, 1.3),
    Qj31 = Uniform(-0.5, 0.5)
)

# ----------------------------------------------------
# Observables
# ----------------------------------------------------
# Ckk (Qd8, Qj31)
function ckk(params)
    coeffs = [ -0.0174332, -0.0152784, -0.00952972, -0.276063, 6.01721e-05, -0.198477 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Crr (Qd8, Qj31)
function crr(params)
    coeffs = [ -0.0110969, 0.000677951, -0.00340285, 0.0621208, 0.00010976, -0.193549 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Cnn (Qd8, Qj31)
function cnn(params)
    coeffs = [ -0.000202056, -0.00282145, -0.0110085, -0.0954381, -0.000729606, 0.219988 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end



# ----------------------------------------------------
# myfunc
# ----------------------------------------------------
function myfunc(params, c)
    return c[1]*params.Qd8 +
           c[2]*params.Qj31 +
           c[3]*params.Qd8^2 +
           c[4]*params.Qj31^2 +
           c[5]*params.Qd8*params.Qj31 +
           c[6]
end

# ----------------------------------------------------
# Measurements for CMS 
# ----------------------------------------------------
measurements = (
    Ckk = Measurement(ckk, 0.040, uncertainties = (stat = 0.050, ), active=true),
    Crr = Measurement(crr, -0.202, uncertainties = (stat = 0.044, ), active=true),
    Cnn = Measurement(cnn, 0.175, uncertainties = (stat = 0.028, ), active=true),
)

# ----------------------------------------------------
# Correlations
# ----------------------------------------------------
correlations = (
    stat = NoCorrelation(active=true),
)

# ----------------------------------------------------
# Nuisance Correlations
# ----------------------------------------------------
nuisance_correlations = (

    ρ_Ckk_Crr = NuisanceCorrelation(:syst, :Ckk, :Crr, -1..1),
    ρ_Ckk_Cnn = NuisanceCorrelation(:syst, :Ckk, :Cnn, -1..1),
    ρ_Crr_Cnn = NuisanceCorrelation(:syst, :Crr, :Cnn, -1..1),
)

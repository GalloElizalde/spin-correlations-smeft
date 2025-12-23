using EFTfitter
using Distributions
using IntervalSets
using Plots

# ----------------------------------------------------
# Parameters
# ----------------------------------------------------
parameters = BAT.distprod(
    tu8 = Uniform(-1.5, 0.8),
    td8 = Uniform(-2.0, 1.2)
)

# number_events (tu8, td8)
function number_events(params)
    coeffs = [ 79091.5, 55785.4, 2511.25, 1758.67, 0.000247649, 13329000.0 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckk (tu8, td8)
function ckk(params)
    coeffs = [ -0.00599847, -0.00465597, -0.000167123, -0.00013699, 5.23875e-05, 0.552169 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Crr (tu8, td8)
function crr(params)
    coeffs = [ -0.00570269, -0.00407642, -0.000129017, -0.000100318, 4.77787e-05, 0.305683 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Cnn (tu8, td8)
function cnn(params)
    coeffs = [ -0.00292373, -0.0022066, -8.7024e-05, -6.72912e-05, 2.51518e-05, 0.51283 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckr (tu8, td8)
function ckr(params)
    coeffs = [ -0.00165576, -0.00119106, -3.82288e-05, -3.04598e-05, 1.39153e-05, -0.127937 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# D3 (tu8, td8)
function d3(params)
    coeffs = [ 0.00292581, 0.00217526, 6.97053e-05, 5.66722e-05, -2.50048e-05, -0.115007 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# A_plus (tu8, td8)
function a_plus(params)
    coeffs = [ 0.00229752, 0.0017146, 5.81469e-05, 4.65953e-05, -1.96676e-05, -0.168439 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# ----------------------------------------------------
# myfunc
# ----------------------------------------------------
function myfunc(params, c)
    return c[1]*params.tu8 +
           c[2]*params.td8 +
           c[3]*params.tu8^2 +
           c[4]*params.td8^2 +
           c[5]*params.tu8*params.td8 +
           c[6]
end

# ----------------------------------------------------
# Measurements
# ----------------------------------------------------
measurements = (
    N_events = Measurement(number_events, 13329000.0, uncertainties = (stat = 57381.8, ), active=true),
    Ckk = Measurement(ckk, 0.552168, uncertainties = (stat = 0.0125974, ), active=true),
    Crr = Measurement(crr, 0.305682, uncertainties = (stat = 0.0128093, ), active=true),
    Cnn = Measurement(cnn, 0.51283, uncertainties = (stat = 0.0127149, ), active=true),
    Ckr = Measurement(ckr, -0.127937, uncertainties = (stat = 0.0128734, ), active=true),
    D3 = Measurement(d3, -0.115007, uncertainties = (stat = 0.00742473, ), active=true),
    A_plus = Measurement(a_plus, -0.168438, uncertainties = (stat = 0.00352759, ), active=true),
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
    ρ_Ckk_N = NuisanceCorrelation(:syst, :N_events, :Ckk, -1..1),
    ρ_Crr_N = NuisanceCorrelation(:syst, :N_events, :Crr, -1..1),
    ρ_Cnn_N = NuisanceCorrelation(:syst, :N_events, :Cnn, -1..1),
    ρ_Ckr_N = NuisanceCorrelation(:syst, :N_events, :Ckr, -1..1),
    ρ_D3_N = NuisanceCorrelation(:syst, :N_events, :D3, -1..1),
    ρ_A_plus_N = NuisanceCorrelation(:syst, :N_events, :A_plus, -1..1),
    ρ_N_events_Ckk = NuisanceCorrelation(:syst, :N_events, :Ckk, -1..1),
    ρ_N_events_Crr = NuisanceCorrelation(:syst, :N_events, :Crr, -1..1),
    ρ_N_events_Cnn = NuisanceCorrelation(:syst, :N_events, :Cnn, -1..1),
    ρ_N_events_Ckr = NuisanceCorrelation(:syst, :N_events, :Ckr, -1..1),
    ρ_N_events_D3 = NuisanceCorrelation(:syst, :N_events, :D3, -1..1),
    ρ_N_events_A_plus = NuisanceCorrelation(:syst, :N_events, :A_plus, -1..1),
    ρ_Ckk_Crr = NuisanceCorrelation(:syst, :Ckk, :Crr, -1..1),
    ρ_Ckk_Cnn = NuisanceCorrelation(:syst, :Ckk, :Cnn, -1..1),
    ρ_Ckk_Ckr = NuisanceCorrelation(:syst, :Ckk, :Ckr, -1..1),
    ρ_Ckk_D3 = NuisanceCorrelation(:syst, :Ckk, :D3, -1..1),
    ρ_Ckk_A_plus = NuisanceCorrelation(:syst, :Ckk, :A_plus, -1..1),
    ρ_Crr_Cnn = NuisanceCorrelation(:syst, :Crr, :Cnn, -1..1),
    ρ_Crr_Ckr = NuisanceCorrelation(:syst, :Crr, :Ckr, -1..1),
    ρ_Crr_D3 = NuisanceCorrelation(:syst, :Crr, :D3, -1..1),
    ρ_Crr_A_plus = NuisanceCorrelation(:syst, :Crr, :A_plus, -1..1),
    ρ_Cnn_Ckr = NuisanceCorrelation(:syst, :Cnn, :Ckr, -1..1),
    ρ_Cnn_D3 = NuisanceCorrelation(:syst, :Cnn, :D3, -1..1),
    ρ_Cnn_A_plus = NuisanceCorrelation(:syst, :Cnn, :A_plus, -1..1),
    ρ_Ckr_D3 = NuisanceCorrelation(:syst, :Ckr, :D3, -1..1),
    ρ_Ckr_A_plus = NuisanceCorrelation(:syst, :Ckr, :A_plus, -1..1),
    ρ_D3_A_plus = NuisanceCorrelation(:syst, :D3, :A_plus, -1..1),
)

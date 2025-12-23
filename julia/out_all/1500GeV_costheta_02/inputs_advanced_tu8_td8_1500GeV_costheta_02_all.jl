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
    coeffs = [ 2580.31, 3687.41, 4126.45, 4998.08, 9.71708e-06, 9137.79 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckk (tu8, td8)
function ckk(params)
    coeffs = [ -0.195926, 0.186944, -0.598713, 0.273169, -0.0254152, -1.29721 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Crr (tu8, td8)
function crr(params)
    coeffs = [ 0.356394, -0.325785, 0.612045, -0.424645, 0.012042, -0.761177 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Cnn (tu8, td8)
function cnn(params)
    coeffs = [ 0.454875, 0.098624, 0.129201, -0.168474, -0.106149, 1.01689 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckr (tu8, td8)
function ckr(params)
    coeffs = [ -0.107064, -0.317125, 0.117538, -0.525817, 0.0131493, -0.00767794 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# D3 (tu8, td8)
function d3(params)
    coeffs = [ 0.0981355, 0.0791552, 0.0386231, -0.00566616, -0.0309254, 1.02509 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# A_plus (tu8, td8)
function a_plus(params)
    coeffs = [ -0.0315079, 0.0272615, -0.00261774, 0.0297422, 0.00262583, 0.404163 ]
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
    N_events = Measurement(number_events, 9137.79, uncertainties = (stat = 1502.25, ), active=true),
    Ckk = Measurement(ckk, -1.33473, uncertainties = (stat = 0.356403, ), active=true),
    Crr = Measurement(crr, -0.678161, uncertainties = (stat = 0.394086, ), active=true),
    Cnn = Measurement(cnn, 1.13978, uncertainties = (stat = 0.569432, ), active=true),
    Ckr = Measurement(ckr, 0.139701, uncertainties = (stat = 0.516097, ), active=true),
    D3 = Measurement(d3, 1.05089, uncertainties = (stat = 0.222918, ), active=true),
    A_plus = Measurement(a_plus, 0.39523, uncertainties = (stat = 0.104329, ), active=true),
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

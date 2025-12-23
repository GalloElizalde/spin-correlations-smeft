using EFTfitter
using Distributions
using IntervalSets
using Plots

# ----------------------------------------------------
# Parameters
# ----------------------------------------------------
parameters = BAT.distprod(
    tu8 = Uniform(-1.5, 0.8),
    tj8 = Uniform(-1.4, 0.5)
)

# number_events (tu8, tj8)
function number_events(params)
    coeffs = [ 79091.5, 131600.0, 2511.25, 4180.41, -0.000667937, 13329000.0 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckk (tu8, tj8)
function ckk(params)
    coeffs = [ -0.0059984, -0.0093924, -0.000167303, -0.000227619, 0.000113938, 0.552169 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Crr (tu8, tj8)
function crr(params)
    coeffs = [ -0.00570265, -0.00950053, -0.000129242, -0.00018819, 0.000111846, 0.305683 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Cnn (tu8, tj8)
function cnn(params)
    coeffs = [ -0.00292371, -0.00474829, -8.7141e-05, -0.000126389, 5.6477e-05, 0.51283 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckr (tu8, tj8)
function ckr(params)
    coeffs = [ -0.00194087, -0.00637735, -0.000496437, -0.00377012, -0.000568928, -0.127682 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# D3 (tu8, tj8)
function d3(params)
    coeffs = [ 0.00292578, 0.00471488, 6.98013e-05, 9.64732e-05, -5.64357e-05, -0.115007 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# A_plus (tu8, tj8)
function a_plus(params)
    coeffs = [ 0.0022975, 0.00370962, 5.82264e-05, 8.16439e-05, -4.43326e-05, -0.168439 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# ----------------------------------------------------
# myfunc
# ----------------------------------------------------
function myfunc(params, c)
    return c[1]*params.tu8 +
           c[2]*params.tj8 +
           c[3]*params.tu8^2 +
           c[4]*params.tj8^2 +
           c[5]*params.tu8*params.tj8 +
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

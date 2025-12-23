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
    coeffs = [ 165951.0, 265052.0, 104352.0, 178085.0, 3.22184e-05, 4768930.0 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckk (tu8, tj8)
function ckk(params)
    coeffs = [ -0.0270817, -0.0478224, -0.0158507, -0.0326888, 0.00157363, -0.197748 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Crr (tu8, tj8)
function crr(params)
    coeffs = [ -0.0101307, -0.00769184, -0.000496585, 0.0186919, 0.00106517, -0.19365 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Cnn (tu8, tj8)
function cnn(params)
    coeffs = [ 0.0098545, 0.0152633, -0.000660699, -0.00440151, -0.000937397, 0.219883 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckr (tu8, tj8)
function ckr(params)
    coeffs = [ -0.00365508, 0.0053654, 0.00165544, 0.00563994, 0.000150595, -0.0952974 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# D3 (tu8, tj8)
function d3(params)
    coeffs = [ 0.0159648, 0.025865, 0.00560675, 0.00539107, -0.000766516, 0.203561 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# A_plus (tu8, tj8)
function a_plus(params)
    coeffs = [ 0.00733786, 0.0111171, 0.00325086, 0.00295911, -0.000474085, 0.0768294 ]
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
    N_events = Measurement(number_events, 4768930.0, uncertainties = (stat = 34323.1, ), active=true),
    Ckk = Measurement(ckk, -0.19816, uncertainties = (stat = 0.021469, ), active=true),
    Crr = Measurement(crr, -0.193724, uncertainties = (stat = 0.0217452, ), active=true),
    Cnn = Measurement(cnn, 0.220023, uncertainties = (stat = 0.0212365, ), active=true),
    Ckr = Measurement(ckr, -0.0952748, uncertainties = (stat = 0.021558, ), active=true),
    D3 = Measurement(d3, 0.203969, uncertainties = (stat = 0.0123823, ), active=true),
    A_plus = Measurement(a_plus, 0.0769464, uncertainties = (stat = 0.00599999, ), active=true),
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

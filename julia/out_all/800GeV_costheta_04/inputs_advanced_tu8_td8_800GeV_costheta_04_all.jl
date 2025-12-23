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
    coeffs = [ 53534.5, 32064.1, 29433.0, 19654.3, -1.67843e-05, 679056.0 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckk (tu8, td8)
function ckk(params)
    coeffs = [ -0.0598342, -0.0211848, -0.063079, -0.00618541, 0.00278053, -0.504886 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Crr (tu8, td8)
function crr(params)
    coeffs = [ 0.00965843, -0.0268254, 0.0303945, -0.00789038, 0.00193933, -0.709492 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Cnn (tu8, td8)
function cnn(params)
    coeffs = [ 0.0295936, 0.00777624, 1.21192e-05, -0.021861, -0.00227514, 0.671103 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# Ckr (tu8, td8)
function ckr(params)
    coeffs = [ -0.0208252, 0.00756577, -0.00788359, -0.0181536, -0.000283896, -0.0628136 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# D3 (tu8, td8)
function d3(params)
    coeffs = [ 0.0269004, 0.0186206, 0.0110927, -0.00244997, -0.00204435, 0.628353 ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end

# A_plus (tu8, td8)
function a_plus(params)
    coeffs = [ 0.00980976, 0.00941792, 0.00639933, 0.00273396, -0.000984453, 0.238466 ]
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
    N_events = Measurement(number_events, 679056.0, uncertainties = (stat = 12951.5, ), active=true),
    Ckk = Measurement(ckk, -0.505595, uncertainties = (stat = 0.0563798, ), active=true),
    Crr = Measurement(crr, -0.709709, uncertainties = (stat = 0.0555165, ), active=true),
    Cnn = Measurement(cnn, 0.671824, uncertainties = (stat = 0.054248, ), active=true),
    Ckr = Measurement(ckr, -0.0625654, uncertainties = (stat = 0.0574207, ), active=true),
    D3 = Measurement(d3, 0.629043, uncertainties = (stat = 0.0301844, ), active=true),
    A_plus = Measurement(a_plus, 0.238624, uncertainties = (stat = 0.0155362, ), active=true),
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

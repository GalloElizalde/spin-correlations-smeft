import pandas as pd
from pathlib import Path

# Load csv with data
df = pd.read_csv("2d_fits_800GeV_cms_normalized.txt")

# Select only Spin correlations
target_observables = ["Ckk", "Crr", "Cnn"]
df = df[df["observable"].isin(target_observables)]

# Select operators
groups = df.groupby(["operator_x", "operator_y"])

# Create output file
output_dir = Path("cms_inputs")
output_dir.mkdir(exist_ok=True)

# Template
header_template = """using EFTfitter
using Distributions
using IntervalSets
using Plots

# ----------------------------------------------------
# Parameters
# ----------------------------------------------------
parameters = BAT.distprod(
    {opx} = Uniform({minx}, {maxx}),
    {opy} = Uniform({miny}, {maxy})
)
"""

observable_template = """
# {observable} ({opx}, {opy})
function {func}(params)
    coeffs = [ {coeffs} ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end
"""

footer_template = """
# ----------------------------------------------------
# myfunc
# ----------------------------------------------------
function myfunc(params, c)
    return c[1]*params.{opx} +
           c[2]*params.{opy} +
           c[3]*params.{opx}^2 +
           c[4]*params.{opy}^2 +
           c[5]*params.{opx}*params.{opy} +
           c[6]
end

# ----------------------------------------------------
# Measurements for CMS Region: m(ttbar) > 800 GeV  (sin corte angular)
# ----------------------------------------------------
measurements = (
    Crr = Measurement(crr, -0.13309, uncertainties = (stat = 0.0079779,), active = true),
    Cnn = Measurement(cnn,  0.16968, uncertainties = (stat = 0.0045698,), active = true),
    Ckk = Measurement(ckk, -0.041208, uncertainties = (stat = 0.012229,), active = true),
)

# ----------------------------------------------------
# Measurements for CMS Region: m(ttbar) > 800 GeV, |cos(theta)| < 0.4
# ----------------------------------------------------
#measurements = (
#    Crr = Measurement(crr, -0.69407, uncertainties = (stat = 0.0059392,), active = true),
#    Cnn = Measurement(cnn,  0.56297, uncertainties = (stat = 0.00415114,), active = true),
#    Ckk = Measurement(ckk, -0.50611, uncertainties = (stat = 0.0041601,), active = true),
#)


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
"""

# Loop over each operator pair
for (opx, opy), subdf in groups:
    minx = subdf["min_x"].iloc[0]
    maxx = subdf["max_x"].iloc[0]
    miny = subdf["min_y"].iloc[0]
    maxy = subdf["max_y"].iloc[0]

    script = header_template.format(opx=opx, opy=opy, minx=minx, maxx=maxx, miny=miny, maxy=maxy)

    for _, row in subdf.iterrows():
        coeffs = ", ".join([str(row["[1]"]), str(row["[2]"]), str(row["[3]"]),
                            str(row["[4]"]), str(row["[5]"]), str(row["[0]"])])
        func_code = observable_template.format(
            observable=row["observable"],
            opx=opx,
            opy=opy,
            func=row["observable"].lower(),
            coeffs=coeffs
        )
        script += func_code

    script += footer_template.format(opx=opx, opy=opy)

    # Save file
    #filename = f"cms_800GeV_costheta_04_{opx}_{opy}_inputs.jl"
    filename = f"cms_800GeV_{opx}_{opy}_inputs.jl"
    with open(output_dir / filename, "w") as f:
        f.write(script)

print("Done! files .jl are located in directory: 'EFTfitter_inputs'.")
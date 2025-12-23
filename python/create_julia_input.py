import pandas as pd
import os

# Load file
df = pd.read_csv("2d_fits_800GeV.txt")

# Show first rows (opcional)
print(df.head())

# Colums names should be like 'operator_x' y 'operator_y'
# Obtain unique pairs
unique_pairs = df[['operator_x', 'operator_y']].drop_duplicates()

# Show unique operator pairs
print(unique_pairs)




# List of input files and output files
files_and_dirs = {
    "2d_fits_400GeV.txt": "400GeV",
    "2d_fits_800GeV.txt": "800GeV",
    "2d_fits_800GeV_costheta_04.txt": "800GeV_costheta_04",
    "2d_fits_1500GeV_costheta_02.txt": "1500GeV_costheta_02",
}

# Define Julia Template
header_template = """using EFTfitter
using Distributions
using IntervalSets
using Plots

# ----------------------------------------------------
# Parameters
# ----------------------------------------------------
parameters = BAT.distprod(
    {op_x} = Uniform({min_x}, {max_x}),
    {op_y} = Uniform({min_y}, {max_y})
)
"""

func_template = """
# {name} ({op_x}, {op_y})
function {func_name}(params)
    coeffs = [ {c1}, {c2}, {c3}, {c4}, {c5}, {c0} ]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction
    return (prediction, uncertainty)
end
"""

myfunc_template = """
# ----------------------------------------------------
# myfunc
# ----------------------------------------------------
function myfunc(params, c)
    return c[1]*params.{op_x} +
           c[2]*params.{op_y} +
           c[3]*params.{op_x}^2 +
           c[4]*params.{op_y}^2 +
           c[5]*params.{op_x}*params.{op_y} +
           c[6]
end
"""

measurement_template = "    {label} = Measurement({func}, {val}, uncertainties = (stat = {err}, ), active={active}),\n"

# Loop over files
for filename, region in files_and_dirs.items():
    df = pd.read_csv(filename)
    grouped = df.groupby(["operator_x", "operator_y"])

    for mode in ["all", "noN", "onlyN"]:
        output_dir = f"./out_{mode}/{region}"
        os.makedirs(output_dir, exist_ok=True)

        for (op_x, op_y), group in grouped:
            row = group.iloc[0]
            min_x, max_x = row["min_x"], row["max_x"]
            min_y, max_y = row["min_y"], row["max_y"]

            base_name = f"inputs_advanced_{op_x}_{op_y}_{region}_{mode}.jl"
            path = os.path.join(output_dir, base_name)
            with open(path, "w") as f:
                f.write(header_template.format(op_x=op_x, op_y=op_y, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y))

                for _, row in group.iterrows():
                    observable = row["observable"]
                    func = observable.lower().replace("+", "_plus").replace("-", "_minus")
                    f.write(func_template.format(
                        name=observable,
                        func_name=func,
                        op_x=op_x,
                        op_y=op_y,
                        c1=row["[1]"],
                        c2=row["[2]"],
                        c3=row["[3]"],
                        c4=row["[4]"],
                        c5=row["[5]"],
                        c0=row["[0]"]
                    ))

                f.write(myfunc_template.format(op_x=op_x, op_y=op_y))

                # Measurements
                f.write("\n# ----------------------------------------------------\n# Measurements\n# ----------------------------------------------------\n")
                f.write("measurements = (\n")
                for _, row in group.iterrows():
                    observable = row["observable"]
                    func = observable.lower().replace("+", "_plus").replace("-", "_minus")
                    label = "N_events" if observable == "number_events" else observable

                    if mode == "all":
                        active = "true"
                    elif mode == "noN":
                        active = "false" if label == "N_events" else "true"
                    elif mode == "onlyN":
                        active = "true" if label == "N_events" else "false"

                    f.write(measurement_template.format(
                        label=label,
                        func=func,
                        val=row["sm_val"],
                        err=row["sm_val_error"],
                        active=active
                    ))
                f.write(")\n")

                # Correlations
                if mode in ["all", "noN"]:
                    f.write("""
# ----------------------------------------------------
# Correlations
# ----------------------------------------------------
correlations = (
    stat = NoCorrelation(active=true),
)
""")

                    f.write("""
# ----------------------------------------------------
# Nuisance Correlations
# ----------------------------------------------------
nuisance_correlations = (
""")
                    observables = [row["observable"] for _, row in group.iterrows()]
                    obs_labels = ["N_events" if o == "number_events" else o for o in observables]

                    for o in obs_labels:
                        if o != "N_events" and not (mode == "noN"):
                            f.write(f"    ρ_{o}_N = NuisanceCorrelation(:syst, :N_events, :{o}, -1..1),\n")

                    for i in range(len(obs_labels)):
                        for j in range(i + 1, len(obs_labels)):
                            oi, oj = obs_labels[i], obs_labels[j]
                            if "N_events" in [oi, oj] and mode == "noN":
                                continue
                            if "N_events" not in [oi, oj] or mode == "all":
                                f.write(f"    ρ_{oi}_{oj} = NuisanceCorrelation(:syst, :{oi}, :{oj}, -1..1),\n")

                    f.write(")\n")
                else:
                    f.write("\ncorrelations = (stat = NoCorrelation(active=true),)\n")
                    f.write("\nnuisance_correlations = NamedTuple()\n")


print("Files created for type: all, noN and onlyN in all regions.")
# =============================================================================
# EE 480 - Course Project 1: The Low Voltage Problem
# A Distribution Planning Study
#
# Installment 1 - Deliverable 3
#
# This script:
# 1. Loads the IEEE 14-bus case14.m file.
# 2. Runs an AC power flow on the original system.
# 3. Increases all loads by 40%.
# 4. Runs an AC power flow on the stressed system.
# 5. Identifies the bus with the lowest voltage.
# 6. Checks for buses below 0.95 p.u.
# 7. Checks generator operating limits.
#
# HOW TO RUN:
# 1. Install Julia and the Julia extension in VS Code.
# 2. Place case14.m in the same folder as this Julia file.
# 3. Run Cells 1 through 7 in order.
# =============================================================================


## Cell 1 - Setup: load the packages

using Pkg

if isfile(joinpath(@__DIR__, "Project.toml"))

    # Use the package versions defined in the project
    Pkg.activate(@__DIR__)
    Pkg.instantiate()

elseif Base.find_package("PowerModels") === nothing ||
       Base.find_package("Ipopt") === nothing

    # Install the required packages if they are not installed
    Pkg.add("PowerModels")
    Pkg.add("Ipopt")

end


## Cell 2 - Load the packages and choose the solver

using PowerModels, Ipopt, LinearAlgebra, Printf

PowerModels.silence()

# Ipopt solves the AC power flow equations.
# Hide the solver's banner and iteration log.

p1_solver = optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level" => 0,
    "sb" => "yes"
)

# PowerModels uses string keys ("1", "2", etc.).
# Sort the keys numerically.

p1_ids(d) = sort(collect(keys(d)), by = k -> parse(Int, k))

# Check whether the solver successfully completed.

p1_solved(r) = string(r["termination_status"]) in (
    "LOCALLY_SOLVED",
    "ALMOST_LOCALLY_SOLVED",
    "OPTIMAL"
)

# Bus type:
# 3 = Slack
# 2 = PV
# 1 = PQ

p1_type(t) = t == 3 ? "slack" : t == 2 ? "PV" : "PQ"

println(
    "PowerModels ", pkgversion(PowerModels),
    " and Ipopt ", pkgversion(Ipopt),
    " loaded."
)


## Cell 3 - Load the IEEE 14-bus case and examine the system

# Load the case14.m file from the same folder as this Julia script.

p1_file = joinpath(@__DIR__, "case14.m")

p1_data = PowerModels.parse_file(p1_file)

# Load the original file without validation for comparison.

p1_raw = PowerModels.parse_file(p1_file; validate = false)

# PowerModels works in per unit.
# Multiply P and Q by baseMVA to convert to MW and MVAr.

p1_MVA = p1_data["baseMVA"]

@printf(
    "IEEE 14-bus case (baseMVA = %.0f): %d buses, %d branches, %d gens, %d loads\n",
    p1_MVA,
    length(p1_data["bus"]),
    length(p1_data["branch"]),
    length(p1_data["gen"]),
    length(p1_data["load"])
)

@printf(
    "Total load: %.1f MW, %.1f MVAr\n\n",
    p1_MVA * sum(l["pd"] for l in values(p1_data["load"])),
    p1_MVA * sum(l["qd"] for l in values(p1_data["load"]))
)

# Display generator information.

@printf(
    "%3s %3s %-5s %6s %6s %6s %6s %6s\n",
    "Gen", "Bus", "type", "P set",
    "Pmax", "Qmin", "Qmax", "V set"
)

@printf(
    "%3s %3s %-5s %6s %6s %6s %6s %6s\n",
    "", "", "", "MW", "MW", "MVAr", "MVAr", "p.u."
)

for k in p1_ids(p1_data["gen"])

    g = p1_data["gen"][k]
    b = p1_data["bus"][string(g["gen_bus"])]

    # The slack generator's P is determined by the power flow.

    pset = b["bus_type"] == 3 ?
        "free" : @sprintf("%.1f", p1_MVA * g["pg"])

    @printf(
        "%3s %3d %-5s %6s %6.1f %6.1f %6.1f %6.3f\n",
        k,
        g["gen_bus"],
        p1_type(b["bus_type"]),
        pset,
        p1_MVA * g["pmax"],
        p1_MVA * g["qmin"],
        p1_MVA * g["qmax"],
        b["vm"]
    )

end

# Display notes about the original system data.

println("\nNotes on the data:")

for k in p1_ids(p1_data["bus"])

    t_file = p1_raw["bus"][k]["bus_type"]
    t_pm = p1_data["bus"][k]["bus_type"]

    if t_file != t_pm

        @printf(
            " * Bus %s is type %d (%s) in the file but has a generator,\n",
            k, t_file, p1_type(t_file)
        )

        @printf(
            " so PowerModels makes it type %d (%s).\n",
            t_pm, p1_type(t_pm)
        )

    end

end

# Display existing shunt capacitor information.

for k in p1_ids(p1_data["shunt"])

    s = p1_data["shunt"][k]

    @printf(
        " * Existing shunt capacitor at Bus %d: bs = %.2f p.u.",
        s["shunt_bus"],
        s["bs"]
    )

    @printf(
        " = %.1f MVAr at 1.0 p.u.\n",
        p1_MVA * s["bs"]
    )

end

# Display branch rating information.

@printf(
    " * rate_a missing on %d of %d branches: rateA = 0 in the file",
    count(!haskey(br, "rate_a") for br in values(p1_data["branch"])),
    length(p1_data["branch"])
)

println(" means \"no limit\"")

println(
    " and PowerModels drops it (thermal ratings are not needed in this project)."
)

if all(b["base_kv"] == 0 for b in values(p1_data["bus"]))

    println(
        " * base_kv = 0 for every bus in the file: voltages are in p.u. only."
    )

end


## Cell 4 - AC power flow on the original case

# Run the original IEEE 14-bus system.

p1_res0 = solve_ac_pf(p1_data, p1_solver)

println(
    "Original case: AC power flow status = ",
    p1_res0["termination_status"]
)

p1_solved(p1_res0) ||
    error("The power flow did not solve - check the data")

# Get the voltage magnitude at every bus.

p1_vm0 = Dict(
    parse(Int, k) => b["vm"]
    for (k, b) in p1_res0["solution"]["bus"]
)

# Find the bus with the lowest voltage.

p1_vmin0, p1_bmin0 = findmin(p1_vm0)

@printf(
    " Lowest voltage: Bus %d = %.4f p.u.\n",
    p1_bmin0,
    p1_vmin0
)

# Display slack generator real and reactive power.

@printf(
    " Slack generator (Bus 1): P = %.1f MW, Q = %.1f MVAr\n",
    p1_MVA * p1_res0["solution"]["gen"]["1"]["pg"],
    p1_MVA * p1_res0["solution"]["gen"]["1"]["qg"]
)


## Cell 5 - Create the stressed case: every load +40%

# Increase all real and reactive power loads by 40%.

p1_stress = 1.40

# Copy the original case without changing the original data.

p1_case = deepcopy(p1_data)

for (k, load) in p1_case["load"]

    # Increase real power demand by 40%.
    load["pd"] *= p1_stress

    # Increase reactive power demand by 40%.
    load["qd"] *= p1_stress

end

# Calculate original total demand.

p1_P0 = p1_MVA * sum(
    l["pd"] for l in values(p1_data["load"])
)

p1_Q0 = p1_MVA * sum(
    l["qd"] for l in values(p1_data["load"])
)

# Calculate stressed total demand.

p1_P1 = p1_MVA * sum(
    l["pd"] for l in values(p1_case["load"])
)

p1_Q1 = p1_MVA * sum(
    l["qd"] for l in values(p1_case["load"])
)

# Display load comparison.

@printf(
    "Stressed case: every load x %.2f",
    p1_stress
)

println(
    " (P and Q, so each load keeps its power factor)"
)

@printf(
    " Total load: %.1f MW -> %.1f MW, %.1f MVAr -> %.1f MVAr\n",
    p1_P0,
    p1_P1,
    p1_Q0,
    p1_Q1
)

println(
    " Generator P setpoints are unchanged: the slack (Bus 1) must pick up"
)

@printf(
    " the extra %.1f MW plus the extra losses.\n",
    p1_P1 - p1_P0
)


## Cell 6 - AC power flow on the stressed case

# Run the power flow with all loads increased by 40%.

p1_res = solve_ac_pf(p1_case, p1_solver)

println(
    "Stressed case: AC power flow status = ",
    p1_res["termination_status"]
)

p1_solved(p1_res) ||
    error("The power flow did not solve - check the data")

# Extract voltage magnitudes for all buses.

p1_vm = Dict(
    parse(Int, k) => b["vm"]
    for (k, b) in p1_res["solution"]["bus"]
)

# Find the lowest-voltage bus.

p1_vmin, p1_bmin = findmin(p1_vm)

# Compare original and stressed voltages.

@printf(
    "\n%3s %-5s %10s %10s %8s\n",
    "Bus",
    "type",
    "V original",
    "V stressed",
    "change"
)

@printf(
    "%3s %-5s %10s %10s %8s\n",
    "",
    "",
    "p.u.",
    "p.u.",
    "p.u."
)

for i in sort(collect(keys(p1_vm)))

    # Flag buses below 0.95 p.u.

    flag = p1_vm[i] < 0.95 ? " <-- below 0.95" : ""

    # Mark the bus with the lowest voltage.

    i == p1_bmin && (flag *= " <-- lowest")

    # Print the voltage comparison.

    @printf(
        "%3d %-5s %10.4f %10.4f %8.4f%s\n",
        i,
        p1_type(p1_case["bus"][string(i)]["bus_type"]),
        p1_vm0[i],
        p1_vm[i],
        p1_vm[i] - p1_vm0[i],
        flag
    )

end

# Find all buses below 0.95 p.u.

p1_lowbuses = sort(
    [i for (i, v) in p1_vm if v < 0.95]
)

# Display the lowest-voltage bus.

@printf(
    "\nLowest voltage in the stressed case: Bus %d = %.4f p.u.\n",
    p1_bmin,
    p1_vmin
)

# Display voltage violations.

println(
    "Buses below 0.95 p.u.: ",
    isempty(p1_lowbuses) ? "none" : join(p1_lowbuses, ", ")
)


## Cell 7 - Check generator operating limits

# IMPORTANT:
# solve_ac_pf holds every PV bus at its voltage setpoint.
# It does not enforce generator reactive power limits (Qmin/Qmax).
#
# The slack generator supplies the remaining P and Q.
#
# This cell checks whether generator outputs exceed their limits.

@printf(
    "%3s %3s %-5s %6s %6s %6s %6s %6s\n",
    "Gen",
    "Bus",
    "type",
    "P",
    "Pmax",
    "Q",
    "Qmin",
    "Qmax"
)

@printf(
    "%3s %3s %-5s %6s %6s %6s %6s %6s\n",
    "",
    "",
    "",
    "MW",
    "MW",
    "MVAr",
    "MVAr",
    "MVAr"
)

# Lists of generators outside reactive power limits.

p1_out_pv = Int[]
p1_out_slack = Int[]

for k in p1_ids(p1_case["gen"])

    g = p1_case["gen"][k]
    s = p1_res["solution"]["gen"][k]

    btype = p1_type(
        p1_case["bus"][string(g["gen_bus"])]["bus_type"]
    )

    # Generator real and reactive power.

    P = p1_MVA * s["pg"]
    Q = p1_MVA * s["qg"]

    # Generator operating limits.

    Pmax = p1_MVA * g["pmax"]
    Qmin = p1_MVA * g["qmin"]
    Qmax = p1_MVA * g["qmax"]

    flags = String[]

    # Check real power limit.

    P > Pmax + 0.05 &&
        push!(flags, "P above Pmax")

    # Check reactive power upper limit.

    Q > Qmax + 0.05 &&
        push!(flags, "Q above Qmax")

    # Check reactive power lower limit.

    Q < Qmin - 0.05 &&
        push!(flags, "Q below Qmin")

    # Identify generators outside reactive power limits.

    if Q > Qmax + 0.05 || Q < Qmin - 0.05

        push!(
            btype == "slack" ? p1_out_slack : p1_out_pv,
            g["gen_bus"]
        )

    end

    flag = isempty(flags) ? "" : " <-- " * join(flags, ", ")

    # Display generator results.

    @printf(
        "%3s %3d %-5s %6.1f %6.1f %6.1f %6.1f %6.1f%s\n",
        k,
        g["gen_bus"],
        btype,
        P,
        Pmax,
        Q,
        Qmin,
        Qmax,
        flag
    )

end

# Display reactive power limit violations.

println(
    "\nGenerators outside their reactive limits: PV buses ",
    isempty(p1_out_pv) ? "none" : join(p1_out_pv, ", "),
    "; slack bus ",
    isempty(p1_out_slack) ? "none" : join(p1_out_slack, ", ")
)

println(
    "solve_ac_pf did not stop them: it does not enforce Qmin/Qmax."
)


# =============================================================================
# END OF INSTALLMENT 1
#
# NEXT STEPS (LATER INSTALLMENTS)
#
# Installment 2:
# - Investigate generator reactive power limit violations.
# - Recalculate voltages while respecting reactive power limits.
# - Determine the reactive power deficit at the problem bus.
#
# Installment 3:
# - Add a new shunt capacitor at the problem bus.
# - Test different capacitor sizes.
# - Determine the minimum compensation needed to bring all buses
#   to at least 0.95 p.u.
# - Compare voltage profiles before and after compensation.
#
# =============================================================================
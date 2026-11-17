# =============================================================================
# EE 480 - PROJECT 1: THE LOW VOLTAGE PROBLEM
# INSTALLMENT 2: PROBLEM ANALYSIS AND REACTIVE POWER DEFICIT

# In Installment 2 you quantify
# the reactive power deficit at the problem bus; in Installment 3 you size a new
# shunt capacitor. See NEXT STEPS at the end of the file.
# =============================================================================
## Cell 1 - Setup: load the packages
using Pkg
if isfile(joinpath(@__DIR__, "Project.toml"))
# this file sits in a folder with its own pinned package versions: use them
Pkg.activate(@__DIR__)
Pkg.instantiate()
elseif Base.find_package("PowerModels") === nothing || Base.find_package("Ipopt") === nothing
# this file on its own: install the two packages once
Pkg.add("PowerModels")
Pkg.add("Ipopt")
end
## Cell 2 - Load the packages and choose the solver
using PowerModels, Ipopt, LinearAlgebra, Printf
PowerModels.silence() # hide PowerModels' warnings; Cell 3 prints what matters
# Ipopt solves the AC power flow equations. print_level = 0 and sb = "yes" hide its
# banner and iteration log so the classroom output stays readable.
p1_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")
# Three small helpers used below
# PowerModels keys every bus/gen/load by a STRING ("1", "2", ...): sort them as numbers
p1_ids(d) = sort(collect(keys(d)), by = k -> parse(Int, k))
# Ipopt reports success as LOCALLY_SOLVED (or ALMOST_LOCALLY_SOLVED / OPTIMAL)
p1_solved(r) = string(r["termination_status"]) in ("LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED", "OPTIMAL")
# bus type code -> name (3 = slack, 2 = PV, 1 = PQ)
p1_type(t) = t == 3 ? "slack" : t == 2 ? "PV" : "PQ"
println("PowerModels ", pkgversion(PowerModels), " and Ipopt ", pkgversion(Ipopt), " loaded.")
## Cell 3 - Load the IEEE 14-bus case and look at what is in it
# Use the provided case14.m stored alongside this Julia script.
p1_file = joinpath(@__DIR__, "case14.m")
p1_data = PowerModels.parse_file(p1_file) # what PowerModels solves
p1_raw = PowerModels.parse_file(p1_file; validate = false) # the file as written
# PowerModels works in per unit on baseMVA = 100: pd = 0.478 means 47.8 MW.
# Multiply P and Q by baseMVA to get MW and MVAr.
p1_MVA = p1_data["baseMVA"]
@printf("IEEE 14-bus case (baseMVA = %.0f): %d buses, %d branches, %d gens, %d loads\n",
p1_MVA, length(p1_data["bus"]), length(p1_data["branch"]),
length(p1_data["gen"]), length(p1_data["load"]))
@printf("Total load: %.1f MW, %.1f MVAr\n\n",
p1_MVA * sum(l["pd"] for l in values(p1_data["load"])),
p1_MVA * sum(l["qd"] for l in values(p1_data["load"])))
# Generator data. In PowerModels the voltage setpoint of a PV or slack bus is the
# BUS "vm" value (changing the generator's "vg" alone does nothing in solve_ac_pf).
@printf("%3s %3s %-5s %6s %6s %6s %6s %6s\n",
"Gen", "Bus", "type", "P set", "Pmax", "Qmin", "Qmax", "V set")
@printf("%3s %3s %-5s %6s %6s %6s %6s %6s\n",
"", "", "", "MW", "MW", "MVAr", "MVAr", "p.u.")
for k in p1_ids(p1_data["gen"])
g = p1_data["gen"][k]
b = p1_data["bus"][string(g["gen_bus"])]
# the slack's P is not a setpoint: the power flow solves for it
pset = b["bus_type"] == 3 ? "free" : @sprintf("%.1f", p1_MVA * g["pg"])
@printf("%3s %3d %-5s %6s %6.1f %6.1f %6.1f %6.3f\n", k, g["gen_bus"],
p1_type(b["bus_type"]), pset, p1_MVA * g["pmax"], p1_MVA * g["qmin"],
p1_MVA * g["qmax"], b["vm"])
end
# What PowerModels changed or dropped while reading the file (silenced warnings)
println("\nNotes on the data:")
for k in p1_ids(p1_data["bus"])
t_file, t_pm = p1_raw["bus"][k]["bus_type"], p1_data["bus"][k]["bus_type"]
if t_file != t_pm
@printf(" * Bus %s is type %d (%s) in the file but has a generator,\n",
k, t_file, p1_type(t_file))
@printf(" so PowerModels makes it type %d (%s).\n", t_pm, p1_type(t_pm))
end
end
for k in p1_ids(p1_data["shunt"])
s = p1_data["shunt"][k]
@printf(" * Existing shunt capacitor at Bus %d: bs = %.2f p.u.",
s["shunt_bus"], s["bs"])
@printf(" = %.1f MVAr at 1.0 p.u.\n", p1_MVA * s["bs"])
end
@printf(" * rate_a missing on %d of %d branches: rateA = 0 in the file",
count(!haskey(br, "rate_a") for br in values(p1_data["branch"])),
length(p1_data["branch"]))
println(" means \"no limit\"")
println(" and PowerModels drops it (thermal ratings are not needed in this project).")
if all(b["base_kv"] == 0 for b in values(p1_data["bus"]))
println(" * base_kv = 0 for every bus in the file: voltages are in p.u. only.")
end
## Cell 4 - AC power flow on the original case (for comparison)
p1_res0 = solve_ac_pf(p1_data, p1_solver)
println("Original case: AC power flow status = ", p1_res0["termination_status"])
p1_solved(p1_res0) || error("the power flow did not solve - check the data")
# solution voltages: vm in p.u., va in radians; gen pg/qg in per unit (x 100 for MW/MVAr)
p1_vm0 = Dict(parse(Int, k) => b["vm"] for (k, b) in p1_res0["solution"]["bus"])
p1_vmin0, p1_bmin0 = findmin(p1_vm0) # findmin on a Dict gives (value, key)
@printf(" Lowest voltage: Bus %d = %.4f p.u.\n", p1_bmin0, p1_vmin0)
@printf(" Slack generator (Bus 1): P = %.1f MW, Q = %.1f MVAr\n",
p1_MVA * p1_res0["solution"]["gen"]["1"]["pg"],
p1_MVA * p1_res0["solution"]["gen"]["1"]["qg"])
## Cell 5 - Create the "stressed" case: every load + 40%
p1_stress = 1.40 # load multiplier: 1.40 means +40%
p1_case = deepcopy(p1_data) # deepcopy: the original case stays untouched
for (k, load) in p1_case["load"]
load["pd"] *= p1_stress # real power
load["qd"] *= p1_stress # reactive power, same factor -> constant power factor
end
p1_P0 = p1_MVA * sum(l["pd"] for l in values(p1_data["load"]))
p1_Q0 = p1_MVA * sum(l["qd"] for l in values(p1_data["load"]))
p1_P1 = p1_MVA * sum(l["pd"] for l in values(p1_case["load"]))
p1_Q1 = p1_MVA * sum(l["qd"] for l in values(p1_case["load"]))
@printf("Stressed case: every load x %.2f", p1_stress)
println(" (P and Q, so each load keeps its power factor)")
@printf(" Total load: %.1f MW -> %.1f MW, %.1f MVAr -> %.1f MVAr\n",
p1_P0, p1_P1, p1_Q0, p1_Q1)
println(" Generator P setpoints are unchanged: the slack (Bus 1) must pick up")
@printf(" the extra %.1f MW plus the extra losses.\n", p1_P1 - p1_P0)
## Cell 6 - AC power flow on the stressed case: bus voltages and the lowest bus
p1_res = solve_ac_pf(p1_case, p1_solver)
println("Stressed case: AC power flow status = ", p1_res["termination_status"])
p1_solved(p1_res) || error("the power flow did not solve - check the data")
p1_vm = Dict(parse(Int, k) => b["vm"] for (k, b) in p1_res["solution"]["bus"])
p1_vmin, p1_bmin = findmin(p1_vm)
# (a) original vs stressed voltage for all 14 buses
@printf("\n%3s %-5s %10s %10s %8s\n",
"Bus", "type", "V original", "V stressed", "change")
@printf("%3s %-5s %10s %10s %8s\n", "", "", "p.u.", "p.u.", "p.u.")
for i in sort(collect(keys(p1_vm)))
flag = p1_vm[i] < 0.95 ? " <-- below 0.95" : ""
i == p1_bmin && (flag *= " <-- lowest")
@printf("%3d %-5s %10.4f %10.4f %8.4f%s\n", i,
p1_type(p1_case["bus"][string(i)]["bus_type"]), p1_vm0[i], p1_vm[i],
p1_vm[i] - p1_vm0[i], flag)
end
# (b) the lowest-voltage bus
p1_lowbuses = sort([i for (i, v) in p1_vm if v < 0.95])
@printf("\nLowest voltage in the stressed case: Bus %d = %.4f p.u.\n", p1_bmin,
p1_vmin)
println("Buses below 0.95 p.u.: ",
isempty(p1_lowbuses) ? "none" : join(p1_lowbuses, ", "))
## Cell 7 - Generators in the stressed case: are they inside their limits?
# IMPORTANT: solve_ac_pf holds every PV bus at its voltage setpoint and does NOT
# enforce the generators' reactive limits (Qmin/Qmax): Q goes wherever it must.
# The slack (Bus 1) supplies whatever P and Q are still missing, also unlimited.
# (The slack's 0..10 MVAr range in case14.m comes from the file's data conversion -
# see the warnings at the end of case14.m - so it is usually not taken literally.)
@printf("%3s %3s %-5s %6s %6s %6s %6s %6s\n",
"Gen", "Bus", "type", "P", "Pmax", "Q", "Qmin", "Qmax")
@printf("%3s %3s %-5s %6s %6s %6s %6s %6s\n",
"", "", "", "MW", "MW", "MVAr", "MVAr", "MVAr")
p1_out_pv, p1_out_slack = Int[], Int[] # buses whose generator is outside its Q range
for k in p1_ids(p1_case["gen"])
g = p1_case["gen"][k]
s = p1_res["solution"]["gen"][k]
btype = p1_type(p1_case["bus"][string(g["gen_bus"])]["bus_type"])
P, Q = p1_MVA * s["pg"], p1_MVA * s["qg"] # per unit -> MW, MVAr
Pmax, Qmin, Qmax = p1_MVA * g["pmax"], p1_MVA * g["qmin"], p1_MVA * g["qmax"]
flags = String[]
P > Pmax + 0.05 && push!(flags, "P above Pmax")
Q > Qmax + 0.05 && push!(flags, "Q above Qmax")
Q < Qmin - 0.05 && push!(flags, "Q below Qmin")
if Q > Qmax + 0.05 || Q < Qmin - 0.05
push!(btype == "slack" ? p1_out_slack : p1_out_pv, g["gen_bus"])
end
flag = isempty(flags) ? "" : " <-- " * join(flags, ", ")
@printf("%3s %3d %-5s %6.1f %6.1f %6.1f %6.1f %6.1f%s\n",
k, g["gen_bus"], btype, P, Pmax, Q, Qmin, Qmax, flag)
end
println("\nGenerators outside their reactive limits: PV buses ",
isempty(p1_out_pv) ? "none" : join(p1_out_pv, ", "), "; slack bus ",
isempty(p1_out_slack) ? "none" : join(p1_out_slack, ", "))
println("solve_ac_pf did not stop them: it does not enforce Qmin/Qmax.")

# =============================================================================
# INSTALLMENT 2
# =============================================================================

## Cell 8 - Respect reactive limits of PV generators and re-solve

# The first stressed power flow may demand more Q than a PV generator can supply.
# If so, hold its Q at the nearest bound and let its voltage float (PV -> PQ).
# This follows the instructor's suggested approach:
#  - turn that generator off in the model;
#  - replace its fixed P and Q contribution with an equivalent negative load;
#  - change the bus to PQ; and
#  - re-run AC power flow, repeating if new PV generators reach their limits.
# The slack generator remains unconstrained in this simplified analysis.

p2_case = deepcopy(p1_case)
p2_res = p1_res
p2_converted = Int[]

for iteration in 1:(length(p2_case["gen"]) + 1)
    new_conversions = Int[]

    for k in p1_ids(p2_case["gen"])
        g = p2_case["gen"][k]

        # Skip generators already converted in a previous iteration.
        get(g, "gen_status", 1) == 1 || continue

        busnum = g["gen_bus"]
        busdata = p2_case["bus"][string(busnum)]

        # Do not convert the reference/slack generator.
        busdata["bus_type"] == 3 && continue

        # Q is in per unit in the PowerModels data and solution.
        Qcalc = p2_res["solution"]["gen"][k]["qg"]
        Qlimit = clamp(Qcalc, g["qmin"], g["qmax"])
        abs(Qcalc - Qlimit) < 1.0e-6 && continue

        # Generator with fixed P and fixed Q can be represented as a
        # negative PQ load at its bus (the injection is -load).
        load_id = maximum(parse(Int, key) for key in keys(p2_case["load"])) + 1
        p2_case["load"][string(load_id)] = Dict{String,Any}(
            "index" => load_id,
            "load_bus" => busnum,
            "pd" => -g["pg"],
            "qd" => -Qlimit,
            "status" => 1,
            "source_id" => ["load", load_id]
        )

        # Remove the PV generator's voltage-control action.
        g["gen_status"] = 0
        busdata["bus_type"] = 1

        push!(new_conversions, busnum)
        @printf("Iteration %d: Bus %d exceeded Q limit (%.2f MVAr). Fixed Q = %.2f MVAr; PV -> PQ\n",
                iteration, busnum, Qcalc * p1_MVA, Qlimit * p1_MVA)
    end

    if isempty(new_conversions)
        println("Reactive-limit iteration complete. No remaining PV Q violations.")
        break
    end

    append!(p2_converted, new_conversions)

    # Re-solve with the Q-limited generators represented as fixed injections.
    global p2_res = solve_ac_pf(p2_case, p1_solver)
    println("Q-limited power flow status = ", p2_res["termination_status"])
    p1_solved(p2_res) || error("Q-limited AC power flow did not converge")
end

println("PV buses converted to PQ: ",
        isempty(p2_converted) ? "none" : join(sort(unique(p2_converted)), ", "))
println("Note: Bus 1's slack generator reactive limit is NOT enforced here.")


## Cell 9 - Voltage table for ALL 14 buses after enforcing PV Q limits

p2_vm = Dict(parse(Int, k) => b["vm"]
             for (k, b) in p2_res["solution"]["bus"])
p2_vmin, p2_problem_bus = findmin(p2_vm)
p2_below = sort([i for (i, v) in p2_vm if v < 0.95])

println("\nINSTALLMENT 2 -- ALL BUS VOLTAGES")
@printf("%3s %-5s %12s %12s %12s %s\n",
        "Bus", "type", "Original", "Stressed", "Q-limited", "Status")
@printf("%3s %-5s %12s %12s %12s %s\n",
        "", "", "(p.u.)", "(p.u.)", "(p.u.)", "")

for i in sort(collect(keys(p2_vm)))
    flag = p2_vm[i] < 0.95 ? "BELOW 0.95" : "OK"
    i == p2_problem_bus && (flag *= " <-- LOWEST")
    @printf("%3d %-5s %12.5f %12.5f %12.5f %s\n",
            i,
            p1_type(p2_case["bus"][string(i)]["bus_type"]),
            p1_vm0[i], p1_vm[i], p2_vm[i], flag)
end

@printf("\nLowest Q-limited voltage: Bus %d = %.5f p.u.\n",
        p2_problem_bus, p2_vmin)
println("Buses below 0.95 p.u.: ",
        isempty(p2_below) ? "none" : join(p2_below, ", "))

# Optional: save the voltage table for use in a report.
p2_csv = joinpath(@__DIR__, "Installment2_bus_voltages.csv")
open(p2_csv, "w") do io
    println(io, "bus,original_pu,stressed_unlimited_Q_pu,stressed_Q_limited_pu,below_0_95")
    for i in sort(collect(keys(p2_vm)))
        println(io, join((i, round(p1_vm0[i]; digits=7),
                           round(p1_vm[i]; digits=7),
                           round(p2_vm[i]; digits=7), p2_vm[i] < 0.95), ","))
    end
end
println("Voltage table saved to: ", p2_csv)


## Cell 10 - Quantify Q needed to hold the problem bus at 0.95 p.u.

# Define the deficit precisely: extra local reactive power injection (MVAr)
# needed to raise the weakest bus to exactly 0.95 p.u., while keeping
# the other operating assumptions fixed.
# This is NOT the final capacitor rating (Installment 3).

if p2_vmin >= 0.95
    println("\nAll buses are already at or above 0.95 p.u.")
    println("Minimum extra Q injection under this definition = 0.0 MVAr.")
else
    p2_test = deepcopy(p2_case)
    p2_bus_key = string(p2_problem_bus)

    # Create a fictitious synchronous condenser
    # (generator with P = 0 and a wide reactive power capability range).
    p2_test["bus"][p2_bus_key]["bus_type"] == 1 ||
        error("The problem bus is not PQ; inspect the generator model before proceeding")

    gen_id = maximum(parse(Int, k) for k in keys(p2_test["gen"])) + 1
    condenser = deepcopy(p2_test["gen"]["1"])
    condenser["index"] = gen_id
    condenser["source_id"] = ["gen", gen_id]
    condenser["gen_bus"] = p2_problem_bus
    condenser["gen_status"] = 1
    condenser["pg"] = 0.0
    condenser["qg"] = 0.0
    condenser["pmin"] = 0.0
    condenser["pmax"] = 0.0
    condenser["qmin"] = -10.0  # per unit = -1000 MVAr at 100 MVA
    condenser["qmax"] = 10.0   # per unit = +1000 MVAr at 100 MVA
    condenser["vg"] = 0.95
    p2_test["gen"][string(gen_id)] = condenser

    # PowerModels solve_ac_pf uses the BUS voltage setpoint for PV buses.
    p2_test["bus"][p2_bus_key]["bus_type"] = 2
    p2_test["bus"][p2_bus_key]["vm"] = 0.95

    p2_target_res = solve_ac_pf(p2_test, p1_solver)
    println("\n0.95-p.u. target power flow status = ",
            p2_target_res["termination_status"])
    p1_solved(p2_target_res) || error("Target-voltage power flow did not converge")

    p2_Q_deficit = p1_MVA * p2_target_res["solution"]["gen"][string(gen_id)]["qg"]
    p2_target_vm = p2_target_res["solution"]["bus"][p2_bus_key]["vm"]

    @printf("Problem bus: Bus %d\n", p2_problem_bus)
    @printf("Voltage before local Q support: %.5f p.u.\n", p2_vmin)
    @printf("Target voltage: %.5f p.u.\n", p2_target_vm)
    @printf("Extra reactive power injection required: %.4f MVAr\n", p2_Q_deficit)

    if p2_Q_deficit < 0
        println("Negative Q indicates reactive power absorption rather than injection.")
    end
    println("This is a local target-voltage Q requirement, NOT final capacitor sizing.")
    println("The other bus voltages and remaining operating limits need checking")
    println("when choosing a physical capacitor in Installment 3.")
end

# End of Installment 2.  The baseline and stressed 40% load cases are unchanged.
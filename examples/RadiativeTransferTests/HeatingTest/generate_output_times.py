#!/usr/bin/env python3

# ----------------------------------------------
# Generate an output times list to your liking
# ----------------------------------------------

import yaml

with open(r"rt_heating_test.yml") as paramfile:
    params = yaml.load(paramfile, Loader=yaml.FullLoader)

    unit_l = params["InternalUnitSystem"]["UnitLength_in_cgs"]
    unit_l = float(unit_l)
    unit_v = params["InternalUnitSystem"]["UnitVelocity_in_cgs"]
    unit_v = float(unit_v)
    t_begin = params["TimeIntegration"]["time_begin"]
    t_begin = float(t_begin)
    t_end = params["TimeIntegration"]["time_end"]
    t_end = float(t_end)
    dt_min = params["TimeIntegration"]["dt_min"]
    dt_min = float(dt_min)
    dt_max = params["TimeIntegration"]["dt_max"]
    dt_max = float(dt_max)

unit_t = unit_l / unit_v
unit_myr = unit_t / (3600 * 24 * 365 * 1e6)


# If you're rebuilding this output times list:
# I read this 'dt' out from a run, then used it to generate the output list
dt_heat = 1.001182e-03

current_t = t_begin
outputtimes = []
# first 16 snapshots
for i in range(16):
    current_t += dt_heat
    outputtimes.append(current_t)
for i in range(16):
    current_t += 4 * dt_heat
    outputtimes.append(current_t)
for i in range(16):
    current_t += 8 * dt_heat
    outputtimes.append(current_t)
for i in range(16):
    current_t += 16 * dt_heat
    outputtimes.append(current_t)
# fill up until heating stops
while current_t + 32 * dt_heat < t_end:
    current_t += 32 * dt_heat
    outputtimes.append(current_t)
outputtimes.append(t_end)

print("preparing", len(outputtimes), "snapshots")

with open(r"snaplist.txt", "w") as outfile:
    # THIS HEADER IS IMPORTANT
    outfile.write("# Time\n")
    for t in outputtimes:
        outfile.write("{0:.6g}\n".format(t))

print("written snaplist.txt")

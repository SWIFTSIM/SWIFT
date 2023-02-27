#!/bin/bash
../../../swift_van --hydro --threads=16 BrioWu.yml
julia plotSolution.jl BrioWu_0004.hdf5 Solution.png
../../../swift_whr --hydro --threads=16 BrioWu_whr.yml
julia plotSolution.jl BrioWu_whr_0004.hdf5 Solution_whr.png
../../../swift_wr --hydro --threads=16 BrioWu_wr.yml
julia plotSolution.jl BrioWu_wr_0004.hdf5 Solution_wr.png
../../../swift_wd --hydro --threads=16 BrioWu_wd.yml
julia plotSolution.jl BrioWu_wd_0004.hdf5 Solution_wd.png
../../../swift_wd_wr --hydro --threads=16 BrioWu_wd_wr.yml
julia plotSolution.jl BrioWu_wd_wr_0004.hdf5 Solution_wd_wr.png

The 'run_test_program.py' script automatically configures and runs with parameters selected in the test_run_parameters.csv file
The results of the run are stored in 'test_results' folder
After a run is complete the script puts all snapshots and run data into corresponding folder in 'test_results' and fills the status column of the corresponding run with 'done'. If there is a need to skip some run just write 'done' in the status column of .csv

Run parameter file 
test_run_parameters.csv has the following parameters:
1) v0 - rms velocity of the flow
2) Vz_factor - prefactor before z component of the velocity field
3) eta - physical resistivity in the run
4) kv - wavenumber of the velocity field
5) kb - wavenumber of the magnetic field
6) Lbox - simulation box size
7) Flow_kind - the kind of Roberts Flow. 0 - as in A.Brandenburg papers, 1 to 4 - the corresponding flows from Roberts paper
8) Scheme - MHD scheme to configure for a run. ODI - direct induction, VP - vector potential, FDI - simpler direct induction.
9) IAfile - initial arrangement file. g32 means glassCube32.hdf5
10) Npar, Nper - how many boxes to stack along xy or z direction
11) LparMul, LperMul - deform box length without adding particles in xy or z direction 
12) monopole_subtraction, artificial_diffusion, hyperbolic_dedner, hyperbolic_dedner_divv, parabolic_dedner - parameters for Direct induction schemes
13) t_max - maximal simulation time
14) max_dt - maximal timestep size
15) Forcing_kind - forcing regime for RobertsFlow dynamo. If 'v' was picked, RobertsFlow is in kinematic regime (that is B field does not affect the flow). If 'a' was picked, RobertsFlow is driven by external force proportional to artificial viscosity.
16) viscosity_alpha - determines the magnitude of the external force
17) Status - the script employs this column to see what runs to ignore (ignores the runs that are filled with 'done')
18) Comment - stores additional information about a run

To run the test program type the following:
'python run_test_program.py'
We recommend running the script in separate screen
Note: number of threads to run the script can be specified by modifying 'threads' parameter in the script

Plotting slices:
'plot_snapshot_slice.py' script allows to plot an SPH slice of the simulation volume. To run the script type the following:
'python plot_snapshot_slice.py #addr_to_snapshot #addr_to_output_png #slice_height #slice_plane'
where #slice_height - z_slice/L_z_box - position of the slice in the box, is in the interval [0.0,1.0]
#slice_plane - in what plane to take a slice (xy,yz,zx)




The 'run_test_program.py' script automatically configures and runs with parameters selected in the .csv file
The results of the run are stored in 'test_results' folder
After a run is complete the script puts all snapshots and run data into corresponding folder in 'test_results' and fills the status column of the corresponding run with 'done'. If there is a need to skip some run just write 'done' in the status column of .csv

To run the test program type the following:
'python run_test_program.py'
We recommend running the script in separate screen

'plot_snapshot_slice.py' script allows to plot an SPH slice of the simulation volume. To run the script type the following:
'python plot_snapshot_slice.py #addr_to_snapshot #addr_to_output_png #slice_height #slice_plane'
where #slice_height - z_slice/L_z_box - position of the slice in the box, is in the interval [0.0,1.0]
#slice_plane - in what plane to take a slice (xy,yz,zx)




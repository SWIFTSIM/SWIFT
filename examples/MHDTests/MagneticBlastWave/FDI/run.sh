# Run SWIFT
../../../../sw_FDI --hydro --threads=16 ../BW_schemes.yml 2>&1 > out.log 

# Plot the temperature evolution
python3 ../plot_all.py 0 60

#Notes:
#Load swift module and gcc compiler to run sanitizer.

#Configure
./configure --enable-sanitizer --enable-debug --enable-optimization=no  --enable-mpi=no


# -d: minimum time step
# -e: maximum time step
# -c: end of simulation
 ./swift_fixdt -m .1 -s "50 50 50"  -t 1 -d 0.001 -e 0.001 -c 0.2 -f ExternalGravity/Sphere.hdf5

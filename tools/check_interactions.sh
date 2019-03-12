#!/bin/bash -l
#
# Runs the SedovBlast_3D and EAGLE_12 examples using naive, serial and vectorised particle interactions. 
# Then compares the output between versions using check_ngbs.py. Test is performed with and without MPI.

echo
echo "# Running SedovBlast_3D and EAGLE_12 with naive interactions and neighbour logging, 16 thread"
echo

cd ../

./autogen.sh

# Naive interactions run
./configure --disable-mpi --enable-debug-interactions=1024 --disable-vec --enable-naive-interactions
make clean; make -j 6

cd examples/SedovBlast_3D/

./getGlass.sh
python makeIC.py

../swift --hydro --threads=16 --steps=5 sedov.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7 

mv sedov_0000.hdf5 sedov_naive.hdf5

cd ../EAGLE_12/

# Link to ICs
ln -s /gpfs/data/Swift/web-storage/ICs/EAGLE_ICs_12.hdf5 EAGLE_ICs_12.hdf5

../swift --hydro --threads=16 --steps=5 eagle_12.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv eagle_0000.hdf5 eagle_12_naive.hdf5

cd ../../

echo
echo "# Running SedovBlast_3D and EAGLE_12 with serial interactions and neighbour logging, 16 thread"
echo

# Serial interactions run
./configure --disable-mpi --enable-debug-interactions=1024 --disable-vec
make clean; make -j 6

cd examples/SedovBlast_3D/

../swift --hydro --threads=16 --steps=5 sedov.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv sedov_0000.hdf5 sedov_serial.hdf5

cd ../EAGLE_12/

../swift --hydro --threads=16 --steps=5 eagle_12.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7 

mv eagle_0000.hdf5 eagle_12_serial.hdf5

cd ../../

echo
echo "# Running SedovBlast_3D and EAGLE_12 with vectorised interactions and neighbour logging, 16 thread"
echo

# Vectorised interactions run
./configure --disable-mpi --enable-debug-interactions=1024
make clean; make -j 6

cd examples/SedovBlast_3D/

../swift --hydro --threads=16 --steps=5 sedov.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv sedov_0000.hdf5 sedov_vec.hdf5

# Compare outputs
if python ../check_ngbs.py sedov_naive.hdf5 sedov_serial.hdf5 
then
  echo "SedovBlast_3D comparison between naive and serial passed"
else
  echo "SedovBlast_3D comparison between naive and serial failed"
  exit 1
fi

if python ../check_ngbs.py sedov_naive.hdf5 sedov_vec.hdf5 
then
  echo "SedovBlast_3D comparison between naive and vectorised passed"
else
  echo "SedovBlast_3D comparison between naive and vectorised failed"
  exit 1
fi

if python ../check_ngbs.py sedov_serial.hdf5 sedov_vec.hdf5 
then
  echo "SedovBlast_3D comparison between serial and vectorised passed"
else
  echo "SedovBlast_3D comparison between serial and vectorised failed"
  exit 1
fi

cd ../EAGLE_12/

../swift --hydro --threads=16 --steps=5 eagle_12.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv eagle_0000.hdf5 eagle_12_vec.hdf5

# Compare outputs
if python ../check_ngbs.py eagle_12_naive.hdf5 eagle_12_serial.hdf5 
then
  echo "EAGLE_12 comparison between naive and serial passed"
else
  echo "EAGLE_12 comparison between naive and serial failed"
  exit 1
fi

if python ../check_ngbs.py eagle_12_naive.hdf5 eagle_12_vec.hdf5 
then
  echo "EAGLE_12 comparison between naive and vectorised passed"
else
  echo "EAGLE_12 comparison between naive and vectorised failed"
  exit 1
fi

if python ../check_ngbs.py eagle_12_serial.hdf5 eagle_12_vec.hdf5 
then
  echo "EAGLE_12 comparison between serial and vectorised passed"
else
  echo "EAGLE_12 comparison between serial and vectorised failed"
  exit 1
fi

# Now run the same test using MPI

# Runs the SedovBlast_3D and EAGLE_12 examples using naive, serial and vectorised particle interactions. Then compares the output between versions using check_ngbs.py

unset I_MPI_HYDRA_BOOTSTRAP

echo
echo "# Running SedovBlast_3D and EAGLE_12 with naive interactions and neighbour logging, 16 thread, 4 MPI ranks"
echo

cd ../../

# Naive interactions run
./configure --with-metis --enable-debug-interactions=1024 --disable-vec --enable-naive-interactions
make clean; make -j 6

cd examples/SedovBlast_3D/

mpirun -np 4 ../swift_mpi --hydro --threads=16 --steps=5 sedov.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv sedov_0000.hdf5 sedov_naive.hdf5

cd ../EAGLE_12/

mpirun -np 4 ../swift_mpi --hydro --threads=16 --steps=5 eagle_12.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv eagle_0000.hdf5 eagle_12_naive.hdf5

cd ../../

echo
echo "# Running SedovBlast_3D and EAGLE_12 with serial interactions and neighbour logging, 16 thread, 4 MPI ranks"
echo

# Serial interactions run
./configure --with-metis --enable-debug-interactions=1024 --disable-vec
make clean; make -j 6

cd examples/SedovBlast_3D/

mpirun -np 4 ../swift_mpi --hydro --threads=16 --steps=5 sedov.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv sedov_0000.hdf5 sedov_serial.hdf5

cd ../EAGLE_12/

mpirun -np 4 ../swift_mpi --hydro --threads=16 --steps=5 eagle_12.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv eagle_0000.hdf5 eagle_12_serial.hdf5

cd ../../

echo
echo "# Running SedovBlast_3D and EAGLE_12 with vectorised interactions and neighbour logging, 16 thread, 4 MPI ranks"
echo

# Vectorised interactions run
./configure --with-metis --enable-debug-interactions=1024
make clean; make -j 6

cd examples/SedovBlast_3D/

mpirun -np 4 ../swift_mpi --hydro --threads=16 --steps=5 sedov.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv sedov_0000.hdf5 sedov_vec.hdf5

# Compare outputs
if python ../check_ngbs.py sedov_naive.hdf5 sedov_serial.hdf5 
then
  echo "SedovBlast_3D comparison between naive and serial passed (MPI)"
else
  echo "SedovBlast_3D comparison between naive and serial failed (MPI)"
  exit 1
fi

if python ../check_ngbs.py sedov_naive.hdf5 sedov_vec.hdf5 
then
  echo "SedovBlast_3D comparison between naive and vectorised passed (MPI)"
else
  echo "SedovBlast_3D comparison between naive and vectorised failed (MPI)"
  exit 1
fi

if python ../check_ngbs.py sedov_serial.hdf5 sedov_vec.hdf5 
then
  echo "SedovBlast_3D comparison between serial and vectorised passed (MPI)"
else
  echo "SedovBlast_3D comparison between serial and vectorised failed (MPI)"
  exit 1
fi

cd ../EAGLE_12/

mpirun -np 4 ../swift_mpi --hydro --threads=16 --steps=5 eagle_12.yml -P SPH:h_tolerance:10 -P Snapshots:compression:7

mv eagle_0000.hdf5 eagle_12_vec.hdf5

# Compare outputs
if python ../check_ngbs.py eagle_12_naive.hdf5 eagle_12_serial.hdf5 
then
  echo "EAGLE_12 comparison between naive and serial passed (MPI)"
else
  echo "EAGLE_12 comparison between naive and serial failed (MPI)"
  exit 1
fi

if python ../check_ngbs.py eagle_12_naive.hdf5 eagle_12_vec.hdf5 
then
  echo "EAGLE_12 comparison between naive and vectorised passed (MPI)"
else
  echo "EAGLE_12 comparison between naive and vectorised failed (MPI)"
  exit 1
fi

if python ../check_ngbs.py eagle_12_serial.hdf5 eagle_12_vec.hdf5 
then
  echo "EAGLE_12 comparison between serial and vectorised passed (MPI)"
  exit 0
else
  echo "EAGLE_12 comparison between serial and vectorised failed (MPI)"
  exit 1
fi

#export CPPFLAGS="-I/opt/ohpc/pub/libs/intel/hdf5/1.10.8/include -I/data/contrib/pisca_79/grackle_intel2023.2.1/include"
#export LDFLAGS="-L/opt/ohpc/pub/libs/intel/hdf5/1.10.8/lib -L/data/contrib/pisca_79/grackle_intel2023.2.1/lib"
module purge
module load cmake autotools 
module load intel/2024.0.0
#module load openmpi5
module load impi
module load gsl hdf5 fftw
module load ucx
module load metis
module list

echo $FFTW_DIR
echo $METIS_DIR

echo " -----------------------------"

./configure  --with-fftw=$FFTW_DIR \
	     --with-metis=$METIS_DIR \
	     --with-tbbmalloc \
             --with-kernel=quintic-spline  --with-spmhd=direct-induction \
             --with-hydro=sphenix \
	     --disable-compiler-warnings --disable-hand-vec \
             --enable-fof
 ##            CC=icc MPICC=mpiicc FC=ifort

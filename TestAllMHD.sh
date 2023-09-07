#/bin/bash

echo "Testing MHD schemes"

rm MHD_Pass.log

case $1 in
   [01])
	rm sw_vp*  MHD_VP.log
   	echo "Compiling Vector Potential version" 

   	./configure --with-spmhd=vector-potential > MHD_VP.log
	make -j 32 >> MHD_VP.log
	mv swift sw_vp
	mv swift_mpi sw_vp_mpi
	
	echo "MHD: VP OK!" >> MHD_Pass.log
   ;;&
   [02])
	rm sw_di*  MHD_DI.log
        echo "Compiling Direct Induction *Orestis* version" 

	./configure --with-spmhd=direct-induction > MHD_DI.log
	make -j 32 >> MHD_DI.log
	mv swift sw_di
	mv swift_mpi sw_di_mpi
	
	echo "MHD: DI Orestis OK!" >> MHD_Pass.log
   ;;&
   [03])
	rm sw_diF*  MHD_DIF.log
	echo "Compiling Direct Induction *Fede* version" 
	
	./configure --with-spmhd=direct-induction-fede > MHD_DIF.log
	make -j 32 >> MHD_DIF.log
	mv swift sw_diF
	mv swift_mpi sw_diF_mpi
	
	echo "MHD: DI Fede OK!" >> MHD_Pass.log
   ;;
   *)
	echo "Usage $0 [0123]"
	echo "0 all schemes"
	echo "1 vector potenial schemes"
	echo "2 Direct Induction Orestis schemes"
	echo "3 Direct Induction Federico schemes"
   ;;
esac

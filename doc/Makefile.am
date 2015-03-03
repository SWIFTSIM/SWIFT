doxyfile.stamp: 
	@DX_ENV@ @DX_DOXYGEN@ Doxyfile
	echo `date` > doxyfile.stamp

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp

clean-local:
	rm -rf html/ latex/ man/

distclean: clean

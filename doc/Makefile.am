
doxyfile.stamp: 
if HAVE_DOXYGEN
	@DX_ENV@ @DX_DOXYGEN@ Doxyfile
	echo `date` > doxyfile.stamp
endif

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp

clean-local:
	rm -rf html/ latex/ man/

distclean: clean

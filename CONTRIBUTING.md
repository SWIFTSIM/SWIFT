The SWIFT source code is using a variation of the 'Google' formatting style. 
The script 'format.sh' in the root directory applies the clang-format-18
tool with our style choices to all the SWIFT C source file. Please apply 
the formatting script to the files before submitting a merge request.

The Python scripts in the repository use default formatting enforced by 
'black'. Since 'black' formatting still strongly depends on the version, 
we currently enforce use of 'black==19.3b0' using a virtual environment. 
The setup of the virtual environment and formatting of all PYthon 
scripts is done automatically by the Python formatting script 
'format_python.sh'. Please apply this script before submitting a merge 
request containing new or updated Python scripts.

The SWIFT code comes with a series of unit tests that are run automatically 
when a push to the master branch occurs. The suite can be run by doing a `make 
check` in the root directory. Please check that the test suite still
runs with your changes applied before submitting a merge request and add 
relevant unit tests probing the correctness of new modules. An example of how
to add a test to the suite can be found by considering the tests/testGreeting 
case.

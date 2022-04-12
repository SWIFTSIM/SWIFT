#!/bin/bash

# Check if we can run pip
# This also serves as a check for python3
python3 -m pip --version > /dev/null
if [[ $? -ne 0 ]]
then
  echo "ERROR: cannot run 'python3 -m pip'"
  exit 1
fi

# Check if the virtual environment with black exists
if [ ! -d black_formatting_env ]
then
  echo "Formatting environment not found, installing it..."
  python3 -m venv black_formatting_env
  ./black_formatting_env/bin/python3 -m pip install click==8.0.4 black==19.3b0
fi
# Now we know exactly which black to use
black="./black_formatting_env/bin/python3 -m black"

# Formatting command
cmd="$black -t py38 $(git ls-files | grep '\.py$')"

# Print the help
function show_help {
    echo -e "This script formats SWIFT's Python scripts using black"
    echo -e "  -h, --help \t Show this help"
    echo -e "  -t, --test \t Test if SWIFT's Python scripts are well formatted"
}

# Parse arguments
TEST=0
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
	# print the help and exit
	-h|--help)
	    show_help
	    exit
	    ;;
	# check if the code is well formatted
	-t|--test)
	    TEST=1
	    shift
	    ;;
	# unknown option
	*)
	    echo "Argument '$1' not implemented"
	    show_help
	    exit
	    ;;
    esac
done

# Run the required commands
if [[ $TEST -eq 1 ]]
then
    # Note trapping the exit status from both commands in the pipe. Also note
    # do not use -q in grep as that closes the pipe on first match and we get
    # a SIGPIPE error.
    echo "Testing if SWIFT's Python is correctly formatted"
    $cmd --check
    status=$?

    # Check formatting
    if [[ ! ${status} -eq 0 ]]
    then
 	echo "ERROR: needs formatting"
 	exit 1
    else
        echo "Everything is correctly formatted"
    fi
else
    echo "Formatting SWIFT Python scripts"
    $cmd
fi

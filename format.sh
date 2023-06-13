#!/bin/bash

CLANG_FORMAT_VERSION="16.0.5"

# Check if we can run pip
# This also serves as a check for python3
python3 -m pip --version > /dev/null
if [[ $? -ne 0 ]]
then
  echo "ERROR: cannot run 'python3 -m pip'"
  exit 1
fi


# Clang format command, can be overridden using CLANG_FORMAT_CMD.
# is CLANG_FORMAT_CMD provided? Then use that.
if [ ! -z ${CLANG_FORMAT_CMD+x} ]
then
    echo GOT THE COMMAND
    echo $CLANG_FORMAT_CMD
    clang="$CLANG_FORMAT_CMD"
else
    # Check if the virtual environment exists
    if [ ! -d .formatting_python_env ]
    then
      echo "Formatting environment not found, installing it..."
      python3 -m venv .formatting_python_env
    fi

    # Check if clang-format executable exists
    clang=".formatting_python_env/bin/clang-format"
    if [ ! -f "$clang" ]
    then
        echo "Installing clang-format"
        ./.formatting_python_env/bin/python3 -m pip install clang-format=="$CLANG_FORMAT_VERSION"
    fi
fi

# Test if `clang-format` works
command -v $clang > /dev/null
if [[ $? -ne 0 ]]
then
    echo "ERROR: cannot find $clang"
    exit 1
fi

# Check that we have the correct version
$clang --version | /usr/bin/grep "$CLANG_FORMAT_VERSION" >> /dev/null
if [ "$?" -eq 1 ]
then
    echo "Wrong version of clang-format installed. I need" "$CLANG_FORMAT_VERSION"
    echo "You've got"
    $clang --version

    if [ -z ${CLANG_FORMAT_CMD+x} ]
    then
        # warning for pip-installed clang-format only
        echo "remove the contents of directory '.formatting_python_env' and re-run this script"
    fi
    exit 1
fi

# Formatting command
cmd="$clang -style=file $(git ls-files | grep '\.[ch]$')"

# Print the help
function show_help {
    echo -e "This script formats SWIFT according to Google style"
    echo -e "  -h, --help \t Show this help"
    echo -e "  -t, --test \t Test if SWIFT is well formatted"
}

# Parse arguments (based on https://stackoverflow.com/questions/192249)
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
    echo "Testing if SWIFT is correctly formatted"
    $cmd -output-replacements-xml | grep "<replacement " > /dev/null
    status=("${PIPESTATUS[@]}")

    #  Trap if first command failed. Note 141 is SIGPIPE, that happens when no
    #  output
    if [[ ${status[0]} -ne 0 ]]
    then
       echo "ERROR: $clang command failed"
       exit 1
    fi

    # Check formatting
    if [[ ${status[1]} -eq 0 ]]
    then
 	echo "ERROR: needs formatting"
 	exit 1
    else
        echo "...is correctly formatted"
    fi
else
    echo "Formatting SWIFT"
    $cmd -i
fi

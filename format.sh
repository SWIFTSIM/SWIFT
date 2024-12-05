#!/bin/bash

# The clang-format command can be overridden using CLANG_FORMAT_CMD.
# We currrently use version 18.0 so any overrides should use that version
# or one known to be compatible with it for instance if your standard
# command is version 18 use:
#    CLANG_FORMAT_CMD=clang-format ./format.sh
clang=${CLANG_FORMAT_CMD:="clang-format-18"}

# Only works in a git checkout. Eat error and just make a report
# so we don't break CI.
git status > /dev/null 2>&1
if test "$?" != "0"; then
   echo "Not operating in a git checkout. Cannot procede."
   exit
fi

# Formatting command
cmd="$clang -style=file $(git ls-files | grep '\.[ch]$' | grep -v csds_io\.h | grep -v argparse\.h | grep -v vector\.h)"

# Test if `clang-format-18` works
command -v $clang > /dev/null
if [[ $? -ne 0 ]]
then
    echo "ERROR: cannot find the command $clang."
    echo
    head -8 "$0" | grep -v "/bin/bash" | grep '^#'
    echo
    exit 1
fi

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

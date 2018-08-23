#!/bin/bash

# Clang software name
clang="clang-format-5.0"

# Formatting command
cmd="$clang -style=file src/*.[ch] src/*/*.[ch] src/*/*/*.[ch] examples/main.c tests/*.[ch]"

# Test if `clang-format-5.0` works
command -v $clang > /dev/null
if [[ $? -ne 0 ]]
then
    echo "Cannot find $clang"
    exit 1
fi

set -e


# Print the help
function show_help {
    echo -e "This script formats Swift according to Google style"
    echo -e "  -h, --help \t Show this help"
    echo -e "  -t, --test \t Test if Swift is well formatted"
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
    echo "Testing Swift formatting"
    $cmd -output-replacements-xml | grep -c "<replacement " >/dev/null
    res=$?

    # Check formatting
    if [[ $res -ne 1 ]]
    then
	echo "Swift needs to be formatted"
	exit 1
    fi
else
    echo "Formatting Swift"
    $cmd -i
fi

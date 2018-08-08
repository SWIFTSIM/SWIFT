#!/bin/bash

clang-format-5.0 -style=file -i src/*.[ch] src/*/*.[ch] src/*/*/*.[ch] examples/main.c tests/*.[ch]

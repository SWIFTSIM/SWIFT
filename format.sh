#!/bin/bash

clang-format-3.8 -style=file -i src/*.[ch] src/*/*.[ch] src/*/*/*.[ch] examples/main.c tests/*.[ch]

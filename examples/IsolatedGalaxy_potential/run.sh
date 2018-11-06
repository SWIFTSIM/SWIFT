#!/bin/bash

../swift -g -G -S -t 8 isolated_galaxy.yml 2>&1 | tee output.log


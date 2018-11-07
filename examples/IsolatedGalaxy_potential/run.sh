#!/bin/bash

../swift -g -G -S -t 16 isolated_galaxy.yml 2>&1 | tee output.log


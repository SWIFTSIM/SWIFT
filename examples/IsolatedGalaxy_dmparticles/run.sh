#!/bin/bash

../swift -G -S -t 64 isolated_galaxy.yml 2>&1 | tee output.log


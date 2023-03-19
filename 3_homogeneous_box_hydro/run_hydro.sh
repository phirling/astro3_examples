#!/bin/bash

../../swiftsim/swift --hydro --threads=4 ../params_hydro.yml 2>&1 | tee output.log

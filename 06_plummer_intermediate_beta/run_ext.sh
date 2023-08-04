#!/bin/bash
SWIFT=../swiftsim/swift     # Location of the SWIFT executable file
NTHREADS=4                  # Number of MPI ranks to use (e.g. number of cores)

$SWIFT --external-gravity --threads=$NTHREADS params.yml 2>&1 | tee output.log

#!/bin/bash

../7_miyamoto_nagai_circular_extgrav/swiftsim/swift --self-gravity --threads=4 params.yml 2>&1 | tee output.log

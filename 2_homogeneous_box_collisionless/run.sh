#!/bin/bash

../../../swiftsim/swift --self-gravity --threads=4 params.yml 2>&1 | tee output.log

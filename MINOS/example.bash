#!/bin/bash
#
# Calculate modes for all five end-member models.
#
for N in 1 2 3 4 5; do
bash run_minos.bash sharp/Model${N}wDiscontinuities
done

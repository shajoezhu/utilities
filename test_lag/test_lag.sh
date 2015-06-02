#!/bin/bash

for lag in 10 50 100 500 1000
do
output="lag${lag}.txt"
~/github/my_fork/smcsmc/smcsmc -seg sim-1Samples2msdata1.seg -p "1*3+15*4+1" -tmax 4 -Np 100 -EM 20 -xc 2-17 -xr 2-17 -lag $lag > $output
done

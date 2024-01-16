#!/bin/bash

size=64
gens=1
method=$1

mkdir ${method}_small
dir=${method}_small
echo $dir

mpirun --oversubscribe -np 2 ./Ex2 $gens $size $size 0 50 1 1 2 $method 1 > $dir/testout_2_1_2.txt
mpirun --oversubscribe -np 2 ./Ex2 $gens $size $size 0 50 1 2 1 $method 1 > $dir/testout_2_2_1.txt

mpirun --oversubscribe -np 4 ./Ex2 $gens $size $size 0 50 1 2 2 $method 1 > $dir/testout_4_2_2.txt
mpirun --oversubscribe -np 4 ./Ex2 $gens $size $size 0 50 1 1 4 $method 1 > $dir/testout_4_1_4.txt
mpirun --oversubscribe -np 4 ./Ex2 $gens $size $size 0 50 1 4 1 $method 1 > $dir/testout_4_4_1.txt



mpirun --oversubscribe -np 1 Ex1 $gens $size $size 0 50 1  > $dir/testout_serial.txt





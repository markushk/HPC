#!/bin/bash

size=128
gens=10
method=$1

mkdir $method

mpirun --oversubscribe -np 2 ./main $gens $size $size 0 50 1 1 2 $method > ${method}/testout_2_1_2.txt
mpirun --oversubscribe -np 2 ./main $gens $size $size 0 50 1 2 1 $method > $method/testout_2_2_1.txt

mpirun --oversubscribe -np 4 ./main $gens $size $size 0 50 1 2 2 $method > $method/testout_4_2_2.txt
mpirun --oversubscribe -np 4 ./main $gens $size $size 0 50 1 1 4 $method > $method/testout_4_1_4.txt
mpirun --oversubscribe -np 4 ./main $gens $size $size 0 50 1 4 1 $method > $method/testout_4_4_1.txt

mpirun --oversubscribe -np 8 ./main $gens $size $size 0 50 1 1 8 $method > $method/testout_8_1_8.txt
mpirun --oversubscribe -np 8 ./main $gens $size $size 0 50 1 8 1 $method > $method/testout_8_8_1.txt
mpirun --oversubscribe -np 8 ./main $gens $size $size 0 50 1 2 4 $method > $method/testout_8_2_4.txt
mpirun --oversubscribe -np 8 ./main $gens $size $size 0 50 1 4 2 $method > $method/testout_8_4_2.txt

mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 1 16 $method > $method/testout_16_1_16.txt
mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 2 8 $method > $method/testout_16_2_8.txt
mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 4 4 $method > $method/testout_16_4_4.txt
mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 16 1 $method > $method/testout_16_16_1.txt
mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 8 2 $method > $method/testout_16_8_2.txt


mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 1 32 $method > $method/testout_32_1_32.txt
mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 2 16 $method > $method/testout_32_2_16.txt
mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 4 8 $method > $method/testout_32_4_8.txt
mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 8 4 $method > $method/testout_32_8_4.txt
mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 16 2 $method > $method/testout_32_16_2.txt
mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 32 1 $method > $method/testout_32_32_1.txt

#mpirun --oversubscribe -np 256 ./main $gens $size $size 0 50 1 1 256 $method > $method/testout_256_1_256.txt
#mpirun --oversubscribe -np 256 ./main $gens $size $size 0 50 1 1 256 $method > $method/testout_256_256_1.txt

#mpirun --oversubscribe -np 256 ./main $gens $size $size 0 50 1 8 32 $method > $method/testout_256_8_32.txt
#mpirun --oversubscribe -np 256 ./main $gens $size $size 0 50 1 32 8 $method > $method/testout_256_32_8.txt


#mpirun --oversubscribe -np 512 ./main $gens $size $size 0 50 1 1 512 $method > $method/testout_512_1_512.txt
#mpirun --oversubscribe -np 512 ./main $gens $size $size 0 50 1 512 1 $method > $method/testout_512_512_1.txt

#mpirun --oversubscribe -np 512 ./main $gens $size $size 0 50 1 16 32 $method > $method/testout_512_16_32.txt
#mpirun --oversubscribe -np 512 ./main $gens $size $size 0 50 1 32 16 $method > $method/testout_512_32_16.txt

#mpirun --oversubscribe -np 1024 ./main $gens $size $size 0 50 1 1024 1 $method > $method/testout_1024_1024_1.txt
#mpirun --oversubscribe -np 1024 ./main $gens $size $size 0 50 1 1 1024 $method > $method/testout_1024_1_1024.txt
#mpirun --oversubscribe -np 1024 ./main $gens $size $size 0 50 1 32 32 $method > $method/testout_1024_32_32.txt



mpirun --oversubscribe -np 1 Ex1 $gens $size $size 0 50 1  > $method/testout_serial.txt





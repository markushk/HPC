#!/bin/bash

size=64
gens=1
method=$1

mkdir ${method}_small
dir=${method}_small
echo $dir

mpirun --oversubscribe -np 2 ./main $gens $size $size 0 50 1 1 2 $method > $dir/testout_2_1_2.txt
mpirun --oversubscribe -np 2 ./main $gens $size $size 0 50 1 2 1 $method > $dir/testout_2_2_1.txt

mpirun --oversubscribe -np 4 ./main $gens $size $size 0 50 1 2 2 $method > $dir/testout_4_2_2.txt
mpirun --oversubscribe -np 4 ./main $gens $size $size 0 50 1 1 4 $method > $dir/testout_4_1_4.txt
mpirun --oversubscribe -np 4 ./main $gens $size $size 0 50 1 4 1 $method > $dir/testout_4_4_1.txt

#mpirun --oversubscribe -np 9 ./main $gens $size $size 0 50 1 1 9 $method > $dir/testout_9_1_1.txt
#mpirun --oversubscribe -np 9 ./main $gens $size $size 0 50 1 9 1 $method > $dir/testout_9_9_1.txt
#mpirun --oversubscribe -np 9 ./main $gens $size $size 0 50 1 3 3 $method > $dir/testout_9_3_3.txt

mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 1 16 $method > $dir/testout_16_1_16.txt
mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 2 8 $method > $dir/testout_16_2_8.txt
mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 4 4 $method > $dir/testout_16_4_4.txt
mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 16 1 $method > $dir/testout_16_16_1.txt
mpirun --oversubscribe -np 16 ./main $gens $size $size 0 50 1 8 2 $method > $dir/testout_16_8_2.txt


#mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 1 32 $method > $dir/testout_32_1_32.txt
#mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 2 16 $method > $dir/testout_32_2_16.txt
#mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 4 8 $method > $dir/testout_32_4_8.txt
#mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 8 4 $method > $dir/testout_32_8_4.txt
#mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 16 2 $method > $dir/testout_32_16_2.txt
#mpirun --oversubscribe -np 32 ./main $gens $size $size 0 50 1 32 1 $method > $dir/testout_32_32_1.txt

#mpirun --oversubscribe -np 256 ./main $gens $size $size 0 50 1 1 256 $method > $dir/testout_256_1_256.txt
#mpirun --oversubscribe -np 256 ./main $gens $size $size 0 50 1 1 256 $method > $dir/testout_256_256_1.txt

#mpirun --oversubscribe -np 256 ./main $gens $size $size 0 50 1 8 32 $method > $dir/testout_256_8_32.txt
#mpirun --oversubscribe -np 256 ./main $gens $size $size 0 50 1 32 8 $method > $dir/testout_256_32_8.txt


#mpirun --oversubscribe -np 512 ./main $gens $size $size 0 50 1 1 512 $method > $dir/testout_512_1_512.txt
#mpirun --oversubscribe -np 512 ./main $gens $size $size 0 50 1 512 1 $method > $dir/testout_512_512_1.txt

#mpirun --oversubscribe -np 512 ./main $gens $size $size 0 50 1 16 32 $method > $dir/testout_512_16_32.txt
#mpirun --oversubscribe -np 512 ./main $gens $size $size 0 50 1 32 16 $method > $dir/testout_512_32_16.txt

#mpirun --oversubscribe -np 1024 ./main $gens $size $size 0 50 1 1024 1 $method > $dir/testout_1024_1024_1.txt
#mpirun --oversubscribe -np 1024 ./main $gens $size $size 0 50 1 1 1024 $method > $dir/testout_1024_1_1024.txt
#mpirun --oversubscribe -np 1024 ./main $gens $size $size 0 50 1 32 32 $method > $dir/testout_1024_32_32.txt



mpirun --oversubscribe -np 1 Ex1 $gens $size $size 0 50 1  > $dir/testout_serial.txt





#!/bin/bash 

method=$1

for file in ${method}/testout_*
do
  correct=`diff -Naur $file ${method}/testout_serial.txt | grep -A 1 testout_`
  if [ ! -z "$correct" ]
  then
    echo "$file is wrong"
	fi 
	unset correct
done

#!/bin/bash
for i in 0 1 2 3
do
   mpicc "mpi"$i".c" -o "mpi"$i".out" -lm
done 
if [ -e output.txt ]
then rm output.txt
fi
if [ $# -eq 1 ]
then n=$1;
else
n=100000000;
fi
for i in 0 1 2 3
do
   echo "mpi"$i" result: " >> output.txt
   for q in 1 2 4 8 16 32
   do
      mpiexec -oversubscribe -np $q ./"mpi"$i".out" $n | tail -n 2 >> output.txt
   done
   echo "" >> output.txt
done
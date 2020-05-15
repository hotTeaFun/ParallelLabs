#!/bin/bash
for i in 0 1 2 3
do
   mpicc "mpi"$i".c" -o "mpi"$i".out" -lm
done 
if [ -e output.txt ]
then rm output.txt
fi
for i in 0 1 2 3
do
   echo "mpi"$i" result: \n" >> output.txt
   for n in 1 2 4 8 16
   do
      mpiexec -n $n ./"mpi"$i".out" 100000000 | tail -n 1 >> output.txt
   done
   echo "\n" >> output.txt
done
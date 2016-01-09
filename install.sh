#!/bin/bash

cd src

make clean && make 

cd ..
cd build
ndirec="$(pwd)"
cd ..
cd bin
bdirc="$(pwd)"
cd ..
cd src


for i in *; do
   if [ "${i}" != "${i%.o}" ] || [ "${i}" != "${i%.mod}" ];then
      mv "${i}" "$ndirec"
   fi
done


mv mcgrid "$bdirc" 

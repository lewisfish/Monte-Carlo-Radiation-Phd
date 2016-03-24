#!/bin/bash


cd src

make clean && make 

cd ..
if [ ! -d "build" ]; then
   mkdir "build"
fi
cd build
ndirec="$(pwd)"
cd ..
if [ ! -d "bin" ]; then
   mkdir "bin"
fi
cd bin
rm mcgrid
bdirc="$(pwd)"
cd ..
cd src


for i in *; do
   if [ "${i}" != "${i%.o}" ] || [ "${i}" != "${i%.mod}" ];then
      mv "${i}" "$ndirec"
   fi
done


mv mcgrid "$bdirc" && echo " " && echo "*****Install complete*****" && echo " "

#!/bin/bash

OUT=hardware_info.txt

echo "===== CPU =====" > $OUT
lscpu >> $OUT

echo -e "\n===== CPU FLAGS =====" >> $OUT
lscpu | grep Flags >> $OUT

echo -e "\n===== MEMORIA =====" >> $OUT
free -h >> $OUT

echo -e "\n===== CACHE =====" >> $OUT
lscpu | grep -E "L1d|L1i|L2|L3" >> $OUT

echo -e "\n===== NUMA =====" >> $OUT
lscpu | grep -E "NUMA|node" >> $OUT

echo -e "\n===== FRECUENCIAS =====" >> $OUT
cpupower frequency-info 2>/dev/null || echo "cpupower no disponible" >> $OUT

echo -e "\n===== KERNEL =====" >> $OUT
uname -a >> $OUT

echo -e "\n===== GCC / GFORTRAN =====" >> $OUT
gfortran --version >> $OUT
gcc --version >> $OUT

echo -e "\n===== OPENMP =====" >> $OUT
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS" >> $OUT

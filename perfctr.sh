#!/bin/bash

METRICA="L2CACHE FLOPS_DP"
TAMANHOS="32 64 128 256 512 1000 2000 4000 8000"

echo -n > op1_TIME_Resultados
echo -n > op1_L2CACHE_Resultados
echo -n > op1_FLOPS_DP_Resultados
echo -n > op1_FLOPS_AVX_Resultados

echo -n > op2_TIME_Resultados
echo -n > op2_L2CACHE_Resultados
echo -n > op2_FLOPS_DP_Resultados
echo -n > op2_FLOPS_AVX_Resultados

echo -n > _temp

make purge
make
for n in $TAMANHOS
do
        echo -n "${n} " >> op1_TIME_Resultados
        echo -n "${n} " >> op2_TIME_Resultados
        ./cgSolver -n ${n} -k 7 -p 1 -i 150 -o _temp
        cat _temp | grep "Tempo iter total:" | cut -d ":" -f 2 | tr -d " " >> op1_TIME_Resultados
        cat _temp | grep "Tempo residuo:" | cut -d ":" -f 2 | tr -d " " >> op2_TIME_Resultados
        echo >> op1_TIME_Resultados 
        echo >> op2_TIME_Resultados 
        rm _temp

        echo -n "${n} " >> op1_L2CACHE_Resultados
        echo -n "${n} " >> op2_L2CACHE_Resultados
        likwid-perfctr -C 3 -g L2CACHE -m ./cgSolver -n ${n} -k 7 -p 1 -i 150 -o /dev/null > _temp 
        cat _temp | grep "L2 miss ratio" | cut -d "|" -f 3 | tr -d ' ' | sed -n '1p' >> op1_L2CACHE_Resultados
        cat _temp | grep "L2 miss ratio" | cut -d "|" -f 3 | tr -d ' ' | sed -n '2p' >> op2_L2CACHE_Resultados
        echo >> op1_L2CACHE_Resultados
        echo >> op2_L2CACHE_Resultados

        echo -n "${n} " >> op1_FLOPS_DP_Resultados
        echo -n "${n} " >> op2_FLOPS_DP_Resultados
        likwid-perfctr -C 3 -g FLOPS_DP -m ./cgSolver -n ${n} -k 7 -p 1 -i 150 -o /dev/null > _temp 
        cat _temp | grep "  DP MFLOP/s" | cut -d "|" -f 3 | tr -d ' ' | sed -n '1p' >> op1_FLOPS_DP_Resultados
        cat _temp | grep "  DP MFLOP/s" | cut -d "|" -f 3 | tr -d ' ' | sed -n '2p' >> op2_FLOPS_DP_Resultados
        echo >> op1_FLOPS_DP_Resultados
        echo >> op2_FLOPS_DP_Resultados

        echo -n "${n} " >> op1_FLOPS_AVX_Resultados
        echo -n "${n} " >> op2_FLOPS_AVX_Resultados
        likwid-perfctr -C 3 -g FLOPS_AVX -m ./cgSolver -n ${n} -k 7 -p 1 -i 150 -o /dev/null > _temp 
        cat _temp | grep "Packed DP MFLOP/s" | cut -d "|" -f 3 | tr -d ' ' | sed -n '1p' >> op1_FLOPS_AVX_Resultados
        cat _temp | grep "Packed DP MFLOP/s" | cut -d "|" -f 3 | tr -d ' ' | sed -n '2p' >> op2_FLOPS_AVX_Resultados
        echo >> op1_FLOPS_AVX_Resultados
        echo >> op2_FLOPS_AVX_Resultados

done
make purge

rm _temp
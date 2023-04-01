#!/usr/bin/gnuplot

set terminal png size 800,600
set grid
set logscale x
set xtics ("32" 32, "64" 64, "128" 128, "256" 256, "512" 512, "1000" 1000, "2000" 2000, "4000" 4000, "8000" 8000)


set output "l2miss.png"
set title "L2 MISS"
set xlabel "N"
set ylabel "MISS RATIO"
plot "./velho/op1_L2CACHE_Resultados" title 'OP_1 TRAB_1' with lines linecolor "red", \
    "./velho/op2_L2CACHE_Resultados" title 'OP_2 TRAB_1' with lines linecolor "blue", \
    "./novo/op1_L2CACHE_Resultados" title 'OP_1 TRAB_2' with lines linecolor "green", \
    "./novo/op2_L2CACHE_Resultados" title 'OP_2 TRAB_2' with lines linecolor "magenta"


set output "time.png"
set title "TEMPO"
set xlabel "N"
set ylabel "Milissegundo"
plot "./velho/op1_TIME_Resultados" title 'OP_1 TRAB_1' with lines linecolor "red", \
    "./velho/op2_TIME_Resultados" title 'OP_2 TRAB_1' with lines linecolor "blue", \
    "./novo/op1_TIME_Resultados" title 'OP_1 TRAB_2' with lines linecolor "green", \
    "./novo/op2_TIME_Resultados" title 'OP_2 TRAB_2' with lines linecolor "magenta"


set output "scaled_time.png"
set title "TEMPO"
set xlabel "N"
set ylabel "Milissegundo" 
plot "./velho/op2_TIME_Resultados" title 'OP_2 TRAB_1' with lines linecolor "blue", \
    "./novo/op1_TIME_Resultados" title 'OP_1 TRAB_2' with lines linecolor "green", \
    "./novo/op2_TIME_Resultados" title 'OP_2 TRAB_2' with lines linecolor "magenta"

set output "avx.png"
set title "FLOPs AVX"
set xlabel "N"
set ylabel "MFLOP/s"
plot "./velho/op1_FLOPS_AVX_Resultados" title 'OP_1 TRAB_1' with lines linecolor "red", \
    "./velho/op2_FLOPS_AVX_Resultados" title 'OP_2 TRAB_1' with lines linecolor "blue", \
    "./novo/op1_FLOPS_AVX_Resultados" title 'OP_1 TRAB_2' with lines linecolor "green", \
    "./novo/op2_FLOPS_AVX_Resultados" title 'OP_2 TRAB_2' with lines linecolor "magenta"

set output "dp.png"
set title "FLOPs DP"
set xlabel "N"
set ylabel "MFLOP/s"
plot "./velho/op1_FLOPS_DP_Resultados" title 'OP_1 TRAB_1' with lines linecolor "red", \
    "./velho/op2_FLOPS_DP_Resultados" title 'OP_2 TRAB_1' with lines linecolor "blue", \
    "./novo/op1_FLOPS_DP_Resultados" title 'OP_1 TRAB_2' with lines linecolor "green", \
    "./novo/op2_FLOPS_DP_Resultados" title 'OP_2 TRAB_2' with lines linecolor "magenta"






set output "unroll_l2miss.png"
set title "L2 MISS"
set xlabel "N"
set ylabel "MISS RATIO"
plot "./velho/op1_L2CACHE_Resultados" title 'OP_1 TRAB_1' with lines linecolor "red", \
    "./velho/op2_L2CACHE_Resultados" title 'OP_2 TRAB_1' with lines linecolor "blue", \
    "./novo/UNROLL_op1_L2CACHE_Resultados" title 'OP_1 TRAB_2' with lines linecolor "green", \
    "./novo/UNROLL_op2_L2CACHE_Resultados" title 'OP_2 TRAB_2' with lines linecolor "magenta"


set output "unroll_time.png"
set title "TEMPO"
set xlabel "N"
set ylabel "Milissegundo"
plot "./velho/op1_TIME_Resultados" title 'OP_1 TRAB_1' with lines linecolor "red", \
    "./velho/op2_TIME_Resultados" title 'OP_2 TRAB_1' with lines linecolor "blue", \
    "./novo/UNROLL_op1_TIME_Resultados" title 'OP_1 TRAB_2' with lines linecolor "green", \
    "./novo/UNROLL_op2_TIME_Resultados" title 'OP_2 TRAB_2' with lines linecolor "magenta"

set output "unroll_avx.png"
set title "FLOPs AVX"
set xlabel "N"
set ylabel "MFLOP/s"
plot "./velho/op1_FLOPS_AVX_Resultados" title 'OP_1 TRAB_1' with lines linecolor "red", \
    "./velho/op2_FLOPS_AVX_Resultados" title 'OP_2 TRAB_1' with lines linecolor "blue", \
    "./novo/UNROLL_op1_FLOPS_AVX_Resultados" title 'OP_1 TRAB_2' with lines linecolor "green", \
    "./novo/UNROLL_op2_FLOPS_AVX_Resultados" title 'OP_2 TRAB_2' with lines linecolor "magenta"

set output "unroll_dp.png"
set title "FLOPs DP"
set xlabel "N"
set ylabel "MFLOP/s"
plot "./velho/op1_FLOPS_DP_Resultados" title 'OP_1 TRAB_1' with lines linecolor "red", \
    "./velho/op2_FLOPS_DP_Resultados" title 'OP_2 TRAB_1' with lines linecolor "blue", \
    "./novo/UNROLL_op1_FLOPS_DP_Resultados" title 'OP_1 TRAB_2' with lines linecolor "green", \
    "./novo/UNROLL_op2_FLOPS_DP_Resultados" title 'OP_2 TRAB_2' with lines linecolor "magenta"
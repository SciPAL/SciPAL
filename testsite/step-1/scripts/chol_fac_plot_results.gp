set term postscript landscape enhanced color solid linewidth 2.0 "Helvetica" 20
set xlabel "n_{rows}"
set ylabel "execution time"
set logscale xy
set grid 
set output "./chol_fac_times_double.ps"
set key inside left Left
plot "./chol_fac_times_16_double.dat" using 1:2 title "CPU" w lp, \
    "" using 1:($3) title "GPU netto" w lp, "" using 1:4 title "GPU brutto" w lp

set ylabel "CPU time/GPU time"
set output "./chol_fac_speedup.ps"
plot "./chol_fac_times_16_double.dat" using 1:($2/$3) title "netto" w lp, \
    "" using 1:($2/$4) title "brutto" w lp

set output "./chol_fac_gpu_individual_components.ps"
set ylabel "execution time"
plot "./chol_fac_times_16_double.dat" using 1:(abs($4-$3)) title "memory transfer" w lp, \
    "" using 1:6 title "factorize diagonal block" w lp, "" using 1:8 title "strip update" w lp, \
    "" using 1:($5) title "diag update" w lp, ""using 1:7 title "lo update" w lp
set term x11

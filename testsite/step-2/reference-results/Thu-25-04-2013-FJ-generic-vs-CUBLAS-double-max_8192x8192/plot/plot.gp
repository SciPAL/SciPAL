set term postscript landscape enhanced color solid  linewidth 2.0 "Helvetica" 20
set xlabel "matrix entries"
set ylabel "execution time"
set logscale xy
set grid
set output "runtimes.ps"
set key inside left Left box lw 0.5
plot "../MVTest_results.out" using ($1*$2):3 title "Fujimoto 1 double" w p,"../MVTest_results.out" using ($1*$2):4 title "CUBLAS double" w p
set ylabel "speedup"
set output "speedup.ps"
set key inside left Left
unset logscale y
plot"../MVTest_results.out" using ($1*$2):($3 / $4) title "Fujimoto vs CUBLAS (double)" w p
!ps2pdf runtimes.ps runtimes.pdf
!rm runtimes.ps
!ps2pdf speedup.ps speedup.pdf
!rm speedup.ps
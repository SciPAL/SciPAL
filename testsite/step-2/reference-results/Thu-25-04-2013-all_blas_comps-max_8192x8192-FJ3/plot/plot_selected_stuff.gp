set term postscript landscape enhanced color solid  linewidth 2.0 "Helvetica" 20
set xlabel "matrix entries"
set ylabel "execution time"
set logscale xy
set grid
set yrange[*:1]
set output "runtimes.ps"
set key inside left Left samplen 0 spacing 1.5 box lw 0.5 
plot "../MVTest_results.out" using ($1*$2):3 title "Fujimoto 3 float" w p, \
"../MVTest_results.out" using ($1*$2):4 title " CUBLAS float" w p, \
"../MVTest_results.out" using ($1*$2):5 title " ATLAS float" w p, \
"../MVTest_results.out" using ($1*$2):6 title " Fujimoto 3 double" w p, \
"../MVTest_results.out" using ($1*$2):7 title " CUBLAS double" w p, \
"../MVTest_results.out" using ($1*$2):8 title " ATLAS double" w p


set yrange[*:.005]
set output "runtimes-FJ-orig-vs-generic.ps"
#set key inside left Left samplen 0 spacing 1.5 box lw 0.5 
plot "../../Thu-25-04-2013-FJ-orig-vs-CUBLAS-float-max_8192x8192/MVTest_results.out" using ($1*$2):3 title " original Fujimoto" w p, \
"../MVTest_results.out" using ($1*$2):3 title " templated, no bitshift" w p
!ps2pdf runtimes-FJ-orig-vs-generic.ps ../../../doc/images/runtimes-FJ-orig-vs-generic.pdf
!rm runtimes-FJ-orig-vs-generic.ps


set output "runtimes-FJ-generic-vs-optimized-double.ps"
#set key inside left Left samplen 0 spacing 1.5 box lw 0.5 
plot "../../Thu-25-04-2013-FJ-generic-vs-CUBLAS-double-max_8192x8192/MVTest_results.out" using ($1*$2):3 title " generic Fujimoto (double)" w p, \
"../MVTest_results.out" using ($1*$2):6 title " optimized (double)" w p
!ps2pdf runtimes-FJ-generic-vs-optimized-double.ps ../../../doc/images/runtimes-FJ-generic-vs-optimized-double.pdf
!rm runtimes-FJ-generic-vs-optimized-double.ps



# rows | columns | Fujimoto 0 float | CUBLAS float | Fujimoto 3 float |Fujimoto 3 double | Fujimoto 1 double | CUBLAS double
#  1      2                3                 4            5                 6                 7                 8
set ylabel "speedup"
set output "speedup-FJ-bitshift.ps"
set key inside left Left
unset logscale y
set yrange[*:*]
plot "../MVTest_results-Fujimoto-comparisons.out" using ($1*$2):($5 / $3) title " Fujimoto original vs. no-bitshift (float)" w p



set output "speedup-FJ-double-opt.ps"
plot "../MVTest_results-Fujimoto-comparisons.out" using ($1*$2):($7 / $6) title " Fujimoto optimized vs generic (double)" w p

!ps2pdf speedup-FJ-bitshift.ps ../../../doc/images/speedup-FJ-bitshift.pdf
!rm speedup-FJ-bitshift.ps
!ps2pdf speedup-FJ-double-opt.ps ../../../doc/images/speedup-FJ-double-opt.pdf
!rm speedup-FJ-double-opt.ps





gbyte = 1024*1024*1024
set ylabel "achieved memory throughput [GByte/s]"
set output "bandwidths.ps"
plot "../MVTest_results-Fujimoto-comparisons.out" using ($1*$2):(4*$1*$2/(gbyte*$3)) title " Fujimoto original (float)" w p, \
""  using ($1*$2):(8*$1*$2/(gbyte*$6)) title " Fujimoto optimized (double)" w p, \
"" using ($1*$2):(4*$1*$2/(gbyte*$4)) title " CUBLAS (float)" w p, \
"" using ($1*$2):(8*$1*$2/(gbyte*$8)) title " CUBLAS (double)" w p



!ps2pdf bandwidths.ps ../../../doc/images/bandwidths.pdf
!rm bandwidths.ps



return




set ylabel "speedup"
set output "speedup.ps"
set key inside left Left
unset logscale y
set yrange[*:30]
plot"../MVTest_results.out" using ($1*$2):($5 / $4) title " ATLAS vs CUBLAS (float)" w p, \
"../MVTest_results.out" using ($1*$2):($8 / $7) title " ATLAS vs CUBLAS (double)" w p, \
"../MVTest_results.out" using ($1*$2):($5 / $3) title " ATLAS vs Fujimoto (float)" w p, \
"../MVTest_results.out" using ($1*$2):($8 / $6) title " ATLAS vs Fujimoto (double)" w p
!ps2pdf runtimes.ps ../../../doc/images/runtimes.pdf
!rm runtimes.ps
!ps2pdf speedup.ps ../../../doc/images/speedup.pdf
!rm speedup.ps

set output "speedup-logplot.ps"
set key inside right bottom Left samplen 0 spacing 1.5 box lw 0.5 
set logscale y
set yrange[*:100]
plot"../MVTest_results.out" using ($1*$2):($5 / $4) title " ATLAS vs CUBLAS (float)" w p, \
"../MVTest_results.out" using ($1*$2):($8 / $7) title " ATLAS vs CUBLAS (double)" w p, \
"../MVTest_results.out" using ($1*$2):($5 / $3) title " ATLAS vs Fujimoto (float)" w p, \
"../MVTest_results.out" using ($1*$2):($8 / $6) title " ATLAS vs Fujimoto (double)" w p
#, "../MVTest_results.out" using ($1*$2):($3 / $4) title "Fujimoto vs CUBLAS (float)" w p, "../MVTest_results.out" using ($1*$2):($6 / $7) title "Fujimoto vs CUBLAS (double)" w p

!ps2pdf speedup-logplot.ps ../../../doc/images/speedup-logplot.pdf
!rm speedup-logplot.ps


set output "speedup-FJ-vs-CUBLAS.ps"
set key inside left top Left samplen 0 spacing 1.5 box lw 0.5 
#set key inside left Left
unset logscale y
set yrange[*:6]
plot "../MVTest_results.out" using ($1*$2):($3 / $4) title " Fujimoto vs CUBLAS (float)" w p, \
"../MVTest_results.out" using ($1*$2):($6 / $7) title " Fujimoto vs CUBLAS (double)" w p

!ps2pdf speedup-FJ-vs-CUBLAS.ps ../../../doc/images/speedup-FJ-vs-CUBLAS.pdf
!rm speedup-FJ-vs-CUBLAS.ps


set output "speedup-logplot-FJ-vs-CUBLAS.ps"
#set key inside left Left
set logscale y
set yrange[*:10]
plot "../MVTest_results.out" using ($1*$2):($3 / $4) title " Fujimoto vs CUBLAS (float)" w p, \
"../MVTest_results.out" using ($1*$2):($6 / $7) title " Fujimoto vs CUBLAS (double)" w p

!ps2pdf speedup-logplot-FJ-vs-CUBLAS.ps ../../../doc/images/speedup-logplot-FJ-vs-CUBLAS.pdf
!rm speedup-logplot-FJ-vs-CUBLAS.ps

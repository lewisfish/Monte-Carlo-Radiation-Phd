set terminal pngcairo size 700,700 enhanced

set output 'validation.png'
set multiplot layout 2,1 title "{/:Bold=18 Validation against S.Jacques Fluro model (1993)}"

set fit errorvariables
set fit quiet
fit(x) = c1*exp(-k1*x/.261)-c2*exp(-k2*x/.261)

c1 = 6.27
c2 = 1.18
k1 = 1.0
k2 = 14.4

set xlabel "Depth/cm" font ",9"
set ylabel "Fraction of escaped fluro per depth" font ",9"
set title "Depth resolved Fluence" font ",14"

fit[0:200] fit(x) 'jmean430.dat' via c1,c2,k1,k2
set print 'fit_params.txt'
print c1, c2, k1, k2
set print

ti = sprintf("%.3f*exp(-%.3f*x/.261)-%.3f*exp(-%.3f*x/.261)",c1,k1,c2,k2)
plot fit(x) t ti, 'jmean430.dat' w p pointtype 2 t "Monte Carlo data", (6.27*exp(-x/.261))-(1.18*exp((-14.4*x)/.261)) t "S.Jacques model"
#plot fit(x) t ti, 'jmean430.dat' w p pointtype 2 t "Monte Carlo data", (5.76*exp(-x/.047))-(1.31*exp((-10.2*x)/.047)) t "S.Jacques model"
set title 'Relative error'
set ylabel 'relative error'
unset key
plot 'error.dat' w l

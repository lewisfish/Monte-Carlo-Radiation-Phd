set terminal pngcairo size 700,700
set output 'validation.png'
set multiplot layout 2, 1 title "Validation against S.Jacques Fluro model (1993)" font ",14"

set fit errorvariables
set fit quiet
fit(x) = a*exp(-b*x/.29)


a = .771
b = .994

set xlabel "Depth/cm"
set ylabel "Fraction of escaped fluro per depth"
set title "Fluro escape function"

fit[0:1] fit(x) 'testff.dat' via a,b
set print "fit_params.txt"
print a, b
set print
set xrange[0:1]
ti = sprintf("%.3f*exp(-%.3f*x/.29)",a,b)
plot fit(x) t ti, 'testff.dat' w p t "Monte Carlo data", .771*exp(-.994*x/.29) t "S.Jacques model"
set title 'Relative Error'
set ylabel 'Relative error'
plot 'error.dat' w l
unset multiplot

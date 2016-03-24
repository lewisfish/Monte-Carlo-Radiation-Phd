set terminal pngcairo size 700,700

set xrange[500:800]

set output 'fluro.png'
plot 'resfluro8.dat' w l,'resfluro5.dat' w l,'resfluro2.dat' w l

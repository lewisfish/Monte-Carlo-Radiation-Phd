set terminal pngcairo enhanced size 1366,768
cd 'data'
set output 'image.png'

unset key
unset border
unset xtics
unset ytics
unset colorbox
set size square

plot 'image.dat' matrix w image

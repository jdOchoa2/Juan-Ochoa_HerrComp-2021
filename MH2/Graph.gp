set xlabel "size"; set ylabel "time[s]"
set title "Strong scaling"
set term png; set out "scaling.png"
set log x; set log y
set key right top
plot "data.txt" u 1:2 w p

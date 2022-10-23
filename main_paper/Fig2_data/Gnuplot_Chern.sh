set term epslatex color colortext size 14cm,12cm standalone
#set term postscript color size 14,12
set output "Chern.tex"

s = 1.7
w = 2
key_sa = 1
key_sp = 0.85

set style line 1 pt 64 lc rgb "blue" lw w ps s
set style line 2 pt 65 lc rgb "red" lw w ps s
set style line 3 pt 66 lc rgb "dark-violet" lw w ps s
set style line 4 pt 11 lc rgb "grey0" lw w ps 2.2
set style line 5 pt 68 lc rgb "dark-green" lw w ps s
#set style line 6 pt 69 lc rgb "dark-orange" lw w ps s

set style line 10 lc 0

set multiplot layout 2,3 rowsfirst

#1
set ylabel "$\\Delta Q_{s}$" offset 3
set xlabel "$\\theta_{F} /\\pi $" offset 0,1
set xr [0:2]
set yr [-0.5:2]
set xtics 0,0.5 offset 0,0.3
set ytics ("$-0.5$" -0.5, "$0$" 0, "$0.5$" 0.5, "$1$" 1, "$1.5$" 1.5, "$2$" 2) offset 0.4
set key width -4 samplen key_sa spacing key_sp at 1.45,1.93
set border lw 4
set lmargin at screen 0.08
set rmargin at screen 0.49
set tmargin at screen 0.96
set bmargin at screen 0.58
set label "(a)" at 0.1,2.1
set arrow 11 from 0,0 to 2,2 nohead dt 3 lc 0 lw 5
set arrow 9 from 0,1 to 2,1 nohead dt 3 lc rgb "gray0" lw 2
set arrow 10 from 0,0 to 2,0 nohead ls 10 dt 3 lw 2 

plot "chern_J2_0.05_Jchi_0.05_Ly_6.dat" index 0 using ($1/pi/2.):3 w lp title "$J_{2}=0.05, J_{\\chi }=0.05$" ls 4,\
"chern_J2_0.1_Jchi_0.05_Ly_6.dat" index 0 using ($1/pi/2.):3 w lp title "$J_{2}=0.1, J_{\\chi }=0.05$" ls 2,\
"chern_J2_0.1_Jchi_0.1_Ly_6.dat" index 0 using ($1/pi/2.):3 w lp title "$J_{2}=0.1, J_{\\chi }=0.1$" ls 1


unset arrow
unset label 
unset xlabel
unset ylabel
unset xtics
unset ytics
set key width -3 samplen key_sa spacing key_sp at 1.43,-0.1
plot "chern_J2_0.01_Jchi_0.05_Ly_6.dat" index 0 using ($1/pi/2.):3 w lp title "$J_{2}=0.01, J_{\\chi }=0.05$" ls 3,\
"chern_J2_0.2_Jchi_0.2_Ly_6.dat" index 0 using ($1/pi/2.):3 w lp title "$J_{2}=0.2, J_{\\chi }=0.2$" ls 5



#2
set format x '%g'
set format y '%g'
set ylabel "$\\Delta Q_{s}$" offset 3.2
set xlabel "$\\theta_{F} /\\pi $" offset 0,1
set xr [0:2]
set yr [0:2]
set xtics 0,0.5 offset 0,0.3
set ytics 0,0.5 offset 0.4
set key width 1 at 0.95,1.8
set key samplen key_sa
set key spacing key_sp
set border lw 4
set lmargin at screen 0.59
set rmargin at screen 0.99
set tmargin at screen 0.96
set bmargin at screen 0.58
set label "(b)" at 0.1,2.1
set label "$J_{2}=0.05, J_{\\chi }=0.05$" at 0.75,0.15
set arrow 10 from 0,1 to 1,1 nohead dt 3 lw 2 lc 0
set arrow 8 from 0,0 to 2,2 nohead dt 3 lw 5 lc rgb "gray0"

plot "chern_J2_0.05_Jchi_0.05_Ly_4.dat" index 0 using ($1/pi/2.):3 w lp title "$L_{y}=4$" ls 1,\
"chern_J2_0.05_Jchi_0.05_Ly_5.dat" index 0 using ($1/pi/2.):3 w lp title "$L_{y}=5$" ls 2,\
"chern_J2_0.05_Jchi_0.05_Ly_6.dat" index 0 using ($1/pi/2.):3 w lp title "$L_{y}=6$" ls 4



#3
unset arrow
unset label
set ylabel "$C$" offset 2.5
set xlabel "$J_{2}$" offset 0,0.4
set xr [0:0.2]
set xtics ("$0$" 0, "$0.05$" 0.05, "$0.1$" 0.1, "$0.15$" 0.15, "$0.2$" 0.2)
set yr [-0.02:2.1]
set ytics 0,0.5 offset 0.4
set key samplen key_sa
set key spacing key_sp
set key width 0 at 0.12,0.2
set border lw 4
set lmargin at screen 0.08
set rmargin at screen 0.49
set tmargin at screen 0.47
set bmargin at screen 0.09
set label "(c)" at 0.01,2.2
set style arrow 4 head noborder size screen 0.03,15,45 dt "-" lc rgb "red" lw 2
set arrow 1 from 0.0175,-0.02 to 0.0175,1.5 as 4 lc rgb "red"
set arrow 2 from 0.1375,-0.02 to 0.1375,1 as 4 lc rgb "red"
set label "$J_{2c}^{(2)}$" at 0.125,-0.33 textcolor rgb "red"
set label "$J_{2c}^{(1)}$" at 0.005,-0.33 textcolor rgb "red"

plot "chern_num_Jchi_0.05_Ly_6_J2.dat" index 0 using 1:2 w lp title "$J_{\\chi }=0.05$" pt 7 lc rgb "grey10" lw 3 ps s


#4
reset
set style line 1 pt 64 lc rgb "blue" lw 3 ps 1.7
set style line 4 pt 11 lc rgb "grey0" lw 3 ps 2.5
set style line 6 pt 69 lc rgb "dark-orange" lw 3 ps 3

set ylabel "$C$" offset 3.2
set xlabel "$\\delta $" offset 0,0.4
set xr [0:0.125]
set yr [0:2.2]
set xtics ("$0$" 0, "$\\frac{1}{36} $" 1/36., "$\\frac{1}{24}$" 1/24., "$\\frac{1}{12}$" 1/12., "$\\frac{1}{8}$" 1/8.) 
set ytics 0,0.5 offset 0.4
set key at 0.12,0.6
set key samplen 1.0
set key spacing 0.8
set border lw 4
set lmargin at screen 0.59
set rmargin at screen 0.99
set tmargin at screen 0.47
set bmargin at screen 0.09
set label "(d)" at 0.00625,2.3

plot "chern_num_Ly_6_doping.dat" index 0 using 1:4 w lp title "$J_{2}=0.1, J_{\\chi }=0.1$" ls 1,\
"chern_num_Ly_6_doping.dat" index 0 using 1:3 w lp title "$J_{2}=0.05, J_{\\chi }=0.05$" ls 4,\
"chern_num_Ly_6_doping.dat" index 0 using 1:2 w lp title "$J_{2}=0.04, J_{\\chi }=0.05$" ls 6

unset multiplot




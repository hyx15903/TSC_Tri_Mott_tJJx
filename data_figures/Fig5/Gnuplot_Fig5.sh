set term epslatex color colortext size 14cm,12cm standalone
#set term postscript color size 10,5
set output "Fig5.tex"
#Fig4

s = 1.7
w = 1.2

set style line 1 pt 5 lc rgb "dark-blue" lw w ps s
set style line 2 pt 7 lc rgb "dark-red" lw w ps 2
set style line 3 pt 13 lc rgb "web-green" lw w ps 1.5
set style line 4 pt 4 lc rgb "dark-blue" lw w ps s
set style line 5 pt 6 lc rgb "dark-red" lw w ps 2
set style line 6 pt 12 lc rgb "web-green" lw w ps 1.5
set style line 7 pt 15 lc rgb "grey20" lw w ps s
set style line 8 pt 14 lc rgb "grey20" lw w ps s

set multiplot layout 2,4 rowsfirst

#1
set ylabel "$\\theta _{\\alpha \\beta }(r)$" offset 5.4,-1.1
set xlabel "$r$" offset 0,1
set xr [2:25]
set xtics 0,5 offset 0,0.4
set yr [-pi-0.2:2*pi]
set ytics offset 0.55
set ytics ("$-\\pi $" -pi, "$0$" 0, "$-\\frac{2}{3}\\pi $" -2*pi/3, "$\\frac{2}{3}\\pi $" 2*pi/3, "$\\pi $" pi, "$2\\pi $" 2*pi)

set key samplen 0.8
set key spacing 1
#set grid
set border lw 6
#set tics nomirror
set lmargin at screen 0.08
set rmargin at screen 0.45
set tmargin at screen 0.96
set bmargin at screen 0.57
set label "$(a) N=36\\times 6$" at 1.6,6.67

set arrow 4 from 2,2*pi/3 to 25,2*pi/3 nohead dt 3 lc rgb "black" lw 2
set arrow 6 from 2,-2*pi/3 to 25,-2*pi/3 nohead dt 3 lc 0 lw 2
#illustration
rx = 14
ry = 0
set arrow 1 from rx,ry to rx+3.76,ry head lc 0 lw 2
set arrow 2 from rx,ry to rx-2,ry-1.717 head lc 0 lw 2
set arrow 3 from rx,ry to rx-2,ry+1.717 head lc 0 lw 2
#arrow 1
#ar1x = 0
#ar1y = 0.31
#set arrow 10 from rx+ar1x,ry-1+ar1y to rx+0.8+ar1x,ry-0.8+ar1y nohead lc rgb "black" lw 2
#set arrow 11 from rx+ar1x,ry-0.6+ar1y to rx+0.8+ar1x,ry-0.8+ar1y nohead lc rgb "black" lw 2
#arrow 2
#ar2x = 0
#ar2y = 1.29
#set arrow 12 from rx+ar2x,ry-1+ar2y to rx+0.8+ar2x,ry-0.8+ar2y nohead lc rgb "black" lw 2
#set arrow 13 from rx+ar2x,ry-0.6+ar2y to rx+0.8+ar2x,ry-0.8+ar2y nohead lc rgb "black" lw 2
#arrow 3
#ar3x = -2.1
#ar3y = 1.05
#set arrow 14 from rx+1+ar3x,ry-1.2+ar3y to rx+0.8+ar3x,ry-0.8+ar3y nohead lc rgb "black" lw 2
#set arrow 15 from rx+ar3x,ry-1.0+ar3y to rx+0.8+ar3x,ry-0.8+ar3y nohead lc rgb "black" lw 2
set label "$\\Delta _{a} $" at rx+3.6,ry+0.4
set label "$\\Delta _{c} $" at rx-1.2,ry-1.3
set label "$\\Delta _{b} $" at rx-4.4,ry+1.4
set label "$\\theta _{ba} $" at rx+1,ry+0.8
set label "$\\theta _{ac} $" at rx+1,ry-0.8
set label "$\\theta _{cb} $" at rx-3.6,ry
set object circle at rx,ry size scr 0.02 fs empty border lc rgb "grey0"


#set title "doping=1/12, 24*6, t=3, $J_{2}=0.1, J_{\\chi}=0.1$, aa-bond"
set key vertical maxrows 1
set key width -3.1
set key at 25,6.0
set label "$J_{2}=0.05$" at 5,5.6
plot "Fig5a.dat" index 0 using 1:2 w lp title "$\\theta _{ba }$" ls 1,\
"Fig5a.dat" index 0 using 1:3 w lp title "$\\theta _{bc }$" ls 4

unset ylabel
unset label
unset xlabel
unset xtics
unset ytics
unset arrow
unset object
set key at 25,5.2
set label "$J_{2}=0.1$" at 5,4.8
plot "Fig5a.dat" index 1 using 1:2 w lp title "$\\theta _{ba }$" ls 2,\
"Fig5a.dat" index 1 using 1:3 w lp title "$\\theta _{bc }$" ls 5

unset label
set key at 25,4.4
set label "$J_{2}=0.15$" at 5,4.0
plot "Fig5a.dat" index 2 using 1:2 w lp title "$\\theta _{ba }$" ls 3,\
"Fig5a.dat" index 2 using 1:3 w lp title "$\\theta _{bc }$" ls 6


#2
unset label
set ylabel "$\\theta _{\\alpha \\beta }(r)$" offset 5.4,-1.1
set xlabel "$r$" offset 0,1
set xr [2:25]
set xtics 0,5 offset 0,0.4
set yr [-pi-0.2:2*pi]
set ytics offset 0.55
set ytics ("$-\\pi $" -pi, "$0$" 0, "$-\\frac{2}{3}\\pi $" -2*pi/3, "$\\frac{2}{3}\\pi $" 2*pi/3, "$\\pi $" pi, "$2\\pi $" 2*pi)

set key samplen 0.8
set key spacing 1
#set grid
set border lw 6
#set tics nomirror
set lmargin at screen 0.58
set rmargin at screen 0.94
set tmargin at screen 0.96
set bmargin at screen 0.57
set label "$(b) N=48\\times 4$" at 1.6,6.67

set arrow 4 from 2,2*pi/3 to 25,2*pi/3 nohead dt 3 lc rgb "black" lw 2
set arrow 6 from 2,-2*pi/3 to 25,-2*pi/3 nohead dt 3 lc 0 lw 2


#set title "doping=1/12, 24*6, t=3, $J_{2}=0.1, J_{\\chi}=0.1$, aa-bond"
set key vertical maxrows 1
set key width -3.1
set key at 25,6.0
set label "$J_{2}=0.05$" at 5,5.6
plot "Fig5b.dat" index 0 using 1:2 w lp title "$\\theta _{ba }$" ls 1,\
"Fig5b.dat" index 0 using 1:3 w lp title "$\\theta _{bc }$" ls 4

unset ylabel
unset label
unset xlabel
unset xtics
unset ytics
unset arrow
unset object
set key at 25,5.2
set label "$J_{2}=0.1$" at 5,4.8
plot "Fig5b.dat" index 1 using 1:2 w lp title "$\\theta _{ba }$" ls 2,\
"Fig5b.dat" index 1 using 1:3 w lp title "$\\theta _{bc }$" ls 5

unset label
set key at 25,4.4
set label "$J_{2}=0.15$" at 5,4.0
plot "Fig5b.dat" index 2 using 1:2 w lp title "$\\theta _{ba }$" ls 3,\
"Fig5b.dat" index 2 using 1:3 w lp title "$\\theta _{bc }$" ls 6

unset label
set key at 25,1.4
set label "$J_{2}=0.01$" at 5,1
plot "Fig5b.dat" index 3 using 1:2 w lp title "$\\theta _{ba }$" ls 7,\
"Fig5b.dat" index 3 using 1:3 w lp title "$\\theta _{bc }$" ls 8



#3
reset
set style line 1 pt 64 lc rgb "blue" lw w ps 2
set style line 2 pt 65 lc rgb "dark-red" lw w ps 1.4
set style line 11 pt 5 lc rgb "dark-blue" lw w ps 1.8
set style line 12 pt 7 lc rgb "green" lw w ps 1.1
set ylabel "$\\left | P_{\\alpha \\beta }(r)/ P_{bb }(r)\\right |$" offset 2
set xlabel "$r$" offset 0,1
set xr [2:25]
set yr [0:1.6]
unset logscale y
unset logscale x
set format y "$%g$"
#set xtics ("$1$" 1, "" 2, "" 3, "" 4, "" 5, "" 6, "" 7, "" 8, "" 9, "$10$" 10, "$20$" 20)
set ytics offset 0.5
set xtics 0,5 offset 0,0.3
set key at 23,1.35
set key samplen 0.8
set key spacing 0.4
set key vertical maxrows 1
set key width -3.5
#set grid
set border lw 4
#set tics nomirror
set lmargin at screen 0.08
set rmargin at screen 0.45
set tmargin at screen 0.46
set bmargin at screen 0.07
set label "$(c) N=36\\times 6$" at 2,1.68

#a(x) = m0+m1*x
#fit a(x) "Sc_cor1_J2_0.1_Jchi_0.05_Lx_36_bb_bond_infinite.dat" index 0 u (log($1)):(log($2)) every 4::5::17 via m0,m1

#b(x) = m2+m3*x
#fit b(x) "Sc_cor1_J2_0.1_Jchi_0.05_Lx_36_ba_bond_infinite.dat" index 0 u (log($1)):(log($2)) every 4::5::17 via m2,m3

set label "$J_{2}=0.1$" at 5.5,0.22
set label "$J_{2}=0.05$" at 5.5,1.45

set arrow 7 from 2,1.02 to 25,1.02 nohead dt 3 lc 0 lw 2
set arrow 8 from 2,0.4633 to 25,0.4633 nohead dt 3 lc 0 lw 2

#plot "merge_J2_0.05_Jchi_0.05_Lx_36.dat" index 0 using 1:($6/$2) every ::1 w lp title "$|P_{ba}/P_{bb}|$" ls 1,\
#"merge_J2_0.05_Jchi_0.05_Lx_36.dat" index 0 using 1:($10/$2) every ::1 w lp title "$|P_{bc}/P_{bb}|$" ls 2

plot "Fig5c.dat" index 0 using 1:2 every ::1 w lp title "$|P_{ba}/P_{bb}|$" ls 1,\
"Fig5c.dat" index 0 using 1:3 every ::1 w lp title "$|P_{bc}/P_{bb}|$" ls 2

set key at 23,0.12
unset label
unset xlabel
unset ylabel
unset xtics
unset ytics
unset arrow

#plot "merge_J2_0.1_Jchi_0.05_Lx_36.dat" index 0 using 1:($6/$2) every ::1 w lp title "$|P_{ba}/P_{bb}|$" ls 11,\
#"merge_J2_0.1_Jchi_0.05_Lx_36.dat" index 0 using 1:($10/$2) every ::1 w lp title "$|P_{bc}/P_{bb}|$" ls 12

plot "Fig5c.dat" index 1 using 1:2 every ::1 w lp title "$|P_{ba}/P_{bb}|$" ls 11,\
"Fig5c.dat" index 1 using 1:3 every ::1 w lp title "$|P_{bc}/P_{bb}|$" ls 12

#set key at 20.5,0.0
#plot "jx0.05j20.15/Sc_cor_ph.dat" index 0 using 4:($6/$5) every ::1 w lp title "$|P_{ba}/P_{bb}|$" ls 1,\
#"jx0.05j20.15/Sc_cor_ph.dat" index 0 using 4:($7/$5) every ::1 w lp title "$|P_{bc}/P_{bb}|$" ls 2


#4
reset
w = 0.2
set style line 1 pt 64 lc rgb "dark-blue" lw w ps 2.2
set style line 2 pt 7 lc rgb "web-green" lw w ps 1.5
set style line 3 pt 8 lc rgb "grey0" lw w ps 2.2
set style line 4 pt 11 lc rgb "red" lw w ps 1.7
set style arrow 3 head filled size screen 0.03,15,45 ls 3

set ylabel "$|\\overline{P_{\\alpha \\beta}/P_{bb}}|$" offset 2.2,0
set y2label "$\\overline{\\theta _{\\alpha \\beta}}$" offset -5.3,-0.6 rotate by -90

set xlabel "$J_{2}$" offset 0,1
set xr [0:0.2]
set yr [0:1.25]
set y2r [2*pi/3-0.15:3.2]


#set xtics ("$1$" 1, "" 2, "" 3, "" 4, "" 5, "" 6, "" 7, "" 8, "" 9, "$10$" 10, "$20$" 20)
set y2tics ("$\\frac{2}{3}\\pi $" 2*pi/3, "$\\frac{5}{6}\\pi $" 5*pi/6, "$\\pi $" pi) offset -0.8
set ytics 0,0.2 offset 0.8
set ytics nomirror
set xtics 0,0.05 offset 0,0.4
#set key at 19.5,0.8
set key at 0.065,0.85
set key samplen 0.8
set key spacing 1.9
#set key vertical maxrows 1
#set key width -1.5
#set grid
set border lw 6
#set tics nomirror
set lmargin at screen 0.58
set rmargin at screen 0.94
set tmargin at screen 0.46
set bmargin at screen 0.07
set label "$(d) N=36\\times 6$" at 0.007,1.31
set arrow from 0.03,1.1 to 0.01,1.1 as 3
set arrow from 0.17,1.1 to 0.19,1.1 as 3

set style arrow 4 head noborder size screen 0.03,15,45 dt "-" lc rgb "red" lw 2
set arrow 2 from 0.08,-0 to 0.08,0.4 as 4 lc rgb "red"
set label "$J_{2c}'$" at 0.07,-0.1 textcolor rgb "red"

#plot "SC_corr_ave_jx_0.05_different_j2.dat" index 0 using 1:2:7 w yerr title "$|\\frac{P_{ba}}{P_{bb}}|$" ls 3 axis x1y1,\
#"SC_corr_ave_jx_0.05_different_j2.dat" index 0 using 1:3:8 w yerr title "$|\\frac{P_{bc}}{P_{bb}}|$" ls 4 axis x1y1,\

plot "Fig5d.dat" index 0 using 1:2 w lp title "$|\\overline{\\frac{P_{ba}}{P_{bb}}}|$" ls 3 axis x1y1,\
"Fig5d.dat" index 0 using 1:3 w lp title "$|\\overline{\\frac{P_{bc}}{P_{bb}}}|$" ls 4 axis x1y1


set key at 0.19,0.86
set key spacing 1.2
unset label
unset xlabel
unset ylabel
unset xtics
unset ytics
unset arrow

#"SC_corr_ave_jx_0.05_different_j2.dat" index 0 using 1:(-$4):9 w yerr title "$\\theta _{ba}$" ls 1 axis x1y2,\
#"SC_corr_ave_jx_0.05_different_j2.dat" index 0 using 1:5:10 w yerr title "$-\\theta _{bc}$" ls 2 axis x1y2,\

plot "Fig5d.dat" index 0 using 1:(-$4) w lp title "$\\overline{\\theta _{ba}}$" ls 1 axis x1y2,\
"Fig5d.dat" index 0 using 1:5 w lp title "$-\\overline{\\theta _{bc}}$" ls 2 axis x1y2



unset multiplot




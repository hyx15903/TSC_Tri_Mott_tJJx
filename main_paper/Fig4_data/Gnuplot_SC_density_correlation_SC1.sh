set term epslatex color size 15cm,13cm standalone
#set term postscript color size 10,5
set output "SC_density_correlation_SC1.tex"
#Fig3

s = 2.5
s4 = 2
w = 1
x_range = 24
ks = 1.2
ks_fit = 2.2

fit_w = 4

set style line 1 pt 7 lc rgb "blue" lw w ps s4
set style line 2 pt 65 lc rgb "red" lw w ps s
set style line 3 pt 67 lc rgb "dark-violet" lw w ps s
set style line 4 pt 11 lc rgb "dark-green" lw w ps 1.6
set style line 5 pt 63 lc rgb "black" lw 0.2 ps 1
set style line 6 pt 63 lc rgb "black" lw fit_w dt "." ps s
set style line 7 pt 64 lc rgb "black" lw fit_w dt "_" ps s
set style line 8 pt 65 lc rgb "red" lw fit_w dt "-" ps s

set style line 11 pt 5 lc rgb "dark-violet" lw w ps s4
set style line 12 pt 7 lc rgb "dark-green" lw w ps s4
set style line 13 pt 9 lc rgb "dark-orange" lw w ps s4
set style line 14 pt 11 lc rgb "grey0" lw w ps s4
set style line 15 pt 13 lc rgb "dark-goldenrod" lw w ps s4

set multiplot layout 2,2 rowsfirst

#1
set ylabel "$\\left | P_{bb}(r)\\right |$" offset 2.5
set xlabel "$r$" offset 0,1.2
set xr [1:x_range]
set yr [5e-6:5e-3]
set logscale y
set logscale x
set format y "$10^{%T}$"
set xtics ("$1$" 1, "" 2, "" 3, "" 4, "" 5, "" 6, "" 7, "" 8, "" 9, "$10$" 10, "$20$" 20) offset 0,0.3
set ytics offset 1,-0.2
set key at 21.2,2e-3
set key samplen ks
set key spacing 1
#set grid
set border lw 4
#set tics nomirror
set lmargin at screen 0.08
set rmargin at screen 0.47
set tmargin at screen 0.96
set bmargin at screen 0.59
set label "$(a)$" at 1.1,7e-3
set label "$L_{y} = 4$" at 10,3e-3
#set label "$J_{2}=0.01, J_{\\chi }=0.05$" at 3.2,3e-3

a(x) = m0+m1*x
fit a(x) "jx0.05j20.01/m10000Lx48/Sc_cor1.dat" index 0 u (log($4)):(log($5)) every 5::6::23 via m0,m1

b(x) = n0+n1*x
fit b(x) "jx0.05j20/m10000Lx48/Sc_cor1.dat" index 0 u (log($4)):(log($5)) every 6::8::23 via n0,n1

#set title "doping=1/12, 24*6, t=3, $J_{2}=0.1, J_{\\chi}=0.1$, aa-bond"
plot "jx0.05j20.01/m10000Lx48/Sc_cor1.dat" index 0 using 4:5 w lp title "$J_{2}=0.01$" ls 3,\
"jx0.05j20/m10000Lx48/Sc_cor1.dat" index 0 using 4:5 w lp title "$J_{2}=0$" ls 1

set key samplen ks_fit
set key at 5.9,2.4e-5
unset label
unset xlabel
unset ylabel
unset xtics
unset ytics
plot exp(m0)*x**m1 title sprintf("$K_{SC}=%.2f$", -m1) ls 6,\
exp(n0)*x**n1 title sprintf("$K_{SC}=%.2f$", -n1) ls 7


#2
reset
w = 1.2
del = 0.08333333

set style line 1 pt 7 lc rgb "blue" lw w ps s4
set style line 2 pt 64 lc rgb "red" lw fit_w dt 2 ps s
set style line 11 pt 5 lc rgb "dark-violet" lw w ps 1.4
set style line 6 pt 63 lc rgb "black" lw fit_w dt "." ps s
set style line 7 pt 64 lc rgb "black" lw fit_w dt "_" ps s
set style line 3 pt 67 lc rgb "dark-violet" lw w ps s
set style line 4 pt 67 lc rgb "yellow4" lw fit_w dt 2 ps s
set style line 13 pt 11 lc rgb "blue" lw w ps s

set ylabel "$\\left |D(r)\\right |$" offset 2.5
set xlabel "$r$" offset 0,1.2
set xr [1:x_range]
set yr [1e-5:1e-2]
set xtics ("$1$" 1, "" 2, "" 3, "" 4, "" 5, "" 6, "" 7, "" 8, "" 9, "$10$" 10, "$20$" 20) offset 0,0.3
set ytics offset 1,-0.2
set logscale y
set logscale x
set format y "$10^{%T}$"
set key at 21,6e-3
set key samplen ks
set key spacing 1

set lmargin at screen 0.58
set rmargin at screen 0.98
set tmargin at screen 0.96
set bmargin at screen 0.59

set border lw 5
set label "$(b)$" at 1.1,0.014


d(x) = d0+d1*x
fit d(x) "jx0.05j20.01/m10000Lx48/density_cor.dat" index 0 u (log($3)):(log(abs($4-($5*$6)))) every 6::3::16 via d0,d1

e(x) = e0+e1*x
fit e(x) "jx0.05j20/m10000Lx48/density_cor.dat" index 0 u (log($3)):(log(abs($4-($5*$6)))) every 6::3::23 via e0,e1

plot "jx0.05j20.01/m10000Lx48/density_cor.dat" index 0 using 3:(abs($4-($5*$6))) w lp title "$J_{2}=0.01$" ls 3,\
"jx0.05j20/m10000Lx48/density_cor.dat" index 0 using 3:(abs($4-($5*$6))) w lp title "$J_{2}=0$" ls 1

set key samplen ks_fit
set key at 7,5e-5
unset label
unset xlabel
unset ylabel
unset xtics
unset ytics

plot exp(d0)*x**d1 title sprintf("$K _{CDW}=%.2f$", -d1) ls 6,\
exp(e0)*x**e1 title sprintf("$K _{CDW}=%.2f$", -e1) ls 7



#3
reset
set style line 1 pt 7 lc rgb "blue" lw w ps s4
set style line 2 pt 64 lc rgb "red" lw fit_w dt 2 ps s
set style line 11 pt 5 lc rgb "dark-violet" lw w ps 1.4
set style line 6 pt 63 lc rgb "black" lw fit_w dt "." ps s
set style line 7 pt 64 lc rgb "black" lw fit_w dt "_" ps s
set style line 3 pt 67 lc rgb "dark-violet" lw w ps s
set style line 4 pt 67 lc rgb "yellow4" lw fit_w dt 2 ps s
set style line 13 pt 11 lc rgb "blue" lw w ps s

set ylabel "$\\left |G(r)\\right |$" offset 2.5
set xlabel "$r$" offset 0,0.8
set xr [1:x_range]
set yr [2e-7:1e-1]
set xtics 0,5 offset 0,0.3
set ytics offset 1,-0.2
set logscale y
set format y '$10^{%T}$'
set format x '%g'
set key at 23,3.5e-2
set key samplen ks
set key spacing 1

set lmargin at screen 0.08
set rmargin at screen 0.47
set tmargin at screen 0.49
set bmargin at screen 0.10

set border lw 5
set label "$(c)$" at 1.5,0.17


c(x) = c0+c1*x
fit c(x) "jx0.05j20.01/m10000Lx48/elec_cor.dat" index 0 u 3:(log(abs($4))) every 4::8::18 via c0,c1

d(x) = d0+d1*x
fit d(x) "jx0.05j20/m10000Lx48/elec_cor.dat" index 0 u 3:(log(abs($4))) every 7::7::18 via d0,d1

plot "jx0.05j20.01/m10000Lx48/elec_cor.dat" index 0 using 3:(abs($4)) w lp title "$J_{2}=0.01$" ls 3,\
"jx0.05j20/m10000Lx48/elec_cor.dat" index 0 using 3:(abs($4)) w lp title "$J_{2}=0$" ls 1

set key samplen ks_fit
set key at 13,3e-6
unset label
unset xlabel
unset ylabel
unset xtics
unset ytics

plot exp(c0+c1*x) title sprintf("$\\xi _{G}=%.2f$", -1/c1) ls 6,\
exp(d0+d1*x) title sprintf("$\\xi _{G}=%.2f$", -1/d1) ls 7


#4
reset
set style line 1 pt 7 lc rgb "blue" lw w ps s4
set style line 3 pt 67 lc rgb "dark-violet" lw w ps s
#set style line 1 pt 66 lc rgb "blue" lw w ps s
#set style line 2 pt 67 lc rgb "dark-green" lw w ps s
#set style line 3 pt 2 lc rgb "red" lw w ps s
#set style line 4 pt 4 lc rgb "yellow4" lw w ps s

set style line 6 pt 63 lc rgb "black" lw fit_w dt "." ps s
set style line 7 pt 64 lc rgb "black" lw fit_w dt "_" ps s
#set style line 8 pt 63 lc rgb "black" lw fit_w dt "-" ps s
#set style line 9 pt 64 lc rgb "black" lw fit_w dt 2 ps s

set ylabel "$\\left |S(r)\\right |$" offset 2.5
set xlabel "$r$" offset 0,0.8
set xr [1:x_range]
set yr [8e-7:0.2]
set xtics 0,5 offset 0,0.3
set ytics offset 1,-0.2
set logscale y
set format y '$10^{%T}$'
set format x '%g'
set key at 23,8e-2
set key samplen ks
set key spacing 1

set lmargin at screen 0.58
set rmargin at screen 0.98
set tmargin at screen 0.49
set bmargin at screen 0.10
set border lw 4


k(x) = e0+e1*x
fit k(x) "jx0.05j20.01/m10000Lx48/spin_cor.dat" index 0 u 3:(log(abs($4))) every 6::6::24 via e0,e1

l(x) = f0+f1*x
fit l(x) "jx0.05j20/m10000Lx48/spin_cor.dat" index 0 u 3:(log(abs($4))) every 4::3::24 via f0,f1


set label "$(d)$" at 1.5,0.34
#set label "$L_{x}=24$" at 4,7e-2
#set label "$L_{x}=36$" at 4,1.5e-2
#set label "$J_{2}=0.1$" at 3,7e-6
#set label "$L_{x}=36$" at 4,1.5e-2
set key samplen 1


plot "jx0.05j20.01/m10000Lx48/spin_cor.dat" index 0 using 3:(abs($4)) w lp title "$J_{2}=0.01$" ls 3,\
"jx0.05j20/m10000Lx48/spin_cor.dat" index 0 using 3:(abs($4)) w lp title "$J_{2}=0$" ls 1


unset label
unset xlabel
unset ylabel
unset xtics
unset ytics
set key samplen ks_fit
set key at 12,9e-6


#set label sprintf("$\\xi _{S}=%.2f$", -1./c1) at 12,1e-5

plot (x > 35 ? 1/0 : exp(e0+e1*x)) title sprintf("$\\xi _{S}=%.2f$", -1/e1) ls 6,\
(x > 35 ? 1/0 : exp(f0+f1*x)) title sprintf("$\\xi _{S}=%.2f$", -1/f1) ls 7



unset multiplot




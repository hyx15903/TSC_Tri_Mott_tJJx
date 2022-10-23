set term epslatex color size 15cm,13cm standalone
#set term postscript color size 10,5
set output "SC_density_correlation_fit_new.tex"
#Fig3

s = 1.4
s4 = 1.2
w = 1

fit_w = 4

set style line 1 pt 64 lc rgb "blue" lw w ps s
set style line 2 pt 65 lc rgb "red" lw w ps s
set style line 3 pt 66 lc rgb "magenta" lw w ps s
set style line 4 pt 11 lc rgb "dark-green" lw w ps 1.6
set style line 5 pt 63 lc rgb "black" lw 0.2 ps 1
set style line 6 pt 63 lc rgb "black" lw fit_w dt 2 ps s
set style line 7 pt 64 lc rgb "blue" lw fit_w dt "-" ps s
set style line 8 pt 65 lc rgb "red" lw fit_w dt "-" ps s

set style line 11 pt 5 lc rgb "dark-violet" lw w ps s4
set style line 12 pt 7 lc rgb "dark-green" lw w ps s4
set style line 13 pt 9 lc rgb "dark-orange" lw w ps s4
set style line 14 pt 11 lc rgb "grey0" lw w ps s4
set style line 15 pt 13 lc rgb "dark-goldenrod" lw w ps s4

set multiplot layout 3,4 rowsfirst

#1
set ylabel "$\\left | P_{bb}(r)\\right |$" offset 2.5
set xlabel "$r$" offset 0,1
set xr [1:20]
set yr [1e-5:5e-3]
set logscale y
set logscale x
set format y "$10^{%T}$"
set xtics ("$1$" 1, "" 2, "" 3, "" 4, "" 5, "" 6, "" 7, "" 8, "" 9, "$10$" 10, "$20$" 20) offset 0,0.3
set ytics offset 1
set key at 5.2,1.65e-4
set key samplen 0.8
set key spacing 0.8
#set grid
set border lw 4
#set tics nomirror
set lmargin at screen 0.08
set rmargin at screen 0.47
set tmargin at screen 0.96
set bmargin at screen 0.59
set label "$(a)$" at 1.1,6.4e-3
set label "$L_{x} = 36$" at 7,3e-3
set label "$J_{2}=0.1, J_{\\chi }=0.05$" at 3.2,1.5e-3

a(x) = m0+m1*x
fit a(x) "Sc_cor1_J2_0.1_Jchi_0.05_Lx_36_bb_bond_infinite.dat" index 0 u (log($1)):(log($2)) every 4::5::17 via m0,m1

#set title "doping=1/12, 24*6, t=3, $J_{2}=0.1, J_{\\chi}=0.1$, aa-bond"
plot "jx0.05j20.1/m6000Lx36/Sc_cor1.dat" index 0 using 4:5 w lp title "$M=6000$" ls 1,\
"jx0.05j20.1/m8000Lx36/Sc_cor1.dat" index 0 using 4:5 w lp title "$M=8000$" ls 2,\
"jx0.05j20.1/m10000Lx36/Sc_cor1.dat" index 0 using 4:5 w lp title "$M=10000$" ls 3,\
"jx0.05j20.1/m12000Lx36/Sc_cor1.dat" index 0 using 4:5 w lp title "$M=12000$" ls 4,\
"Sc_cor1_J2_0.1_Jchi_0.05_Lx_36_bb_bond_infinite.dat" index 0 using 1:2 every ::1 w lp title "$M\\rightarrow \\infty $" ls 5

set key samplen 2
set key at 5.4,2e-5
unset label
unset xlabel
unset ylabel
unset xtics
unset ytics
plot exp(m0)*x**m1 title sprintf("$K_{SC}=%.2f$", -m1) ls 6


#2
reset
s = 1.4
s4 = 1.2
w = 1.2
del = 0.08333333

set style line 1 pt 64 lc rgb "red" lw w ps 1.5
set style line 2 pt 64 lc rgb "red" lw fit_w dt 2 ps s
set style line 11 pt 5 lc rgb "dark-violet" lw w ps 1.4

set style line 3 pt 67 lc rgb "yellow4" lw w ps 2
set style line 4 pt 67 lc rgb "yellow4" lw fit_w dt 2 ps s
set style line 13 pt 11 lc rgb "blue" lw w ps 1.8

set ylabel "$\\left |G(r)\\right |,\\left |D(r)\\right |$" offset 2.5
set xlabel "$r$" offset 0,0.8
set xr [1:20]
set yr [1e-5:1e-1]
set xtics ("$1$" 1, "" 2, "" 3, "" 4, "" 5, "" 6, "" 7, "" 8, "" 9, "$10$" 10, "$20$" 20) offset 0,0.3
set ytics offset 1
set logscale y
set logscale x
set format y "$10^{%T}$"
set key at 11.6,2e-5
set key samplen 0.8
set key spacing 0.4
set key vertical maxrows 1
set key width -6.5
set lmargin at screen 0.58
set rmargin at screen 0.98
set tmargin at screen 0.96
set bmargin at screen 0.59

set border lw 5
set label "$(b)$" at 1.1,0.14
set label "$L_{x}=36$" at 1.03,1.5e-5
set label "$L_{x}=24$" at 1.03,4e-5

c(x) = c0+c1*x
fit c(x) "elec_cor_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" index 0 u (log($1)):(log(abs($2))) every 4::6::14 via c0,c1

d(x) = d0+d1*x
fit d(x) "density_cor_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" index 0 u (log($1)):(log(abs($2))) every 4::6::14 via d0,d1

set label sprintf("$K _{G}=%.2f$", -c1) at 4,3e-2
set label sprintf("$K _{CDW}=%.2f$", -d1) at 1.05,2e-4

plot "elec_cor_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" index 0 using 1:(abs($2)) every ::1 w lp title "$\\left |G(r)\\right |$" ls 1,\
exp(c0)*x**c1 notitle ls 2,\
"density_cor_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" index 0 using 1:(abs($2)) every ::1 w lp title "$\\left |D(r)\\right |$" ls 3,\
exp(d0)*x**d1 notitle ls 4

set key at 11.6,5e-5
unset label
unset xlabel
unset ylabel
unset xtics
unset ytics
plot "elec_cor_J2_0.1_Jchi_0.05_Lx_24_infinite.dat" index 0 using 1:(abs($2)) every ::1::14 w lp title "$\\left |G(r)\\right |$" ls 11,\
"density_cor_J2_0.1_Jchi_0.05_Lx_24_infinite.dat" index 0 using 1:(abs($2)) every ::1::14 w lp title "$\\left |D(r)\\right |$" ls 13



#3
reset
set style line 1 pt 64 lc rgb "blue" lw w ps 1.5
set style line 2 pt 64 lc rgb "blue" lw fit_w dt 2 ps s
set style line 3 pt 67 lc rgb "red" lw w ps 2
set style line 4 pt 67 lc rgb "red" lw fit_w dt 2 ps s
set style line 11 pt 5 lc rgb "dark-violet" lw w ps 1.4
set style line 12 pt 11 lc rgb "dark-green" lw w ps 1.8
set ylabel "$\\left |G(r)\\right |$" offset 2.5
set xlabel "$r$" offset 0,1

set logscale y
set format y "$10^{%T}$"
set xr [1:20]
set yr [1e-6:1e-1]
#set xtics ("$1$" 1, "" 2, "" 3, "" 4, "" 5, "" 6, "" 7, "" 8, "" 9, "$10$" 10, "$20$" 20)
set ytics offset 1
set xtics 0,5 offset 0,0.3
set key at 18,3e-6
set key samplen 0.8
set key spacing 0.9
set key vertical maxrows 1
set key width -0.1
#set grid
set border lw 4
#set tics nomirror
set lmargin at screen 0.08
set rmargin at screen 0.47
set tmargin at screen 0.49
set bmargin at screen 0.10
set label "$(c)$" at 1.5,1.7e-1
set label "$J_{2} = 0.05$" at 1.4,5e-6
set label "$J_{2} = 0.1$" at 2.3,5e-2

w(x) = w0+w1*x
fit w(x) "elec_cor_J2_0.05_Jchi_0.05_Lx_36_m10000.dat" index 0 using 3:(log(abs($4))) every 4::9::17 via w0,w1

c(x) = c0+c1*x
fit c(x) "elec_cor_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" index 0 u 1:(log(abs($2))) every 4::6::18 via c0,c1

set label sprintf("$\\xi _{G}=%.2f$", -1./w1) at 13,1.2e-4 rotate by -30

set label sprintf("$\\xi _{G}=%.2f$", -1./c1) at 13,7e-3 rotate by -10.6

plot "elec_cor_J2_0.05_Jchi_0.05_Lx_36_m10000.dat" index 0 using 3:(abs($4)) every ::1 w lp title "$L_{x}=36$" ls 1,\
"elec_cor_J2_0.05_Jchi_0.05_Lx_24_m10000.dat" index 0 using 3:(abs($4)) every ::1::14 w lp title "$L_{x}=24$" ls 11,\
(x < 3.5 ? 1/0 : exp(w0+w1*x)) notitle ls 2,\
(x < 3.2 ? 1/0 : exp(c0+c1*x)) notitle ls 4

set key at 20.5,8e-2
unset label
unset xlabel
unset ylabel
unset xtics
unset ytics

plot "elec_cor_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" index 0 using 1:(abs($2)) every ::1 w lp title "$L_{x}=36$" ls 3,\
"elec_cor_J2_0.1_Jchi_0.05_Lx_24_infinite.dat" index 0 using 1:(abs($2)) every ::1::14 w lp title "$L_{x}=24$" ls 12



#4
reset
set style line 1 pt 65 lc rgb "blue" lw w ps s
set style line 2 pt 65 lc rgb "blue" lw fit_w dt 2 ps s
set style line 11 pt 7 lc rgb "dark-green" lw w ps s4

set style line 3 pt 64 lc rgb "red" lw w ps s
set style line 4 pt 64 lc rgb "red" lw 1 dt "-" ps s
set style line 12 pt 5 lc rgb "dark-violet" lw w ps s4

set ylabel "$\\left |S(r)\\right |$" offset 2.5
set xlabel "$r$" offset 0,0.8
set xr [1:20]
set yr [1e-6:1e-1]
set xtics 0,5 offset 0,0.3
set ytics offset 1
set logscale y
set format y '$10^{%T}$'
set format x '%g'
set key at 14,3e-6
set key samplen 0.8
set key spacing 0.4
set key vertical maxrows 1
set key width -0.2
set lmargin at screen 0.58
set rmargin at screen 0.98
set tmargin at screen 0.49
set bmargin at screen 0.10
set border lw 4
b(x) = b0+b1*x
fit b(x) "spin_cor_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" index 0 u 1:(log(abs($2))) every 2::2::12 via b0,b1

k(x) = k0+k1*x
fit k(x) "spin_cor_J2_0.05_Jchi_0.05_Lx_36_infinite.dat" index 0 u 1:(log(abs($2))) every 6::6::24 via k0,k1


set label "$(d)$" at 1.5,0.15
#set label "$L_{x}=24$" at 4,7e-2
#set label "$L_{x}=36$" at 4,1.5e-2
set label "$J_{2}=0.1$" at 3,7e-6
#set label "$L_{x}=36$" at 4,1.5e-2

set label sprintf("$\\xi _{S}=%.2f$", -1./b1) at 4,3e-2 rotate by -42

plot "spin_cor_J2_0.1_Jchi_0.05_Lx_24_infinite.dat" index 0 using 1:(abs($2)) every ::1::14 w lp title "$L_{x}=24$" ls 11,\
"spin_cor_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" index 0 using 1:(abs($2)) every ::1 w lp title "$L_{x}=36$" ls 1,\
(x > 15 ? 1/0 : exp(b0+b1*x)) notitle ls 2


#4inset
set key at 21,8e-2
set key samplen 0.8
set key spacing 0.7
set key horizontal
unset label
set ylabel "$\\left |S(r)\\right |$" offset 2.4,1
set xlabel "$r$" offset 0,1.3
set ytics 1e-5,100
set xtics 0,5 offset 0,0.4
set format y '$10^{%T}$'
set format x '%g'
set label "$J_{2}=0.05$" at 4,1.1
set label sprintf("$\\xi _{S}=%.2f$", -1./k1) at 3,5e-4 rotate by -42

set border lw 1
set lmargin at screen 0.8
set rmargin at screen 0.96
set tmargin at screen 0.47
set bmargin at screen 0.31

plot "spin_cor_J2_0.05_Jchi_0.05_Lx_24_infinite.dat" index 0 using 1:(abs($2)) every ::1::14 w lp title "$L_{x}=24$" ls 12,\
"spin_cor_J2_0.05_Jchi_0.05_Lx_36_infinite.dat" index 0 using 1:(abs($2)) every ::1 w lp title "$L_{x}=36$" ls 3,\
(x > 20 ? 1/0 : exp(k0+k1*x)) notitle ls 4

unset multiplot




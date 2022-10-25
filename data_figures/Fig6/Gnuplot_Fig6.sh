set term epslatex color colortext size 13cm,14cm standalone
#set term postscript color
set output "Fig6.tex"


s=2
w=3
border_w = 4
line_w = 3

set multiplot layout 2,3 rowsfirst
#1

set label "$(a)$" at 0.005,-1.895

set style line 1 pt 6 lc rgb "blue" lw w ps s
set style line 2 pt 7 lc rgb "dark-goldenrod" lw w ps s
set style line 3 pt 8 lc rgb "black" lw w ps 2.3
set style line 4 pt 11 lc rgb "red" lw w ps s


set ylabel "$ E_{0}$" offset 2.5 tc rgb "blue"
set y2label "$S$" offset -2.5,0 rotate by -90
set xlabel "$ J_{2}$" offset 0,0.7
set lmargin at screen 0.11
set rmargin at screen 0.44
set tmargin at screen 0.97
set bmargin at screen 0.73
set border lw border_w
set key samplen 0.8
set key spacing 0.8
set key right bottom
set yr [-2:-1.9]
set xr [0:0.16]
set y2r [2.5:6.1]
set y2tics 3,0.5 offset -0.8
set ytics -2,0.02 offset 0.8 tc rgb "blue"
set ytics nomirror
set xtics 0,0.04 offset 0,0.4


set style arrow 4 head noborder size screen 0.03,15,45 dt "-" lc rgb "red" lw 1.5
set arrow 1 from 0.02125,-2 to 0.02125,-1.96 as 4 lc rgb "red"
set label "$J_{2c}^{(1)}$" at 0.004,-2.014 textcolor rgb "red"


plot "Fig6.dat" u 1:($2/96) w lp notitle ls 1 axis x1y1,\
"Fig6.dat" u 1:3 w lp notitle ls 3 axis x1y2



#2
reset
set style line 1 pt 6 lc rgb "blue" lw w ps s
set style line 2 pt 7 lc rgb "dark-goldenrod" lw w ps s
set style line 3 pt 8 lc rgb "black" lw w ps 2.3
set style line 4 pt 11 lc rgb "red" lw w ps s

set lmargin at screen 0.625
set rmargin at screen 0.975
set tmargin at screen 0.97
set bmargin at screen 0.73
set border lw border_w
set key samplen 0.8
set key spacing 0.8

set xr [0.002:0.06]

set ytics -4,1 offset 0.8,-0.3 tc rgb "blue"


set xtics 0,0.02  offset 0,0.4
set ylabel "$ -dE_{0}/dJ_{2}$" offset 1.5,-0 tc rgb "blue"

set xlabel "$ J_{2}$" offset 0,0.7

set style arrow 4 head noborder size screen 0.03,15,45 dt "-" lc rgb "red" lw 1.5
set arrow 1 from 0.02125,-4.2 to 0.02125,-1.3 as 4 lc rgb "red"
set label "$J_{2c}^{(1)}$" at 0.019,-5.25 textcolor rgb "red"


set label "$(b)$" at 0.0034,1.26

x0=NaN
y0=NaN
x20=NaN
y20=NaN

plot "Fig6.dat" u (dx=$1-x0,x0=$1,$1-dx/2):(dy=($2/96)-y0,y0=($2/96),-dy/dx) w lp notitle ls 1 axis x1y1





#3
reset
set style line 1 pt 6 lc rgb "blue" lw w ps s
set style line 2 pt 7 lc rgb "dark-goldenrod" lw w ps s
set style line 3 pt 8 lc rgb "black" lw w ps 2.3
set style line 4 pt 11 lc rgb "red" lw w ps s
set style arrow 3 head filled size screen 0.03,15,45 ls 4

set ytics 0,0.4 offset 0.8 tc rgb "blue"
set xtics 0.04,0.04  offset 0,0.4
set ylabel "$ -dE_{0}/dJ_{2}$" offset 1 tc rgb "blue"
set xlabel "$ J_{2}$" offset 0,0.7
set lmargin at screen 0.11
set rmargin at screen 0.44
set tmargin at screen 0.64
set bmargin at screen 0.40
set border lw border_w
set xr [0.04:0.16]

#set style arrow 4 head noborder size screen 0.03,15,45 dt "-" lc rgb "red" lw 1.5
#set arrow 1 from 0.0825,-35 to 0.0825,95 as 4 lc rgb "red"
#set label "$J_{2c}^{(2)}$" at 0.0825,-40 textcolor rgb "red"
set label "$(c)$" at 0.044,1.28

x0=NaN
y0=NaN

plot "Fig6.dat" u (dx=$1-x0,x0=$1,$1-dx/2):(dy=($2/96)-y0,y0=($2/96),-dy/dx) w lp notitle ls 1 axis x1y1


#4
unset arrow
unset label
set style arrow 4 head filled size screen 0.03,15,45 ls 4
set ytics -20,4 offset 0.8 tc rgb "blue"
set label "$(c)$" at 0.043,-100
set ylabel "$ d^{2}E_{0}/dJ_{2}^{2}$" offset 1 tc rgb "blue"
set xlabel "$ J_{2}$" offset 0,0.7
set lmargin at screen 0.625
set rmargin at screen 0.975
set tmargin at screen 0.64
set bmargin at screen 0.40
set border lw border_w
set xr [0.04:0.16]
set arrow 1 from 0.075,-2300 to 0.075,-1900 as 4 lc rgb "red"
set label "$J_{2c}^{(2)}$" at 0.075,-2350 textcolor rgb "red"



set label "$(d)$" at 0.044,1

x0=NaN
y0=NaN
x1=NaN
y1=NaN


plot "Fig6.dat" u (dx=$1-x1,dx1=x1-x0,x0=x1,x1=$1,$1-dx/2-(dx+dx1)/4):(dy=($2/96-y1)/dx-(y1-y0)/dx1,y0=y1,y1=$2/96,dy/((dx+dx1)/2)) w lp notitle ls 1


#5
reset
set style line 1 pt 6 lc rgb "blue" lw w ps s
set style line 2 pt 7 lc rgb "dark-goldenrod" lw w ps s
set style line 3 pt 8 lc rgb "black" lw w ps 2.3
set style line 4 pt 11 lc rgb "red" lw w ps s
set style arrow 3 head filled size screen 0.03,15,45 ls 4

set ytics 0,20 offset 0.8

set xtics 0.04,0.04  offset 0,0.4

set ylabel "$-dS/dJ_{2}$" offset 1 
set xlabel "$ J_{2}$" offset 0,0.7
set lmargin at screen 0.11
set rmargin at screen 0.44
set tmargin at screen 0.31
set bmargin at screen 0.07
set border lw border_w
set xr [0.04:0.16]
set yrange [-3:60]

set style arrow 4 head noborder size screen 0.03,15,45 dt "-" lc rgb "red" lw 1.5
set arrow 1 from 0.0825,-3 to 0.0825,43 as 4 lc rgb "red"

set label "$J_{2c}'$" at 0.078,-15 textcolor rgb "red"
set label "$(e)$" at 0.043,63.5

x20=NaN
y20=NaN

plot "Fig6.dat" u (dx2=$1-x20,x20=$1,$1-dx2/2):(dy2=($3)-y20,y20=($3),-dy2/dx2) w lp notitle ls 3



#6
unset arrow
unset label
set ytics -3000,1000 offset 0.8

set ylabel "$d^{2}S/dJ_{2}^{2}$" offset 1.5
set xlabel "$ J_{2}$" offset 0,0.7
set lmargin at screen 0.625
set rmargin at screen 0.975
set tmargin at screen 0.31
set bmargin at screen 0.07
set border lw border_w
set xr [0.04:0.16]
set yr [-3000:2000]

set arrow 2 from 0.1,-3000 to 0.1,1100 as 4 lc rgb "red"
set label "$J_{2c}^{(2)}$" at 0.101,-2500 textcolor rgb "red"
set label "$(f)$" at 0.043,2300

x20=NaN
y20=NaN
x21=NaN
y21=NaN

plot "Fig6.dat" u (dx2=$1-x21,dx21=x21-x20,x20=x21,x21=$1,$1-dx2/2-(dx2+dx21)/4):(dy2=($3-y21)/dx2-(y21-y20)/dx21,y20=y21,y21=$3,dy2/((dx2+dx21)/2)) w lp notitle ls 3



unset multiplot

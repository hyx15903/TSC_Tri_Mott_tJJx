extended t-J-Jchi triangular model main paper Fig. 2

The file "chern_J2_{x1}_Jchi_{x2}_Ly_{y1}.dat" contains the adiabatical flux insertion data for J2=x1, Jchi=x2, Ly=y1, doping=1/12
The first column stands for the flux, the second column the energy, and the third column the accumulated spin 


The file "chern_num_Jchi_{x2}_Ly_{y1}_J2.dat" contains the extracted Chern number for Jchi=x2, Ly=y1, doping=1/12 and various J2
The first column stands for J2, and the second column the Chern number


The file "chern_num_Ly_{y1}_doping.dat" contains the extracted Chern number for Ly=y1 and various doping
The first column stands for doping, the second column the Chern number at jchi=0.05, J2=0.04, the third column the Chern number at jx=0.05, j2=0.05 , the fourth column the Chern number at jx=0.1, j2=0.1


To compile, use
$gnuplot Gnuplot_Chern.sh
$latex Chern.tex
$dvips Chern.dvi
$ps2pdf Chern.ps

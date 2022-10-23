extended t-J-Jchi triangular model main paper Fig. 3

The files "jx0.05j20.1/m{M}Lx36/Sc_cor1.dat" contains the SC pairing correlation data for J2=0.05, Jchi=0.1, Ly=6, Lx=36, doping=1/12, bond dimension=M
The files "elec_cor_J2_0.05_Jchi_0.05_Lx_36_m10000.dat" contains the SC pairing correlation data for J2=0.05, Jchi=0.05, Ly=6, Lx=36, doping=1/12, bond dimension=10000
The first column stands for the site i=xi*Ly+yi, the second column the site j=xj*Ly+yj, we can use the fourth column for the distance in the x-direction between i and j.
The 5th column stands for the b-b bond SC correlations, 6th column the b-a bond, 7th column the b-c bond.
The 8th column stand for the SC pairing phase between b-b bond, 9th column b-a bond, 10th column b-c bond


The file "Sc_cor1_J2_0.1_Jchi_0.05_Lx_36_bb_bond_infinite.dat" contains the SC pairing correlation data on b-b bond
The file "elec_J2_0.1_Jchi_0.05_Lx_36__infinite.dat" contains the electron hopping correlation data
The file "density_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" contains the density-density correlation data
The file "spin_cor_J2_0.1_Jchi_0.05_Lx_36_infinite.dat" contains the density-density correlation data
all for J2=0.1, Jchi=0.05, Ly=6, doping=1/12 in the infinite M limit.
The first column stands for the distance x between i and j site, and the second column the bond dimension scaled correlation. The 3rd and 4th column are the parameters for the fitting polynomials, and the 5th column is the standard deviation.


To compile, use
$gnuplot Gnuplot_SC_density_correlation_fit_new.sh
$latex SC_density_correlation_fit_new.tex
$dvips SC_density_correlation_fit_new.dvi
$ps2pdf SC_density_correlation_fit_new.ps

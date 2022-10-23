extended t-J-Jchi triangular model main paper Fig. 5

The files "jx{Jx}j2{J2}/Sc_cor_ph.dat" contains the SC pairing correlation data for J2=J2, Jchi=Jx, Ly=6, Lx=36, doping=1/12, bond dimension=10000
The files "jx{Jx}j2{J2}Ny4/Sc_cor_ph.dat.dat" contains the SC pairing correlation data for J2=J2, Jchi=Jx, Ly=4, Lx=48, doping=1/12, bond dimension=8000

The first column stands for the site i=xi*Ly+yi, the second column the site j=xj*Ly+yj, we can use the fourth column for the distance in the x-direction between i and j.
The 5th column stands for the b-b bond SC correlations, 6th column the b-a bond, 7th column the b-c bond.
The 8th column stand for the SC pairing phase between b-b bond, 9th column b-a bond, 10th column b-c bond


The files "SC_corr_ave_jx_0.05_different_j2.dat" contains the SC pairing correlation average over distance of x = 6 - 25, for different J2 at Jchi=0.05, Ly=6, Lx=36, doping=1/12.
The values representing each column is below:
j2 ; Pba/Pbb ; Pbc/Pbb ; \theta_ba ; \theta_bc ; \theta ca ; sta_div(Pba/Pbb) ; sta_div(Pbc/Pbb) ; sta_div(\theta_ba) ; sta_div(\theta_bc) ; sta_div(\theta ca)


To compile, use
$bash run.sh

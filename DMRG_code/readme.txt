
The numerical code calculates the ground state of t1-t2-J1-J2-Jx model see
Yixuan Huang and D. N. Sheng, Phys. Rev. X 12, 031009 (2022) for Hamiltonian.

The lattice size Num_site=Nx*Ny, lengths parameters (Nx, Ny) can be found in common.f90.  The lattice geometry can be determined as square=1,  or triangle=1, we set other lattice parameters as zero.

The code main.f90 has more parameters you can adjust:

From lines 55, 56, kept_min, and kept_max are range of bond dimensions.
Set kept_max as initial disired bond dimensions.   kept_maxf as final bond
dimension. Sweep_num usually needs to be around 6-10+ (more sweeps give better convergence). For example, one can set:

kept_min=1600
kept_max=2000
Kept_max=4000
sweep_num=6
you job will be done with bond dimension 4000, afetr six complete sweeps.

For triangular lattice:
jz(1:6, 1:num_site)=J1     nearest neighbor (NN) Heisenberg coupling
jz(7:12, 1:num_site)=J2   next nearest neighbor (NNN) Heisenberg coupling
jt(1:6, 1:num_site)=-t1*dsqrt(2.0) NN hopping, the  factor sqrt(2.0) is due to 
                                    SU2 code (which adjusts coupling constant). 

jt(7:12, 1:num_site)=-t2*dsqrt(2.0)   NNN hopping,  the factor sqrt(2.0) is due to SU2 code.
jn is the coupling coef. for n_in_j terms,  we have jn=-0.25*jz
jring(1:num_site)= Jx    for chiral interactions
delta is an integer, which gives hole doping level, 1/delta.
e.g., set delta=12 for hole doping level 1/12. 

to start the job, you do 
make 
which generates the compiled program tjmodel
The intel new compiler ifort version 2021.5.0 is used.   For other version, one needs to do adjustment for the makefile.

For initial job,  you need to edit restart.dat to be in the form
0 0
0

The useful files and measurement results
 density_cor.dat   density correlations 
      3rd   4th      5th      6th  (columns)
    r_ji   <n_in_j>  <n_i>  <n_j>
 elec_cor.dat       hopping corr.      use 3rd and 4th columns
 spin_cor.dat       spin correlations,  use 3rd and 4th columns  
 elec_density.dat    onsite electron density,    x, y, n(x,y) format      
 Energy_All.dat      total energy
 entropy1.dat        X  entropy ,    X is subsystem length
 Kept_State.dat   actual bond dimension during the DMRG
 restart.dat       for initial of job

 Sc_cor_x.dat
 Sc_cor_xy.dat
 Sc_cor_y.dat
	above are pairing correlations for measurements of different bonds.
	The first bond direction is indicated by the name of the files as
        x==a,  y==b, xy==c (see the PRX paper mentioned above).
	data format  for SC_cor_x.dat
	4th   5th     6th   7th     8th   9th   10th
         r    |Pxb|, |Pxa| |Pxc|   \phi_xb  \phi_xa  \phi_xc
	r is distance,  P is complex pairing corre, \phi is its phase.

	For other files, replacing x by y, or xy.

common.f90  have commonly used parameters and data structures 

main.f90 have parameters for the job and perform simulations for ground state and measurements.

dmrg_sub.f90  have subroutines for setting up the dmrg, set up neighboring coupings, and chiral interactions. It also contains the dmrg codes for doping infinite process (warm_up), left_to_right and  right_to_left sweepings.

den_mat find the eigenstates of reduced density matrix.

eigen.f90 perform the lanczos calculations for ground state.

get_sys update the dmrg blocks (including operators), perform truncations
 and also operator*vectors for lanczos.

measure.f90 perform the measurements

The DMRG method follows Steve White's original algorithm with including SU2 spin rotational symmetry.    The code was built based on original code of Hong Chen Jiang (2008) and detailed note of Shoushu Gong (2020) for SU2 implimentation.   Many thanks to Hong Chen and Shoushu!

Please contact donna.sheng44@gmail.com for questions. 

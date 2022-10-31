!=================================================
!For the t-J model on the triangle/square Lattice
!=================================================
program tjmodelc
	use pubdata
	implicit none

	real(8),parameter :: trun_err=1.0E-15
	integer :: i,x,y,trun_idx,levels
	real(8) :: bg_time,end_time
  call factrl
         !! for Wigner 3j 6j and 9j coefs^M
 !FUNCTION W3JS(J1,J2,J3,M1,M2,M3)^M

        tjmodel=1
        hubmodel=0
        realcode=.false.
        if(tjmodel==1)N_maX=2
        ss=1  !! for spin 1/2

        !! selected measurements
        density1=.false.
        measbl=.false.
        measli=.true.
        superc=.true.


	!<1>: Set Model parameter
	call Get_Lattice
	call Get_site_operator_su2
	call Get_single_site_basis
	call Get_model_name(sysname,envname)
	call Get_wave_name

!!!!!!! using them for initial or following up calculations with increased bond dimension
        open(1212,file='restart.dat')
        read(1212,*)restart1,sys_len_start
        read(1212,*)restart2
        close(1212)        
        iw6j1=0
        w6j1=0.d0
        iw6j2=0
        w6j2=0.d0
        iw6j3=0
        w6j3=0.d0

        infi_del=0
        cone=1.0d0
        czero=0.0d0
        openy=0  !! periodic condition along y
        ci=dcmplx(0.0d0, 1.0d0)

	!Start program
		kept_min=110
		kept_max=300
        kept_maxf=300
		sweep_num=6

        jz=0.d0
        jt=0.d0
        jd=0.d0
                jz(1:6, 1:num_site)=1.0d0
                jz(7:12, 1:num_site)=0.10d0
        jt(1:6,1:num_site )=-3.0d0*dsqrt(2.0d0)
        jt(7:12, 1:num_site)=-3.0d0*dsqrt(cdabs(jz(7,1)))*dsqrt(2.0d0)
        jd(1:18,1:num_site)=jz(1:18,1:num_site)*(-dsqrt(3.0d0))
        jn=0
        if(tjmodel==1)then !!!  t-J
        hubbard=0.0d0
        jn(1:12, 1:num_site)=-0.250d0*jz(1:12, 1:num_site)
        v123=1
        endif

        lring=0
        jring(1:num_site)=0.0d0
        jring(1:num_site)=jring(1:num_site)*dsqrt(6.0d0)
        if(jring(1).ne.0.0)lring=1
        delta=12    !! inverse of hole doping level

        !! adding operators related to chiral interactions
        call dmrg_setup(num_site/2-1, num_site/2-1)
        call dmrg_setup0()
         if(lring.eq.1)then
        call get_triangle()
        endif


		!Get the Wavefunction
        trun_idx=1

		!<1>: For the ground state

        do i=1, num_site
        qn0(i,1)=i/2-(i+2*delta-1)/(2*delta)
        qn0(i,1)=qn0(i,1)*2
        qn0(i,2)=0
        enddo
		tot_num_up=qn0(num_site,1)
		tot_num_down=0   

		!Get ground state
        if(restart1.ne.4)then
		levels=1
		call Get_Ground_State(trun_idx,trun_err,levels)
                endif

		!<2>: Measurement
		call Measurement(trun_idx)

end program tjmodelc


!Get ground_state wavefunction
subroutine Get_Ground_State(trun_idx,trun_err,levels)
	use pubdata
	implicit none

	real(8),intent(in) :: trun_err
	integer,intent(in) :: trun_idx,levels

	logical :: newflag
	integer :: i,iter,tmp_keep,point, idx1
	real(8) :: bg_time,end_time,tmp_err
	real(8),parameter :: lanerr=1.0e-6

	!<1>: For the Warm_up process
	call cpu_time(bg_time)
	open(7,file="Energy_All.dat",position='append')
	write(7,*) "Warm_up:","point=",trun_idx,kept_min, kept_max
	write(7,*) "jz, jt",jz(1:8,1), jt(1:8,1),tot_num_up, num_site
	write(7,*) "kept",kept_min, kept_max
	close(7)

	spectrum=.false.
	lanzero=lanerr*(10**3)
	tmp_err=trun_err*(10**2)

        if(restart1.eq.0)then
	call warm_up_point(trun_idx,Num_site/2,tmp_err,levels)
                endif


	newflag=.false.
	lanzero=lanerr*(15)
	tmp_err=trun_err*(15)

        if(restart1.ne.0)then
	lanzero=lanerr*2 
	tmp_err=trun_err*2 
                endif

        if(restart2.ne.0)then
	lanzero=lanerr 
	tmp_err=trun_err 
                endif


        idx1=trun_idx
        if(restart1.le.1)then
        if(restart1.eq.0)sys_len_start=num_site/2-1
        restart1=restart2
call left_to_right(sys_len_start,Num_site-idx1-1,tmp_err,levels,.true.)
        restart1=1
        restart2=1
                endif

                if(restart1.le.2)then
        if(restart1.le.1)sys_len_start=num_site-idx1
        restart1=restart2
	call right_to_left(sys_len_start-1,idx1,tmp_err,levels,.true.)
        restart1=2
                        endif

        restart1=1

	!<2>: Sweep
	iter=1
	do while(iter<=sweep_num)
                if(iter==sweep_num)infi_del=1
		if(iter.eq.1)kept_max=(kept_max+kept_maxf)/2
		if(iter.eq.2)kept_max=kept_maxf

	lanzero=lanerr
	tmp_err=trun_err
        if(iter.ge.5) lanzero=lanerr*0.50d0

	open(7,file="Energy_All.dat",position='append')
		write(7,*) "Sweep=",iter,"trun_err=",tmp_err
	write(7,*) "jz, jt",jz(1:8,1), jt(1:8,1),tot_num_up, num_site
	write(7,*) "kept",kept_min, kept_max, iter
	close(7)

		call left_to_right(idx1,Num_site-idx1-1,tmp_err,levels,newflag)

		call right_to_left(Num_site-idx1-1,idx1,tmp_err,levels,newflag)

		call cpu_time(end_time)
		open(7,file="Energy_All.dat",position='append')
		write(7,"(A6,F12.4,2X,A3)") "Time =",end_time-bg_time,"Sec"
		write(7,*)
		close(7)

		iter=iter+1
	end do
	spectrum=.true.
	call left_to_right(idx1,5*num_site/6-1,tmp_err,levels,newflag)
        !!finalstop

end subroutine Get_Ground_State


!Perform measurement
subroutine Measurement(trun_idx)
	use pubdata
	implicit none

	integer,intent(in) :: trun_idx
	real(8) :: bg_time,end_time,value
	integer :: i,j,i1, j1,label,sys_len,env_len,x,y

	integer :: name_len,knum
	character(len=30) :: Filename
	type(Total_Model) :: sys,env
	type(Total_Basis) :: sys_bs,env_bs
	type(Wavefunction) :: eigvec
        character ca


	double precision, dimension(Num_site) :: ave_num_elec
	double complex ,dimension(Num_site,Num_site) :: num_elec_cor,spin_sz_cor,spin_sd_cor,num_hole_cor
	double complex,dimension(Num_site,Num_site) :: elec_cor,elec_down_cor,spin_cor
	double precision,dimension(Num_site,Num_site) :: elec_cor1,spin_cor1

	!Read eigen_state from disk (Ground state only)
	call wave_from_disk(eigvec,1001,wave_name(1),7)
	sys_len=eigvec%sys_len
	env_len=eigvec%env_len
	! Read truncation operators from disk
	do i=trun_idx+1,sys_len
		call truns_from_disk(systruns(i),1001,i,.true.)
	end do
	do i=trun_idx+1,env_len
		call truns_from_disk(envtruns(i),1001,i,.false.)
	end do
	! Read basis from disk
	do i=2,sys_len
		call basis_from_disk(sys_bsm(i),1001,i,.true.)
	end do
	do i=2,env_len
		call basis_from_disk(env_bsm(i),1001,i,.false.)
	end do



        xi1=nx/4*ny+2-ny
        xj1=nx*3/4*ny+2
    ! correlations between xi1 point until xj1 point,  "1,3,5" first NN bond

        if(superc)then
        call Get_mSCOP_Operator_Cor0(eigvec,trun_idx,'S0',"Sc_cor_y.dat", xi1, xj1,0, 1)
        call Get_mSCOP_Operator_Cor0(eigvec,trun_idx,'S0',"Sc_cor_x.dat", xi1, xj1,0, 3)
        call Get_mSCOP_Operator_Cor0(eigvec,trun_idx,'S0',"Sc_cor_xy.dat", xi1, xj1,0, 5)
        endif



	Filename="elec_density.dat"
	call measure_operator_dia(ave_num_elec,num_elec,eigvec,trun_idx,Filename)

        if(measli)then
        spin_cor=0.0d0
        elec_cor=0.0d0
        num_elec_cor=0.0d0
filename='spin_cor11.dat'
        call Get_density_cor(spin_cor, st_sd,st_sd,eigvec,trun_idx,'B',Filename,xi1,xi1)
        open(251,file='spin_cor.dat')
        write(251,*)'Nx=',Nx, 'Ny=',Ny, 'J2=',jz(5,1), 'kept_max=',kept_max
        do i=1, num_site
        do j=1,num_site
        if(abs(spin_cor(i,j)).ne.0.0)then
        write(251,28)i, j,(j-1)/ny-(i-1)/ny,real(spin_cor(i, j)), imag(spin_cor(i,j))
        endif
        enddo
        enddo
28     format(3i8, 3f21.14)

filename='elec_cor11.dat'
        call Get_density_cor(elec_cor, st_elec_up,st_elec_down,eigvec,trun_idx,'F',Filename, xi1, xi1)
        open(251,file='elec_cor.dat')
        write(251,*)'Nx=',Nx, 'Ny=',Ny, 'J2=',jz(5,1), 'kept_max=',kept_max
        do i=1, num_site
        do j=1,num_site
        if(abs(elec_cor(i,j)).ne.0.0)then
        write(251,28)i, j,(j-1)/ny-(i-1)/ny,real(elec_cor(i, j)), imag(elec_cor(i,j))
        endif
        enddo
        enddo
filename='density_cor11.dat'
        call Get_density_cor(num_elec_cor, num_elec,num_elec,eigvec,trun_idx,'B',Filename, xi1, xi1)

        open(251,file='density_cor.dat')
        write(251,*)'Nx=',Nx, 'Ny=',Ny, 'J2=',jz(5,1), 'kept_max=',kept_max
        do i=1, num_site
        do j=1,num_site
        if(abs(num_elec_cor(i,j)).ne.0.0)then
        write(251,28)i, j,(j-1)/ny-(i-1)/ny,real(num_elec_cor(i, j)), real(ave_num_elec(i)),real(ave_num_elec(j))
        endif
        enddo
        enddo
        endif
        stop


        if(measbl)then
        xi1=1
        xj1=1
        Filename="spin_cor1.dat"
        call Get_oper_cor_ndia1(spin_cor,st_sd,st_sd,eigvec,trun_idx,'B',Filename)
        spin_cor1=spin_cor
        Filename="strfac_spin.dat"
        call Get_structure_factor_cut(spin_cor1,Filename, (xi1-1)/ny)

        Filename="density_cor1.dat"
        call Get_oper_cor_ndia1(num_elec_cor,num_elec, num_elec,eigvec,trun_idx,'B',Filename)
        spin_cor1=num_elec_cor
        Filename="strfac_spin.dat"
        call Get_structure_factor_cut(spin_cor1,Filename, (xi1-1)/ny)

        Filename="elec_cor1.dat"
        call Get_oper_cor_ndia1(elec_cor,st_elec_up,st_elec_down,eigvec,trun_idx,'F',Filename)
        spin_cor1=elec_cor
        Filename="strfac_elec.dat"
        call Get_structure_factor_cut(spin_cor1,Filename, (xi1-1)/ny)
        endif

	!<6>: ========== Free space ===========
	do i=trun_idx+1,sys_len
		call deallocate_block(systruns(i))
	end do

	do i=trun_idx+1,env_len
		call deallocate_block(envtruns(i))
	end do

	do i=2,sys_len
		call deallocate_basis(sys_bsm(i))
	end do

	do i=2,env_len
		call deallocate_basis(env_bsm(i))
	end do
	call deallocate_wave(eigvec)

end subroutine Measurement



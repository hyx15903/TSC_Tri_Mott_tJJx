!======================================================================
!Get all <n_i> in system and environment blocks
subroutine measure_operator_dia(ave_dia,st_oper,wave,trun_idx,Filename)
	use pubdata
	implicit none

	character(len=30) :: Filename
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
	double precision,intent(inout) :: ave_dia(Num_site)
        double complex  ave1(num_site)
        double precision a11

	integer :: i,sys_len,env_len,dx,dy
	type(Total_Block) :: new_oper

	!<1>: Get general information
	sys_len=wave%sys_len
	env_len=wave%env_len
	ave_dia=0.0d0
        ave1=dcmplx(0.0, 0.0)

	!<2>: Get <n_i> in system block
	do i=1,sys_len-1
		call Get_sys_oper_dia(st_oper,i,new_oper,sys_len-1,trun_idx,.true.)
		call measure_sys_block(ave1(i),new_oper,sys_bsm(sys_len),wave)
		call deallocate_block(new_oper)
	end do

	!Get <n_i> in the middle of the system
	call measure_sys_site(ave1(sys_len),st_oper,sys_bsm(sys_len),wave)

	!<2>: Get <n_i> in environment block
	do i=1,env_len-1
		call Get_env_oper_dia(st_oper,i,new_oper,env_len-1,trun_idx,.true.)
		call measure_env_block(ave1(Num_site-i+1),new_oper,env_bsm(env_len),wave)
		call deallocate_block(new_oper)
	end do

	!Get <n_i> in the middle of the system
	call measure_env_site(ave1(Num_site-env_len+1),st_oper,env_bsm(env_len),wave)

	!<3>: Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max

        ave_dia=ave1
        a11=0.0d0
	do i=1,Num_site
		dx=Lattice(1,i)
		dy=Lattice(2,i)
        a11=a11+dreal(ave1(i))
		write(10,101) (i-1)/ny+1,mod(i-1, ny)+1,dreal(ave1(i)), a11
	end do
	write(10,*)
	close(10)
111 format(I4,1X,I2,1X,2F16.12)
101     format(2i6, 3f14.8)

end subroutine measure_operator_dia
!========================================================================
subroutine measure_operator_ndia21(ave_dia2,st_oper,wave,trun_idx,Filename)
	use pubdata
	implicit none

	character(len=30) :: Filename
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: ave_dia2(Num_site)

	integer :: i,sys_len,env_len,dx,dy
	type(Total_Block) :: tmp_oper2,new_oper,st_oper2

	!<1>: Get general information
	sys_len=wave%sys_len
	env_len=wave%env_len
	call block_mul_block_ndia(st_oper,'T',st_oper,'N',cone,st_oper2,st_basis)
	!!!ave_dia2=0.0d0

	!<2>: Get <n_i^2> in system block
	do i=1,sys_len-1
		call Get_sys_oper_dia(st_oper2,i,new_oper,sys_len-1,trun_idx,.true.)
		call measure_sys_block(ave_dia2(i),new_oper,sys_bsm(sys_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_sys_site(ave_dia2(sys_len),st_oper2,sys_bsm(sys_len),wave)

	!<2>: Get <n_i^2> in environment block
	do i=1,env_len-1
		call Get_env_oper_dia(st_oper2,i,new_oper,env_len-1,trun_idx,.true.)
		call measure_env_block(ave_dia2(Num_site-i+1),new_oper,env_bsm(env_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_env_site(ave_dia2(Num_site-env_len+1),st_oper2,env_bsm(env_len),wave)
	call deallocate_block(st_oper2)

	!<3>: Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do i=1,Num_site
		dx=Lattice(1,i)
		dy=Lattice(2,i)
		write(10,*) dx,dy,i, ave_dia2(i)
	end do
	write(10,*)
	close(10)
111 format(I4,1X,I2,1X,F16.12)

end subroutine measure_operator_ndia21


!============================================================================
!Get all <n_i^2> in system and environment blocks
subroutine measure_operator_ndia2(ave_dia2,st_oper,wave,trun_idx,Filename)
	use pubdata
	implicit none

	character(len=30) :: Filename
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: ave_dia2(Num_site)

	integer :: i,sys_len,env_len,dx,dy
	type(Total_Block) :: tmp_oper2,new_oper,st_oper2

	!<1>: Get general information
	sys_len=wave%sys_len
	env_len=wave%env_len
	call block_mul_block_ndia(st_oper,'N',st_oper,'T',cone,st_oper2,st_basis)
	!ave_dia2=0.0d0


	!<2>: Get <n_i^2> in system block
	do i=1,sys_len-1
		call Get_sys_oper_dia(st_oper2,i,new_oper,sys_len-1,trun_idx,.true.)
		call measure_sys_block(ave_dia2(i),new_oper,sys_bsm(sys_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_sys_site(ave_dia2(sys_len),st_oper2,sys_bsm(sys_len),wave)

	!<2>: Get <n_i^2> in environment block
	do i=1,env_len-1
		call Get_env_oper_dia(st_oper2,i,new_oper,env_len-1,trun_idx,.true.)
		call measure_env_block(ave_dia2(Num_site-i+1),new_oper,env_bsm(env_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_env_site(ave_dia2(Num_site-env_len+1),st_oper2,env_bsm(env_len),wave)
	call deallocate_block(st_oper2)

	!<3>: Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do i=1,Num_site
		dx=Lattice(1,i)
		dy=Lattice(2,i)
		write(10,*) dx,dy,i, ave_dia2(i)
	end do
	write(10,*)
	close(10)
111 format(I4,1X,I2,1X,F16.12)

end subroutine measure_operator_ndia2


!============================================================================


!========================================================================
!Get all <n_i^2> in system and environment blocks
subroutine measure_operator_dia2(ave_dia2,st_oper,wave,trun_idx,Filename)
	use pubdata
	implicit none

	character(len=30) :: Filename
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: ave_dia2(Num_site)

	integer :: i,sys_len,env_len,dx,dy
	type(Total_Block) :: tmp_oper2,new_oper,st_oper2

	!<1>: Get general information
	sys_len=wave%sys_len
	env_len=wave%env_len
	call block_mul_block_dia(st_oper,'N',st_oper,'N',cone,st_oper2)
	ave_dia2=0.0d0

	!<2>: Get <n_i^2> in system block
	do i=1,sys_len-1
		call Get_sys_oper_dia(st_oper2,i,new_oper,sys_len-1,trun_idx,.true.)
		call measure_sys_block(ave_dia2(i),new_oper,sys_bsm(sys_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_sys_site(ave_dia2(sys_len),st_oper2,sys_bsm(sys_len),wave)

	!<2>: Get <n_i^2> in environment block
	do i=1,env_len-1
		call Get_env_oper_dia(st_oper2,i,new_oper,env_len-1,trun_idx,.true.)
		call measure_env_block(ave_dia2(Num_site-i+1),new_oper,env_bsm(env_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_env_site(ave_dia2(Num_site-env_len+1),st_oper2,env_bsm(env_len),wave)
	call deallocate_block(st_oper2)

	!<3>: Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do i=1,Num_site
		dx=Lattice(1,i)
		dy=Lattice(2,i)
		write(10,111) dx,dy,ave_dia2(i)
	end do
	write(10,*)
	close(10)
111 format(I4,1X,I2,1X,F16.12)

end subroutine measure_operator_dia2


!============================================================================
!Get all diagonal operator correlatioins using memory directly
!============================================================================
subroutine Get_oper_cor_dia(oper_dia,st_oper,wave,trun_idx,Filename)
	use pubdata
	implicit none

	character(len=30) :: Filename
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: oper_dia(Num_site,Num_site)

	complex*16 :: value
	integer :: i,j,sys_len,env_len,ti,tj
	type(Total_Basis) :: sys_bs,env_bs
	type(Total_Block) :: st_oper2,new_oper

	!<1>: For information saving to the disk
	open(10,file="oper_dia_tmp.dat",position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max,Filename
	109	format(A3,I4,1X,A3,I2,1X,A2,F6.3,1X,A2,F6.3,1X,A5,I3,1X,A7,I3,1X,A3,1X,I4,1X,A3,1X,I4,1X,A30)
	close(10)


	!<2>: Get general information
	sys_len=wave%sys_len
	env_len=wave%env_len
	call basis_transfer(sys_bsm(sys_len),sys_bs)
	call basis_transfer(env_bsm(env_len),env_bs)
	call block_mul_block_dia(st_oper,'N',st_oper,'N',cone,st_oper2)
	oper_dia=0.0d0

	!On-site density correlation
	!<2>: Get <n_i^2> in system block
	do i=1,sys_len-1
		call Get_sys_oper_dia(st_oper2,i,new_oper,sys_len-1,trun_idx,.true.)
		call measure_sys_block(oper_dia(i,i),new_oper,sys_bsm(sys_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_sys_site(oper_dia(sys_len,sys_len),st_oper2,sys_bsm(sys_len),wave)

	!<2>: Get <n_i^2> in environment block
	do i=1,env_len-1
		call Get_env_oper_dia(st_oper2,i,new_oper,env_len-1,trun_idx,.true.)
		call measure_env_block(oper_dia(Num_site-i+1,Num_site-i+1),new_oper,env_bsm(env_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_env_site(oper_dia(Num_site-env_len+1,Num_site-env_len+1),st_oper2,env_bsm(env_len),wave)
	call deallocate_block(st_oper2)
	
	!Save to the disk
	do i=1,Num_site
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) i,i,oper_dia(i,i)
		close(10)
	end do


	!<2>: Get sysopers and envopers
	do i=1,sys_len-1
		call block_transfer(st_oper,sysopers(i))
	end do

	do i=1,env_len-1
		call block_transfer(st_oper,envopers(i))
	end do


	!<3>: Get sys_bl_bl and env_bl_bl operator correlation
	call Get_sys_oper_cor_dia(oper_dia,st_oper,wave,trun_idx)
	call Get_env_oper_cor_dia(oper_dia,st_oper,wave,trun_idx)


	!<4>: Get sys_bl_st and sys_bl_env_st
	do i=1,sys_len-1
		call sys_block_site_cor_dia(oper_dia(i,sys_len),sysopers(i),st_oper,sys_bs,wave)
		oper_dia(sys_len,i)=dconjg(oper_dia(i,sys_len))

		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) i,sys_len,oper_dia(i,sys_len)
		write(10,*) sys_len,i,oper_dia(sys_len,i)
		close(10)

		call sys_block_env_site_cor_dia(oper_dia(i,Num_site-env_len+1),sysopers(i),st_oper,sys_bs,env_bs,wave)
		oper_dia(Num_site-env_len+1,i)=oper_dia(i,Num_site-env_len+1)

		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) i,Num_site-env_len+1,oper_dia(i,Num_site-env_len+1)
		write(10,*) Num_site-env_len+1,i,oper_dia(Num_site-env_len+1,i)
		close(10)
	end do

	!<5>: Get env_bl_st and sys_st_env_bl
	do i=1,env_len-1
		call env_block_site_cor_dia(oper_dia(Num_site-i+1,Num_site-env_len+1),envopers(i),st_oper,env_bs,wave)
		oper_dia(Num_site-env_len+1,Num_site-i+1)=oper_dia(Num_site-i+1,Num_site-env_len+1)

		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) Num_site-i+1,Num_site-env_len+1,real(oper_dia(Num_site-i+1,Num_site-env_len+1))
		write(10,*) Num_site-env_len+1,Num_site-i+1,real(oper_dia(Num_site-env_len+1,Num_site-i+1))
		close(10)

		call sys_site_env_block_cor_dia(oper_dia(sys_len,Num_site-i+1),st_oper,envopers(i),sys_bs,env_bs,wave)
		oper_dia(Num_site-i+1,sys_len)=dconjg(oper_dia(sys_len,Num_site-i+1))

		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) Num_site-i+1,sys_len,real(oper_dia(Num_site-i+1,sys_len))
		write(10,*) sys_len,Num_site-i+1,real(oper_dia(sys_len,Num_site-i+1))
		close(10)
	end do

	!<6>: Get sys_bl_env_bl
	do i=1,sys_len-1
	do j=1,env_len-1
		call sys_block_env_block_cor_dia(oper_dia(i,Num_site-j+1),sysopers(i),envopers(j),sys_bs,env_bs,wave)
		oper_dia(Num_site-j+1,i)=dconjg(oper_dia(i,Num_site-j+1))
			
		!Save to disk
		open(10,file="oper_dia_tmp.dat",position='append')
		write(10,*) i,Num_site-j+1,real(oper_dia(i,Num_site-j+1))
		write(10,*) Num_site-j+1,i,real(oper_dia(Num_site-j+1,i))
		close(10)
	end do
	end do

	call sys_site_env_site_cor_dia(oper_dia(sys_len,Num_site-env_len+1),st_oper,st_oper,sys_bs,env_bs,wave)
	oper_dia(Num_site-env_len+1,sys_len)=dconjg(oper_dia(sys_len,Num_site-env_len+1))

	!Save to disk
	open(10,file="oper_dia_tmp.dat",position='append')
	write(10,*) sys_len,Num_site-env_len+1,real(oper_dia(sys_len,Num_site-env_len+1))
	write(10,*) Num_site-env_len+1,sys_len,real(oper_dia(Num_site-env_len+1,sys_len))
	write(10,*)
	close(10)


	!Free space
	do i=1,sys_len-1
		call deallocate_block(sysopers(i))
	end do

	do i=1,env_len-1
		call deallocate_block(envopers(i))
	end do
	call deallocate_basis(sys_bs)
	call deallocate_basis(env_bs)


	!<5>: Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do i=1,Num_site
		do j=1,Num_site
			write(10,111) i,j,real(oper_dia(i,j)),imag(oper_dia(I,j))
		end do
	end do
	write(10,*)
	close(10)
111 format(I4,1X,I4,1X,2F16.12)

end subroutine Get_oper_cor_dia
!==============================================================
!Get all density correlcations in system block using memory
!==============================================================
subroutine Get_sys_oper_cor_dia(oper_dia,st_oper,wave,trun_idx)
	use pubdata
	implicit none

	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
        complex*16,intent(inout) :: oper_dia(Num_site,Num_site)

	integer :: i,j,x,y,sys_len
	type(Total_Block) :: new_oper1,new_oper2
	type(Total_Block) :: tmp_oper,mid_oper,sys_oper


	!<1>: For general information
	sys_len=wave%sys_len

	!Get g(i,j) for (i<j)
	do j=2,sys_len-1
		call update_site_dia(st_oper,sys_bsm(j),new_oper2)

		!Update sysopers(sys_len-1)
		if(j==sys_len-1) then
			call deallocate_block(sysopers(j))
			if(j<=trun_idx) then
				call update_site_dia(st_oper,sys_bsm(j),sysopers(j))
			else
				call update_site_dia(st_oper,sys_bsm(j),mid_oper)
				call update_trun_dia(mid_oper,sysopers(j),systruns(j))
				call deallocate_block(mid_oper)
			endif
		endif

		do i=1,j-1
			if(i==j-1) then
				if(i>1) then !(i>1)
					call update_site_dia(st_oper,sys_bsm(i),mid_oper)
					if(i<=trun_idx) then
						call block_transfer(mid_oper,tmp_oper)
					else
						call update_trun_dia(mid_oper,tmp_oper,systruns(i))
					endif
					call update_block_dia(tmp_oper,sys_bsm(j),new_oper1)
					call deallocate_block(mid_oper)
					call deallocate_block(tmp_oper)
				else !(i==1)
					call update_block_dia(sysopers(i),sys_bsm(j),new_oper1)
				endif
			else
				call update_block_dia(sysopers(i),sys_bsm(j),new_oper1)
			endif

			!Update sysopers(i)
			call deallocate_block(sysopers(i))
			if(j<=trun_idx) then
				call block_transfer(new_oper1,sysopers(i))
			else
				call update_trun_dia(new_oper1,sysopers(i),systruns(j))
			endif
			
			!Multiply oper(i) and oper(j)
			call block_mul_block_dia(new_oper1,'N',new_oper2,'N',cone,tmp_oper)
			call deallocate_block(new_oper1)

			!Update tmp_oper(i,j) to new configuration
			!For (j<sys_len-1)
			if(j<sys_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,mid_oper)
				else
					call update_trun_dia(tmp_oper,mid_oper,systruns(j))
				endif
				call change_sys_oper_dia(mid_oper,j,sys_oper,sys_len-1,trun_idx,.true.)
				call deallocate_block(mid_oper)
			
			!For (j==sys_len-1)
			else if(j==sys_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,sys_oper)
				else
					call update_trun_dia(tmp_oper,sys_oper,systruns(sys_len-1))
				endif
			endif
			call deallocate_block(tmp_oper)

			!Measurement
			call measure_sys_block(oper_dia(i,j),sys_oper,sys_bsm(sys_len),wave)
			oper_dia(j,i)=dconjg(oper_dia(i,j))

			!Save to disk
			open(10,file="oper_dia_tmp.dat",position='append')
			write(10,111) i,j,real(oper_dia(i,j))
			write(10,111) j,i,real(oper_dia(j,i))
			close(10)
			
			call deallocate_block(sys_oper)
		end do

		call deallocate_block(new_oper2)
	end do
111 format(I4,1X,I4,1X,F16.12)
end subroutine Get_sys_oper_cor_dia

!====================================================================
!Get all density correlcations in environment block using memory
!====================================================================
subroutine Get_env_oper_cor_dia(oper_dia,st_oper,wave,trun_idx)
	use pubdata
	implicit none

	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: oper_dia(Num_site,Num_site)

	integer :: i,j,x,y,env_len
	type(Total_Block) :: new_oper1,new_oper2
	type(Total_Block) :: tmp_oper,mid_oper,env_oper


	!<1>: For general information
	env_len=wave%env_len

	!Get g(i,j) for (i<j)
	do j=2,env_len-1
		call update_site_dia(st_oper,env_bsm(j),new_oper2)

		!Update envopers(env_len-1)
		if(j==env_len-1) then
			call deallocate_block(envopers(j))
			if(j<=trun_idx) then
				call update_site_dia(st_oper,env_bsm(j),envopers(j))
			else
				call update_site_dia(st_oper,env_bsm(j),mid_oper)
				call update_trun_dia(mid_oper,envopers(j),envtruns(j))
				call deallocate_block(mid_oper)
			endif
		endif

		do i=1,j-1
			if(i==j-1) then
				if(i>1) then !(i>1)
					call update_site_dia(st_oper,env_bsm(i),mid_oper)
					if(i<=trun_idx) then
						call block_transfer(mid_oper,tmp_oper)
					else
						call update_trun_dia(mid_oper,tmp_oper,envtruns(i))
					endif
					call update_block_dia(tmp_oper,env_bsm(j),new_oper1)
					call deallocate_block(mid_oper)
					call deallocate_block(tmp_oper)
				else !(i==1)
					call update_block_dia(envopers(i),env_bsm(j),new_oper1)
				endif
			else
				call update_block_dia(envopers(i),env_bsm(j),new_oper1)
			endif

			!Update envopers(i)
			call deallocate_block(envopers(i))
			if(j<=trun_idx) then
				call block_transfer(new_oper1,envopers(i))
			else
				call update_trun_dia(new_oper1,envopers(i),envtruns(j))
			endif
			
			!Multiply oper(i) and oper(j)
			call block_mul_block_dia(new_oper1,'N',new_oper2,'N',cone,tmp_oper)
			call deallocate_block(new_oper1)

			!Update tmp_oper(i,j) to new configuration
			!For (j<env_len-1)
			if(j<env_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,mid_oper)
				else
					call update_trun_dia(tmp_oper,mid_oper,envtruns(j))
				endif
				call change_env_oper_dia(mid_oper,j,env_oper,env_len-1,trun_idx,.true.)
				call deallocate_block(mid_oper)
			
			!For (j==env_len-1)
			else if(j==env_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,env_oper)
				else
					call update_trun_dia(tmp_oper,env_oper,envtruns(env_len-1))
				endif
			endif
			call deallocate_block(tmp_oper)

			!Measurement
			call measure_env_block(oper_dia(Num_site-i+1,Num_site-j+1),env_oper,env_bsm(env_len),wave)
			oper_dia(Num_site-j+1,Num_site-i+1)=dconjg(oper_dia(Num_site-i+1,Num_site-j+1))

			!Save to disk
			open(10,file="oper_dia_tmp.dat",position='append')
			write(10,111) Num_site-i+1,Num_site-j+1,real(oper_dia(Num_site-i+1,Num_site-j+1))
			write(10,111) Num_site-j+1,Num_site-i+1,real(oper_dia(Num_site-j+1,Num_site-i+1))
			close(10)
			
			call deallocate_block(env_oper)
		end do

		call deallocate_block(new_oper2)
	end do
111 format(I4,1X,I4,1X,F16.12)

end subroutine Get_env_oper_cor_dia
!==========================================================================
!Get the correlation functions between arbitrary two sites: (idx1<=idx2)
!==========================================================================
subroutine Get_oper_cor_bond_dia(value,oper1,idx1,oper2,idx2,trun_idx,wave)
	use pubdata
	implicit none

	integer,intent(in) :: idx1,idx2,trun_idx
	type(Total_Block),intent(in) :: oper1,oper2
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,j
	complex*16 :: value_sd,value_sz
	integer :: idx_one,idx_two,one,two,sys_len,env_len
	type(Total_Block) :: st_oper1,st_oper2,tmp_oper
	type(Total_Block) :: mid_oper,oper,new_oper
	type(Total_Block) :: sys_sz,sys_sd,env_sz,env_sd
	logical :: Sys_Flag,Env_Flag


	!<1>: For general information
	sys_len=wave%sys_len
	env_len=wave%env_len

	!Get proper operator for specific site
	if(idx1<=idx2) then
		idx_one=idx1
		idx_two=idx2
		call block_transfer(oper1,st_oper1)
		call block_transfer(oper2,st_oper2)
	else !For (idx1>idx2)
		idx_one=idx2
		idx_two=idx1
		call block_transfer(oper1,st_oper2)
		call block_transfer(oper2,st_oper1)
	endif

	value=0.0d0
	!Get truncation flag
	Sys_Flag=.true.
	if(trun_idx==(sys_len-1)) then
		Sys_Flag=.false.
	endif
	Env_Flag=.true.
	if(trun_idx==(env_len-1)) then
		Env_Flag=.false.
	endif


	!<2>: For (idx_one<sys_len)
	if(idx_one<sys_len) then
		!<2-1>: For (idx_two<sys_len)
		if(idx_two<sys_len) then
			call Get_sys_oper_dia(st_oper1,idx_one,tmp_oper,idx_two,trun_idx,.false.)
			call Get_sys_oper_dia(st_oper2,idx_two,mid_oper,idx_two,trun_idx,.false.)
			call block_mul_block_dia(tmp_oper,'N',mid_oper,'N',cone,new_oper)
			call deallocate_block(tmp_oper)
			call deallocate_block(mid_oper)

			!Change to (sys_len-1)'s configuration
			if(idx_two<=trun_idx) then
				call block_transfer(new_oper,tmp_oper)
			else
				call update_trun_dia(new_oper,tmp_oper,systruns(idx_two))
			endif
			call deallocate_block(new_oper)
			if(idx_two<sys_len-1) then
				call Change_sys_oper_dia(tmp_oper,idx_two,new_oper,sys_len-1,trun_idx,Sys_Flag)
				call deallocate_block(tmp_oper)
				call block_transfer(new_oper,tmp_oper)
				call deallocate_block(new_oper)
			endif

			call measure_sys_block(value,tmp_oper,sys_bsm(sys_len),wave)
			call deallocate_block(tmp_oper)

		!<2-2>: For (idx_two>=sys_len)
		else if(idx_two>=sys_len) then
			call Get_sys_oper_dia(st_oper1,idx_one,new_oper,sys_len-1,trun_idx,Sys_Flag)

			!For (idx_two=sys_len)
			if(idx_two==sys_len) then
				call sys_block_site_cor_dia(value,new_oper,st_oper2,sys_bsm(sys_len),wave)

			!For (idx_two=Num_site-env_len+1)
			else if(idx_two==Num_site-env_len+1) then
				call sys_block_env_site_cor_dia(value,new_oper,st_oper2,sys_bsm(sys_len),env_bsm(env_len),wave)

			!For (idx_two>Num_site-env_len+1)
			else if(idx_two>Num_site-env_len+1) then
				call Get_env_oper_dia(st_oper2,Num_site-idx_two+1,tmp_oper,env_len-1,trun_idx,Env_Flag)
				call sys_block_env_block_cor_dia(value,new_oper,tmp_oper,sys_bsm(sys_len),env_bsm(env_len),wave)
				call deallocate_block(tmp_oper)
			endif
			call deallocate_block(new_oper)
		endif


	!<3>: For (idx_one=sys_len)
	else if(idx_one==sys_len) then
		!For (idx_two=Num_site-env_len+1)
		if(idx_two==Num_site-env_len+1) then
			call sys_site_env_site_cor_dia(value,st_oper1,st_oper2,sys_bsm(sys_len),env_bsm(env_len),wave)

		!For (idx_two>Num_site-env_len+1)
		else if(idx_two>Num_site-env_len+1) then
			call Get_env_oper_dia(st_oper2,Num_site-idx_two+1,new_oper,env_len-1,trun_idx,Env_Flag)
			call sys_site_env_block_cor_dia(value,st_oper1,new_oper,sys_bsm(sys_len),env_bsm(env_len),wave)
			call deallocate_block(new_oper)
		endif


	!<4>: For (idx_one>sys_len)
	else if(idx_one>sys_len) then
		one=Num_site-idx_two+1
		two=Num_site-idx_one+1

		call block_transfer(st_oper1,tmp_oper)
		call deallocate_block(st_oper1)
		call block_transfer(st_oper2,st_oper1)
		call deallocate_block(st_oper2)
		call block_transfer(tmp_oper,st_oper2)
		call deallocate_block(tmp_oper)

		!For (two<env_len)
		if(two<env_len) then
			call Get_env_oper_dia(st_oper1,one,tmp_oper,two,trun_idx,.false.)
			call Get_env_oper_dia(st_oper2,two,mid_oper,two,trun_idx,.false.)
			call block_mul_block_dia(tmp_oper,'N',mid_oper,'N',cone,new_oper)
			call deallocate_block(tmp_oper)
			call deallocate_block(mid_oper)

			!Change to (env_len-1)'s configuration
			if(two<=trun_idx) then
				call block_transfer(new_oper,tmp_oper)
			else
				call update_trun_dia(new_oper,tmp_oper,envtruns(two))
			endif
			call deallocate_block(new_oper)

			if(two<env_len-1) then
				call Change_env_oper_dia(tmp_oper,two,new_oper,env_len-1,trun_idx,Env_Flag)
				call deallocate_block(tmp_oper)
				call block_transfer(new_oper,tmp_oper)
				call deallocate_block(new_oper)
			endif

			call measure_env_block(value,tmp_oper,env_bsm(env_len),wave)
			call deallocate_block(tmp_oper)

		!For (two=env_len)
		else if(two==env_len) then
			call Get_env_oper_dia(st_oper1,one,new_oper,env_len-1,trun_idx,Env_Flag)
			call env_block_site_cor_dia(value,new_oper,st_oper2,env_bsm(env_len),wave)
			call deallocate_block(new_oper)
		endif
	endif

	call deallocate_block(st_oper1)
	call deallocate_block(st_oper2)

end subroutine Get_oper_cor_bond_dia



!==========================================================
!Measure <n_i> and <n_i^2> in system block
!==========================================================
!Get <n_i> in the system block
subroutine measure_sys_block(value,sys_bl,sys_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	integer :: sys_num_up,sys_num_down,env_dim
	integer :: bl_num_up,bl_num_down,spos,sdim
	complex*16,allocatable :: mid(:,:)
        double complex coef1

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_dim=wave%sub(i)%env_dim

		do k=1,sys_bs%num
			if(sys_bs%sub(k)%new_num_up==sys_num_up) then
			if(sys_bs%sub(k)%new_num_down==sys_num_down) then
				do x=1,sys_bs%sub(k)%num
					bl_num_up=sys_bs%sub(k)%sub(x)%bl_num_up
					bl_num_down=sys_bs%sub(k)%sub(x)%bl_num_down
					spos=sys_bs%sub(k)%sub(x)%spos
					sdim=sys_bs%sub(k)%sub(x)%sdim
                 coef1=1.0d0/dsqrt(1.0d0+bl_num_down)


					do y=1,sys_bl%num
						if(sys_bl%sub(y)%num_up==bl_num_up) then
						if(sys_bl%sub(y)%num_down==bl_num_down) then
							allocate(mid(sdim,env_dim))

                                                                if(realcode)then
                                   call DGEMM('N','N',sdim,env_dim,sdim,coef1,sys_bl%sub(y)%mat,sdim&
                                                                        &,wave%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim,0.0d0,mid,sdim)


                                                                else
							call ZGEMM('N','N',sdim,env_dim,sdim,coef1,sys_bl%sub(y)%mat,sdim&
									&,wave%sub(i)%vec(spos+1:spos+sdim,1:env_dim),sdim,czero,mid,sdim)
                                                                endif
							do z=1,env_dim
								value=value+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,z),mid(1:sdim,z))
							end do

							deallocate(mid)
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do

end subroutine measure_sys_block
!Get <n_i> in the system block
subroutine measure_sys_site(value,sys_st,sys_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs
	type(Total_Block),intent(in) :: sys_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	integer :: sys_num_up,sys_num_down,env_dim
	integer :: st_num_up,st_num_down,spos,sdim
	complex*16 :: yita
        complex*16 coef1

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_dim=wave%sub(i)%env_dim

		do k=1,sys_bs%num
			if(sys_bs%sub(k)%new_num_up==sys_num_up) then
			if(sys_bs%sub(k)%new_num_down==sys_num_down) then
				do x=1,sys_bs%sub(k)%num
					st_num_up=sys_bs%sub(k)%sub(x)%st_num_up
					st_num_down=sys_bs%sub(k)%sub(x)%st_num_down
					spos=sys_bs%sub(k)%sub(x)%spos
					sdim=sys_bs%sub(k)%sub(x)%sdim

					do y=1,sys_st%num
						if(sys_st%sub(y)%num_up==st_num_up) then
						if(sys_st%sub(y)%num_down==st_num_down) then

                 coef1=1.0d0/dsqrt(1.0d0+st_num_down)
							yita=0.0d0
							do z=1,env_dim
								yita=yita+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,z)&
													&,wave%sub(i)%vec(spos+1:spos+sdim,z))
							end do

							value=value+yita*sys_st%sub(y)%mat(1,1)*coef1
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do
end subroutine measure_sys_site


!==========================================================
!Measure <n_i> and <n_i^2> in environment block
!==========================================================

!==========================================================
!Measure <n_i> and <n_i^2> in environment block
!==========================================================
!Get <n_i> in the environment block
subroutine measure_env_block(value,env_bl,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	integer :: env_num_up,env_num_down,sys_dim
	integer :: bl_num_up,bl_num_down,spos,sdim
	complex*16,allocatable :: mid(:,:)
        complex*16 coef1

	value=0.0d0
	do i=1,wave%num
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim

		do k=1,env_bs%num
			if(env_bs%sub(k)%new_num_up==env_num_up) then
			if(env_bs%sub(k)%new_num_down==env_num_down) then
				do x=1,env_bs%sub(k)%num
					bl_num_up=env_bs%sub(k)%sub(x)%bl_num_up
					bl_num_down=env_bs%sub(k)%sub(x)%bl_num_down
					spos=env_bs%sub(k)%sub(x)%spos
					sdim=env_bs%sub(k)%sub(x)%sdim

					do y=1,env_bl%num
						if(env_bl%sub(y)%num_up==bl_num_up) then
						if(env_bl%sub(y)%num_down==bl_num_down) then
							allocate(mid(sys_dim,sdim))
         coef1=1.0d0/dsqrt(1.0d0+bl_num_down)



                if(realcode)then
							call DGEMM('N','N',sys_dim,sdim,sdim,coef1&
									&,wave%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim&
									&,env_bl%sub(y)%mat,sdim,czero,mid,sys_dim)

                else
							call ZGEMM('N','N',sys_dim,sdim,sdim,coef1&
									&,wave%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim&
									&,env_bl%sub(y)%mat,sdim,czero,mid,sys_dim)
                        endif
							do z=1,sdim
								value=value+dot_product(mid(1:sys_dim,z),wave%sub(i)%vec(1:sys_dim,spos+z))
							end do

							deallocate(mid)
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do

end subroutine measure_env_block
!Get <n_i> in the environment block
subroutine measure_env_site(value,env_st,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: env_bs
	type(Total_Block),intent(in) :: env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	integer :: env_num_up,env_num_down,sys_dim
	integer :: st_num_up,st_num_down,spos,sdim
	complex*16 :: yita
        complex*16 coef1

	value=0.0d0
	do i=1,wave%num
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim

		do k=1,env_bs%num
			if(env_bs%sub(k)%new_num_up==env_num_up) then
			if(env_bs%sub(k)%new_num_down==env_num_down) then
				do x=1,env_bs%sub(k)%num
					st_num_up=env_bs%sub(k)%sub(x)%st_num_up
					st_num_down=env_bs%sub(k)%sub(x)%st_num_down
					spos=env_bs%sub(k)%sub(x)%spos
					sdim=env_bs%sub(k)%sub(x)%sdim

					do y=1,env_st%num
						if(env_st%sub(y)%num_up==st_num_up) then
						if(env_st%sub(y)%num_down==st_num_down) then

        
                                                                coef1=1.0d0/dsqrt(1.0d0+st_num_down)

							yita=0.0d0
							do z=1,sdim
								yita=yita+dot_product(wave%sub(i)%vec(1:sys_dim,spos+z)&
													&,wave%sub(i)%vec(1:sys_dim,spos+z))
							end do

							value=value+yita*env_st%sub(y)%mat(1,1)*coef1
							goto 101
						endif
						endif
					end do
					101 continue
				end do
			endif
			endif
		end do
	end do

end subroutine measure_env_site




!Note: value=<psi|n_1^+.n_2|psi> in system block
subroutine sys_block_site_cor_dia(value,sys_bl,sys_st,sys_bs,wave)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: sys_bl,sys_st
	type(Total_Basis),intent(in) :: sys_bs
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: value

	integer :: i,k,x,y,z
	logical :: bl_flag,st_flag
	integer :: bl_id,st_id,spos,sdim
	integer :: sys_num_up,sys_num_down,env_dim
	integer :: bl_num_up,bl_num_down,st_num_up,st_num_down
	double complex,allocatable :: mid(:,:)
	double complex :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim
					
					bl_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==bl_num_up) then
						if(sys_bl%sub(k)%num_down==bl_num_down) then
							bl_flag=.true.
							bl_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==st_num_up) then
						if(sys_st%sub(k)%num_down==st_num_down) then
							st_flag=.true.
							st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					if(bl_flag.and.st_flag) then
						coef=sys_st%sub(st_id)%mat(1,1)
                           coef=coef*dsqrt((1.0d0+bl_num_down)*(1.0d0+st_num_down))

						allocate(mid(sdim,env_dim))
                        if(realcode)then
						call DGEMM('N','N',sdim,env_dim,sdim,coef,sys_bl%sub(bl_id)%mat&
								&,sdim,wave%sub(i)%vec(spos+1:spos+sdim,1:env_dim)&
								&,sdim,0.0d0,mid,sdim)
                                else



						call zGEMM('N','N',sdim,env_dim,sdim,coef,sys_bl%sub(bl_id)%mat&
								&,sdim,wave%sub(i)%vec(spos+1:spos+sdim,1:env_dim)&
								&,sdim,czero,mid,sdim)


                                endif

						do k=1,env_dim
							value=value+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,k),mid(1:sdim,k))
						end do
						deallocate(mid)
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_site_cor_dia



!Note: value=<psi|n_1^+.n_2|psi> in environment block
subroutine env_block_site_cor_dia(value,env_bl,env_st,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Block),intent(in) :: env_bl,env_st
	type(Total_Basis),intent(in) :: env_bs
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z
	logical :: bl_flag,st_flag
	integer :: bl_id,st_id,spos,sdim
	integer :: env_num_up,env_num_down,sys_dim
	integer :: bl_num_up,bl_num_down,st_num_up,st_num_down
	complex*16,allocatable :: mid(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim

		do x=1,env_bs%num
			if(env_bs%sub(x)%new_num_up==env_num_up) then
			if(env_bs%sub(x)%new_num_down==env_num_down) then
				do y=1,env_bs%sub(x)%num
					bl_num_up=env_bs%sub(x)%sub(y)%bl_num_up
					bl_num_down=env_bs%sub(x)%sub(y)%bl_num_down
					st_num_up=env_bs%sub(x)%sub(y)%st_num_up
					st_num_down=env_bs%sub(x)%sub(y)%st_num_down
					spos=env_bs%sub(x)%sub(y)%spos
					sdim=env_bs%sub(x)%sub(y)%sdim
					
					bl_flag=.false.
					do k=1,env_bl%num
						if(env_bl%sub(k)%num_up==bl_num_up) then
						if(env_bl%sub(k)%num_down==bl_num_down) then
							bl_flag=.true.
							bl_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,env_st%num
						if(env_st%sub(k)%num_up==st_num_up) then
						if(env_st%sub(k)%num_down==st_num_down) then
							st_flag=.true.
							st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					if(bl_flag.and.st_flag) then
						coef=env_st%sub(st_id)%mat(1,1)
                           coef=coef*dsqrt((1.0d0+bl_num_down)*(1.0d0+st_num_down))

						allocate(mid(sys_dim,sdim))

                                                        if(realcode)then
						call DGEMM('N','N',sys_dim,sdim,sdim,coef&
								&,wave%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim&
								&,env_bl%sub(bl_id)%mat,sdim,0.0d0,mid,sys_dim)
                                                                else
						call ZGEMM('N','N',sys_dim,sdim,sdim,coef&
								&,wave%sub(i)%vec(1:sys_dim,spos+1:spos+sdim),sys_dim&
								&,env_bl%sub(bl_id)%mat,sdim,czero,mid,sys_dim)
                                                        endif
						do k=1,sdim
							value=value+dot_product(mid(1:sys_dim,k),wave%sub(i)%vec(1:sys_dim,spos+k))
						end do
						deallocate(mid)
					endif
				end do
			endif
			endif
		end do
	end do
end subroutine env_block_site_cor_dia




!Note: value=<psi|n_s^+.n_e|psi>
subroutine sys_block_env_block_cor_dia(value,sys_bl,env_bl,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16, allocatable :: mat1(:,:),mat2(:,:)
        complex*16 coef1

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								if(sys_flag.and.env_flag) then
									allocate(mat1(sdim,edim),mat2(sdim,edim))

         coef1=1.0d0/dsqrt((1.0d0+env_bl_num_down)*(1.0d0+sys_bl_num_down))
                                                        if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,cone&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat1,sdim)
									call DGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat2,sdim)
                                                        else 
									call ZGEMM('N','N',sdim,edim,edim,cone&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat1,sdim)
									call ZGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat2,sdim)
                                
                                                        endif

									do z=1,edim
										value=value+dot_product(mat1(:,z),mat2(:,z))*coef1
									end do
									deallocate(mat1,mat2)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_env_block_cor_dia



!Note: value=<psi|n_s^+.n_e|psi>
subroutine sys_block_env_site_cor_dia(value,sys_bl,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: env_st_num_up,env_st_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								if(sys_flag.and.env_flag) then
									allocate(mat(sdim,edim))
									coef=env_st%sub(env_id)%mat(1,1)

                                  coef=coef/dsqrt((1.0d0+sys_bl_num_down)*(1.0d0+env_st_num_down))
                                                                                if(realcode)then
									call DGEMM('N','N',sdim,edim,sdim,coef,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat,sdim)


                                                                                        else

									call ZGEMM('N','N',sdim,edim,sdim,coef,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat,sdim)

                                                                                endif

									do z=1,edim
										value=value+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,epos+z),mat(1:sdim,z))
									end do
									deallocate(mat)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_env_site_cor_dia
!Note: value=<psi|n_s^+.n_e|psi>


subroutine sys_site_env_block_cor_dia(value,sys_st,env_bl,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag
	integer :: sys_st_num_up,sys_st_num_down,sys_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								if(sys_flag.and.env_flag) then
									coef=sys_st%sub(sys_id)%mat(1,1)
									allocate(mat(sdim,edim))
                                                                        if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat,sdim)


                        else
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat,sdim)

                        endif

									do z=1,edim
										value=value+dot_product(mat(1:sdim,z),wave%sub(i)%vec(spos+1:spos+sdim,epos+z))
									end do
									deallocate(mat)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_site_env_block_cor_dia

!Note: value=<psi|n_s^+.n_e|psi>
subroutine sys_block_site_env_site_cor_dia(value,sys_bl,sys_st,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st,env_st
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,st_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: sys_st_num_up,sys_st_num_down,st_id
	integer :: env_st_num_up,env_st_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	double complex,allocatable :: mat(:,:)
	double complex :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							st_flag=.true.
							st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										env_flag=.true.
										env_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								if(sys_flag.and.st_flag.and.env_flag) then
									allocate(mat(sdim,edim))
									coef=sys_st%sub(st_id)%mat(1,1)*env_st%sub(env_id)%mat(1,1)
                                !!  coef=coef/dsqrt((1.0d0+env_bl_num_down)*(1.0d0+sys_st_num_down))

                        if(realcode)then

									call DGEMM('N','N',sdim,edim,sdim,coef,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat,sdim)
                        else

			call ZGEMM('N','N',sdim,edim,sdim,coef,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat,sdim)


                        endif

									do z=1,edim
										value=value+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,epos+z),mat(1:sdim,z))
									end do
									deallocate(mat)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_site_env_site_cor_dia


!Note: value=<psi|n_s^+.n_e|psi>
subroutine sys_site_env_block_site_cor_dia(value,sys_st,env_bl,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_bl,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,st_flag
	integer :: sys_st_num_up,sys_st_num_down,sys_id
	integer :: env_st_num_up,env_st_num_down,st_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								st_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										st_flag=.true.
										st_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								if(sys_flag.and.env_flag.and.st_flag) then
									coef=sys_st%sub(sys_id)%mat(1,1)*env_st%sub(st_id)%mat(1,1)
                                coef=sys_st%sub(st_id)%mat(1,1)*env_st%sub(env_id)%mat(1,1)

									allocate(mat(sdim,edim))

        if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat,sdim)

                                else
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat,sdim)

                                                                        endif
									do z=1,edim
										value=value+dot_product(mat(1:sdim,z),wave%sub(i)%vec(spos+1:spos+sdim,epos+z))
									end do
									deallocate(mat)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do
end subroutine sys_site_env_block_site_cor_dia

!Note: value=<psi|n_s.n_e|psi>



!Note: value=<psi|n_s.n_e|psi>
subroutine sys_site_env_site_cor_dia(value,sys_st,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag
	integer :: sys_st_num_up,sys_st_num_down,sys_id
	integer :: env_st_num_up,env_st_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16 :: coef,tmp_gl

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										env_flag=.true.
										env_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								if(sys_flag.and.env_flag) then
									coef=sys_st%sub(sys_id)%mat(1,1)*env_st%sub(env_id)%mat(1,1)

									tmp_gl=0.0d0
									do z=1,edim
										tmp_gl=tmp_gl+dot_product(wave%sub(i)%vec(spos+1:spos+sdim,epos+z),wave%sub(i)%vec(spos+1:spos+sdim,epos+z))
									end do
									value=value+tmp_gl*coef
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_site_env_site_cor_dia


!Note: value=<psi|n_s.n_st.n_e|psi>
subroutine sys_block_site_env_block_cor_dia(value,sys_bl,sys_st,env_bl,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st,env_bl
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,st_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: sys_st_num_up,sys_st_num_down,st_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat1(:,:),mat2(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							st_flag=.true.
							st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								if(sys_flag.and.st_flag.and.env_flag) then
									allocate(mat1(sdim,edim),mat2(sdim,edim))

									coef=sys_st%sub(st_id)%mat(1,1)
                                                if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat1,sdim)
									call DGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat2,sdim)
                                                        else
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat1,sdim)
									call ZGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat2,sdim)

                                                endif
									do z=1,edim
										value=value+dot_product(mat1(:,z),mat2(:,z))
									end do
									deallocate(mat1,mat2)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_site_env_block_cor_dia


!Note: value=<psi|n_s.n_st.n_e|psi>
subroutine sys_block_env_block_site_cor_dia(value,sys_bl,env_bl,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,st_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: env_st_num_up,env_st_num_down,st_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat1(:,:),mat2(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								st_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										st_flag=.true.
										st_id=k
										goto 102
									endif
									endif
								end do
								102 continue

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								if(sys_flag.and.st_flag.and.env_flag) then
									allocate(mat1(sdim,edim),mat2(sdim,edim))

									coef=env_st%sub(st_id)%mat(1,1)

                                if(realcode)then
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat1,sdim)
									call DGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat2,sdim)
                                        else
        
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat1,sdim)
									call ZGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat2,sdim)
                                        endif

									do z=1,edim
										value=value+dot_product(mat1(:,z),mat2(:,z))
									end do
									deallocate(mat1,mat2)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_env_block_site_cor_dia


!Note: value=<psi|n_s.n_st.n_e|psi>
subroutine sys_block_site_env_block_site_cor_dia(value,sys_bl,sys_st,env_bl,env_st,sys_bs,env_bs,wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,sys_st,env_bl,env_st
	type(Wavefunction),intent(in) :: wave
	complex*16,intent(inout) :: value

	integer :: i,k,x,y,z,m,n
	logical :: sys_flag,env_flag,sys_st_flag,env_st_flag
	integer :: sys_bl_num_up,sys_bl_num_down,sys_id
	integer :: sys_st_num_up,sys_st_num_down,sys_st_id
	integer :: env_st_num_up,env_st_num_down,env_st_id
	integer :: env_bl_num_up,env_bl_num_down,env_id
	integer :: sys_num_up,sys_num_down,sys_dim,spos,sdim
	integer :: env_num_up,env_num_down,env_dim,epos,edim
	complex*16,allocatable :: mat1(:,:),mat2(:,:)
	complex*16 :: coef

	value=0.0d0
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do x=1,sys_bs%num
			if(sys_bs%sub(x)%new_num_up==sys_num_up) then
			if(sys_bs%sub(x)%new_num_down==sys_num_down) then
				do y=1,sys_bs%sub(x)%num
					sys_bl_num_up=sys_bs%sub(x)%sub(y)%bl_num_up
					sys_bl_num_down=sys_bs%sub(x)%sub(y)%bl_num_down
					sys_st_num_up=sys_bs%sub(x)%sub(y)%st_num_up
					sys_st_num_down=sys_bs%sub(x)%sub(y)%st_num_down
					spos=sys_bs%sub(x)%sub(y)%spos
					sdim=sys_bs%sub(x)%sub(y)%sdim

					sys_flag=.false.
					do k=1,sys_bl%num
						if(sys_bl%sub(k)%num_up==sys_bl_num_up) then
						if(sys_bl%sub(k)%num_down==sys_bl_num_down) then
							sys_flag=.true.
							sys_id=k
							goto 101
						endif
						endif
					end do
					101 continue

					sys_st_flag=.false.
					do k=1,sys_st%num
						if(sys_st%sub(k)%num_up==sys_st_num_up) then
						if(sys_st%sub(k)%num_down==sys_st_num_down) then
							sys_st_flag=.true.
							sys_st_id=k
							goto 102
						endif
						endif
					end do
					102 continue

					do m=1,env_bs%num
						if(env_bs%sub(m)%new_num_up==env_num_up) then
						if(env_bs%sub(m)%new_num_down==env_num_down) then
							do n=1,env_bs%sub(m)%num
								env_bl_num_up=env_bs%sub(m)%sub(n)%bl_num_up
								env_bl_num_down=env_bs%sub(m)%sub(n)%bl_num_down
								env_st_num_up=env_bs%sub(m)%sub(n)%st_num_up
								env_st_num_down=env_bs%sub(m)%sub(n)%st_num_down
								epos=env_bs%sub(m)%sub(n)%spos
								edim=env_bs%sub(m)%sub(n)%sdim

								env_st_flag=.false.
								do k=1,env_st%num
									if(env_st%sub(k)%num_up==env_st_num_up) then
									if(env_st%sub(k)%num_down==env_st_num_down) then
										env_st_flag=.true.
										env_st_id=k
										goto 103
									endif
									endif
								end do
								103 continue

								env_flag=.false.
								do k=1,env_bl%num
									if(env_bl%sub(k)%num_up==env_bl_num_up) then
									if(env_bl%sub(k)%num_down==env_bl_num_down) then
										env_flag=.true.
										env_id=k
										goto 104
									endif
									endif
								end do
								104 continue

								if(sys_flag.and.sys_st_flag.and.env_st_flag.and.env_flag) then
									allocate(mat1(sdim,edim),mat2(sdim,edim))

									coef=sys_st%sub(sys_st_id)%mat(1,1)*env_st%sub(env_st_id)%mat(1,1)
                                                        if(realcode)then 
									call DGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,0.0d0,mat1,sdim)
									call DGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,0.0d0,mat2,sdim)

                                                                else
									call ZGEMM('N','N',sdim,edim,edim,coef&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,env_bl%sub(env_id)%mat,edim,czero,mat1,sdim)
									call ZGEMM('N','N',sdim,edim,sdim,cone,sys_bl%sub(sys_id)%mat,sdim&
											&,wave%sub(i)%vec(spos+1:spos+sdim,epos+1:epos+edim),sdim&
											&,czero,mat2,sdim)
                                                                endif
									do z=1,edim
										value=value+dot_product(mat1(:,z),mat2(:,z))
									end do
									deallocate(mat1,mat2)
								endif
							end do
						endif
						endif
					end do
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_site_env_block_site_cor_dia


!===========================================================
!Get structure factor for square/triangleLattice
!===========================================================
subroutine Get_structure_factor(oper_cor,Filename)
	use pubdata
	implicit none

	character(len=30) :: Filename
	real(8),intent(in) :: oper_cor(Num_site,Num_site)

	integer :: i,j,x,y,dx1,dy1,dx2,dy2
	integer :: x1p,x2p,kxn,kyn,xidx,yidx
	real(8) :: mkx(0:Nx),mky(0:Ny),Stri,Strinew,Spi0
	complex(8) :: nk(0:Nx,0:Ny),nk_mid

	kxn=Nx
	kyn=Ny

	mkx=0.0d0
	do i=0,kxn
		mkx(i)=2.0d0*pi*i/kxn
	end do

	mky=0.0d0
	do i=0,kyn
		mky(i)=2.0d0*pi*i/kyn
	end do

	!Get structure factor
	nk=0.0d0
	do i=0,kxn
	do j=0,kyn
		do x=1,Num_site
			nk(i,j)=nk(i,j)+dcmplx(oper_cor(x,x),0.0d0)
		end do

		nk_mid=0.0d0
		do x=1,Num_site
		do y=1,Num_site
			if(x/=y) then
				dx1=Lattice(1,x)
				dy1=Lattice(2,x)
				dx2=Lattice(1,y)
				dy2=Lattice(2,y)
				nk_mid=nk_mid+cdexp(dcmplx(0.0d0,mkx(i)*(dx1-dx2)+mky(j)*(dy1-dy2)))*dcmplx(oper_cor(x,y),0.0d0)
			endif
		end do
		end do

		nk(i,j)=nk(i,j)+nk_mid
	end do
	end do
	nk=nk/(Num_site)

	!Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do j=0,kyn
		do i=0,kxn
			write(10,111) mkx(i),mky(j),dreal(nk(i,j))
		end do
	end do
	write(10,*)
	close(10)
111 format(F16.12,1X,F16.12,1X,F16.12)

end subroutine Get_structure_factor

subroutine Get_structure_factor_cut(oper_cor,Filename,NxCut)
	use pubdata
	implicit none

	character(len=30) :: Filename
	real(8),intent(in) :: oper_cor(Num_site,Num_site)
    integer,intent(in) :: NxCut

	integer :: i,j,x,y,dx1,dy1,dx2,dy2
	integer :: x1p,x2p,kxn,kyn,xidx,yidx
	real(8) :: mkx(0:(Nx-NxCut*2)),mky(0:Ny),Stri,Strinew,Spi0
	complex(8) :: nk(0:(Nx-NxCut*2),0:Ny),nk_mid

	kxn=Nx-NxCut*2
	kyn=Ny

	mkx=0.0d0
	do i=0,kxn
		mkx(i)=2.0d0*pi*i/kxn
	end do

	mky=0.0d0
	do i=0,kyn
		mky(i)=2.0d0*pi*i/kyn
	end do

	!Get structure factor
	nk=0.0d0
	do i=0,kxn
	do j=0,kyn
		do x=NxCut*Ny+1,Num_site-NxCut*Ny
			nk(i,j)=nk(i,j)+dcmplx(oper_cor(x,x),0.0d0)
		end do

		nk_mid=0.0d0
		do x=NxCut*Ny+1,Num_site-NxCut*Ny
		do y=NxCut*Ny+1,Num_site-NxCut*Ny
			if(x/=y) then
				dx1=Lattice(1,x)
				dy1=Lattice(2,x)
				dx2=Lattice(1,y)
				dy2=Lattice(2,y)
				nk_mid=nk_mid+cdexp(dcmplx(0.0d0,mkx(i)*(dx1-dx2)+mky(j)*(dy1-dy2)))*dcmplx(oper_cor(x,y),0.0d0)
			endif
		end do
		end do

		nk(i,j)=nk(i,j)+nk_mid
	end do
	end do
	nk=nk/(Num_site-(NxCut*2)*Ny)

	!Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do j=0,kyn
		do i=0,kxn
			write(10,111) mkx(i)/PI,mky(j)/PI,dreal(nk(i,j))
		end do
	end do
	write(10,*)
	close(10)
    open(11,file="-PI_to_PI",position='append')
    	    do j=0,kyn
		    do i=0,kxn
                  if(mkx(i)<=1.0001)then
			      write(11,111) mkx(i)/PI,mky(j)/PI,dreal(nk(i,j))
                  endif
                  if(mkx(i)>=0.99999)then
			      write(11,111) (mkx(i)/PI)-2,mky(j)/PI,dreal(nk(i,j))
                  endif
		    end do
	       end do
    close(11)
111 format(F16.12,1X,F16.12,1X,F16.12)

end subroutine Get_structure_factor_cut
!==================================================================================
!<1>: Get bond-bond correlation function (idx1<idx2): (n_i.n_j)
!==================================================================================
subroutine Get_Operator_Bond_Sys(st_oper,idx1,idx2,new_oper,new_idx,trun_idx,truns)
	use pubdata
	implicit none

	logical,intent(in) :: truns
	integer,intent(in) :: idx1,idx2,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: i,j,x,y,sys_len
	type(Total_Block) :: new_oper1,new_oper2,sys_oper

	!<0>: Check the relation of idx1 and idx2
	if(idx1>=idx2) then
		write(*,*) "idx1>=idx2 in Get_Operator_Bond_Sys"
		return
	endif

	!<1>: For bond operator: (n_i.n_j)
	call Get_sys_oper_dia(st_oper,idx1,new_oper1,idx2,trun_idx,.false.)
	call Get_sys_oper_dia(st_oper,idx2,new_oper2,idx2,trun_idx,.false.)
	call block_mul_block_dia(new_oper1,'N',new_oper2,'N',cone,new_oper)
	call deallocate_block(new_oper1)
	call deallocate_block(new_oper2)

	!<2>: Update the operator in idx2 to new_idx
	if(idx2<=trun_idx) then
		call block_transfer(new_oper,sys_oper)
	else
		call update_trun_dia(new_oper,sys_oper,systruns(idx2))
	endif
	call deallocate_block(new_oper)

	call Change_sys_oper_dia(sys_oper,idx2,new_oper,new_idx,trun_idx,Truns)
	call deallocate_block(sys_oper)

end subroutine Get_Operator_Bond_Sys


!==================================================================================
!<1>: Get bond-bond correlation function (idx1<idx2): (n_i.n_j)
!==================================================================================
subroutine Get_Operator_Bond_Env(st_oper,idx1,idx2,new_oper,new_idx,trun_idx,truns)
	use pubdata
	implicit none

	logical,intent(in) :: truns
	integer,intent(in) :: idx1,idx2,new_idx,trun_idx
	type(Total_Block),intent(in) :: st_oper
	type(Total_Block),intent(inout) :: new_oper

	integer :: i,j,x,y,env_len
	type(Total_Block) :: new_oper1,new_oper2,env_oper

	!<0>: Check the relation of idx1 and idx2
	if(idx1>=idx2) then
		write(*,*) "idx1>=idx2 in Get_Operator_Bond_Env"
		return
	endif

	!<1>: For bond operator: (n_i.n_j)
	call Get_env_oper_dia(st_oper,idx1,new_oper1,idx2,trun_idx,.false.)
	call Get_env_oper_dia(st_oper,idx2,new_oper2,idx2,trun_idx,.false.)
	call block_mul_block_dia(new_oper1,'N',new_oper2,'N',cone,new_oper)
	call deallocate_block(new_oper1)
	call deallocate_block(new_oper2)

	!<2>: Update the operator in idx2 to new_idx
	if(idx2<=trun_idx) then
		call block_transfer(new_oper,env_oper)
	else
		call update_trun_dia(new_oper,env_oper,envtruns(idx2))
	endif
	call deallocate_block(new_oper)

	call Change_env_oper_dia(env_oper,idx2,new_oper,new_idx,trun_idx,Truns)
	call deallocate_block(env_oper)

end subroutine Get_Operator_Bond_Env
!call Get_Chiral_cor (eigvec,trun_idx,Filename)

subroutine Get_sys2b(idx1,idx2,oper1, oper2, sys_oper,trun_idx,wave,flags)

!! for superconductivity
    use pubdata
        implicit none

   integer,intent(in) :: idx1,idx2
        type(Total_Block),intent(inout) :: sys_oper
   integer,intent(in) :: trun_idx
   type(Wavefunction),intent(in) :: wave
        character(len=1),intent(in) :: Flags

   real (8) :: one, coefs
   integer :: x,y,i,j,sys_len,env_len
        type(Total_Block) :: new_oper1,new_oper2
    type(Total_Block) ::  new_oper,oper1, oper2

     one=1.0d0
    sys_len=wave%sys_len
        env_len=wave%env_len
   call Get_sys_oper_ndia(oper1,idx1,new_oper1,idx2,trun_idx,.false.,flags)
   call Get_sys_oper_ndia(oper2,idx2,new_oper2,idx2,trun_idx,.false.,flags)

                new_oper%down_dif=0
    call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',cone,new_oper,sys_bsm(idx2))
                call deallocate_block(new_oper1)
                call deallocate_block(new_oper2)

            if(idx2<=trun_idx) then
              call block_transfer(new_oper,new_oper1)
              else
           call update_trun_dia(new_oper,new_oper1,systruns(idx2))
                      endif
                    call deallocate_block(new_oper)

   call Change_sys_oper_dia(new_oper1,idx2,sys_oper,sys_len-1,trun_idx,.true.)
     call deallocate_block(new_oper1)
end subroutine Get_sys2b

subroutine Get_density_cor(oper_ndia,oper1, oper2, wave,trun_idx,Flags,Filename, xi, xj)
        use pubdata
        implicit none

        character(len=1),intent(in) :: Flags
    type(Wavefunction),intent(in) :: wave
    integer,intent(inout) :: trun_idx
    character(len=30) :: Filename
           double complex,intent(inout) :: oper_ndia(Num_site,Num_site)
        integer xi, xj
    integer :: sid1,sid2,sid3, sid4,  sid12
    integer :: eid1,eid2,eid3, eid4
    integer :: tri_x1,tri_y1,tri_x2,tri_y2
    integer :: x,y,i,ij,ij1,ik,j,i1, j1,sys_len,env_len,sx,sy,ex,ey
    type(Total_Basis) :: sys_bs,env_bs
    type(Total_Block) ::sys_oper, env_oper, oper1, oper2,tmp1, new_oper, new_oper1, new_oper2
    logical :: sys_flag1,env_flag1,sys_flag2,env_flag2
    double complex  :: value
        integer  chis1, chis2

       sys_len=wave%sys_len
	env_len=wave%env_len
        trun_idx=1
        go to 11

         !<1-2>: Read truncation operators from disk
        do i=trun_idx+1,sys_len
                call truns_from_disk(systruns(i),1001,i,.true.)
        end do

        do i=trun_idx+1,env_len
                call truns_from_disk(envtruns(i),1001,i,.false.)
        end do

        !<1-3>: Read basis from disk
        do i=2,sys_len
                call basis_from_disk(sys_bsm(i),1001,i,.true.)
        end do

        do i=2,env_len
                call basis_from_disk(env_bsm(i),1001,i,.false.)
        end do


11              continue

	call basis_transfer(sys_bsm(sys_len),sys_bs)
	call basis_transfer(env_bsm(env_len),env_bs)

            !<2>: Open file and save to disk
	        open(10,file=Filename,position='append')
	        write(10,*) "Nx=",Nx,"Ny=",Ny,"Jd=",jd(7,1),"t2=",jt(1, 1),"N_up=",tot_num_up&
				        &,"N_down=",tot_num_down, kept_min, kept_max
	        !!!write(10,*) "hubbard", hubbard(1)


      do sid1=xi, xj    !!!!!Lattice(1,sys_len)
        do sid2=sid1+ny, sys_len-1, ny
        call Get_sys2b(min(sid1,sid2),max(sid1,sid2),oper1,oper2, sys_oper,trun_idx,wave,flags)
		call measure_sys_block(value,sys_oper,sys_bs,wave)

        value=value*dsqrt(1.0d0+oper1%down_dif)  !!!*2.0d0
        if(oper1%down_dif==2) value=-value !!!*2.0d0
        oper_ndia(sid1,sid2)=value
        
41      format(2i8, 3f18.9)

                call deallocate_block(sys_oper)
21      continue
        enddo

        if(sys_len.le.num_site-num_site/6+ny)then
   call Get_sys_oper_ndia(oper1,sid1,new_oper,sid1,trun_idx,.false.,flags)
            if(sid1<=trun_idx) then
              call block_transfer(new_oper,new_oper1)
              else
           call update_trun_ndia(new_oper,new_oper1,systruns(sid1))
                      endif
                    call deallocate_block(new_oper)

   call Change_sys_oper_ndia(new_oper1,sid1,sys_oper,sys_len-1,trun_idx,.true.,flags)
     call deallocate_block(new_oper1)
                endif

        do sid3=sys_len, num_site-num_site/6+ny

        if(mod(sid1-sid3, ny).ge.0)then

        value=0.0d0
        if(sid3==sys_len)then
        
        if(oper2%down_dif.eq.1)then
        call block_transfer_trans(oper2, tmp1)
         call sys_block_site_cor_ndia(value,sys_oper,tmp1,sys_bs,wave,FlagS)
        call deallocate_block(tmp1)
                else
         call sys_block_site_cor_ndia(value,sys_oper,oper2,sys_bs,wave,FlagS)
                endif
        endif

        if(sid3==sys_len+1)then
        if(oper2%down_dif.eq.1)then
        call block_transfer_trans(oper2, tmp1)
 call sys_block_env_site_cor_ndia(value,sys_oper,tmp1,sys_bs,env_bs,wave,FlagS)
        call deallocate_block(tmp1)
                else
 call sys_block_env_site_cor_ndia(value,sys_oper,oper2,sys_bs,env_bs,wave,FlagS)
                endif
        endif

        if(sid3.ge.sys_len+2)then
        
        eid1=num_site-sid3+1

   call Get_env_oper_ndia(oper2,eid1,new_oper,eid1,trun_idx,.false.,flags)
            if(eid1<=trun_idx) then
              call block_transfer(new_oper,new_oper2)
              else
           call update_trun_ndia(new_oper,new_oper2,envtruns(eid1))
!!update_trun_ndia update_trun_dia
                      endif
                    call deallocate_block(new_oper)

     call Change_env_oper_ndia(new_oper2,eid1,env_oper,env_len-1,trun_idx,.true.,flags)
     call deallocate_block(new_oper2)

        if(oper2%down_dif.ne.0)then
        call block_transfer_trans(env_oper, tmp1)
          call sys_block_env_block_cor_ndia(value,sys_oper,tmp1,sys_bs,env_bs,wave,flags)
        call deallocate_block(tmp1)
                        else
          call sys_block_env_block_cor_ndia(value,sys_oper,env_oper,sys_bs,env_bs,wave,flags)

                        endif
                    call deallocate_block(env_oper)
                endif

        value=value*dsqrt(1.0d0+oper1%down_dif)  !!!*2.0d0
        if(oper1%down_dif==2) value=-value !!!*2.0d0
        oper_ndia(sid1,sid3)=value
        
        write(10,41)sid1, sid2, real(value), imag(value)

       endif
        enddo
                    if(allocated(sys_oper%sub))call deallocate_block(sys_oper)
                    if(allocated(new_oper%sub)) call deallocate_block(new_oper)

        enddo




22      continue
close(10)

end subroutine Get_density_cor
!=============================================================
subroutine Get_mSCOP_Operator_Cor0(wave,trun_idx,Flags,Filename,xc1, xc2, ij, ij1)
        use pubdata
        implicit none

        character(len=2),intent(in) :: Flags
    type(Wavefunction),intent(in) :: wave
    integer,intent(inout) :: trun_idx, xc1, xc2
    character(len=30) :: Filename

    integer :: sid1,sid2,sid3, sid4,  sid12
    integer :: eid1,eid2,eid3, eid4
    integer :: tri_x1,tri_y1,tri_x2,tri_y2
    integer :: x,y,i,ij,ij1,ik,j,i1, j1,sys_len,env_len,sx,sy,ex,ey
    type(Total_Basis) :: sys_bs,env_bs
    type(Total_Block) ::sys_oper, sys_oper2, sys_oper1, env_oper
    logical :: sys_flag1,env_flag1,sys_flag2,env_flag2
    double complex :: value, val1(3), val11(num_site, -num_site:num_site, 3)
    integer  :: kv(num_site, -num_site:num_site, 3)


    sys_len=wave%sys_len
	env_len=wave%env_len
        trun_idx=1

	call basis_transfer(sys_bsm(sys_len),sys_bs)
	call basis_transfer(env_bsm(env_len),env_bs)

            !<2>: Open file and save to disk
	        open(14,file=Filename,position='append')
	        write(14,*) "Nx=",Nx,"Ny=",Ny,"Jd=",jd(1,1),"t=",jt(1, 1),"Jchi=", jring(1),"N_up=",tot_num_up&
				        &,"N_down=",tot_num_down, kept_min, kept_max, "bond", ij1

        val11=0.0d0
        kv=0

        go to 11

        do sid1=1,  xc1-1
        val1=0.0d0
        do i=1,5,2
        j=i
        sid2=neib(sid1,j)
        call Get_sys_two(min(sid1,sid2),max(sid1,sid2),sys_oper2,trun_idx,wave)

      do sid3=xc1, xc1    !!!!!Lattice(1,sys_len)
        sid4= neib(sid3, 1)  !! consider x direction

        if(max(sid1,sid2).ge.min(sid3,sid4))go to 211

                if(max(sid3,sid4)<=(sys_len-1)) then
        call Get_sys_two1 (max(sid1,sid2),min(sid3,sid4),max(sid3,sid4),sys_oper2,sys_oper,trun_idx,wave)
		call measure_sys_block(value,sys_oper,sys_bsm(sys_len),wave)
        val1((i+1)/2)=value
        val11(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)=val11(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)+val1((i+1)/2)
        kv(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)=kv(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)+1
                call deallocate_block(sys_oper)
                 endif

211      continue
        enddo
                call deallocate_block(sys_oper2)
        enddo

        sid3=xc1
        write(14,151)sid3,sid1, mod(sid1-sid3+num_site, ny),(sid1-1)/ny-(sid3-1)/ny,abs(val1(1:3)), imag(log(val1(1:3)))
151      format(4i8, 3f18.10)
        enddo



11      continue

      do sid1=xc1, xc1     !!!!!Lattice(1,sys_len)
        sid2= neib(sid1, ij1)  !! consider x direction
        call Get_sys_two(min(sid1,sid2),max(sid1,sid2),sys_oper2,trun_idx,wave)


        do sid3=sid1+ny+ij, sys_len-1, ny
        val1=0.0d0
        if(sid3.le.max(sid1,sid2))go to 23
        if(sid3.gt.xc2)go to 23
        do i=1,5,2
        j=i
        sid4=neib(sid3,j)
        if(min(sid4,sid3).le.max(sid1,sid2))go to 21

             if(max(sid3,sid4)<=(sys_len-1)) then
        call Get_sys_two1 (max(sid1,sid2),min(sid3,sid4),max(sid3,sid4),sys_oper2,sys_oper,trun_idx,wave)
		call measure_sys_block(value,sys_oper,sys_bsm(sys_len),wave)
        val1((i+1)/2)=value
        val11(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)=val11(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)+val1((i+1)/2)
        kv(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)=kv(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)+1
                call deallocate_block(sys_oper)
                 endif

21      continue
        enddo
        write(14,15)sid1,sid3, mod(sid3-sid1+num_site, ny),(sid3-1)/ny-(sid1-1)/ny, abs(val1(1:3)), imag(log(val1(1:3))),sid2, sid4

15      format(4i8, 6f21.12, 2i8)
23      continue

        write(*,*)sid1, sid2, sid3, sid4, i1
        enddo


12      continue

     call Change_sys_oper_ndia(sys_oper2,max(sid1,sid2),sys_oper,sys_len-1,trun_idx,.true.,'B') 
                call deallocate_block(sys_oper2)

        do i1=sys_len+2-mod(sys_len+2,ny)+mod(sid1-1,ny)+1+ij, xc2, ny

        sid3=i1
        eid1=num_site-sid3+1

        if(eid1.gt.env_len-1)go to 16


        val1=0.0d0

        do i=1,5,2
        j=i
        sid4=neib(sid3,j)
        eid2=num_site-sid4+1

        
        if(max(eid1, eid2).le.env_len-1)then
	call Get_SCOP_Operator_Env(min(eid1,eid2),max(eid1,eid2),env_oper,env_len-1,trun_idx,.true.,Flags)
	call sys_block_env_block_cor_ndia(value,sys_oper,env_oper,sys_bs,env_bs,wave,'B')
        call deallocate_block(env_oper)

        val1((i+1)/2)=dconjg(value)
        val11(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)=val11(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)+val1((i+1)/2)
        kv(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)=kv(sid1, (sid3-1)/ny-(sid1-1)/ny, (i+1)/2)+1
        endif
        enddo

        write(14,15)sid1,sid3, mod(sid3-sid1+num_site, ny),(sid3-1)/ny-(sid1-1)/ny, abs(val1(1:3)), imag(log(val1(1:3))),sid2, sid4
16      continue
                enddo

        call deallocate_block(sys_oper)

        enddo
22      continue

        close(10)

19      format(2i8, 2f21.12)

end subroutine Get_mSCOP_Operator_cor0
!=============================================================

subroutine Get_sys_two (idx1,idx2,sys_oper,trun_idx,wave)
!! for superconductivity
    use pubdata
	implicit none

   integer,intent(in) :: idx1,idx2
   	type(Total_Block),intent(inout) :: sys_oper
   integer,intent(in) :: trun_idx
   type(Wavefunction),intent(in) :: wave
   
   real (8) :: one, coefs
   integer :: x,y,i,j,sys_len,env_len
   	type(Total_Block) :: new_oper1,new_oper2,new_oper3,new_oper12 
    type(Total_Block) ::  new_oper

     one=cone
    sys_len=wave%sys_len
	env_len=wave%env_len


   call Get_sys_oper_ndia(st_elec_up,idx1,new_oper1,idx2,trun_idx,.false.,'F')

    call Get_sys_oper_ndia(st_elec_up,idx2,new_oper2,idx2,trun_idx,.false.,'F')

                new_oper%down_dif=0
    call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',cone,new_oper,sys_bsm(idx2))
                call deallocate_block(new_oper1)
                call deallocate_block(new_oper2)

            if(idx2<=trun_idx) then
              call block_transfer(new_oper,sys_oper)
              else
           call update_trun_ndia(new_oper,sys_oper,systruns(idx2))
!!update_trun_ndia update_trun_dia
                      endif
                    call deallocate_block(new_oper)

end subroutine Get_sys_two


subroutine Get_sys_two1(idx2,idx3,idx4,sys_oper,sys_four,trun_idx,wave)

    use pubdata
	implicit none

   integer,intent(in) ::idx2,idx3, idx4
   	type(Total_Block),intent(inout) :: sys_four
   	type(Total_Block),intent(in) :: sys_oper
   integer,intent(in) :: trun_idx
   type(Wavefunction),intent(in) :: wave
   
   real (8) :: one, coefs
   integer :: x,y,i,j,sys_len,env_len
   	type(Total_Block) :: new_oper1,new_oper2,new_oper3,new_oper12,new_oper123 
    type(Total_Block) ::  new_oper

     one=cone
    sys_len=wave%sys_len
	env_len=wave%env_len
                
   call Change_sys_oper_ndia(sys_oper,idx2,new_oper,idx3,trun_idx,.false.,'B') !!! sys_op2
   call Get_sys_oper_ndia(st_elec_down,idx3,new_oper1,idx3,trun_idx,.false.,'F')
                        new_oper2%down_dif=1
    call block_mul_block_ndia(new_oper,'N',new_oper1,'N',cone,new_oper2,sys_bsm(idx3))
        call deallocate_block(new_oper1)

            if(idx3<=trun_idx) then
              call block_transfer(new_oper2,new_oper1)
              else
           call update_trun_ndia(new_oper2,new_oper1,systruns(idx3))
                      endif
                    call deallocate_block(new_oper2)

    call Change_sys_oper_ndia(new_oper1,idx3,new_oper12,idx4,trun_idx,.false.,'F') !!! sys_op2
    call Get_sys_oper_ndia(st_elec_down,idx4,new_oper2,idx4,trun_idx,.false.,'F')
                        new_oper3%down_dif=0
        call block_mul_block_ndia(new_oper12,'N',new_oper2,'N',cone,new_oper3,sys_bsm(idx4))
     call deallocate_block(new_oper1)
     call deallocate_block(new_oper2)
     call deallocate_block(new_oper12)

                call block_transfer(new_oper3,new_oper)
                call deallocate_block(new_oper3)


         if(idx4<=trun_idx) then
		              call block_transfer (new_oper,new_oper1)
 	              else
 		              call update_trun_dia(new_Oper,new_oper1,systruns(idx4))
 	              endif
	           call deallocate_block(new_oper) 

     call Change_sys_oper_dia(new_oper1,idx4,sys_four,sys_len-1,trun_idx,.true.)  ! true ±£Ö¤ÁËtrun
     call deallocate_block(new_oper1)

end subroutine Get_sys_two1



subroutine Get_sys_four (idx1,idx2,idx3,idx4,sys_four,trun_idx,wave)

    use pubdata
	implicit none

   integer,intent(in) :: idx1,idx2,idx3, idx4
   	type(Total_Block),intent(inout) :: sys_four
   integer,intent(in) :: trun_idx
   type(Wavefunction),intent(in) :: wave
   
   real (8) :: one, coefs
   integer :: x,y,i,j,sys_len,env_len
   	type(Total_Block) :: new_oper1,new_oper2,new_oper3,new_oper12,new_oper123 
    type(Total_Block) ::  new_oper

     one=cone
    sys_len=wave%sys_len
	env_len=wave%env_len


   call Get_sys_oper_ndia(st_elec_up,idx1,new_oper1,idx2,trun_idx,.false.,'F')
    call Get_sys_oper_ndia(st_elec_up,idx2,new_oper2,idx2,trun_idx,.false.,'F')
                new_oper12%down_dif=0
    call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',cone,new_oper12,sys_bsm(idx2))
                call deallocate_block(new_oper1)
                call deallocate_block(new_oper2)

                call block_transfer(new_oper12,new_oper)
                call deallocate_block(new_oper12)

            if(idx2<=trun_idx) then
              call block_transfer(new_oper,new_oper12)
              else
           call update_trun_ndia(new_oper,new_oper12,systruns(idx2))
!!update_trun_ndia update_trun_dia
                      endif
                    call deallocate_block(new_oper)

     call Change_sys_oper_ndia(new_oper12,idx2,new_oper,idx3,trun_idx,.false.,'B') !!! sys_op2
     call deallocate_block(new_oper12)

   call Get_sys_oper_ndia(st_elec_up,idx3,new_oper1,idx3,trun_idx,.false.,'F')
                        new_Oper2%down_dif=1
    call block_mul_block_ndia(new_oper,'N',new_oper1,'T',cone,new_oper2,sys_bsm(idx3))
        call deallocate_block(new_oper1)

            if(idx3<=trun_idx) then
              call block_transfer(new_oper2,new_oper1)
              else
           call update_trun_ndia(new_oper2,new_oper1,systruns(idx3))
                      endif
                    call deallocate_block(new_oper2)


     call Change_sys_oper_ndia(new_oper1,idx3,new_oper12,idx4,trun_idx,.false.,'F') !!! sys_op2
    call Get_sys_oper_ndia(st_elec_up,idx4,new_oper2,idx4,trun_idx,.false.,'F')
                        new_oper3%down_dif=0
        call block_mul_block_ndia(new_oper12,'N',new_oper2,'T',cone,new_oper3,sys_bsm(idx4))
     call deallocate_block(new_oper1)
     call deallocate_block(new_oper2)
     call deallocate_block(new_oper12)

                call block_pass_info(new_oper3,new_oper)
                call deallocate_block(new_oper3)

!!!!!!!!!!!xxxxxxxxxxxxxx

         if(idx4<=trun_idx) then
		              call block_transfer (new_oper,new_oper1)
 	              else
 		              call update_trun_dia(new_Oper,new_oper1,systruns(idx4))
 	              endif
	           call deallocate_block(new_oper) 

     call Change_sys_oper_dia(new_oper1,idx4,sys_four,sys_len-1,trun_idx,.true.)  ! true ±£Ö¤ÁËtrun
     call deallocate_block(new_oper1)

end subroutine Get_sys_four



subroutine Get_sys_triangular (id1,id2,id3,sys_tri,trun_idx,wave)

    use pubdata
	implicit none

   integer,intent(in) :: id1,id2,id3
   	type(Total_Block),intent(inout) :: sys_tri
   integer,intent(in) :: trun_idx
   type(Wavefunction),intent(in) :: wave
   
   real (8) :: one
   integer :: x,y,i,j,sys_len,env_len
   	type(Total_Block) :: new_oper1,new_oper2,new_oper3,new_oper12,new_oper123 
    type(Total_Block) ::  new_oper,new_oper_temp


   one=cone
    sys_len=wave%sys_len
	env_len=wave%env_len

     call Get_sys_oper_ndia(st_sd,id1,new_oper1,id2,trun_idx,.false.,'B')
     call Get_sys_oper_ndia(st_sd,id2,new_oper2,id2,trun_idx,.false.,'B')
                        new_oper%down_dif=su
     call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',one,new_oper,sys_bsm(id2))
	 call deallocate_block(new_oper1)
	 call deallocate_block(new_oper2)
                 if(id2<=trun_idx) then
		              call block_transfer(new_oper,new_oper12)
	              else
 		              call update_trun_ndia(new_oper,new_oper12,systruns(id2))        !!update_trun_ndia update_trun_dia
 	              endif
	            call deallocate_block(new_oper)
 
     call Change_sys_oper_ndia(new_oper12,id2,new_oper_temp,id3,trun_idx,.false.,'B')
     call deallocate_block(new_oper12)

    call Get_sys_oper_ndia(st_sd,id3,new_oper3,id3,trun_idx,.false.,'B') ! false ±£Ö¤ÁË²»ÏÈtrun
                new_oper%down_dif=0
    call block_mul_block_ndia(new_oper_temp,'N',new_oper3,'N',one,new_oper,sys_bsm(id3))
    call deallocate_block(new_oper_temp)
    call deallocate_block(new_oper3)

	call block_transfer(new_oper,sys_tri)

	call deallocate_block(new_oper)

   !=================
         if(id3<=trun_idx) then
		              call block_transfer (sys_tri,new_oper123)
 	              else
 		              call update_trun_dia(sys_tri,new_oper123,systruns(id3))
 	              endif
	           call deallocate_block(sys_tri) 

     call Change_sys_oper_dia(new_oper123,id3,sys_tri,sys_len-1,trun_idx,.true.)  ! true ±£Ö¤ÁËtrun
     call deallocate_block(new_oper123)

end subroutine Get_sys_triangular

  


 subroutine Get_env_triangular (id1,id2,id3,env_tri,trun_idx,wave)

    use pubdata
	implicit none

   integer,intent(in) :: id1,id2,id3
   	type(Total_Block),intent(inout) :: env_tri
   integer,intent(in) :: trun_idx
   type(Wavefunction),intent(in) :: wave
   
   real(8) :: one
   integer :: x,y,i,j,sys_len,env_len
   	type(Total_Block) :: new_oper1,new_oper2,new_oper3,new_oper12,new_oper123 
    type(Total_Block) ::  new_oper,new_oper_temp


   one=cone
    sys_len=wave%sys_len
	env_len=wave%env_len
   !=================
     call Get_env_oper_ndia(st_sd,id1,new_oper1,id2,trun_idx,.false.,'B')
     call Get_env_oper_ndia(st_sd,id2,new_oper2,id2,trun_idx,.false.,'B')
                        new_oper%down_dif=su
     call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',one,new_oper,env_bsm(id2))
	 call deallocate_block(new_oper1)
	 call deallocate_block(new_oper2)
                 if(id2<=trun_idx) then
		              call block_transfer(new_oper,new_oper12)
	              else
 	         call update_trun_ndia(new_oper,new_oper12,envtruns(id2)) 
 	              endif
	            call deallocate_block(new_oper)

     call Change_env_oper_ndia(new_oper12,id2,new_oper_temp,id3,trun_idx,.false.,'B')
     call deallocate_block(new_oper12)

    call Get_env_oper_ndia(st_sd,id3,new_oper3,id3,trun_idx,.false.,'B') ! false ±£Ö¤ÁË²»ÏÈtrun
                new_oper%down_dif=0
    call block_mul_block_ndia(new_oper_temp,'N',new_oper3,'N',one,new_oper,env_bsm(id3))
    call deallocate_block(new_oper_temp)
    call deallocate_block(new_oper3)

	call block_transfer(new_oper,env_tri)
	call deallocate_block(new_oper)
         if(id3<=trun_idx) then
		              call block_transfer (env_tri,new_oper123)
 	              else
 		              call update_trun_dia(env_tri,new_oper123,envtruns(id3))
 	              endif
	           call deallocate_block(env_tri) 

     call Change_env_oper_dia(new_oper123,id3,env_tri,env_len-1,trun_idx,.true.)  ! true ±£Ö¤ÁËtrun
     call deallocate_block(new_oper123)

end subroutine Get_env_triangular

!===============================================================================
!Get all non-diagonal operator correlatioins using memory directly
!===============================================================================
subroutine Get_oper_cor_ndia1(oper_ndia,st_oper, st_oper1, wave,trun_idx,FlagSign,Filename)
	use pubdata
	implicit none

	character(len=1) :: FlagSign
	character(len=30) :: Filename
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper, st_oper1
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: oper_ndia(Num_site,Num_site)
	real*8 :: value
	integer :: i,j,sys_len,env_len,ti,tj
	type(Total_Basis) :: sys_bs,env_bs
	type(Total_Block) :: st_oper2,new_oper, tmp1, tmp2, st_oper11

	!<1>: For information saving to the disk
	open(10,file="oper_ndia_tmp.dat",position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max,Filename
	109	format(A3,I2,1X,A3,I2,1X,A2,F6.3,1X,A2,F6.3,1X,A5,I3,1X,A7,I3,1X,A3,1X,I4,1X,A3,1X,I4,1X,A30)
	close(10)


	!<2>: Get general information
	sys_len=wave%sys_len
	env_len=wave%env_len
	call basis_transfer(sys_bsm(sys_len),sys_bs)
	call basis_transfer(env_bsm(env_len),env_bs)

                st_oper2%down_dif=0 
                call block_transfer_trans(st_oper1, st_oper11)
	call block_mul_block_ndia(st_oper,'N',st_oper11,'T',cone,st_oper2,st_basis)
	oper_ndia=czero

	!On-site density correlation
	!<2>: Get <n_i^2> in system block
	do i=xi1,sys_len-1
		call Get_sys_oper_dia(st_oper2,i,new_oper,sys_len-1,trun_idx,.true.)
		call measure_sys_block(oper_ndia(i,i),new_oper,sys_bsm(sys_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_sys_site(oper_ndia(sys_len,sys_len),st_oper2,sys_bsm(sys_len),wave)

	!<2>: Get <n_i^2> in environment block
	do i=xj1,env_len-1
		call Get_env_oper_dia(st_oper2,i,new_oper,env_len-1,trun_idx,.true.)
		call measure_env_block(oper_ndia(Num_site-i+1,Num_site-i+1),new_oper,env_bsm(env_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_env_site(oper_ndia(Num_site-env_len+1,Num_site-env_len+1),st_oper2,env_bsm(env_len),wave)
	call deallocate_block(st_oper2)
	
	!Save to the disk
	do i=1,Num_site
		open(10,file="oper_ndia_tmp.dat",position='append')
		!oper_ndia(i,i)=ave_dia(i)
		write(10,*) i,i,oper_ndia(i,i)
		close(10)
	end do


	!<2>: Get sysopers and envopers
	do i=1,sys_len-1
		call block_transfer(st_oper,sysopers(i))
	end do

	do i=1,env_len-1
		call block_transfer(st_oper1,envopers(i))
	end do

	!<3>: Get sys_bl_bl and env_bl_bl operator correlation
	call Get_sys_oper_cor_ndia1(oper_ndia,st_oper,st_oper1,wave,trun_idx,FlagSign)
	call Get_env_oper_cor_ndia1(oper_ndia,st_oper,st_oper1,wave,trun_idx,FlagSign)


	!<4>: Get sys_bl_st and sys_bl_env_st
	do i=xi1,sys_len-1
        
                        call block_transfer_trans(st_oper1, tmp1)
		call sys_block_site_cor_ndia(oper_ndia(i,sys_len),sysopers(i),tmp1,sys_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(sys_len,i)=oper_ndia(i,sys_len)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(sys_len,i)=-oper_ndia(i,sys_len)
			oper_ndia(sys_len,i)=oper_ndia(i,sys_len) !(check)
		endif

		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) i,sys_len,oper_ndia(i,sys_len)
		write(10,*) sys_len,i,oper_ndia(sys_len,i)
		close(10)

                        call block_transfer_trans(st_oper1, tmp1)
		call sys_block_env_site_cor_ndia(oper_ndia(i,Num_site-env_len+1),sysopers(i),tmp1,sys_bs,env_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(Num_site-env_len+1,i)=oper_ndia(i,Num_site-env_len+1)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(Num_site-env_len+1,i)=-oper_ndia(i,Num_site-env_len+1)
			oper_ndia(Num_site-env_len+1,i)=oper_ndia(i,Num_site-env_len+1) !(check)
		endif

		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) i,Num_site-env_len+1,oper_ndia(i,Num_site-env_len+1)
		write(10,*) Num_site-env_len+1,i,oper_ndia(Num_site-env_len+1,i)
		close(10)
	end do

	!<5>: Get env_bl_st and sys_st_env_bl
	do i=xj1,env_len-1
                        call block_transfer_trans(envopers(i), tmp1)
		call env_block_site_cor_ndia(oper_ndia(Num_site-i+1,Num_site-env_len+1),tmp1,st_oper,env_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(Num_site-env_len+1,Num_site-i+1)=oper_ndia(Num_site-i+1,Num_site-env_len+1)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(Num_site-env_len+1,Num_site-i+1)=-oper_ndia(Num_site-i+1,Num_site-env_len+1)
			oper_ndia(Num_site-env_len+1,Num_site-i+1)=-oper_ndia(Num_site-i+1,Num_site-env_len+1) !(check)
			oper_ndia(Num_site-i+1,Num_site-env_len+1)=-oper_ndia(Num_site-i+1,Num_site-env_len+1) !(check)
		endif

		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) Num_site-i+1,Num_site-env_len+1,oper_ndia(Num_site-i+1,Num_site-env_len+1)
		write(10,*) Num_site-env_len+1,Num_site-i+1,oper_ndia(Num_site-env_len+1,Num_site-i+1)
		close(10)

                        call block_transfer_trans(envopers(i), tmp1)
		call sys_site_env_block_cor_ndia(oper_ndia(sys_len,Num_site-i+1),st_oper,tmp1,sys_bs,env_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(Num_site-i+1,sys_len)=oper_ndia(sys_len,Num_site-i+1)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(Num_site-i+1,sys_len)=-oper_ndia(sys_len,Num_site-i+1)
			oper_ndia(Num_site-i+1,sys_len)=oper_ndia(sys_len,Num_site-i+1) !(check)
		endif

		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) Num_site-i+1,sys_len,oper_ndia(Num_site-i+1,sys_len)
		write(10,*) sys_len,Num_site-i+1,oper_ndia(sys_len,Num_site-i+1)
		close(10)
	end do

	!<6>: Get sys_bl_env_bl
	do j=xj1,env_len-1
                        call block_transfer_trans(envopers(j), tmp1)
	do i=xi1,sys_len-1
		call sys_block_env_block_cor_ndia(oper_ndia(i,Num_site-j+1),sysopers(i),tmp1,sys_bs,env_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(Num_site-j+1,i)=oper_ndia(i,Num_site-j+1)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(Num_site-j+1,i)=-oper_ndia(i,Num_site-j+1)
			oper_ndia(Num_site-j+1,i)=oper_ndia(i,Num_site-j+1) !(check)
		endif
			
		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) i,Num_site-j+1,oper_ndia(i,Num_site-j+1)
		write(10,*) Num_site-j+1,i,oper_ndia(Num_site-j+1,i)
		close(10)
	end do
	end do

                        call block_transfer_trans(st_oper1, tmp1)
	call sys_site_env_site_cor_ndia(oper_ndia(sys_len,Num_site-env_len+1),st_oper,tmp1,sys_bs,env_bs,wave,FlagSign)
	if(FlagSign=='B') then !For Boson
		oper_ndia(Num_site-env_len+1,sys_len)=oper_ndia(sys_len,Num_site-env_len+1)
	else if(FlagSign=='F') then !For Fermion
		!oper_ndia(Num_site-env_len+1,sys_len)=-oper_ndia(sys_len,Num_site-env_len+1)
		oper_ndia(Num_site-env_len+1,sys_len)=oper_ndia(sys_len,Num_site-env_len+1) !(check)
	endif

	!Save to disk
	open(10,file="oper_ndia_tmp.dat",position='append')
	write(10,*) sys_len,Num_site-env_len+1,oper_ndia(sys_len,Num_site-env_len+1)
	write(10,*) Num_site-env_len+1,sys_len,oper_ndia(Num_site-env_len+1,sys_len)
	write(10,*)
	close(10)


	!Free space
	do i=1,sys_len-1
		call deallocate_block(sysopers(i))
	end do

	do i=1,env_len-1
		call deallocate_block(envopers(i))
	end do
	call deallocate_basis(sys_bs)
	call deallocate_basis(env_bs)


	!<5>: Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do i=xi1,Num_site-xj1+1
		do j=xi1,Num_site-xj1+1
        oper_ndia(i,j)=oper_ndia(i,j)*dsqrt(1.0d0+st_oper%down_dif)  !!!*2.0d0
        if(st_oper%down_dif==2) oper_ndia(i,j)=-oper_ndia(i,j)  !!!*2.0d0
			write(10,*) i,j,oper_ndia(i,j)
		end do
	end do
	write(10,*)
	close(10)

	write(10,*)
	close(10)



111 format(I4,1X,I4,1X,F16.12)



end subroutine Get_oper_cor_ndia1



!===============================================================================
!Get all non-diagonal operator correlatioins using memory directly
!===============================================================================
subroutine Get_oper_cor_ndia(oper_ndia,st_oper, st_oper1, wave,trun_idx,FlagSign,Filename)
	use pubdata
	implicit none

	character(len=1) :: FlagSign
	character(len=30) :: Filename
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper, st_oper1
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: oper_ndia(Num_site,Num_site)
	double complex :: value
	integer :: i,j,sys_len,env_len,ti,tj
	type(Total_Basis) :: sys_bs,env_bs
	type(Total_Block) :: st_oper2,new_oper, tmp1, tmp2, st_oper11

	!<1>: For information saving to the disk
	open(10,file="oper_ndia_tmp.dat",position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max,Filename
	close(10)


	!<2>: Get general information
	sys_len=wave%sys_len
	env_len=wave%env_len
	call basis_transfer(sys_bsm(sys_len),sys_bs)
	call basis_transfer(env_bsm(env_len),env_bs)

                st_oper2%down_dif=0 
                call block_transfer_trans(st_oper1, st_oper11)
	call block_mul_block_ndia(st_oper,'N',st_oper11,'T',cone,st_oper2,st_basis)
	oper_ndia=czero

	!On-site density correlation
	!<2>: Get <n_i^2> in system block
	do i=1,sys_len-1
		call Get_sys_oper_dia(st_oper2,i,new_oper,sys_len-1,trun_idx,.true.)
		call measure_sys_block(oper_ndia(i,i),new_oper,sys_bsm(sys_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_sys_site(oper_ndia(sys_len,sys_len),st_oper2,sys_bsm(sys_len),wave)

	!<2>: Get <n_i^2> in environment block
	do i=1,env_len-1
		call Get_env_oper_dia(st_oper2,i,new_oper,env_len-1,trun_idx,.true.)
		call measure_env_block(oper_ndia(Num_site-i+1,Num_site-i+1),new_oper,env_bsm(env_len),wave)
		call deallocate_block(new_oper)
	end do
	call measure_env_site(oper_ndia(Num_site-env_len+1,Num_site-env_len+1),st_oper2,env_bsm(env_len),wave)
	call deallocate_block(st_oper2)
	
	!Save to the disk
	do i=1,Num_site
		open(10,file="oper_ndia_tmp.dat",position='append')
		!oper_ndia(i,i)=ave_dia(i)
		write(10,*) i,i,oper_ndia(i,i)
		close(10)
	end do


	!<2>: Get sysopers and envopers
	do i=1,sys_len-1
		call block_transfer(st_oper,sysopers(i))
	end do

	do i=1,env_len-1
		call block_transfer(st_oper1,envopers(i))
	end do

	!<3>: Get sys_bl_bl and env_bl_bl operator correlation
	call Get_sys_oper_cor_ndia(oper_ndia,st_oper,st_oper1,wave,trun_idx,FlagSign)
	call Get_env_oper_cor_ndia(oper_ndia,st_oper,st_oper1,wave,trun_idx,FlagSign)


	!<4>: Get sys_bl_st and sys_bl_env_st
	do i=1,sys_len-1
        
                        call block_transfer_trans(st_oper1, tmp1)
		call sys_block_site_cor_ndia(oper_ndia(i,sys_len),sysopers(i),tmp1,sys_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(sys_len,i)=oper_ndia(i,sys_len)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(sys_len,i)=-oper_ndia(i,sys_len)
			oper_ndia(sys_len,i)=oper_ndia(i,sys_len) !(check)
		endif

		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) i,sys_len,oper_ndia(i,sys_len)
		write(10,*) sys_len,i,oper_ndia(sys_len,i)
		close(10)

                        call block_transfer_trans(st_oper1, tmp1)
		call sys_block_env_site_cor_ndia(oper_ndia(i,Num_site-env_len+1),sysopers(i),tmp1,sys_bs,env_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(Num_site-env_len+1,i)=oper_ndia(i,Num_site-env_len+1)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(Num_site-env_len+1,i)=-oper_ndia(i,Num_site-env_len+1)
			oper_ndia(Num_site-env_len+1,i)=oper_ndia(i,Num_site-env_len+1) !(check)
		endif

		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) i,Num_site-env_len+1,oper_ndia(i,Num_site-env_len+1)
		write(10,*) Num_site-env_len+1,i,oper_ndia(Num_site-env_len+1,i)
		close(10)
	end do

	!<5>: Get env_bl_st and sys_st_env_bl
	do i=1,env_len-1
                        call block_transfer_trans(envopers(i), tmp1)

		call env_block_site_cor_ndia(oper_ndia(Num_site-i+1,Num_site-env_len+1),tmp1,st_oper,env_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(Num_site-env_len+1,Num_site-i+1)=oper_ndia(Num_site-i+1,Num_site-env_len+1)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(Num_site-env_len+1,Num_site-i+1)=-oper_ndia(Num_site-i+1,Num_site-env_len+1)
			oper_ndia(Num_site-env_len+1,Num_site-i+1)=-oper_ndia(Num_site-i+1,Num_site-env_len+1) !(check)
			oper_ndia(Num_site-i+1,Num_site-env_len+1)=-oper_ndia(Num_site-i+1,Num_site-env_len+1) !(check)
		endif

		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) Num_site-i+1,Num_site-env_len+1,oper_ndia(Num_site-i+1,Num_site-env_len+1)
		write(10,*) Num_site-env_len+1,Num_site-i+1,oper_ndia(Num_site-env_len+1,Num_site-i+1)
		close(10)

                        call block_transfer_trans(envopers(i), tmp1)
		call sys_site_env_block_cor_ndia(oper_ndia(sys_len,Num_site-i+1),st_oper,tmp1,sys_bs,env_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(Num_site-i+1,sys_len)=oper_ndia(sys_len,Num_site-i+1)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(Num_site-i+1,sys_len)=-oper_ndia(sys_len,Num_site-i+1)
			oper_ndia(Num_site-i+1,sys_len)=oper_ndia(sys_len,Num_site-i+1) !(check)
		endif

		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) Num_site-i+1,sys_len,oper_ndia(Num_site-i+1,sys_len)
		write(10,*) sys_len,Num_site-i+1,oper_ndia(sys_len,Num_site-i+1)
		close(10)
	end do

	!<6>: Get sys_bl_env_bl
	do j=1,env_len-1
                        call block_transfer_trans(envopers(j), tmp1)
	do i=1,sys_len-1
		call sys_block_env_block_cor_ndia(oper_ndia(i,Num_site-j+1),sysopers(i),tmp1,sys_bs,env_bs,wave,FlagSign)
		if(FlagSign=='B') then !For Boson
			oper_ndia(Num_site-j+1,i)=oper_ndia(i,Num_site-j+1)
		else if(FlagSign=='F') then !For Fermion
			!oper_ndia(Num_site-j+1,i)=-oper_ndia(i,Num_site-j+1)
			oper_ndia(Num_site-j+1,i)=oper_ndia(i,Num_site-j+1) !(check)
		endif
			
		!Save to disk
		open(10,file="oper_ndia_tmp.dat",position='append')
		write(10,*) i,Num_site-j+1,oper_ndia(i,Num_site-j+1)
		write(10,*) Num_site-j+1,i,oper_ndia(Num_site-j+1,i)
		close(10)
	end do
	end do

                        call block_transfer_trans(st_oper1, tmp1)
	call sys_site_env_site_cor_ndia(oper_ndia(sys_len,Num_site-env_len+1),st_oper,tmp1,sys_bs,env_bs,wave,FlagSign)
	if(FlagSign=='B') then !For Boson
		oper_ndia(Num_site-env_len+1,sys_len)=oper_ndia(sys_len,Num_site-env_len+1)
	else if(FlagSign=='F') then !For Fermion
		!oper_ndia(Num_site-env_len+1,sys_len)=-oper_ndia(sys_len,Num_site-env_len+1)
		oper_ndia(Num_site-env_len+1,sys_len)=oper_ndia(sys_len,Num_site-env_len+1) !(check)
	endif

	!Save to disk
	open(10,file="oper_ndia_tmp.dat",position='append')
	write(10,*) sys_len,Num_site-env_len+1,oper_ndia(sys_len,Num_site-env_len+1)
	write(10,*) Num_site-env_len+1,sys_len,oper_ndia(Num_site-env_len+1,sys_len)
	write(10,*)
	close(10)


	!Free space
	do i=1,sys_len-1
		call deallocate_block(sysopers(i))
	end do

	do i=1,env_len-1
		call deallocate_block(envopers(i))
	end do
	call deallocate_basis(sys_bs)
	call deallocate_basis(env_bs)


	!<5>: Save to the disk
	open(10,file=Filename,position='append')
	write(10,*) "Nx=",Nx,"Ny=",Ny,"N_up=",tot_num_up&
				&,"N_down=",tot_num_down,"M1=",kept_min,"M2=",kept_max
	do i=1,Num_site
		do j=1,Num_site
        oper_ndia(i,j)=oper_ndia(i,j)/dsqrt(1.0d0+st_oper%down_dif)
			write(10,*) i,j,oper_ndia(i,j)
		end do
	end do
	write(10,*)
	close(10)
111 format(I4,1X,I4,1X,F16.12)

end subroutine Get_oper_cor_ndia

!=========================================================================
subroutine Get_sys_oper_cor_ndia1(oper_ndia,st_oper,st_oper1,wave,trun_idx,FlagSign)
	use pubdata
	implicit none

	character(len=1) :: FlagSign
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper, st_oper1
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: oper_ndia(Num_site,Num_site)

	integer :: i,j,x,y,sys_len
	type(Total_Block) :: new_oper1,new_oper2
	type(Total_Block) :: tmp_oper,mid_oper,sys_oper,tmp1


	!<1>: For general information
	sys_len=wave%sys_len

	!Get g(i,j) for (i<j)
	do j=2,sys_len-1
		call update_site_ndia_sys(st_oper1,sys_bsm(j),new_oper2,FlagSign)

		!Update sysopers(sys_len-1)
		if(j==sys_len-1) then   !! different from new_oper2 
			call deallocate_block(sysopers(j))
			if(j<=trun_idx) then
				call update_site_ndia_sys(st_oper,sys_bsm(j),sysopers(j),FlagSign)
			else
				call update_site_ndia_sys(st_oper,sys_bsm(j),mid_oper,FlagSign)
				call update_trun_ndia(mid_oper,sysopers(j),systruns(j))
				call deallocate_block(mid_oper)
			endif
		endif

		do i=1,j-1
			if(i==j-1) then
				if(i>1) then !(i>1)
					call update_site_ndia_sys(st_oper,sys_bsm(i),mid_oper,FlagSign)
					if(i<=trun_idx) then
						call block_transfer(mid_oper,tmp_oper)
					else
						call update_trun_ndia(mid_oper,tmp_oper,systruns(i))
					endif
					call update_block_ndia_sys(tmp_oper,sys_bsm(j),new_oper1,FlagSign)
					call deallocate_block(mid_oper)
					call deallocate_block(tmp_oper)
				else !(i==1)
					call update_block_ndia_sys(sysopers(i),sys_bsm(j),new_oper1,FlagSign)
				endif
			else
                                                !! general case
				call update_block_ndia_sys(sysopers(i),sys_bsm(j),new_oper1,FlagSign)
			endif
			!Update sysopers(i)
			call deallocate_block(sysopers(i))
			if(j<=trun_idx) then
				call block_transfer(new_oper1,sysopers(i))
			else
				call update_trun_ndia(new_oper1,sysopers(i),systruns(j))
			endif
			
			!Multiply oper(i) and oper(j)

   if(i.lt.xi1)then    
                        call deallocate_block(new_oper1)
        go to 11
        endif
                                tmp_oper%down_dif=0
			!call block_mul_block_ndia(new_oper1,'N',new_oper2,'T',cone,tmp_oper,sys_bsm(j))

                                call block_transfer_trans(new_oper2, tmp1)
			call block_mul_block_ndia(new_oper1,'N',tmp1,'T',cone,tmp_oper,sys_bsm(j))
			call deallocate_block(new_oper1)
			call deallocate_block(tmp1)

			!Update tmp_oper(i,j) to new configuration
			!For (j<sys_len-1)
			if(j<sys_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,mid_oper)
				else
					call update_trun_dia(tmp_oper,mid_oper,systruns(j))
				endif
				call change_sys_oper_dia(mid_oper,j,sys_oper,sys_len-1,trun_idx,.true.)
				call deallocate_block(mid_oper)
			
			!For (j==sys_len-1)
			else if(j==sys_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,sys_oper)
				else
					call update_trun_dia(tmp_oper,sys_oper,systruns(sys_len-1))
				endif
			endif
			call deallocate_block(tmp_oper)

			!Measurement
			call measure_sys_block(oper_ndia(i,j),sys_oper,sys_bsm(sys_len),wave)
			if(FlagSign=='B') then !For Boson
				oper_ndia(j,i)=oper_ndia(i,j)
			else if(FlagSign=='F') then !For Fermion
				!oper_ndia(j,i)=-oper_ndia(i,j) !anti-commute
				oper_ndia(j,i)=oper_ndia(i,j) !anti-commute (check)
			endif

			!Save to disk
			open(10,file="oper_ndia_tmp.dat",position='append')
			write(10,*) i,j,oper_ndia(i,j)
			write(10,*) j,i,oper_ndia(j,i)
			close(10)
			
			call deallocate_block(sys_oper)

!			if(allocated(new_oper1%sub))call deallocate_block(new_oper1)

11      continue

		end do

		call deallocate_block(new_oper2)
	end do
111 format(I4,1X,I4,1X,F16.12)

end subroutine Get_sys_oper_cor_ndia1


!=========================================================================
!Get all density correlcations in system block using memory
!=========================================================================
subroutine Get_sys_oper_cor_ndia(oper_ndia,st_oper,st_oper1,wave,trun_idx,FlagSign)
	use pubdata
	implicit none

	character(len=1) :: FlagSign
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper, st_oper1
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: oper_ndia(Num_site,Num_site)

	integer :: i,j,x,y,sys_len
	type(Total_Block) :: new_oper1,new_oper2
	type(Total_Block) :: tmp_oper,mid_oper,sys_oper,tmp1


	!<1>: For general information
	sys_len=wave%sys_len

	!Get g(i,j) for (i<j)
	do j=1,sys_len-1
		call update_site_ndia_sys(st_oper1,sys_bsm(j),new_oper2,FlagSign)

		!Update sysopers(sys_len-1)
		if(j==sys_len-1) then   !! different from new_oper2 
			call deallocate_block(sysopers(j))
			if(j<=trun_idx) then
				call update_site_ndia_sys(st_oper,sys_bsm(j),sysopers(j),FlagSign)
			else
				call update_site_ndia_sys(st_oper,sys_bsm(j),mid_oper,FlagSign)
				call update_trun_ndia(mid_oper,sysopers(j),systruns(j))
				call deallocate_block(mid_oper)
			endif
		endif

		do i=1,j-1
			if(i==j-1) then
				if(i>1) then !(i>1)
					call update_site_ndia_sys(st_oper,sys_bsm(i),mid_oper,FlagSign)
					if(i<=trun_idx) then
						call block_transfer(mid_oper,tmp_oper)
					else
						call update_trun_ndia(mid_oper,tmp_oper,systruns(i))
					endif
					call update_block_ndia_sys(tmp_oper,sys_bsm(j),new_oper1,FlagSign)
					call deallocate_block(mid_oper)
					call deallocate_block(tmp_oper)
				else !(i==1)
					call update_block_ndia_sys(sysopers(i),sys_bsm(j),new_oper1,FlagSign)
				endif
			else
                                                !! general case
				call update_block_ndia_sys(sysopers(i),sys_bsm(j),new_oper1,FlagSign)
			endif
			!Update sysopers(i)
			call deallocate_block(sysopers(i))
			if(j<=trun_idx) then
				call block_transfer(new_oper1,sysopers(i))
			else
				call update_trun_ndia(new_oper1,sysopers(i),systruns(j))
			endif
			
			!Multiply oper(i) and oper(j)
                                tmp_oper%down_dif=0
			!call block_mul_block_ndia(new_oper1,'N',new_oper2,'T',cone,tmp_oper,sys_bsm(j))

                                call block_transfer_trans(new_oper2, tmp1)
			call block_mul_block_ndia(new_oper1,'N',tmp1,'T',cone,tmp_oper,sys_bsm(j))
			call deallocate_block(new_oper1)
			call deallocate_block(tmp1)

			!Update tmp_oper(i,j) to new configuration
			!For (j<sys_len-1)
			if(j<sys_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,mid_oper)
				else
					call update_trun_dia(tmp_oper,mid_oper,systruns(j))
				endif
				call change_sys_oper_dia(mid_oper,j,sys_oper,sys_len-1,trun_idx,.true.)
				call deallocate_block(mid_oper)
			
			!For (j==sys_len-1)
			else if(j==sys_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,sys_oper)
				else
					call update_trun_dia(tmp_oper,sys_oper,systruns(sys_len-1))
				endif
			endif
			call deallocate_block(tmp_oper)

			!Measurement
			call measure_sys_block(oper_ndia(i,j),sys_oper,sys_bsm(sys_len),wave)
			if(FlagSign=='B') then !For Boson
				oper_ndia(j,i)=oper_ndia(i,j)
			else if(FlagSign=='F') then !For Fermion
				!oper_ndia(j,i)=-oper_ndia(i,j) !anti-commute
				oper_ndia(j,i)=oper_ndia(i,j) !anti-commute (check)
			endif

			!Save to disk
			open(10,file="oper_ndia_tmp.dat",position='append')
			write(10,*) i,j,oper_ndia(i,j)
			write(10,*) j,i,oper_ndia(j,i)
			close(10)
			
			call deallocate_block(sys_oper)
		end do

		call deallocate_block(new_oper2)
	end do
111 format(I4,1X,I4,1X,F16.12)

end subroutine Get_sys_oper_cor_ndia


!=========================================================================
!Get all density correlcations in block block using memory
!=========================================================================

!=========================================================================
subroutine Get_env_oper_cor_ndia1(oper_ndia,st_oper,st_oper1,wave,trun_idx,FlagSign)
	use pubdata
	implicit none

	character(len=1) :: FlagSign
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper, st_oper1
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: oper_ndia(Num_site,Num_site)

	integer :: i,j,x,y,env_len
	type(Total_Block) :: new_oper1,new_oper2
	type(Total_Block) :: tmp_oper,mid_oper,env_oper


	!<1>: For general information
	env_len=wave%env_len

	!Get g(i,j) for (i<j)
	do j=2,env_len-1
		call update_site_ndia_env(st_oper,env_bsm(j),new_oper2,FlagSign)

		!Update envopers(env_len-1)
		if(j==env_len-1) then
			call deallocate_block(envopers(j))
			if(j<=trun_idx) then
				call update_site_ndia_env(st_oper1,env_bsm(j),envopers(j),FlagSign)
			else
				call update_site_ndia_env(st_oper1,env_bsm(j),mid_oper,FlagSign)
				call update_trun_ndia(mid_oper,envopers(j),envtruns(j))
				call deallocate_block(mid_oper)
			endif
		endif

		do i=1,j-1
			if(i==j-1) then
				if(i>1) then !(i>1)
					call update_site_ndia_env(st_oper1,env_bsm(i),mid_oper,FlagSign)
					if(i<=trun_idx) then
						call block_transfer(mid_oper,tmp_oper)
					else
						call update_trun_ndia(mid_oper,tmp_oper,envtruns(i))
					endif
					call update_block_ndia_env(tmp_oper,env_bsm(j),new_oper1,FlagSign)
					call deallocate_block(mid_oper)
					call deallocate_block(tmp_oper)
				else !(i==1)
					call update_block_ndia_env(envopers(i),env_bsm(j),new_oper1,FlagSign)
				endif
			else
				call update_block_ndia_env(envopers(i),env_bsm(j),new_oper1,FlagSign)
			endif

			!Update envopers(i)
			call deallocate_block(envopers(i))
			if(j<=trun_idx) then
				call block_transfer(new_oper1,envopers(i))
			else
				call update_trun_ndia(new_oper1,envopers(i),envtruns(j))
			endif
			
			!Multiply oper(i) and oper(j)

   if(i.lt.xj1)then
                        call deallocate_block(new_oper1)
        go to 11
        endif
                                tmp_oper%down_dif=0
					call block_transfer_trans(new_oper2,mid_oper)
			call block_mul_block_ndia(new_oper1,'N',mid_oper,'T',cone,tmp_oper,env_bsm(j))
			call deallocate_block(new_oper1)
			call deallocate_block(mid_oper)

			!Update tmp_oper(i,j) to new configuration
			!For (j<env_len-1)
			if(j<env_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,mid_oper)
				else
					call update_trun_dia(tmp_oper,mid_oper,envtruns(j))
				endif
				call change_env_oper_dia(mid_oper,j,env_oper,env_len-1,trun_idx,.true.)
				call deallocate_block(mid_oper)
			
			!For (j==env_len-1)
			else if(j==env_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,env_oper)
				else
					call update_trun_dia(tmp_oper,env_oper,envtruns(env_len-1))
				endif
			endif
			call deallocate_block(tmp_oper)

			!Measurement
			call measure_env_block(oper_ndia(Num_site-i+1,Num_site-j+1),env_oper,env_bsm(env_len),wave)
			if(FlagSign=='B') then !For Boson
				oper_ndia(Num_site-j+1,Num_site-i+1)=oper_ndia(Num_site-i+1,Num_site-j+1)
			else if(FlagSign=='F') then !For Fermion
	oper_ndia(Num_site-j+1,Num_site-i+1)=oper_ndia(Num_site-i+1,Num_site-j+1) 
			endif

			!Save to disk
			open(10,file="oper_ndia_tmp.dat",position='append')
			write(10,*) Num_site-i+1,Num_site-j+1,oper_ndia(Num_site-i+1,Num_site-j+1)
			write(10,*) Num_site-j+1,Num_site-i+1,oper_ndia(Num_site-j+1,Num_site-i+1)
			close(10)
			
			call deallocate_block(env_oper)

       11       continue

		end do

		call deallocate_block(new_oper2)
	end do
111 format(I4,1X,I4,1X,F16.12)

end subroutine Get_env_oper_cor_ndia1



!=========================================================================
!Get all density correlcations in block block using memory
!=========================================================================
subroutine Get_env_oper_cor_ndia(oper_ndia,st_oper,st_oper1,wave,trun_idx,FlagSign)
	use pubdata
	implicit none

	character(len=1) :: FlagSign
	integer,intent(in) :: trun_idx
	type(Total_Block),intent(in) :: st_oper, st_oper1
	type(Wavefunction),intent(in) :: wave
	double complex,intent(inout) :: oper_ndia(Num_site,Num_site)

	integer :: i,j,x,y,env_len
	type(Total_Block) :: new_oper1,new_oper2
	type(Total_Block) :: tmp_oper,mid_oper,env_oper


	!<1>: For general information
	env_len=wave%env_len

	!Get g(i,j) for (i<j)
	do j=2,env_len-1
		call update_site_ndia_env(st_oper,env_bsm(j),new_oper2,FlagSign)

		!Update envopers(env_len-1)
		if(j==env_len-1) then
			call deallocate_block(envopers(j))
			if(j<=trun_idx) then
				call update_site_ndia_env(st_oper1,env_bsm(j),envopers(j),FlagSign)
			else
				call update_site_ndia_env(st_oper1,env_bsm(j),mid_oper,FlagSign)
				call update_trun_ndia(mid_oper,envopers(j),envtruns(j))
				call deallocate_block(mid_oper)
			endif
		endif

		do i=1,j-1
			if(i==j-1) then
				if(i>1) then !(i>1)
					call update_site_ndia_env(st_oper1,env_bsm(i),mid_oper,FlagSign)
					if(i<=trun_idx) then
						call block_transfer(mid_oper,tmp_oper)
					else
						call update_trun_ndia(mid_oper,tmp_oper,envtruns(i))
					endif
					call update_block_ndia_env(tmp_oper,env_bsm(j),new_oper1,FlagSign)
					call deallocate_block(mid_oper)
					call deallocate_block(tmp_oper)
				else !(i==1)
					call update_block_ndia_env(envopers(i),env_bsm(j),new_oper1,FlagSign)
				endif
			else
				call update_block_ndia_env(envopers(i),env_bsm(j),new_oper1,FlagSign)
			endif

			!Update envopers(i)
			call deallocate_block(envopers(i))
			if(j<=trun_idx) then
				call block_transfer(new_oper1,envopers(i))
			else
				call update_trun_ndia(new_oper1,envopers(i),envtruns(j))
			endif
			
			!Multiply oper(i) and oper(j)
                                tmp_oper%down_dif=0
					call block_transfer_trans(new_oper2,mid_oper)
			call block_mul_block_ndia(new_oper1,'N',mid_oper,'T',cone,tmp_oper,env_bsm(j))
			call deallocate_block(new_oper1)
			call deallocate_block(mid_oper)

			!Update tmp_oper(i,j) to new configuration
			!For (j<env_len-1)
			if(j<env_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,mid_oper)
				else
					call update_trun_dia(tmp_oper,mid_oper,envtruns(j))
				endif
				call change_env_oper_dia(mid_oper,j,env_oper,env_len-1,trun_idx,.true.)
				call deallocate_block(mid_oper)
			
			!For (j==env_len-1)
			else if(j==env_len-1) then
				if(j<=trun_idx) then
					call block_transfer(tmp_oper,env_oper)
				else
					call update_trun_dia(tmp_oper,env_oper,envtruns(env_len-1))
				endif
			endif
			call deallocate_block(tmp_oper)

			!Measurement
			call measure_env_block(oper_ndia(Num_site-i+1,Num_site-j+1),env_oper,env_bsm(env_len),wave)
			if(FlagSign=='B') then !For Boson
				oper_ndia(Num_site-j+1,Num_site-i+1)=oper_ndia(Num_site-i+1,Num_site-j+1)
			else if(FlagSign=='F') then !For Fermion
				!oper_ndia(Num_site-j+1,Num_site-i+1)=-oper_ndia(Num_site-i+1,Num_site-j+1) !anti-commute
				oper_ndia(Num_site-j+1,Num_site-i+1)=oper_ndia(Num_site-i+1,Num_site-j+1) !anti-commute(check)
			endif

			!Save to disk
			open(10,file="oper_ndia_tmp.dat",position='append')
			write(10,111) Num_site-i+1,Num_site-j+1,oper_ndia(Num_site-i+1,Num_site-j+1)
			write(10,111) Num_site-j+1,Num_site-i+1,oper_ndia(Num_site-j+1,Num_site-i+1)
			close(10)
			
			call deallocate_block(env_oper)
		end do

		call deallocate_block(new_oper2)
	end do
111 format(I4,1X,I4,1X,F16.12)

end subroutine Get_env_oper_cor_ndia



!==================================================================================
!The following subroutines are used to calculate the electron-electron correlation
!function, which don't change the number of up spin and down spin electron
!==================================================================================
!Site-site corrrelation in system block [C^+_i(block 2nd).C_j(site 1st)] with (i<j)
subroutine sys_block_site_cor_ndia(value,sys_bl,sys_st,sys_bs,wave,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Wavefunction),intent(in) :: wave
	type(Total_Block),intent(in) :: sys_bl,sys_st
	type(Total_Basis),intent(in) :: sys_bs
        double complex,intent(inout) :: value

	integer :: i,j,k,x,y,up_dif,down_dif
	logical :: bl_flag,st_flag,bs_flag
	integer :: sys_num_up,sys_num_down,sys_dim
	integer :: env_num_up,env_num_down,env_dim
	integer :: bl_num_up,bl_num_down,st_num_up,st_num_down,tot_num
	integer :: new_bl_num_up,new_bl_num_down,new_st_num_up,new_st_num_down
	integer :: bl_id,st_id,bs_id,old_pos,new_pos,old_dim,new_dim

        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif
          real(8),external :: w6js, w3js, w9js
        real(8) ::  coef1,coef11, coef12, coef13, coef14

	double complex,allocatable :: mid(:,:)
	double complex :: coef,signs

	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif

	value=czero
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do j=1,sys_bs%num
			if(sys_bs%sub(j)%new_num_up==sys_num_up) then
			if(sys_bs%sub(j)%new_num_down==sys_num_down) then
				do k=1,sys_bs%sub(j)%num
					bl_num_up=sys_bs%sub(j)%sub(k)%bl_num_up
					bl_num_down=sys_bs%sub(j)%sub(k)%bl_num_down
					st_num_up=sys_bs%sub(j)%sub(k)%st_num_up
					st_num_down=sys_bs%sub(j)%sub(k)%st_num_down
					old_pos=sys_bs%sub(j)%sub(k)%spos
					old_dim=sys_bs%sub(j)%sub(k)%sdim

   do down_dif=-down_dif1, down_dif1, su  !! for sys_bl only
        do st_down_dif=-down_dif1, down_dif1, su  !! for sys_st

					new_bl_num_up=bl_num_up+up_dif
					new_bl_num_down=bl_num_down+down_dif
					new_st_num_up=st_num_up-up_dif
					new_st_num_down=st_num_down-st_down_dif
					
					bl_flag=.false.
					do x=1,sys_bl%num
						if(sys_bl%sub(x)%num_up==bl_num_up) then
						if(sys_bl%sub(x)%num_down==bl_num_down) then
						if(sys_bl%sub(x)%down_dif==down_dif) then
							bl_flag=.true.
							bl_id=x
							goto 101
						endif
						endif
						endif
					end do
					101 continue

					st_flag=.false.
					do x=1,sys_st%num
						if(sys_st%sub(x)%num_up==new_st_num_up) then
						if(sys_st%sub(x)%num_down==new_st_num_down) then
						if(sys_st%sub(x)%down_dif==st_down_dif) then
							st_flag=.true.
							st_id=x
							goto 102
						endif
						endif
						endif
					end do
					102 continue

					bs_flag=.false.
					do x=1,sys_bs%sub(j)%num
			if(sys_bs%sub(j)%sub(x)%bl_num_up==new_bl_num_up) then
			if(sys_bs%sub(j)%sub(x)%bl_num_down==new_bl_num_down) then
			if(sys_bs%sub(j)%sub(x)%st_num_up==new_st_num_up) then
			if(sys_bs%sub(j)%sub(x)%st_num_down==new_st_num_down) then
							new_pos=sys_bs%sub(j)%sub(x)%spos
							new_dim=sys_bs%sub(j)%sub(x)%sdim
							bs_flag=.true.
							goto 103
						endif
						endif
						endif
						endif
					end do
					103 continue

					if(bl_flag.and.st_flag.and.bs_flag) then
						tot_num=bl_num_up
						signs=czero
						if(FlagSign=='B') then !For Boson
							signs=cone
						else if(FlagSign=='F') then !For Fermion
							if(mod(tot_num,2)==0) then
								signs=cone
							else
								signs=-cone
							endif
						endif


!!block
        j1=new_bl_num_down
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down
        j5=new_st_num_down !! diag new_sys_st_num
        j6=sys_num_down
        coef1=1./dsqrt(1.0d0+down_dif1)*w6js(j1, j2, j3, j4,j5, j6)

        if(coef1.ne.0.0)then

        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)

						coef=signs*sys_st%sub(st_id)%mat(1,1)*coef1
						allocate(mid(new_dim,env_dim))
                if(realcode)then
						call DGEMM('N','N',new_dim,env_dim,old_dim,coef,sys_bl%sub(bl_id)%mat,new_dim&
								&,wave%sub(i)%vec(old_pos+1:old_pos+old_dim,1:env_dim),old_dim,czero,mid,new_dim)

                        else
	call ZGEMM('N','N',new_dim,env_dim,old_dim,coef,sys_bl%sub(bl_id)%mat,new_dim&
				&,wave%sub(i)%vec(old_pos+1:old_pos+old_dim,1:env_dim),old_dim,czero,mid,new_dim)

                        endif
						do x=1,env_dim
							value=value+dot_product(wave%sub(i)%vec(new_pos+1:new_pos+new_dim,x),mid(1:new_dim,x))
						end do
						deallocate(mid)
					endif
					endif
				end do
                        enddo
                        enddo
			endif
			endif
		end do
	end do

end subroutine sys_block_site_cor_ndia


!Site-site corrrelation in environment [C^+_i(site 2nd).C_j(block 1st)] with (i<j)
subroutine env_block_site_cor_ndia(value,env_bl,env_st,env_bs,wave,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Wavefunction),intent(in) :: wave
	type(Total_Block),intent(in) :: env_st,env_bl
	type(Total_Basis),intent(in) :: env_bs
        double complex,intent(inout) :: value

	integer :: i,j,k,x,y,up_dif,down_dif
	logical :: bl_flag,st_flag,bs_flag
	integer :: bl_num,st_num,sys_num,tot_num
	integer :: sys_num_up,sys_num_down,sys_dim
	integer :: env_num_up,env_num_down,env_dim
	integer :: bl_num_up,bl_num_down,st_num_up,st_num_down
	integer :: new_bl_num_up,new_bl_num_down,new_st_num_up,new_st_num_down
	integer :: bl_id,st_id,bs_id,old_pos,new_pos,old_dim,new_dim
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif, bl_down_dif
          real(8),external :: w6js, w3js, w9js
        real(8) ::  coef1,coef11, coef12, coef13, coef14
	double complex,allocatable :: mid(:,:)
	double complex :: coef
        real*8 signs

	up_dif=env_bl%up_dif
	down_dif1=env_bl%down_dif

	value=czero
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down
		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim

		do j=1,env_bs%num
			if(env_bs%sub(j)%new_num_up==env_num_up) then
			if(env_bs%sub(j)%new_num_down==env_num_down) then
				do k=1,env_bs%sub(j)%num
					bl_num_up=env_bs%sub(j)%sub(k)%bl_num_up
					bl_num_down=env_bs%sub(j)%sub(k)%bl_num_down
					st_num_up=env_bs%sub(j)%sub(k)%st_num_up
					st_num_down=env_bs%sub(j)%sub(k)%st_num_down
					old_pos=env_bs%sub(j)%sub(k)%spos
					old_dim=env_bs%sub(j)%sub(k)%sdim


   do bl_down_dif=-down_dif1, down_dif1, su
					!For C^+_site.C_block
					new_bl_num_up=bl_num_up-up_dif
					new_bl_num_down=bl_num_down-bl_down_dif

					
					bl_flag=.false.
					do x=1,env_bl%num
						if(env_bl%sub(x)%num_up==new_bl_num_up) then
						if(env_bl%sub(x)%num_down==new_bl_num_down) then
						if(env_bl%sub(x)%down_dif==bl_down_dif) then
							bl_flag=.true.
							bl_id=x
							goto 101
						endif
						endif
						endif
					end do
					101 continue



 do st_down_dif=-down_dif1, down_dif1, su
					new_st_num_up=st_num_up+up_dif
					new_st_num_down=st_num_down+st_down_dif

					st_flag=.false.
					do x=1,env_st%num
						if(env_st%sub(x)%num_up==st_num_up) then
						if(env_st%sub(x)%num_down==st_num_down) then
						if(env_st%sub(x)%down_dif==st_down_dif) then
							st_flag=.true.
							st_id=x
							goto 102
						endif
						endif
						endif
					end do
					102 continue

					bs_flag=.false.
					do x=1,env_bs%sub(j)%num
						if(env_bs%sub(j)%sub(x)%bl_num_up==new_bl_num_up) then
						if(env_bs%sub(j)%sub(x)%bl_num_down==new_bl_num_down) then
						if(env_bs%sub(j)%sub(x)%st_num_up==new_st_num_up) then
						if(env_bs%sub(j)%sub(x)%st_num_down==new_st_num_down) then
							new_pos=env_bs%sub(j)%sub(x)%spos
							new_dim=env_bs%sub(j)%sub(x)%sdim

							bs_flag=.true.
							bs_id=x
							goto 103
						endif
						endif
						endif
						endif
					end do
					103 continue

					if(bl_flag.and.st_flag.and.bs_flag) then
						!====== For electron sign =====================
						sys_num=sys_num_up !For c^_i(site)
						st_num=st_num_up
						bl_num=sys_num+st_num !For c_j(block)
						tot_num=bl_num+sys_num

						signs=czero
						if(FlagSign=='B') then !For Boson
							signs=cone
						else if(FlagSign=='F') then !For Fermion
							if(mod(tot_num,2)==0) then
								signs=cone
							else
								signs=-cone
							endif
						endif


  j1=new_bl_num_down
        j2=bl_num_down
        j3=down_dif1
        j4=st_num_down
        j5=new_st_num_down !! diag new_sys_st_num
        j6=env_num_down
        coef1=1./dsqrt(1.0d0+down_dif1)*w6js(j1, j2, j3, j4,j5, j6)
        coef1=coef1*(-1)**((j2+j5+j3+j6)/2)



                if(coef1.ne.0.0)then


						coef=signs*env_st%sub(st_id)%mat(1,1)*coef1
						allocate(mid(sys_dim,new_dim))

                if(realcode)then
		call DGEMM('N','N',sys_dim,new_dim,old_dim,coef&
		&,wave%sub(i)%vec(1:sys_dim,old_pos+1:old_pos+old_dim)&
			&,sys_dim,env_bl%sub(bl_id)%mat,old_dim,czero,mid,sys_dim)
                        else

           call ZGEMM('N','N',sys_dim,new_dim,old_dim,coef&
                                 &,wave%sub(i)%vec(1:sys_dim,old_pos+1:old_pos+old_dim)&
                                            &,sys_dim,env_bl%sub(bl_id)%mat,old_dim,czero,mid,sys_dim)
                        endif


						do x=1,new_dim
							value=value+dot_product(wave%sub(i)%vec(1:sys_dim,new_pos+x),mid(1:sys_dim,x))
						end do
						deallocate(mid)
					endif
                                endif
				end do
                        enddo
                        enddo
			endif
			endif
		end do
	end do

end subroutine env_block_site_cor_ndia


!Single site-site hopping corrrelation in system and environment block
!Note: [C^+_i(sys-block 2nd).C_j(env-block 1st)] with (i<j)
subroutine sys_block_env_block_cor_ndia(value,sys_bl,env_bl,sys_bs,env_bs,wave,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_bl
	type(Wavefunction),intent(in) :: wave
        double complex,intent(inout) :: value

	integer :: i,j,k,x,y,m,n,up_dif,down_dif
	integer :: sys_pos,sys_dim,new_sys_pos,new_sys_dim
	integer :: env_pos,env_dim,new_env_pos,new_env_dim
	integer :: sys_bs_id,env_bs_id,sys_id,env_id,wave_id
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,new_sys_bl_num_up,new_sys_bl_num_down
	integer :: env_bl_num_up,env_bl_num_down,new_env_bl_num_up,new_env_bl_num_down
	integer :: env_st_num_up,env_st_num_down,env_st_num,bl_num,sys_num,tot_num

        integer sys_st_num_up, sys_st_num_down

	logical :: sys_bs_flag,env_bs_flag,sys_flag,env_flag,wave_flag
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif, env_down_dif, bl_down_dif, bl_down_dif1
          real(8),external :: w6js, w3js, w9js
        real(8) ::  coef1,coef11, coef12, coef13, coef14
          complex*16,allocatable :: mat1(:,:),mat2(:,:)
        complex*16 :: coef
        real(8) :: signs

	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif

	value=czero
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down

		do j=1,sys_bs%num
			if(sys_bs%sub(j)%new_num_up==sys_num_up) then
			if(sys_bs%sub(j)%new_num_down==sys_num_down) then
				do k=1,env_bs%num
					if(env_bs%sub(k)%new_num_up==env_num_up) then
					if(env_bs%sub(k)%new_num_down==env_num_down) then
						do x=1,sys_bs%sub(j)%num
						do y=1,env_bs%sub(k)%num

							sys_bl_num_up=sys_bs%sub(j)%sub(x)%bl_num_up
							sys_bl_num_down=sys_bs%sub(j)%sub(x)%bl_num_down
							sys_st_num_up=sys_bs%sub(j)%sub(x)%st_num_up
							sys_st_num_down=sys_bs%sub(j)%sub(x)%st_num_down
							sys_pos=sys_bs%sub(j)%sub(x)%spos
							sys_dim=sys_bs%sub(j)%sub(x)%sdim

							env_bl_num_up=env_bs%sub(k)%sub(y)%bl_num_up
							env_bl_num_down=env_bs%sub(k)%sub(y)%bl_num_down
							env_pos=env_bs%sub(k)%sub(y)%spos
							env_dim=env_bs%sub(k)%sub(y)%sdim

							!For environment site
							env_st_num_up=env_bs%sub(k)%sub(y)%st_num_up
							env_st_num_down=env_bs%sub(k)%sub(y)%st_num_down
							env_st_num=env_st_num_up

                 do down_dif=-down_dif1, down_dif1, su

							new_sys_bl_num_up=sys_bl_num_up+up_dif
							new_sys_bl_num_down=sys_bl_num_down+down_dif

							sys_flag=.false.
							do m=1,sys_bl%num
				if(sys_bl%sub(m)%num_up==sys_bl_num_up) then
				if(sys_bl%sub(m)%num_down==sys_bl_num_down) then
				if(sys_bl%sub(m)%down_dif==down_dif) then
									sys_flag=.true.
									sys_id=m
									goto 101
								endif
								endif
								endif
							end do
							101 continue


 do env_down_dif=-down_dif1, down_dif1, su
                                new_env_bl_num_up=env_bl_num_up-up_dif
                                new_env_bl_num_down=env_bl_num_down-env_down_dif

							env_flag=.false.
						do m=1,env_bl%num
					if(env_bl%sub(m)%num_up==new_env_bl_num_up) then
					if(env_bl%sub(m)%num_down==new_env_bl_num_down) then
					if(env_bl%sub(m)%down_dif==env_down_dif) then
									env_flag=.true.
									env_id=m
									goto 102
								endif
								endif
								endif
							end do
							102 continue

 do new_sys_num_down= abs(sys_num_down-down_dif1), sys_num_down+down_dif1, su
					new_sys_num_up=sys_num_up+up_dif
					new_env_num_up=env_num_up-up_dif
					new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
							do m=1,sys_bs%num
								if(sys_bs%sub(m)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(m)%new_num_down==new_sys_num_down) then
									do n=1,sys_bs%sub(m)%num
								if(sys_bs%sub(m)%sub(n)%bl_num_up==new_sys_bl_num_up) then
								if(sys_bs%sub(m)%sub(n)%bl_num_down==new_sys_bl_num_down) then
									if(sys_bs%sub(m)%sub(n)%st_num_up==sys_st_num_up) then
									if(sys_bs%sub(m)%sub(n)%st_num_down==sys_st_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(m)%sub(n)%spos
											new_sys_dim=sys_bs%sub(m)%sub(n)%sdim
											goto 103
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do m=1,env_bs%num
								if(env_bs%sub(m)%new_num_up==new_env_num_up) then
								if(env_bs%sub(m)%new_num_down==new_env_num_down) then
									do n=1,env_bs%sub(m)%num
										if(env_bs%sub(m)%sub(n)%bl_num_up==new_env_bl_num_up) then
										if(env_bs%sub(m)%sub(n)%bl_num_down==new_env_bl_num_down) then
										if(env_bs%sub(m)%sub(n)%st_num_up==env_st_num_up) then
										if(env_bs%sub(m)%sub(n)%st_num_down==env_st_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(m)%sub(n)%spos
											new_env_dim=env_bs%sub(m)%sub(n)%sdim
											goto 104
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							wave_flag=.false.
							do m=1,wave%num
								if((wave%sub(m)%sys_num_up==new_sys_num_up).and.(wave%sub(m)%sys_num_down==new_sys_num_down)) then
								if((wave%sub(m)%env_num_up==new_env_num_up).and.(wave%sub(m)%env_num_down==new_env_num_down)) then
									wave_flag=.true.
									wave_id=m
									goto 105
								endif
								endif
							end do
							105 continue

							if(sys_flag.and.env_flag.and.sys_bs_flag.and.env_bs_flag.and.wave_flag) then
								!====== For electron sign ======
								sys_num=sys_num_up
								tot_num=sys_num+env_st_num_up

								signs=czero
								if(FlagSign=='B') then !For Boson
									signs=cone
								else if(FlagSign=='F') then !For Fermion
									if(mod(tot_num,2)==0) then
										signs=cone
									else
										signs=-cone
									endif
								endif


 coef1=czero

        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down !! new-diag
        j6=new_sys_bl_num_down

        j11=new_env_num_down  !! new-env
        j22=down_dif1
        j33=env_num_down       !! new-env-
        j44= env_bl_num_down
        j55=env_st_num_down
        j66=new_env_bl_num_down

        coef1=w6js(j1, j2, j3, j4,j5, j6)
        if(coef1.eq.0.0)go to 107

        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)
        if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))

        if(mod((j6+j5+j3+j2)/2+(j66+j55+j33+j22)/2, 2)==1)coef1=-coef1





					coef=signs*coef1
			allocate(mat1(new_sys_dim,env_dim),mat2(new_sys_dim,env_dim))

                        if(realcode)then
								call DGEMM('N','N',new_sys_dim,env_dim,sys_dim,cone,sys_bl%sub(sys_id)%mat,new_sys_dim&
										&,wave%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
										&,sys_dim,czero,mat1,new_sys_dim)
								call DGEMM('N','T',new_sys_dim,env_dim,new_env_dim,coef&
										&,wave%sub(wave_id)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim&
										&,env_bl%sub(env_id)%mat,env_dim,czero,mat2,new_sys_dim)

                else
        call ZGEMM('N','N',new_sys_dim,env_dim,sys_dim,cone,sys_bl%sub(sys_id)%mat,new_sys_dim&
                                         &,wave%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
                                                                         &,sys_dim,czero,mat1,new_sys_dim)
                                                 call ZGEMM('N','C',new_sys_dim,env_dim,new_env_dim,coef&
                                      &,wave%sub(wave_id)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim,new_env_pos+1:new_env_pos+new_env_dim),new_sys_dim&
                                                                                &,env_bl%sub(env_id)%mat,env_dim,czero,mat2,new_sys_dim)







                endif
								do m=1,env_dim
									value=value+dot_product(mat1(:,m),mat2(:,m))
								end do
								deallocate(mat1,mat2)


							endif
                        107 continue
                                endif
						end do
						end do
                                        enddo
                                        enddo
                                        enddo
                                        
					endif
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine sys_block_env_block_cor_ndia


!Single site-site hopping corrrelation in system and environment block
!Note: [C^+_i(sys-block 2nd).C_j(env-site 1st)] with (i<j)
subroutine sys_block_env_site_cor_ndia(value,sys_bl,env_st,sys_bs,env_bs,wave,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_bl,env_st
	type(Wavefunction),intent(in) :: wave
        double complex,intent(inout) :: value

	integer :: i,j,k,x,y,m,n,up_dif,down_dif
	integer :: sys_pos,sys_dim,new_sys_pos,new_sys_dim
	integer :: sys_bs_id,env_bs_id,sys_id,env_id,wave_id
	integer :: env_pos,env_dim,new_env_pos,sys_num,tot_num
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,new_sys_bl_num_up,new_sys_bl_num_down
	integer :: env_st_num_up,env_st_num_down,new_env_st_num_up,new_env_st_num_down
	logical :: sys_bs_flag,env_bs_flag,sys_flag,env_flag,wave_flag
	double complex,allocatable :: mat(:,:)
	double complex :: coef
        real*8 tmp_gl,signs
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer :: sys_st_num_down, bl_down_dif, env_bl_num_down
        integer j5, j6,j55, j66, st_down_dif, sys_st_num_up
          real(8),external :: w6js, w3js, w9js
        real(8) ::  coef1,coef11, coef12, coef13, coef14

	up_dif=sys_bl%up_dif
	down_dif1=sys_bl%down_dif

	value=czero
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down

		do j=1,sys_bs%num
			if(sys_bs%sub(j)%new_num_up==sys_num_up) then
			if(sys_bs%sub(j)%new_num_down==sys_num_down) then
				do k=1,env_bs%num
					if(env_bs%sub(k)%new_num_up==env_num_up) then
					if(env_bs%sub(k)%new_num_down==env_num_down) then
						do x=1,sys_bs%sub(j)%num
						do y=1,env_bs%sub(k)%num

							sys_bl_num_up=sys_bs%sub(j)%sub(x)%bl_num_up
							sys_bl_num_down=sys_bs%sub(j)%sub(x)%bl_num_down
							sys_st_num_down=sys_bs%sub(j)%sub(x)%st_num_down
							sys_st_num_up=sys_bs%sub(j)%sub(x)%st_num_up

							sys_pos=sys_bs%sub(j)%sub(x)%spos
							sys_dim=sys_bs%sub(j)%sub(x)%sdim

							env_st_num_up=env_bs%sub(k)%sub(y)%st_num_up
							env_st_num_down=env_bs%sub(k)%sub(y)%st_num_down
							env_bl_num_down=env_bs%sub(k)%sub(y)%bl_num_down
							env_pos=env_bs%sub(k)%sub(y)%spos
							env_dim=env_bs%sub(k)%sub(y)%sdim

                
                 do bl_down_dif=-sys_bl%down_dif, sys_bl%down_dif, su

							new_sys_bl_num_up=sys_bl_num_up+up_dif
							new_sys_bl_num_down=sys_bl_num_down+bl_down_dif


							sys_flag=.false.
							do m=1,sys_bl%num
								if(sys_bl%sub(m)%num_up==sys_bl_num_up) then
								if(sys_bl%sub(m)%num_down==sys_bl_num_down) then
								if(sys_bl%sub(m)%down_dif==bl_down_dif) then
									sys_flag=.true.
									sys_id=m
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                                         do st_down_dif=-env_st%down_dif, env_st%down_dif, su
							new_env_st_num_up=env_st_num_up-up_dif
							new_env_st_num_down=env_st_num_down-st_down_dif

							env_flag=.false.
							do m=1,env_st%num
								if(env_st%sub(m)%num_up==new_env_st_num_up) then
								if(env_st%sub(m)%num_down==new_env_st_num_down) then
								if(env_st%sub(m)%down_dif==st_down_dif) then
									env_flag=.true.
									env_id=m
									goto 102
								endif
								endif
								endif
							end do
							102 continue
                  if(sys_flag.and.env_flag) then
                   do new_sys_num_down=abs(sys_num_down-sys_bl%down_dif), sys_num_down+sys_bl%down_dif, su

							new_sys_num_up=sys_num_up+up_dif
							new_env_num_up=env_num_up-up_dif
							new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
							do m=1,sys_bs%num
								if(sys_bs%sub(m)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(m)%new_num_down==new_sys_num_down) then
									do n=1,sys_bs%sub(m)%num
										if(sys_bs%sub(m)%sub(n)%bl_num_up==new_sys_bl_num_up) then
										if(sys_bs%sub(m)%sub(n)%bl_num_down==new_sys_bl_num_down) then
										if(sys_bs%sub(m)%sub(n)%st_num_up==sys_st_num_up) then
										if(sys_bs%sub(m)%sub(n)%st_num_down==sys_st_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(m)%sub(n)%spos
											new_sys_dim=sys_bs%sub(m)%sub(n)%sdim
											goto 103
										endif
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do m=1,env_bs%num
								if(env_bs%sub(m)%new_num_up==new_env_num_up) then
								if(env_bs%sub(m)%new_num_down==new_env_num_down) then
									do n=1,env_bs%sub(m)%num
										if(env_bs%sub(m)%sub(n)%st_num_up==new_env_st_num_up) then
										if(env_bs%sub(m)%sub(n)%st_num_down==new_env_st_num_down) then
										if(env_bs%sub(m)%sub(n)%bl_num_down==env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(m)%sub(n)%spos
											goto 104
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							wave_flag=.false.
							do m=1,wave%num
								if((wave%sub(m)%sys_num_up==new_sys_num_up).and.(wave%sub(m)%sys_num_down==new_sys_num_down)) then
								if((wave%sub(m)%env_num_up==new_env_num_up).and.(wave%sub(m)%env_num_down==new_env_num_down)) then
									wave_flag=.true.
									wave_id=m
									goto 105
								endif
								endif
							end do
							105 continue

							if(sys_flag.and.env_flag.and.sys_bs_flag.and.env_bs_flag.and.wave_flag) then
								!====== For electron sign ======
								sys_num=sys_num_up
								tot_num=sys_num

								signs=czero
								if(FlagSign=='B') then !For Boson
									signs=cone
								else if(FlagSign=='F') then !For Fermion
									if(mod(tot_num,2)==0) then
										signs=cone
									else
										signs=-cone
									endif
								endif


        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_bl_num_down
        j5=sys_st_num_down !! diag
        j6=new_sys_bl_num_down
        j11=new_env_num_down  !! jenew    site
        j22=down_dif1         !! k
        j33=env_num_down       !! jenew'
        j44= env_st_num_down   !! je'
        j55=env_bl_num_down    !! jenv=jenv'
        j66=new_env_st_num_down    !! je
        coef1=w6js(j1, j2, j3, j4,j5, j6)
        if(coef1.eq.0.0)go to 110

        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)
                if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.d0+j11)*(1.d0+j33)/(1.0d0+down_dif1))
        if(mod((j6+j5+j3+j2)/2+(j55++j44+j11+j22)/2,2).ne.0) coef1=-coef1



				coef=signs*env_st%sub(env_id)%mat(1,1)*coef1
					allocate(mat(new_sys_dim,env_dim))

                if(realcode)then
				call DGEMM('N','N',new_sys_dim,env_dim,sys_dim,coef,sys_bl%sub(sys_id)%mat,new_sys_dim&
						&,wave%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
										&,sys_dim,czero,mat,new_sys_dim)
                        else

         call ZGEMM('N','N',new_sys_dim,env_dim,sys_dim,coef,sys_bl%sub(sys_id)%mat,new_sys_dim&
                  &,wave%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+1:env_pos+env_dim)&
                                                                                &,sys_dim,czero,mat,new_sys_dim)

                       endif

								tmp_gl=czero
								do m=1,env_dim
	value=value+dot_product(wave%sub(wave_id)%vec(new_sys_pos+1:new_sys_pos+new_sys_dim,new_env_pos+m),mat(1:new_sys_dim,m))
								end do
								deallocate(mat)
							endif
                                                endif
                110             continue
						end do
                                                endif
						end do
                                        enddo
                                        enddo
                                        enddo
                                        endif
                                        endif
                                        enddo
					endif
					endif
				end do
                        enddo

end subroutine sys_block_env_site_cor_ndia


!Single site-site hopping corrrelation in system and environment block
!Note: [C^+_i(sys-site 2nd).C_j(env-site 1st)] with (i<j)
subroutine sys_site_env_site_cor_ndia(value,sys_st,env_st,sys_bs,env_bs,wave,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_st
	type(Wavefunction),intent(in) :: wave
        double complex,intent(inout) :: value

	integer :: i,j,k,x,y,m,n,up_dif,down_dif
	integer :: sys_pos,sys_dim,new_sys_pos
	integer :: env_pos,env_dim,new_env_pos
	integer :: sys_bs_id,env_bs_id,sys_id,env_id,wave_id
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_bl_num_up,sys_bl_num_down,sys_bl_num,sys_num,tot_num
	integer :: sys_st_num_up,sys_st_num_down,new_sys_st_num_up,new_sys_st_num_down
	integer :: env_st_num_up,env_st_num_down,new_env_st_num_up,new_env_st_num_down
	logical :: sys_bs_flag,env_bs_flag,sys_flag,env_flag,wave_flag
	real(8) :: coef,tmp_gl,signs
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66,   st_down_dif1, st_down_dif2, env_bl_num_down, env_bl_num_up

          real(8),external :: w6js, w3js, w9js
        real(8) ::  coef1,coef11, coef12, coef13, coef14

	up_dif=sys_st%up_dif
	down_dif1=sys_st%down_dif

	value=czero
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down

		do j=1,sys_bs%num
			if(sys_bs%sub(j)%new_num_up==sys_num_up) then
			if(sys_bs%sub(j)%new_num_down==sys_num_down) then
				do k=1,env_bs%num
					if(env_bs%sub(k)%new_num_up==env_num_up) then
					if(env_bs%sub(k)%new_num_down==env_num_down) then
						do x=1,sys_bs%sub(j)%num
						do y=1,env_bs%sub(k)%num

							sys_st_num_up=sys_bs%sub(j)%sub(x)%st_num_up
							sys_st_num_down=sys_bs%sub(j)%sub(x)%st_num_down
							sys_pos=sys_bs%sub(j)%sub(x)%spos
							sys_dim=sys_bs%sub(j)%sub(x)%sdim
							sys_bl_num_up=sys_bs%sub(j)%sub(x)%bl_num_up
							sys_bl_num_down=sys_bs%sub(j)%sub(x)%bl_num_down
							sys_bl_num=sys_bl_num_up

							env_st_num_up=env_bs%sub(k)%sub(y)%st_num_up
							env_st_num_down=env_bs%sub(k)%sub(y)%st_num_down
							env_bl_num_up=env_bs%sub(k)%sub(y)%bl_num_up
							env_bl_num_down=env_bs%sub(k)%sub(y)%bl_num_down
							env_pos=env_bs%sub(k)%sub(y)%spos
							env_dim=env_bs%sub(k)%sub(y)%sdim



                          do st_down_dif1=-sys_st%down_dif, sys_st%down_dif, su
							new_sys_st_num_up=sys_st_num_up+up_dif
							new_sys_st_num_down=sys_st_num_down+st_down_dif1

							sys_flag=.false.
							do m=1,sys_st%num
								if(sys_st%sub(m)%num_up==sys_st_num_up) then
								if(sys_st%sub(m)%num_down==sys_st_num_down) then
								if(sys_st%sub(m)%down_dif==st_down_dif1) then
									sys_flag=.true.
									sys_id=m
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                          do st_down_dif2=-sys_st%down_dif, sys_st%down_dif, su

							new_env_st_num_up=env_st_num_up-up_dif
							new_env_st_num_down=env_st_num_down-st_down_dif2
							env_flag=.false.
							do m=1,env_st%num
								if(env_st%sub(m)%num_up==new_env_st_num_up) then
								if(env_st%sub(m)%num_down==new_env_st_num_down) then
								if(env_st%sub(m)%down_dif==st_down_dif2) then
									env_flag=.true.
									env_id=m
									goto 102
								endif
								endif
								endif
							end do
							102 continue

                                 if(sys_flag.and.env_flag) then
                do new_sys_num_down=abs(sys_num_down-sys_st%down_dif), sys_num_down+sys_st%down_dif, su
							new_sys_num_up=sys_num_up+up_dif
							new_env_num_up=env_num_up-up_dif
							new_env_num_down=new_sys_num_down

							sys_bs_flag=.false.
							do m=1,sys_bs%num
								if(sys_bs%sub(m)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(m)%new_num_down==new_sys_num_down) then
									do n=1,sys_bs%sub(m)%num
										if(sys_bs%sub(m)%sub(n)%st_num_up==new_sys_st_num_up) then
										if(sys_bs%sub(m)%sub(n)%st_num_down==new_sys_st_num_down) then
										if(sys_bs%sub(m)%sub(n)%bl_num_down==sys_bl_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(m)%sub(n)%spos
											goto 103
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do m=1,env_bs%num
								if(env_bs%sub(m)%new_num_up==new_env_num_up) then
								if(env_bs%sub(m)%new_num_down==new_env_num_down) then
									do n=1,env_bs%sub(m)%num
										if(env_bs%sub(m)%sub(n)%st_num_up==new_env_st_num_up) then
										if(env_bs%sub(m)%sub(n)%st_num_down==new_env_st_num_down) then
										if(env_bs%sub(m)%sub(n)%bl_num_down==env_bl_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(m)%sub(n)%spos
											goto 104
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							wave_flag=.false.
							do m=1,wave%num
								if((wave%sub(m)%sys_num_up==new_sys_num_up).and.(wave%sub(m)%sys_num_down==new_sys_num_down)) then
								if((wave%sub(m)%env_num_up==new_env_num_up).and.(wave%sub(m)%env_num_down==new_env_num_down)) then
									wave_flag=.true.
									wave_id=m
									goto 105
								endif
								endif
							end do
							105 continue

							if(sys_flag.and.env_flag.and.sys_bs_flag.and.env_bs_flag.and.wave_flag) then
								!====== For electron sign ========================
								sys_num=sys_num_up !For c_j(env_site)
								tot_num=sys_num+sys_bl_num !For c^+_i(sys_site)

								signs=czero
								if(FlagSign=='B') then !For Boson
									signs=cone
								else if(FlagSign=='F') then !For Fermion
									if(mod(tot_num,2)==0) then
										signs=cone
									else
										signs=-cone
									endif
								endif
								coef=signs*sys_st%sub(sys_id)%mat(1,1)*env_st%sub(env_id)%mat(1,1)
        if(coef.ne.0.0)then
        j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_st_num_down
        j5=sys_bl_num_down !! diag new_sys_bl
        j6=new_sys_st_num_down
!! site
        j11=new_env_num_down  !! jenew    site
        j22=down_dif1         !! k
        j33=env_num_down       !! jenew'
        j44= env_st_num_down   !! je'
        j55=env_bl_num_down    !! jenv=jenv'
        j66=new_env_st_num_down    !! je

        coef1=w6js(j1, j2, j3, j4,j5, j6)
        if(coef1.eq.0.0)go to 110
        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)

        if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(1.0d0+j33)/(1.0d0+down_dif1))
        if(mod((j55+j44+j11+j22)/2+(j5+j4+j1+j2)/2, 2).ne.0)coef1=-coef1
       coef=coef*coef1


								tmp_gl=czero
								do m=1,env_dim
									tmp_gl=tmp_gl+dot_product(wave%sub(wave_id)%vec(new_sys_pos+1:new_sys_pos+sys_dim&
											&,new_env_pos+m),wave%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+m))
								end do
								value=value+tmp_gl*coef
							endif

                110             continue
							endif
							endif
						end do
                                                endif
						end do
                                        enddo
                                        enddo
                                        enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine sys_site_env_site_cor_ndia


!Single site-site hopping corrrelation in system and environment block
!Note: [C^+_i(sys-site 2nd).C_j(env-block 1st)] with (i<j)
subroutine sys_site_env_block_cor_ndia(value,sys_st,env_bl,sys_bs,env_bs,wave,FlagSign)
	use pubdata
	implicit none

	character(len=1),intent(in) :: FlagSign
	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_st,env_bl
	type(Wavefunction),intent(in) :: wave
        double complex,intent(inout) :: value

	integer :: i,j,k,x,y,m,n,up_dif,down_dif
	integer :: sys_pos,sys_dim,new_sys_pos
	integer :: env_pos,env_dim,new_env_pos,new_env_dim
	integer :: sys_bs_id,env_bs_id,sys_id,env_id,wave_id
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down
	integer :: sys_st_num_up,sys_st_num_down,new_sys_st_num_up,new_sys_st_num_down
	integer :: env_bl_num_up,env_bl_num_down,new_env_bl_num_up,new_env_bl_num_down
	integer :: env_st_num_up,env_st_num_down,env_st_num,sys_num,tot_num
	integer :: sys_bl_num_up,sys_bl_num_down,sys_bl_num
	logical :: sys_bs_flag,env_bs_flag,sys_flag,env_flag,wave_flag
	double complex,allocatable :: mat(:,:)
	double complex  :: coef
        real*8 signs
        integer :: down_dif1, j1,j2,j3,j4,j11,j22,j33,j44, j12, j34, j12n, j34n, j24
        integer j5, j6,j55, j66, st_down_dif,   bl_down_dif
          real(8),external :: w6js, w3js, w9js
        real(8) ::  coef1,coef11, coef12, coef13, coef14

	up_dif=sys_st%up_dif
	down_dif1=sys_st%down_dif

	value=czero
	do i=1,wave%num
		sys_num_up=wave%sub(i)%sys_num_up
		sys_num_down=wave%sub(i)%sys_num_down
		env_num_up=wave%sub(i)%env_num_up
		env_num_down=wave%sub(i)%env_num_down

		do j=1,sys_bs%num
			if(sys_bs%sub(j)%new_num_up==sys_num_up) then
			if(sys_bs%sub(j)%new_num_down==sys_num_down) then
				do k=1,env_bs%num
					if(env_bs%sub(k)%new_num_up==env_num_up) then
					if(env_bs%sub(k)%new_num_down==env_num_down) then
						do x=1,sys_bs%sub(j)%num
						do y=1,env_bs%sub(k)%num
							sys_st_num_up=sys_bs%sub(j)%sub(x)%st_num_up
							sys_st_num_down=sys_bs%sub(j)%sub(x)%st_num_down
							sys_pos=sys_bs%sub(j)%sub(x)%spos
							sys_dim=sys_bs%sub(j)%sub(x)%sdim
							sys_bl_num_up=sys_bs%sub(j)%sub(x)%bl_num_up
							sys_bl_num_down=sys_bs%sub(j)%sub(x)%bl_num_down
							sys_bl_num=sys_bl_num_up

							env_bl_num_up=env_bs%sub(k)%sub(y)%bl_num_up
							env_bl_num_down=env_bs%sub(k)%sub(y)%bl_num_down
							env_pos=env_bs%sub(k)%sub(y)%spos
							env_dim=env_bs%sub(k)%sub(y)%sdim
							env_st_num_up=env_bs%sub(k)%sub(y)%st_num_up
							env_st_num_down=env_bs%sub(k)%sub(y)%st_num_down
							env_st_num=env_st_num_up


                                        do st_down_dif=-sys_st%down_dif, sys_st%down_dif, su
							new_sys_st_num_up=sys_st_num_up+up_dif
							new_sys_st_num_down=sys_st_num_down+st_down_dif


							sys_flag=.false.
							do m=1,sys_st%num
								if(sys_st%sub(m)%num_up==sys_st_num_up) then
								if(sys_st%sub(m)%num_down==sys_st_num_down) then
								if(sys_st%sub(m)%down_dif==st_down_dif) then
									sys_flag=.true.
									sys_id=m
									goto 101
								endif
								endif
								endif
							end do
							101 continue

                                        do bl_down_dif=-env_bl%down_dif, env_bl%down_dif, su
							env_flag=.false.
							new_env_bl_num_up=env_bl_num_up-up_dif
							new_env_bl_num_down=env_bl_num_down-bl_down_dif
							do m=1,env_bl%num
								if(env_bl%sub(m)%num_up==new_env_bl_num_up) then
								if(env_bl%sub(m)%num_down==new_env_bl_num_down) then
								if(env_bl%sub(m)%down_dif==bl_down_dif) then
									env_flag=.true.
									env_id=m
									goto 102
								endif
								endif
								endif
							end do
							102 continue

  if(sys_flag.and.env_flag) then
           do new_sys_num_down=abs(sys_num_down-sys_st%down_dif),sys_num_down+sys_st%down_dif, su
							new_sys_num_up=sys_num_up+up_dif
							new_env_num_up=env_num_up-up_dif
							new_env_num_down=new_sys_num_down
							sys_bs_flag=.false.
							do m=1,sys_bs%num
								if(sys_bs%sub(m)%new_num_up==new_sys_num_up) then
								if(sys_bs%sub(m)%new_num_down==new_sys_num_down) then
									do n=1,sys_bs%sub(m)%num
										if(sys_bs%sub(m)%sub(n)%st_num_up==new_sys_st_num_up) then
										if(sys_bs%sub(m)%sub(n)%st_num_down==new_sys_st_num_down) then
										if(sys_bs%sub(m)%sub(n)%bl_num_down==sys_bl_num_down) then
											sys_bs_flag=.true.
											new_sys_pos=sys_bs%sub(m)%sub(n)%spos
											goto 103
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							103 continue

							env_bs_flag=.false.
							do m=1,env_bs%num
								if(env_bs%sub(m)%new_num_up==new_env_num_up) then
								if(env_bs%sub(m)%new_num_down==new_env_num_down) then
									do n=1,env_bs%sub(m)%num
										if(env_bs%sub(m)%sub(n)%bl_num_up==new_env_bl_num_up) then
										if(env_bs%sub(m)%sub(n)%bl_num_down==new_env_bl_num_down) then
										if(env_bs%sub(m)%sub(n)%st_num_down==env_st_num_down) then
											env_bs_flag=.true.
											new_env_pos=env_bs%sub(m)%sub(n)%spos
											new_env_dim=env_bs%sub(m)%sub(n)%sdim
											goto 104
										endif
										endif
										endif
									end do
								endif
								endif
							end do
							104 continue

							wave_flag=.false.
							do m=1,wave%num
								if((wave%sub(m)%sys_num_up==new_sys_num_up).and.(wave%sub(m)%sys_num_down==new_sys_num_down)) then
								if((wave%sub(m)%env_num_up==new_env_num_up).and.(wave%sub(m)%env_num_down==new_env_num_down)) then
									wave_flag=.true.
									wave_id=m
									goto 105
								endif
								endif
							end do
							105 continue

							if(sys_flag.and.env_flag.and.sys_bs_flag.and.env_bs_flag.and.wave_flag) then
								!====== For electron sign =====================
								sys_num=sys_num_up 
								tot_num=sys_num+env_st_num_up !For c_j(env_block)
								tot_num=tot_num+sys_bl_num_up !For c^+_i(sys_site)

								signs=czero
								if(FlagSign=='B') then !For Boson
									signs=cone
								else if(FlagSign=='F') then !For Fermion
									if(mod(tot_num,2)==0) then
										signs=cone
									else
										signs=-cone
									endif
								endif

								coef=signs*sys_st%sub(sys_id)%mat(1,1)



        if(coef.ne.0.0)then

 j1=new_sys_num_down
        j2=down_dif1
        j3=sys_num_down
        j4= sys_st_num_down
        j5=sys_bl_num_down !! diag new_sys_bl
        j6=new_sys_st_num_down
!! block
        j11=new_env_num_down  !! new-env
        j22=down_dif1
        j33=env_num_down       !! new-env-
        j44= env_bl_num_down
        j55=env_st_num_down !! diag new_env_st_num_down
        j66=new_env_bl_num_down

        coef1=w6js(j1, j2, j3, j4,j5, j6)
        if(coef1.eq.0.0)go to 110
        coef1=coef1*w6js(j11, j22, j33, j44,j55, j66)

        if(coef1.ne.0.0)then
        coef1=coef1*dsqrt((1.0d0+j11)*(j33+1.d0)/(1.0d0+down_dif1))
        if(mod((j66+j55+j33+j22)/2+(j5+j4+j1+j2)/2,2).ne.0)coef1=-coef1
        coef=coef*coef1

								allocate(mat(sys_dim,env_dim))
                if(realcode)then

								call DGEMM('N','T',sys_dim,env_dim,new_env_dim,coef&
										&,wave%sub(wave_id)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+new_env_dim)&
										&,sys_dim,env_bl%sub(env_id)%mat,env_dim,czero,mat,sys_dim)

                else
     call ZGEMM('N','C',sys_dim,env_dim,new_env_dim,coef&
                 &,wave%sub(wave_id)%vec(new_sys_pos+1:new_sys_pos+sys_dim,new_env_pos+1:new_env_pos+new_env_dim)&
                                  &,sys_dim,env_bl%sub(env_id)%mat,env_dim,czero,mat,sys_dim)
                endif




								do m=1,env_dim
									value=value+dot_product(mat(1:sys_dim,m),wave%sub(i)%vec(sys_pos+1:sys_pos+sys_dim,env_pos+m))
								end do
								deallocate(mat)

							endif
                        110             continue
							endif
							endif
						end do
                                                endif
						end do
                                        enddo
                                        enddo
                                        enddo
					endif
					endif
				end do
			endif
			endif
		end do
	end do

end subroutine sys_site_env_block_cor_ndia


!=========================================================================================
subroutine Get_SCOP_Operator_Sys(idx1,idx2,new_oper,new_idx,trun_idx,truns,Flags)
	use pubdata
	implicit none

	logical,intent(in) :: truns
	character(len=2),intent(in) :: Flags
	integer,intent(in) :: idx1,idx2,new_idx,trun_idx
	type(Total_Block),intent(inout) :: new_oper

	real(8) :: coefs
	integer :: i,j,x,y,sys_len
	type(Total_Block) :: new_oper1,new_oper2,sys_oper
	type(Total_Block) :: mid_oper1,mid_oper2,tmp_oper

	!<0>: Check the relation of idx1 and idx2
	if(idx1>=idx2) then
		write(*,*) "idx1>=idx2 in Get_SCOP_Operator_Sys"
		return
	endif

	!<1>: For singlet SCOP: (Flags='S0')
	if(Flags=='S0') then
		!<1a>: Get O^+_12=(C^+_{i,up}.C^+_{j,down})
		call Get_sys_oper_ndia(st_elec_up,idx1,new_oper1,idx2,trun_idx,.false.,'F')
		call Get_sys_oper_ndia(st_elec_up,idx2,new_oper2,idx2,trun_idx,.false.,'F')
                        new_oper%down_dif=0 !! rank 0 for SC
		call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',cone,new_oper,sys_bsm(idx2))
		call deallocate_block(new_oper1)
		call deallocate_block(new_oper2)
	endif


	!<2>: For triplet SCOP: Flags='T0'
	if(Flags=='T0') then !! using X^1 later

		!<2c>: Get O^+=(O^+_12 + O^+_21)/sqrt(2)
		call Get_sys_oper_ndia(st_elec_up,idx1,new_oper1,idx2,trun_idx,.false.,'F')
		call Get_sys_oper_ndia(st_elec_up,idx2,new_oper2,idx2,trun_idx,.false.,'F')
                        new_oper%down_dif=2 !! rank 0 for SC
		call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',cone,new_oper,sys_bsm(idx2))
		call deallocate_block(new_oper1)
		call deallocate_block(new_oper2)
	endif

	!<5>: Update the operator in idx2 to new_idx
	if(idx2<=trun_idx) then
		call block_transfer(new_oper,sys_oper)
	else
		call update_trun_ndia(new_oper,sys_oper,systruns(idx2))
	endif
	call deallocate_block(new_oper)

	call Change_sys_oper_ndia(sys_oper,idx2,new_oper,new_idx,trun_idx,Truns,'B')
	call deallocate_block(sys_oper)

end subroutine Get_SCOP_Operator_Sys



!=========================================================================================
!<1>: Get singlet superconducting order parameter: (i<j)
!<S0>: O^+_{ij}=(C^+_{i,S1}.C^+_{j,S2}-C^+_{i,S2).C^+_{j,S1})/sqrt(2) with (S1=up,S2=down)
!<2>: Get triplet superconducting order parameter: (i<j)
!<T0>: O^+_{ij}=(C^+_{i,S1}.C^+_{j,S2}+C^+_{i,S2).C^+_{j,S1})/sqrt(2) with (S1=up,S2=down)
!<TU>: O^+_{ij}=C^+_{i,S1}.C^+_{j,S2} with (S1=S2=up, for S_tot=1)
!<TD>: O^+_{ij}=C^+_{i,S1}.C^+_{j,S2} with (S1=S2=down, for S_tot=-1)
!=========================================================================================
subroutine Get_SCOP_Operator_Env(idx1,idx2,new_oper,new_idx,trun_idx,truns,Flags)
	use pubdata
	implicit none

	logical,intent(in) :: truns
	character(len=2),intent(in) :: Flags
	integer,intent(in) :: idx1,idx2,new_idx,trun_idx
	type(Total_Block),intent(inout) :: new_oper

	real(8) :: coefs
	integer :: i,j,x,y,env_len
	type(Total_Block) :: new_oper1,new_oper2,env_oper
	type(Total_Block) :: mid_oper1,mid_oper2,tmp_oper

	!<0>: Check the relation of idx1 and idx2
	if(idx1>=idx2) then
		write(*,*) "idx1>=idx2 in Get_SCOP_Operator_Env"
		return
	endif

	!<1>: For singlet SCOP: (Flags='S0')
	if(Flags=='S0') then
		!<1a>: Get O^+_12=(C^+_{i,up}.C^+_{j,down})
		call Get_env_oper_ndia(st_elec_up,idx1,new_oper1,idx2,trun_idx,.false.,'F')
		call Get_env_oper_ndia(st_elec_up,idx2,new_oper2,idx2,trun_idx,.false.,'F')
                        new_oper%down_dif=0
		call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',cone,new_oper,env_bsm(idx2))
		call deallocate_block(new_oper1)
		call deallocate_block(new_oper2)
	endif


	!<2>: For triplet SCOP: Flags='T0'
	if(Flags=='T0') then
		!<2a>: Get O^+_12=(C^+_{i,up}.C^+_{j,down})
		call Get_env_oper_ndia(st_elec_up,idx1,new_oper1,idx2,trun_idx,.false.,'F')
		call Get_env_oper_ndia(st_elec_up,idx2,new_oper2,idx2,trun_idx,.false.,'F')
                        new_oper%down_dif=0
		call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',cone,new_oper,env_bsm(idx2))
		call deallocate_block(new_oper1)
		call deallocate_block(new_oper2)
	endif


	!<5>: Update the operator in idx2 to new_idx
	if(idx2<=trun_idx) then
		call block_transfer(new_oper,env_oper)
	else
		call update_trun_ndia(new_oper,env_oper,envtruns(idx2))
	endif
	call deallocate_block(new_oper)

	call Change_env_oper_ndia(env_oper,idx2,new_oper,new_idx,trun_idx,Truns,'B')
	call deallocate_block(env_oper)

end subroutine Get_SCOP_Operator_Env



!call Get_Chiral_cor (eigvec,trun_idx,Filename)
subroutine Get_Chiral_cor(wave,trun_idx,Filename)

    use pubdata
	implicit none
!    complex*16,intent(inout) : chiral_cor
    type(Wavefunction),intent(in) :: wave
    integer,intent(inout) :: trun_idx
    character(len=30) :: Filename

    integer :: sid1,sid2,sid3
    integer :: eid1,eid2,eid3
    integer :: tri_x1,tri_y1,tri_x2,tri_y2
    integer :: x,y,i,ij,ik,j,i1, j1,sys_len,env_len,sx,sy,ex,ey
    type(Total_Basis) :: sys_bs,env_bs
    type(Total_Block) ::sys_oper,env_oper
    logical :: sys_flag1,env_flag1,sys_flag2,env_flag2
    real(8) :: value
        integer  chis1, chis2

        chis1=1
        chis2=1

    sys_len=wave%sys_len
	env_len=wave%env_len
        trun_idx=1
        go to 11

         !<1-2>: Read truncation operators from disk
        do i=trun_idx+1,sys_len
                call truns_from_disk(systruns(i),1001,i,.true.)
        end do

        do i=trun_idx+1,env_len
                call truns_from_disk(envtruns(i),1001,i,.false.)
        end do

        !<1-3>: Read basis from disk
        do i=2,sys_len
                call basis_from_disk(sys_bsm(i),1001,i,.true.)
        end do

        do i=2,env_len
                call basis_from_disk(env_bsm(i),1001,i,.false.)
        end do


11              continue

	call basis_transfer(sys_bsm(sys_len),sys_bs)
	call basis_transfer(env_bsm(env_len),env_bs)

            !<2>: Open file and save to disk
	        open(10,file=Filename,position='append')
	        write(10,110) "Nx=",Nx,"Ny=",Ny,"Jd=",jd(1,1),"t2=",jt(1,1),"N_up=",tot_num_up&
				        &,"N_down=",tot_num_down
	        110	format(A3,I2,1X,A3,I2,1X,A2,F6.3,1X,A3,F6.3,1X,A5,I3,1X,A7,I3)
	        !!!write(10,*) "hubbard", hubbard(1)


        if(kagome==1)then
do sid1=1,sys_len-1,3   !!!!!Lattice(1,sys_len)
        sid2= sid1+1
         sid3=sid2+1
                 sys_flag1 =.false.
                if(max(sid2,sid3)<=(sys_len-1)) then
                sys_flag1 =.true.
		        call Get_sys_triangular (sid1,sid2,sid3,sys_oper,trun_idx,wave)
                         !   write(*,*) "Get_sys_triangular "
                 endif

         do eid1=1, env_len-1, 3
              eid2=eid1+1
              eid3=eid2+1
 
          env_flag1 =.false.
           if(max(eid2,eid3)<=(env_len-1)) then
            env_flag1 =.true.
		     call Get_env_triangular (eid1,eid2,eid3,env_oper,trun_idx,wave)
             !write(*,*) "Get_env_triangular "
           endif

        	        if(sys_flag1) then
                    if(env_flag1) then
					     call sys_block_env_block_cor_dia(value,sys_oper,env_oper,sys_bs,env_bs,wave)
 					     write(101,111) "i=",sx,sy, "j=",ex,ey,"r=",ex-sx,-value*(-0.25d0)
 					     write(10,15) sid1, sid2, num_site-eid1+1, num_site-eid2+1,-value*(-0.25d0)
 				   endif
                   endif
                   111 format(A2,1X,I2,I2,2X,A2,1X,I2,I2,2X,A2,1X,I2,2X,F16.10)
                   15 format(4i8,12F16.10)

             if(env_flag1) call deallocate_block(env_oper)

	     enddo

     		if(sys_flag1) call deallocate_block(sys_oper)
	 
enddo
                endif


        if(triangle==1)then
        do ik=1,1
        do ij=1,2
        write(10,*)'triangle=', triangle,'ij=', ij, 'kagome=',kagome

do i1=2*ny+1,sys_len-1,ny   !!!!!Lattice(1,sys_len)
        if(ik==1)then
        sid1=i1
        sid2= sid1+1
        chis1=1
         sid3=sid1+nleg
         if(mod(sid1,nleg)==0)then
         chis1=-chis1
        sid1=sid2-nleg
        sid2=i1
         sid3=sid2+nleg
                endif
                else
        chis1=-1
        sid1=i1
        sid2= sid1-1+nleg
         sid3=sid1+nleg

         if(mod(sid1,nleg)==1)then
         chis1=-chis1
        sid2=sid1+nleg
         sid3=sid2-1+nleg
                endif
                endif




                 sys_flag1 =.false.
                if(max(sid2,sid3)<=(sys_len-1)) then
                sys_flag1 =.true.
		        call Get_sys_triangular (sid1,sid2,sid3,sys_oper,trun_idx,wave)
                 endif
        
         do j1=1, env_len-1 
                if(ij==1)then
        eid1=j1
        chis2=-1
              eid2=eid1+nleg-1
              eid3=eid1+nleg
         if(mod(eid1,nleg)==1)then
                i=eid3
                eid3=eid2+nleg
                eid2=i
                chis2=-chis2
                endif
                endif

                if(ij==2)then
        chis2=1
              eid1=j1  
              eid2=eid1+1
              eid3=eid1+nleg
         if(mod(eid1,nleg)==0)then
                eid1=eid2-nleg
                eid2=j1
                chis2=-chis2
              eid3=j1+nleg
                endif
                endif



 
          env_flag1 =.false.
           if(max(eid2,eid3)<=(env_len-1)) then
            env_flag1 =.true.
		     call Get_env_triangular (eid1,eid2,eid3,env_oper,trun_idx,wave)
           endif

        	        if(sys_flag1) then
                    if(env_flag1) then
	     call sys_block_env_block_cor_dia(value,sys_oper,env_oper,sys_bs,env_bs,wave)
 	     write(101,112) "i=",sx,sy, "j=",ex,ey,"r=",ex-sx,-value*(-0.25d0)
     write(10,16) sid1, sid2, sid3,num_site-eid3+1, num_site-eid2+1,num_site-eid1+1,(num_site-eid3)/Ny-(sid1-1)/ny,value*(-6.0d0)*chis1*chis2,  chis1, chis2
 				   endif
                   endif
                   112 format(A2,1X,I2,I2,2X,A2,1X,I2,I2,2X,A2,1X,I2,2X,F16.10)
                   16 format(7i8,F16.10, 2i6)

             if(env_flag1) call deallocate_block(env_oper)

	     enddo

     		if(sys_flag1) call deallocate_block(sys_oper)
enddo
enddo
enddo
	 
                endif


call deallocate_basis(sys_bs)
call deallocate_basis(env_bs)


write(10,*)
close(10)

end subroutine Get_Chiral_cor


  


 subroutine Get_env_four (id1,id2,id3,env_tri,trun_idx,wave)

    use pubdata
	implicit none

   integer,intent(in) :: id1,id2,id3
   	type(Total_Block),intent(inout) :: env_tri
   integer,intent(in) :: trun_idx
   type(Wavefunction),intent(in) :: wave
   
   real(8) :: one
   integer :: x,y,i,j,sys_len,env_len
   	type(Total_Block) :: new_oper1,new_oper2,new_oper3,new_oper12,new_oper123 
    type(Total_Block) ::  new_oper,new_oper_temp


   one=cone
    sys_len=wave%sys_len
	env_len=wave%env_len
   !=================
     call Get_env_oper_ndia(st_sd,id1,new_oper1,id2,trun_idx,.false.,'B')
     call Get_env_oper_dia(st_sz,id2,new_oper2,id2,trun_idx,.false.)
        
     call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',one,new_oper,env_bsm(id2))
	 call deallocate_block(new_oper1)
	 call deallocate_block(new_oper2)
                 if(id2<=trun_idx) then
		              call block_transfer(new_oper,new_oper12)
	              else
 		              call update_trun_ndia(new_oper,new_oper12,envtruns(id2))        !!update_trun_ndia update_trun_dia
 	              endif
	            call deallocate_block(new_oper)

     call Change_env_oper_ndia(new_oper12,id2,new_oper_temp,id3,trun_idx,.false.,'B')
     call deallocate_block(new_oper12)

    call Get_env_oper_ndia(st_sd,id3,new_oper3,id3,trun_idx,.false.,'B') ! false ±£Ö¤ÁË²»ÏÈtrun
    call block_mul_block_ndia(new_oper_temp,'N',new_oper3,'T',one,new_oper,env_bsm(id3))
    call deallocate_block(new_oper_temp)
    call deallocate_block(new_oper3)

	call block_transfer(new_oper,env_tri)

	do i=1,env_tri%num
		env_tri%sub(i)%mat=env_tri%sub(i)%mat-transpose(new_oper%sub(i)%mat)
	end do
	call deallocate_block(new_oper)

   !=================
     call Get_env_oper_dia(st_sz,id1,new_oper1,id2,trun_idx,.false.)
     call Get_env_oper_ndia(st_sd,id2,new_oper2,id2,trun_idx,.false.,'B')
     call block_mul_block_ndia(new_oper1,'N',new_oper2,'N',one,new_oper,env_bsm(id2))
	 call deallocate_block(new_oper1)
	 call deallocate_block(new_oper2)
                 if(id2<=trun_idx) then
		              call block_transfer(new_oper,new_oper12)
	              else
 		              call update_trun_ndia(new_oper,new_oper12,envtruns(id2))        !!update_trun_ndia update_trun_dia
 	              endif
	            call deallocate_block(new_oper)

     call Change_env_oper_ndia(new_oper12,id2,new_oper_temp,id3,trun_idx,.false.,'B')
     call deallocate_block(new_oper12)

    call Get_env_oper_ndia(st_sd,id3,new_oper3,id3,trun_idx,.false.,'B') ! false ±£Ö¤ÁË²»ÏÈtrun
    call block_mul_block_ndia(new_oper_temp,'N',new_oper3,'T',one,new_oper,env_bsm(id3))
    call deallocate_block(new_oper_temp)
    call deallocate_block(new_oper3)


	do i=1,env_tri%num
		env_tri%sub(i)%mat=env_tri%sub(i)%mat-new_oper%sub(i)%mat+transpose(new_oper%sub(i)%mat)
	end do
	call deallocate_block(new_oper)

    !=================
     call Get_env_oper_ndia(st_sd,id1,new_oper1,id2,trun_idx,.false.,'B')
     call Get_env_oper_ndia(st_sd,id2,new_oper2,id2,trun_idx,.false.,'B')
     call block_mul_block_ndia(new_oper1,'N',new_oper2,'T',one,new_oper,env_bsm(id2))
	 call deallocate_block(new_oper1)
	 call deallocate_block(new_oper2)
                 if(id2<=trun_idx) then
		              call block_transfer(new_oper,new_oper12)
	              else
 		              call update_trun_dia(new_oper,new_oper12,envtruns(id2))        !!update_trun_ndia update_trun_dia
 	              endif
	            call deallocate_block(new_oper)

     call Change_env_oper_dia(new_oper12,id2,new_oper_temp,id3,trun_idx,.false.)
     call deallocate_block(new_oper12)

    call Get_env_oper_dia(st_sz,id3,new_oper3,id3,trun_idx,.false.) ! false ±£Ö¤ÁË²»ÏÈtrun
    call block_mul_block_dia(new_oper_temp,'N',new_oper3,'N',one,new_oper)
    call deallocate_block(new_oper_temp)
    call deallocate_block(new_oper3)

	do i=1,env_tri%num
		env_tri%sub(i)%mat=env_tri%sub(i)%mat-new_oper%sub(i)%mat+transpose(new_oper%sub(i)%mat)
	end do
	call deallocate_block(new_oper)


         if(id3<=trun_idx) then
		              call block_transfer (env_tri,new_oper123)
 	              else
 		              call update_trun_dia(env_tri,new_oper123,envtruns(id3))
 	              endif
	           call deallocate_block(env_tri) 

     call Change_env_oper_dia(new_oper123,id3,env_tri,env_len-1,trun_idx,.true.)  ! true ±£Ö¤ÁËtrun
     call deallocate_block(new_oper123)

end subroutine Get_env_four



!==========================================================
!Get reduced density matrix for system block
!==========================================================
subroutine Get_density_matrix_sys(sys_den,wave)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: wave
	type(Total_Block),intent(inout) :: sys_den

	integer :: i,sys_len,sys_dim,env_dim

	!Get general information for sys_den
	sys_den%len=wave%sys_len
	sys_den%num=wave%num
	sys_den%dim=0
	sys_den%up_dif=0
	sys_den%down_dif=0

	allocate(sys_den%sub(sys_den%num))
	do i=1,sys_den%num
		sys_den%sub(i)%num_up=wave%sub(i)%sys_num_up
		sys_den%sub(i)%num_down=wave%sub(i)%sys_num_down
		sys_den%sub(i)%down_dif=0

		sys_den%sub(i)%row_dim=wave%sub(i)%sys_dim
		sys_den%sub(i)%sdim=wave%sub(i)%sys_dim
		sys_den%dim=sys_den%dim+sys_den%sub(i)%sdim

		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim
		allocate(sys_den%sub(i)%mat(sys_dim,sys_dim))
		
                if(realcode)then
		call DGEMM('N','T',sys_dim,sys_dim,env_dim,1.0d0&
				&,wave%sub(i)%vec,sys_dim,wave%sub(i)%vec&
				&,sys_dim,0.0d0,sys_den%sub(i)%mat,sys_dim)
                        else


                 call ZGEMM('N','C',sys_dim,sys_dim,env_dim, cone&
                                &,wave%sub(i)%vec,sys_dim,wave%sub(i)%vec&
                                &,sys_dim,czero,sys_den%sub(i)%mat,sys_dim)


                        endif
	end do

end subroutine Get_density_matrix_sys


!==========================================================
!Get reduced density matrix for environment block
!==========================================================
subroutine Get_density_matrix_env(env_den,wave)
	use pubdata
	implicit none

	type(Wavefunction),intent(in) :: wave
	type(Total_Block),intent(inout) :: env_den

	integer :: i,env_len,sys_dim,env_dim

	!Get general information for env_den
	env_den%len=wave%env_len
	env_den%num=wave%num
	env_den%dim=0
	env_den%up_dif=0
	env_den%down_dif=0

	allocate(env_den%sub(env_den%num))
	do i=1,env_den%num
		env_den%sub(i)%num_up=wave%sub(i)%env_num_up
		env_den%sub(i)%num_down=wave%sub(i)%env_num_down
		env_den%sub(i)%down_dif=0
		env_den%sub(i)%row_dim=wave%sub(i)%env_dim
		env_den%sub(i)%sdim=wave%sub(i)%env_dim
		env_den%dim=env_den%dim+env_den%sub(i)%sdim

		sys_dim=wave%sub(i)%sys_dim
		env_dim=wave%sub(i)%env_dim
		allocate(env_den%sub(i)%mat(env_dim,env_dim))

                        if(realcode)then
		call DGEMM('T','N',env_dim,env_dim,sys_dim,1.0d0&
				&,wave%sub(i)%vec,sys_dim,wave%sub(i)%vec&
				&,sys_dim,0.0d0,env_den%sub(i)%mat,env_dim)
                        else
   call ZGEMM('C','N',env_dim,env_dim,sys_dim,cone&
                                &,wave%sub(i)%vec,sys_dim,wave%sub(i)%vec&
                                &,sys_dim,czero,env_den%sub(i)%mat,env_dim)



                        endif

	end do

end subroutine Get_density_matrix_env



!===============================================================
!Including ground and excited states for system block
!===============================================================
subroutine Get_density_matrix_sys_new(sys_den,wave,levels)
	use pubdata
	implicit none

	integer,intent(in) :: levels
	type(Wavefunction),intent(in) :: wave(levels)
	type(Total_Block),intent(inout) :: sys_den
	type(Total_Block) :: tmp_den

	integer :: i,x
	type(Wavefunction) :: eigvecs

	call wave_transfer(wave(1),eigvecs)
	do i=2,levels
	do x=1,eigvecs%num
		eigvecs%sub(x)%vec=eigvecs%sub(x)%vec+wave(i)%sub(x)%vec
	end do
	end do
	call wave_normalize(eigvecs)

                if(levels.eq.1)then
	call Get_density_matrix_sys(sys_den,eigvecs)
	call deallocate_wave(eigvecs)
                        else
	call deallocate_wave(eigvecs)
        call get_density_matrix_sys(sys_den,wave(1))
        call get_density_matrix_sys(tmp_den,wave(2))
        do i=1,sys_den%num
                sys_den%sub(i)%mat=(sys_den%sub(i)%mat+tmp_den%sub(i)%mat)/2.0d0
        end do
        call deallocate_block(tmp_den)
                endif


end subroutine Get_density_matrix_sys_new


!===============================================================
!Including ground and excited states for environment block
!===============================================================
subroutine Get_density_matrix_env_new(env_den,wave,levels)
	use pubdata
	implicit none

	integer,intent(in) :: levels
	type(Wavefunction),intent(in) :: wave(levels)
	type(Total_Block),intent(inout) :: env_den
	type(Total_Block) :: tmp_den

	integer :: i,x
	type(Wavefunction) :: eigvecs

	call wave_transfer(wave(1),eigvecs)
	do i=2,levels
	do x=1,eigvecs%num
		eigvecs%sub(x)%vec=eigvecs%sub(x)%vec+wave(i)%sub(x)%vec
	end do
	end do
	call wave_normalize(eigvecs)
        if(levels.eq.1)then
	call Get_density_matrix_env(env_den,eigvecs)
                else
	call deallocate_wave(eigvecs)
        call get_density_matrix_env(env_den,wave(1))
        call get_density_matrix_env(tmp_den,wave(2))
        do i=1,env_den%num
                env_den%sub(i)%mat=(env_den%sub(i)%mat+tmp_den%sub(i)%mat)/2.0d0
        end do
        call deallocate_block(tmp_den)
                endif

end subroutine Get_density_matrix_env_new


!=====================================================================
!Get the truncation operator, return the number of
!kept states, according to the discarded error
!=====================================================================
subroutine Get_truncation_operator(den_mat,trun_op,trun_err,num_state)
	use pubdata
	implicit none

	real(8),intent(in) :: trun_err
	type(Total_Block),intent(in) :: den_mat
	type(Total_Block),intent(inout) :: trun_op
	integer,intent(inout) :: num_state

	integer :: i,j,k,x,y,sub_dim,sys_len
	integer,allocatable :: idx(:),eigidx(:)
	real(8) :: pm1,tot_pm,tmp_err,  entropy2, entropy1,val1
           integer, allocatable ::  eigs(:)

	integer :: lda,lwork,info !,eigidx(keep)
	double complex,allocatable :: work(:)
        double precision, allocatable :: rwork(:),  eigvalue(:)
	integer,dimension(den_mat%num) :: tmpidx,tmp_up,tmp_down,row,col
	type(Density) :: tmp

	!<1>: Get tmp%num,%dim,tmp%sub(:)%dim,%qn,%eigvec
	tmp%len=den_mat%len
	tmp%dim=den_mat%dim
	tmp%num=den_mat%num
	tmp%up_dif=den_mat%up_dif
	tmp%down_dif=den_mat%down_dif
	allocate(tmp%sub(tmp%num))
	do i=1,tmp%num
		tmp%sub(i)%num_up=den_mat%sub(i)%num_up
		tmp%sub(i)%num_down=den_mat%sub(i)%num_down

		tmp%sub(i)%row_dim=den_mat%sub(i)%row_dim
		tmp%sub(i)%sdim=den_mat%sub(i)%sdim

		sub_dim=tmp%sub(i)%sdim
		allocate(tmp%sub(i)%idx(sub_dim),tmp%sub(i)%eigval(sub_dim))
		allocate(tmp%sub(i)%eigvec(sub_dim,sub_dim))

		tmp%sub(i)%idx=0
		tmp%sub(i)%eigval=0.0d0
		tmp%sub(i)%eigvec=den_mat%sub(i)%mat
	end do

	!<2>: Diagonalize den_mat to get tmp%sub(:)%eigval,~eigvec,~idx
	tmp%dim=0
	do k=1,tmp%num
		sub_dim=tmp%sub(k)%sdim
		lda=sub_dim
		lwork=3*sub_dim
		allocate(rwork(lwork))

		!Diagonalize tmp%sub(k)%eigvec in ascending order

                        if(realcode)then
		call DSYEV('V','U',sub_dim,tmp%sub(k)%eigvec,lda,tmp%sub(k)%eigval,rwork,lwork,info)
		deallocate(rwork)
                                else

                                allocate(work(lwork))
                        call zheev('V','U',sub_dim,tmp%sub(k)%eigvec,lda,tmp%sub(k)%eigval,work,3*lda,rwork,info)
                deallocate(work,rwork)

                                endif


   do i=1,tmp%sub(k)%sdim
                        tmp%sub(k)%idx(i)=tmp%dim+i
		end do
		tmp%dim=tmp%dim+tmp%sub(k)%sdim
	end do

	allocate(idx(tmp%dim),eigvalue(tmp%dim))
	allocate(eigs(tmp%dim))
	idx=0
	eigvalue=0.0d0
	tmp%dim=0
	do k=1,tmp%num
		idx(tmp%dim+1:tmp%dim+tmp%sub(k)%sdim)=tmp%sub(k)%idx(1:tmp%sub(k)%sdim)
		eigvalue(tmp%dim+1:tmp%dim+tmp%sub(k)%sdim)=tmp%sub(k)%eigval(1:tmp%sub(k)%sdim)
            eigs(tmp%dim+1:tmp%dim+tmp%sub(k)%sdim)=tmp%sub(k)%num_down+1
		tmp%dim=tmp%dim+tmp%sub(k)%sdim
	end do

	!<3>: Sort the eigvalue in descending order to choose the kept states

	!Renormalize density-matrxi spectrum
	tot_pm=0.0d0
	do i=1,tmp%dim
		tot_pm=tot_pm+eigvalue(i)
	end do

         Entropy1=0.0d0
        do i=1, tmp%dim
                !Get Entropy
                if(eigvalue(i)>=1.0E-15) then
                val1=eigvalue(i)/eigs(i)
                        entropy1=entropy1-val1*dlog(val1)*eigs(i)
                endif
                enddo


	call QuickSort_Double(eigvalue,idx,1,tmp%dim,tmp%dim)

!Renormalize density-matrxi spectrum
        tot_pm=0.0d0
        do i=1,tmp%dim
                tot_pm=tot_pm+eigvalue(i)
        end do
        do i=1,tmp%dim
                eigvalue(i)=eigvalue(i)/tot_pm
        end do


	!Get Num_state and tmp_err
	tmp_err=0.0d0
	if(tmp%dim<=kept_min) then
	    Num_state=tmp%dim
	    tmp_err=0.0d0
	else
	    tmp_err=0.0d0
	    do i=kept_min+1,tmp%dim
	        tmp_err=tmp_err+eigvalue(i)
	    end do
	    
	    Num_state=kept_min
	    !do while(Num_state<kept_max)
	    do while(tmp_err>trun_err)
	        Num_state=Num_state+1
	        tmp_err=tmp_err-eigvalue(Num_state)
	        if((tmp_err<=trun_err).or.(Num_state>=kept_max)) then
	            goto 1111
	        endif
	    end do
	    1111 continue
	endif
	allocate(eigidx(num_state))
	eigidx(1:num_state)=idx(1:num_state)


	!Get truncation error and Entanglement entropy
	Entropys=0.0d0
	tot_pm=0.0d0
	do i=1,tmp%dim
		tot_pm=tot_pm+eigvalue(i)

		!Get Entropy
		if(eigvalue(i)>=1.0E-15) then
			Entropys=Entropys-eigvalue(i)*dlog(eigvalue(i))
		endif
	end do

	!For Entanglement Entropy
	if(system_block) then
		sys_len=den_mat%len
	else
		sys_len=Num_site-den_mat%len
	endif

open(60,file="entropy1.dat",position='append')
        write(60,22)sys_len,entropy1, Num_site,tot_num_up,tot_num_down
        close(60)
22      format(i8, f20.10, 3i9)
	!Get Entanglement spectrum
	open(60,file="Kept_State.dat",position='append')
	write(60,112) "Sys_len=",sys_len,"Num_State=",Num_State,"Err=",tmp_err !trun_err
	close(60)
112 format(A8,I3,2X,A10,I5,2X,A4,E18.8)

	!<4>: Get trun%num,%qn,%row,%col
	tmpidx=0
	do k=1,num_state
		do i=1,tmp%num
			do x=1,tmp%sub(i)%sdim
				if(tmp%sub(i)%idx(x)==eigidx(k)) then
					tmpidx(i)=tmpidx(i)+1
					goto 601
				endif
			end do
		end do
		601 continue
	end do

	trun_op%num=0
	tmp_up=0
	tmp_down=0
	col=0
	row=0
	do i=1,tmp%num
		if(tmpidx(i)>0) then
			trun_op%num=trun_op%num+1
			tmp_up(trun_op%num)=tmp%sub(i)%num_up
			tmp_down(trun_op%num)=tmp%sub(i)%num_down
			col(trun_op%num)=tmpidx(i)
			row(trun_op%num)=tmp%sub(i)%sdim
		endif
	end do

	!<5>: Get trun%sub(:)%mat
	allocate(trun_op%sub(trun_op%num))
	do i=1,trun_op%num
		trun_op%sub(i)%num_up=tmp_up(i)
		trun_op%sub(i)%num_down=tmp_down(i)
                 trun_op%sub(i)%down_dif=0

		trun_op%sub(i)%row_dim=row(i)
		trun_op%sub(i)%sdim=col(i)
		allocate(trun_op%sub(i)%mat(row(i),col(i)))
		trun_op%sub(i)%mat=0.0d0
		do x=1,tmp%num
			if(tmp%sub(x)%num_up==trun_op%sub(i)%num_up) then
			if(tmp%sub(x)%num_down==trun_op%sub(i)%num_down) then
				trun_op%sub(i)%mat(1:row(i),1:col(i))&
					&=tmp%sub(x)%eigvec(1:row(i),row(i)-col(i)+1:row(i))
				goto 602
			endif
			endif
		end do
		602 continue
	end do

	!check the relation between trun%dim and keep
	trun_op%len=den_mat%len
	trun_op%dim=0
	trun_op%up_dif=den_mat%up_dif
	trun_op%down_dif=den_mat%down_dif
	do i=1,trun_op%num
		trun_op%dim=trun_op%dim+trun_op%sub(i)%sdim
	end do

	do i=1,tmp%num
		deallocate(tmp%sub(i)%idx,tmp%sub(i)%eigval,tmp%sub(i)%eigvec)
	end do
	deallocate(tmp%sub,idx,eigvalue)

end subroutine Get_truncation_operator

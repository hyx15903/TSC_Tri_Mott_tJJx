

subroutine get_triangle()! get different small triangles for chiral interactions
        use pubdata
	integer kk,  js, j00, i0, i, j, i0a(4)
	integer jkk(10), jw(10), ji, ji1, i1,i2,i3,i4,ii, ns21, ns1
	integer k1,k2,k3,k4,k0,  tri(8*num_site,4)

        tri=0
        tri3=0
        if(kagome.eq.1)then
        kk1=2
        do i=1,num_site/3  !! all up and down triangles
        ii=2*i-1
        tri(ii,1)=3*(i-1)+1
        tri(ii,2)=3*(i-1)+2
        tri(ii,3)=3*(i-1)+3
        ii=ii+1
        tri(ii,1)=tri(ii-1,3)
        tri(ii,2)=tri(ii,1)+1
        if(mod(tri(ii,1),nleg).eq.0)tri(ii,2)=tri(ii,2)-nleg
        tri(ii,3)=tri(ii,2)+1-nleg
        if(open.eq.0.and.tri(ii,3).le.0)tri(ii,3)=tri(ii,3)+num_site !! used PBC 
        if(tri(ii,3).le.0)tri(ii,3)=0  !! non PBC 
        enddo
        endif


        if(triangle.eq.1)then!!ringexchange
        kk1=2
        tri=0
        do i=1,num_site  !! all up and down triangles
        ii=i
        tri(ii,1)=i
        tri(ii,2)=i+nleg
        tri(ii,4)=-1
        if(i+nleg.gt.num_site)then
        if(open.eq.0)tri(ii,2)=i+nleg-num_site
        if(open.ne.0)tri(ii,2)=0  !! non PBC 
        endif
        tri(ii,3)=i+1
        if(mod(i,nleg).eq.0)then
        tri(ii,3)=i+1-nleg
        tri(ii,4)=1
        endif


        ii=ii+num_site
        tri(ii,1)=tri(ii-num_site,2)
        tri(ii,2)=tri(ii-num_site,3)
        tri(ii,3)=tri(ii,2)+nleg
        tri(ii,4)=1
        if(mod(i,nleg).eq.0)then
        tri(ii,4)=-1
        endif
        if(tri(ii,3).gt.num_site)then
        if(open.eq.0)tri(ii,3)=tri(ii,3)-num_site
        if(open.ne.0)tri(ii,3)=0  !! non PBC 
        endif
        enddo
        endif


        if(square.eq.1)then
        kk1=4
        tri=0
        do i=1, num_site  !! all up and down triangles
        ii=i
        tri(ii,1)=i
        tri(ii,2)=i+1
        tri(ii,3)=i+nleg+1
        tri(ii,4)=-1
        if(mod(i,nleg).eq.0)then
        tri(ii,2)=tri(ii,2)-nleg
        tri(ii,3)=tri(ii,3)-nleg
        tri(ii,4)=1
        endif
        if(tri(ii,3).gt.num_site)then
        if(open.eq.0)tri(ii,3)=tri(ii,3)-num_site
        if(open.ne.0)tri(ii,3)=0  !! non PBC 
        endif


        ii=i+num_site
        tri(ii,1)=i
        tri(ii,3)=i+nleg
        tri(ii,2)=i+nleg+1
        tri(ii,4)=1
        if(mod(i,nleg).eq.0)then
        tri(ii,2)=tri(ii,2)-nleg
        tri(ii,4)=-1
        endif
        if(tri(ii,2).gt.num_site)then
        if(open.eq.0)tri(ii,2)=tri(ii,2)-num_site
        if(open.eq.0)tri(ii,3)=tri(ii,3)-num_site
        if(open.ne.0)tri(ii,2)=0  !! non PBC 
        if(open.ne.0)tri(ii,3)=0  !! non PBC 
        endif

        if(kk1.eq.4)then

        ii=i+2*num_site
        tri(ii,1)=i
        tri(ii,3)=i+nleg
        tri(ii,2)=i+nleg-1
        tri(ii,4)=1
        if(mod(i,nleg).eq.1)then
        tri(ii,2)=tri(ii,2)+nleg
        tri(ii,4)=-1
        endif
        if(tri(ii,2).gt.num_site)then
        if(open.eq.0)tri(ii,2)=tri(ii,2)-num_site
        if(open.ne.0)tri(ii,2)=0  !! non PBC 
        if(open.eq.0)tri(ii,3)=tri(ii,3)-num_site
        if(open.ne.0)tri(ii,3)=0  !! non PBC 
        endif


        ii=i+num_site*3
        tri(ii,1)=i
        tri(ii,3)=i-1!!!!nleg
        tri(ii,2)=i+nleg-1
        tri(ii,4)=-1
        if(mod(i,nleg).eq.1)then
        tri(ii,3)=tri(ii,3)+nleg
        tri(ii,2)=tri(ii,2)+nleg
        tri(ii,4)=1
        endif
        if(tri(ii,2).gt.num_site)then
        if(open.eq.0)tri(ii,2)=tri(ii,2)-num_site
        if(open.ne.0)tri(ii,2)=0  !! non PBC 
        endif

        endif

        enddo

        endif




        do ij=1, kk1
        do ii=1,num_site

        i=ii+(ij-1)*num_site
        i1=tri(i,1)
        i2=tri(i,2)
        i3=tri(i,3)

        if(i1*i2*i3.ne.0)then
        tri3(i1,i2,ij)=i3
        tri3(i2,i1,ij)=i3
        tri3(i1,i3,ij)=i2
        tri3(i3,i1,ij)=i2
        tri3(i2,i3,ij)=i1
        tri3(i3,i2,ij)=i1

        tri3(i1,i2,ij+kk1)=tri(i, 4)
        tri3(i2,i1,ij+kk1)=tri(i,4)
        tri3(i1,i3,ij+kk1)=tri(i,4)
        tri3(i3,i1,ij+kk1)=tri(i,4)
        tri3(i2,i3,ij+kk1)=tri(i,4)
        tri3(i3,i2,ij+kk1)=tri(i,4)

        tri3(i1,i2,ij+kk1*3)=i
        tri3(i2,i1,ij+kk1*3)=i
        tri3(i1,i3,ij+kk1*3)=i
        tri3(i3,i1,ij+kk1*3)=i
        tri3(i2,i3,ij+kk1*3)=i
        tri3(i3,i2,ij+kk1*3)=i
        
        endif

        enddo
        enddo

         kk=0
        num_spin2=0
        do i=1, num_site-1
        do j1=i+1, num_site !from two sites, we get the index of pair operators

        do ij=1, kk1
        kk=kk+1
        num_spin2(i,j1,ij)=kk
        num_spin2(j1,i,ij)=kk

        enddo
        enddo
        enddo

        ind2e=0
                 !!!!counting operators we need for a given system size
                !! there are always less than num_site**2/4 operators
        nleg2=0
        tri3b=0
        tri3s=1.0d0
        do ns1=2,num_site
        k0=0
        do j1=1, ns1-1
        do k=1,neibt1  !! NN 
        j2=neib(j1,k)

        if(j2.gt.j1.and.j2.le.ns1)then


        do ij=1,kk1

        ji=num_spin2(j1, j2,ij) !!! ji th bond
        i=tri3(j1,j2, ij)         !! i forms the triangle with j1, j2

        if(i.gt.ns1.and.i.le.num_site)then !!!neibors, unclosed inside sys
        k0=k0+1
        ind2e(ji, ns1, 1)=k0
        tri3b(k0,ns1,1)=i

        if(kagome.eq.1)then
        if(mod(i, nleg).eq.0)then
        if(mod(j2,nleg).le.3)tri3s(k0, ns1,1)=-1
        endif
        if(mod(j2, nleg).eq.0)then
        if(mod(i,nleg).le.3)tri3s(k0, ns1,1)=-1
        endif
        endif


        if(triangle.eq.1)then
        tri3s(k0, ns1, 1)= tri3(j1, j2, ij+kk1)
        
        nrr=0
        if(nrr==0)go to 13
        if(tri3(j1,j2,ij+3*kk1).ne.nrr)then
        tri3s(k0, ns1, 1)= tri3(j1, j2, ij+kk1)*0
        tri3(j1, j2, ij+kk1)= tri3(j1, j2, ij+kk1)*0
        tri3(j2, j1, ij+kk1)= tri3(j1, j2, ij+kk1)*0
        tri3(j1, i, ij+kk1)= tri3(j1, j2, ij+kk1)*0
        tri3(i, j1, ij+kk1)= tri3(j1, j2, ij+kk1)*0
        tri3(i, j2, ij+kk1)= tri3(j1, j2, ij+kk1)*0
        tri3(j2, i, ij+kk1)= tri3(j1, j2, ij+kk1)*0
        endif
13      continue
        
        tri3s(k0, ns1, 3)= i
        endif
        endif
        enddo
        endif
        enddo
        enddo
        nleg2(ns1, 1)=k0  !!! for each system size ns1, there are nleg2 (number) of C_iCj pairs
        enddo
!!!! above is for bond operators
        tri3s(1:num_site, 1:num_site,2)=-1
        tri3s(1:num_site, 1:num_site,2)=1

        do ns1=2,num_site
        k0=0
        do j1=1, ns1-1       !! env site
        j11=num_site-j1+1   !! count from sys
        do k=1,neibt1
        j22=neib(j11,k) 
        j2=num_site-j22+1     !! in env
        if(j2.gt.j1.and.j2.le.ns1)then
        do ij=1,kk1
        ji=num_spin2(j1, j2, ij) !!! ji th ring
        i=tri3(j11,j22,ij)  !! sys ocount
        ii=num_site-i+1   !! env
        if(ii.gt.ns1.and.ii.le.num_site)then !!!neibors, unclosed
        k0=k0+1
        ind2e(ji, ns1, 2)=k0
        tri3b(k0,ns1,2)=ii

        if(kagome.eq.1)then
        if(mod(ii,nleg).eq.nleg-1)then
        tri3s(k0, ns1,2)=-tri3s(k0, ns1, 2)
        endif

        if(mod(ii,nleg).eq.0.and.i.eq.1)then
        if(mod(j22, nleg).eq.0)then
        tri3s(k0, ns1,2)=-tri3s(k0, ns1, 2)
        endif
        endif
        endif


        if(triangle.eq.1)then
        tri3s(k0, ns1, 2)= tri3(j11, j22, ij+kk1)
        if(tri3(j11,j22,ij+3*kk1).gt.2)then
        endif
        tri3s(k0, ns1, 4)= i+num_site
        endif
        endif

        enddo
        endif
        enddo
        enddo
        nleg2(ns1, 2)=k0  !!! for each system size ns1, there nleg2 of C_iCj pairs
        enddo

end subroutine get_triangle


subroutine dmrg_setup0()
	use pubdata
        integer ns, i, syssite, kmax, i1, k0
	nleg11=0
	do ns=1,num_site-2
	kmax=0
	k0=0
	do i1=1,ns
	INB1(i1,ns)=10000
	k00=0

        do j=1, neibt
        if(neib(i1, j).ge.ns)k00=k00+1
        enddo

        if(k00.ne.0)then
	k0=k0+1
	INB1(i1,ns)=k0
	endif

		enddo

		nleg11(ns)=k0
		kmax=max(kmax,k0)
	enddo
	nleg1s=nleg11
        inb1s=inb1
end subroutine dmrg_setup0


subroutine dmrg_setup(ns, ns21)  !! sys and env lengths
	use pubdata
        implicit none
        integer  ns,  ns21
        integer i, j1, j2,nn, nt, totsite, ii, ix, iy, j, jx, jy

	nn=ns+ns21+2
	nt=num_site
	nxc=(nn-1)/nleg+1
        if(nxc==1)nxc=2
	nposi=0
        nposi_lat=0
        nposi_dm=0
        neib=0
        neiba=0
	totsite=nleg*nxc   !!! last point of system during "infinite process".
        do i=1,ns+1
        nposi(i)=1
        nposi_lat(i)=i
        nposi_dm(i)=i
        enddo
        nposi(ns+1)=2
        
        do i=1,ns21+1
        nposi(ns+1+i)=3
        nposi_dm(totsite+1-i)=i+(ns+1)
        nposi_lat(i+ns+1)=totsite+1-i
        enddo
        nposi(ns+1+ns21+1)=4

!! the following part is lattice dependent
        neib=0
        neiba=0

        if(square.eq.1)then
	do ix=1,nxc !!! nxc is the given system total length for any warm_up or sweep
	do iy=1,nleg
	ii=iy+(ix-1)*nleg !! 

        if(nleg.eq.1)then
	jy=ii+1  !! square lattice,  y neighbor 
        if(ix.eq.nxc)then
        if(open.ne.0)jy=0
        if(open.eq.0)jy=jy-nxc
        endif

        if(jy.ne.0)then
	neib(ii,1)=jy
	neib(jy,2)=ii
        neiba(ii,jy)=1
        neiba(jy,ii)=2
        endif

	jy=ii+3  !! square lattice,  y neighbor 
        if(jy.gt.nxc)then
        if(open.ne.0)jy=0
        if(open.eq.0)jy=jy-nxc
        endif

        if(jy.ne.0)then
        if(ii.eq.num_site/2-1)then
	neib(ii,3)=jy
	neib(jy,4)=ii
        neiba(ii,jy)=3
        neiba(jy,ii)=4
        endif
        endif
        endif

!!!! 2D
	ii=iy+(ix-1)*nleg !! 
        jy=ii-iy+mod(iy,nleg)+1  !! square lattice,  y neighbor 
        if(iy.ne.nleg.or.nleg.ne.2)then
        neib(ii,1)=jy
        neib(jy,2)=ii
        neiba(ii,jy)=1
        neiba(jy,ii)=2
                endif

        jx=ii+nleg*nly
        if(jx.gt.nleg*nxc)then
        if(open.eq.0.and.nxc.ge.3)then
        jx=jx-nleg*nxc
        else
        jx=0
        endif
        endif
        if(jx.ne.0)then
        neib(ii,3)=jx
        neib(jx,4)=ii
        neiba(ii,jx)=3
        neiba(jx,ii)=4

        jy= jx-1
        if(iy.eq.1)jy=jy+nleg        !! going to same layer
        if(jy.ne.0)then
        neib(ii,5)=jy
        neib(jy,6)=ii
        neiba(ii,jy)=5
        neiba(jy,ii)=6
        endif
        jy= jx+1
        if(iy.eq.nleg)jy=jy-nleg        !! going to same layer
        if(jy.ne.0)then
        neib(ii,7)=jy
        neib(jy,8)=ii
        neiba(ii,jy)=7
        neiba(jy,ii)=8
        endif


        endif

        enddo
        enddo


        endif

        if(triangle.eq.1)then
	do ix=1,nxc !!! nxc is the given system total length for any warm_up or sweep
	do iy=1,nleg
!!!! 2D
	ii=iy+(ix-1)*nleg !! 
        jy=ii-iy+mod(iy,nleg)+1  !!  y neighbor 
        if(openy.ne.0.and.iy.eq.nleg)jy=0
        if(jy.ne.0)then
        if(iy.ne.nleg.or.nleg.ne.2)then
        neib(ii,1)=jy
        neib(jy,2)=ii
        neiba(ii,jy)=1
        neiba(jy,ii)=2
                endif
                endif

        jx=ii+nleg*nly         !! going to same layer
        if(jx.gt.nleg*nxc)then
        if(open.eq.0.and.nxc.ge.3)then
        jx=jx-nleg*nxc
        else
        jx=0
        endif
        endif

        if(jx.ne.0)then
        neib(ii,3)=jx
        neib(jx,4)=ii
        neiba(ii,jx)=3
        neiba(jx,ii)=4


        jy= jx-1
        if(iy.eq.1)then
        if(openy==0)then
        jy=jy+nleg        !! going to same layer
        else
        jy=0
        endif
        endif

        if(jy.ne.0)then
        neib(ii,5)=jy
        neib(jy,6)=ii
        neiba(ii,jy)=5
        neiba(jy,ii)=6
        endif
        endif



        enddo
        enddo





	do ix=1,nxc !!! nxc is the given system total length for any warm_up or sweep
	do iy=1,nleg
	ii=iy+(ix-1)*nleg !! 

                j1=neib(ii,1)
                j2=neib(j1,3)
                        if(j1*j2.ne.0)then
                neib(ii,7)=j2
                neib(j2, 8)=ii
        neiba(ii,j2)=7
        neiba(j2,ii)=8
                        endif
                
                if(nleg.ne.3)then
                j1=neib(ii,3)
                j2=neib(j1,5)
                        if(j1*j2.ne.0)then
                neib(ii,9)=j2
                neib(j2, 10)=ii
        neiba(ii,j2)=9
        neiba(j2,ii)=10
                        endif
                        endif
                
                j1=neib(ii,2)
                j2=neib(j1,5)
                        if(j1*j2.ne.0)then
                neib(ii,11)=j2
                neib(j2, 12)=ii
        neiba(ii,j2)=11
        neiba(j2,ii)=12
                        endif



 if(nleg==3)go to 13
                j1=neib(ii,1)
                j2=neib(j1,1)
                        if(j1*j2.ne.0)then
                if(neiba(ii,j2).eq.0)then
                neib(ii,13)=j2
                neib(j2,14)=ii
        neiba(ii,j2)=13
        neiba(j2,ii)=14
                        endif
                        endif
13      continue

                j1=neib(ii,3)
                j2=neib(j1,3)
                        if(j1*j2.ne.0)then
                if(neiba(ii,j2).eq.0)then
                neib(ii,15)=j2
                neib(j2, 16)=ii
        neiba(ii,j2)=15
        neiba(j2,ii)=16
                        endif
                        endif

                j1=neib(ii,5)
                j2=neib(j1,5)
                        if(j1*j2.ne.0)then
                if(neiba(ii,j2).eq.0)then
                neib(ii,17)=j2
                neib(j2, 18)=ii
        neiba(ii,j2)=17
        neiba(j2,ii)=18
                        endif
                        endif

                enddo
                enddo
        endif
                
!!!!!!!!!!!!!


        if(honeycomb.eq.1)then
	do ix=1,nxc !!! nxc is the given system total length for any warm_up or sweep
	do iy=1,nleg
!!!! 2D
	ii=iy+(ix-1)*nleg !! 
        jy=ii-iy+mod(iy,nleg)+1  !!  y neighbor 
        if(openy.ne.0.and.iy.eq.nleg)jy=0
        if(jy.ne.0)then
        if(iy.ne.nleg.or.nleg.ne.2)then
        neib(ii,1)=jy
        neib(jy,2)=ii
        neiba(ii,jy)=1
        neiba(jy,ii)=2
                endif
                endif

        if(mod(iy,2)==0)then   !!! B sublattice
        jx=ii+nleg*nly-1         !! going to same layer
        if(jx.gt.nleg*nxc)then
        if(open.eq.0.and.nxc.ge.3)then
        jx=jx-nleg*nxc
        else
        jx=0
        endif
        endif
        if(jx.ne.0)then
        neib(ii,3)=jx
        neib(jx,4)=ii
        neiba(ii,jx)=3
        neiba(jx,ii)=4
        endif
        endif

        if(nly==2.and.mod(ix,2)==1)then   !!! B sublattice
        jx=ii+nleg         !! going to same layer
        if(jx.ne.0)then
        neib(ii,5)=jx
        neib(jx,6)=ii
        neiba(ii,jx)=5
        neiba(jx,ii)=6
        endif
        endif
        enddo
        enddo
        endif

	do i=1,nt
	do j=neibt1+1,neibt
	if(abs(jz(j,1))+abs(jd(j,1))+abs(jt(j,1))+abs(jn(j,1)).eq.0.0)then
        neiba(i, neib(i,j))=0
        neiba(neib(i,j), i)=0
	neib(i,j)=0
	endif
	enddo
	enddo

end subroutine dmrg_setup

! new part for the triangle
!To get the position dxy of site(x,y) 
integer function get_dxy(dx,dy)
	use pubdata
	implicit none
	integer,intent(in) :: dx,dy
	get_dxy=(dy-1)*nx+dx
	return
end function get_dxy

!Get filename for all model data
subroutine Get_model_name(sysmod,envmod)
	use pubdata
	implicit none

	character(len=20) :: sysmod(Num_site,3)
	character(len=20) :: envmod(Num_site,3)

	integer :: i,x,y
	character(len=20) :: systr="systruns",envtr="envtruns"
	character(len=20) :: sysbl="sysblocks",envbl="envblocks"
	character(len=20) :: sysbs="sysbs",envbs="envbs"

	!<1>: Get system model name
	do i=1,Num_site
		sysmod(i,1)=sysbl
		call str_link_str(sysmod(i,1),9,i)

		sysmod(i,2)=systr
		call str_link_str(sysmod(i,2),8,i)

		sysmod(i,3)=sysbs
		call str_link_str(sysmod(i,3),5,i)
	end do


	!<1>: Get environment model name
	do i=1,Num_site
		envmod(i,1)=envbl
		call str_link_str(envmod(i,1),9,i)

		envmod(i,2)=envtr
		call str_link_str(envmod(i,2),8,i)

		envmod(i,3)=envbs
		call str_link_str(envmod(i,3),5,i)
	end do

end subroutine Get_model_name


!Get system and environment's msn name
subroutine Get_wave_name
	use pubdata
	implicit none

	integer :: i
	character(len=6) :: wavename="eigvec"
	do i=1,8
		wave_name(i)=wavename
		call str_link_str(wave_name(i),6,i)
	end do

end subroutine Get_wave_name


!Get Lattice configuration
subroutine Get_Lattice
	use pubdata
	implicit none
	integer :: x,y,idx
	Lattice=0
	Latticever=0
	do x=1,Nx
	do y=1,Ny
		idx=(x-1)*Ny+y
		Lattice(1,idx)=x
		Lattice(2,idx)=y
		Latticever(x,y)=idx
	end do
	end do

end subroutine Get_Lattice


!============================================
!Initiate the site for the t-J model
subroutine Get_site_operator_su2
	use pubdata
	implicit none
	integer :: i,num_ups,num_downs,idx

	!For t-J model
	!<1>: For st_elec_up: c^+_{up spin}

        ! for tjmodel,  s=0, 1state for 0-electron or 1-electron on site

        st_elec_up%len=1
        if(tjmodel.ge.1)then
        st_elec_up%num=1
        st_elec_up%dim=1
                else
        if(hubmodel==1)then
        st_elec_up%num=2
        st_elec_up%dim=2
                endif
                endif

	st_elec_up%up_dif=up_bias !! C_i up^+
	st_elec_up%down_dif=down_bias                 
	allocate(st_elec_up%sub(st_elec_up%num))

	!Matrix element <up,0|c^+_{up}|0,0>
	st_elec_up%sub(1)%num_up=0    !! acting on empty site
	st_elec_up%sub(1)%num_down=0
	st_elec_up%sub(1)%down_dif=1
	st_elec_up%sub(1)%row_dim=1
	st_elec_up%sub(1)%sdim=1
	allocate(st_elec_up%sub(1)%mat(1,1))
	st_elec_up%sub(1)%mat(1,1)=-dsqrt(2.0d0)  


                if(hubmodel==1)then
	st_elec_up%sub(2)%num_up=1   !! acting on s=1
	st_elec_up%sub(2)%num_down=1
	st_elec_up%sub(2)%down_dif=-1
	st_elec_up%sub(2)%row_dim=1
	st_elec_up%sub(2)%sdim=1
	allocate(st_elec_up%sub(2)%mat(1,1))
	st_elec_up%sub(2)%mat(1,1)=dsqrt(2.0d0)
                endif


        call block_transfer(st_elec_up, st_elec_down)
	st_elec_down%up_dif=-up_bias !!  tensor for annihilation
	!Matrix element <up,0|c^+_{up}|0,0>
	st_elec_down%sub(1)%num_up=1    !! acting on single  site
	st_elec_down%sub(1)%num_down=1
	st_elec_down%sub(1)%down_dif=-1

        if(hubmodel==1)then
	st_elec_down%sub(2)%num_up= 2   !! acting on double site
	st_elec_down%sub(2)%num_down=0
	st_elec_down%sub(2)%down_dif=1
	st_elec_down%sub(2)%mat(1,1)=-dsqrt(2.0d0)  
                endif
        

	!<3>: For num_elec_up: n_{up spin}
	num_elec%len=1
	num_elec%num=N_max
	num_elec%dim=N_max
	num_elec%up_dif=0
	num_elec%down_dif=0
	allocate(num_elec%sub(num_elec%num))

	idx=1
		num_elec%sub(idx)%num_up=1
		num_elec%sub(idx)%num_down=ss
		num_elec%sub(idx)%down_dif=0
		num_elec%sub(idx)%row_dim=1
		num_elec%sub(idx)%sdim=1
		allocate(num_elec%sub(idx)%mat(1,1))
		num_elec%sub(idx)%mat(1,1)=dsqrt(1.0d0+ss)

        if(num_elec%num==2)then
	idx=2
		num_elec%sub(idx)%num_up=0
		num_elec%sub(idx)%num_down=0
		num_elec%sub(idx)%down_dif=0   ! down_dif
		num_elec%sub(idx)%row_dim=1
		num_elec%sub(idx)%sdim=1
		allocate(num_elec%sub(idx)%mat(1,1))
		num_elec%sub(idx)%mat(1,1)=0.0d0   !! hole site
                endif

        if(hubmodel==1)then
        idx=3
		num_elec%sub(idx)%num_up=2
		num_elec%sub(idx)%num_down=0
		num_elec%sub(idx)%down_dif=0   ! down_dif
		num_elec%sub(idx)%row_dim=1
		num_elec%sub(idx)%sdim=1
		allocate(num_elec%sub(idx)%mat(1,1))
		num_elec%sub(idx)%mat(1,1)=2.0d0   !! dsqrt(2J+1)=1,  2 electrons
                        endif

        call block_transfer(num_elec, st_double)
	idx=1
		st_double%sub(idx)%mat(1,1)=0.0d0
        st_double%sub(1)%mat(1,1)=0.0d0
        if(hubmodel==1)then
	idx=3
		st_double%sub(idx)%mat(1,1)=1.0d0
                endif



	!<6>: For st_sd: S^+_{up spin}
	st_sd%len=1
	st_sd%num=1
	st_sd%dim=1
	st_sd%up_dif=0
	st_sd%down_dif=su   !! spin change =1,  but unit=su
	allocate(st_sd%sub(st_sd%num))

	!Matrix element <up+1,down-1|S^+|up,down>
	st_sd%sub(1)%num_up=1
	st_sd%sub(1)%num_down=ss
	st_sd%sub(1)%down_dif=0
	st_sd%sub(1)%row_dim=1
	st_sd%sub(1)%sdim=1
	allocate(st_sd%sub(1)%mat(1,1))
	st_sd%sub(1)%mat(1,1)=dsqrt(1.50d0)*ss

    !<7>: For st_sz: Sz_{up spin}
        st_sz%len=1
        st_sz%num=N_max
        st_sz%dim=N_max
        st_sz%up_dif=0
        st_sz%down_dif=0
        allocate(st_sz%sub(st_sz%num))
        idx=1
                st_sz%sub(idx)%num_up=1
                st_sz%sub(idx)%num_down=ss
                st_sz%sub(idx)%down_dif=0
                st_sz%sub(idx)%row_dim=1
                st_sz%sub(idx)%sdim=1
                allocate(st_sz%sub(idx)%mat(1,1))
                st_sz%sub(idx)%mat(1,1)=1.0d0

        if(n_max==2)then
        idx=2
                st_sz%sub(idx)%num_up=0
                st_sz%sub(idx)%num_down=0
                st_sz%sub(idx)%down_dif=0
                st_sz%sub(idx)%row_dim=1
                st_sz%sub(idx)%sdim=1
               allocate(st_sz%sub(idx)%mat(1,1))
                st_sz%sub(idx)%mat(1,1)=0.0d0
                endif

        if(hubmodel==1)then
        idx=3
                st_sz%sub(idx)%num_up=2
                st_sz%sub(idx)%num_down=0
                st_sz%sub(idx)%down_dif=0
                st_sz%sub(idx)%row_dim=1
                st_sz%sub(idx)%sdim=1
                allocate(st_sz%sub(idx)%mat(1,1))
                st_sz%sub(idx)%mat(1,1)=0.0d0
                 endif



	!<7>: For st_si: %%Sz_{up spin}
        call block_transfer(num_elec, st_si)

                do i=1, n_max
		st_si%sub(i)%mat(1,1)=dsqrt(1.0d0+st_si%sub(i)%down_dif) !
                        enddo
21      continue
end subroutine Get_site_operator_su2

!==================================================
!Initiate the single site hamiltonian
subroutine Get_hamiltonian_site(site_ham)
	use pubdata
	implicit none

	type(Total_Model) :: site_ham
	integer :: i,idx,num_ups,num_downs

	!================================================
	!<1>: Get single site Hamiltonian

        site_ham%len=1
		call block_transfer(st_sz,site_ham%ham)
        do idx=1,st_sz%num
		site_ham%ham%sub(idx)%mat(1,1)=0.0d0
                enddo

	!<2>: For operator in the middle
                i=1
		call block_transfer(st_elec_up,site_ham%sub(i)%elec_up)
		call block_transfer(st_elec_down,site_ham%sub(i)%elec_down)
		call block_transfer(st_sd,site_ham%sub(i)%spin_sd)
		call block_transfer(num_elec,site_ham%sub(i)%num_sn)
                if(hubmodel==1)then
		call block_transfer(st_double,site_ham%sub(i)%double)
                endif
12      continue
end subroutine Get_hamiltonian_site


!Get single site basis information
subroutine Get_single_site_basis
	use pubdata
	implicit none

	integer :: i,idx,num_ups,num_downs


     st_basis%len=1
        st_basis%num=N_max
        st_basis%dim=N_max
        allocate(st_basis%sub(st_basis%num))

        do idx=1, st_sz%num
                st_basis%sub(idx)%idx=0
                st_basis%sub(idx)%num=1
                st_basis%sub(idx)%dim=1
                st_basis%sub(idx)%new_num_up=st_sz%sub(idx)%num_up
                st_basis%sub(idx)%new_num_down=st_sz%sub(idx)%num_down
                allocate(st_basis%sub(idx)%sub(1))
                st_basis%sub(idx)%sub(1)%spos=0
                st_basis%sub(idx)%sub(1)%sdim=1
                st_basis%sub(idx)%sub(1)%bl_num_up=st_basis%sub(idx)%new_num_up
                st_basis%sub(idx)%sub(1)%bl_num_down=st_basis%sub(idx)%new_num_down
                st_basis%sub(idx)%sub(1)%st_num_up=0
                st_basis%sub(idx)%sub(1)%st_num_down=0
        end do

end subroutine Get_single_site_basis

!===========================================================
!Get ground_state using two-site algorithm
!===========================================================
subroutine warm_up_point(bg_id,end_id,trun_err,levels)
	use pubdata
	implicit none

	real(8),intent(in) :: trun_err
	integer,intent(in) :: bg_id,end_id,levels

	integer :: i,j,x,y,label,num_idx,num_state,tot_up,tot_down
	type(Super_Basis) :: super
	type(Total_Model) :: sys,new_sys,env,new_env
	type(Total_Block) :: sys_den,sys_trun,env_den,env_trun
	type(Total_Basis) :: sys_bs,env_bs,env_basis
	type(Wavefunction) :: eigvec(levels)


	!<1>: Enlarge system without truncation
	num_idx=1 !Using mat_wave for whole system
	label=1
	system_block=.true.

              call dmrg_setup(1,1) 

	call Get_hamiltonian_site(sys)


	call Get_hamiltonian_site(env)

        
	call model_to_disk(sys,1001,label,.true.)
	call model_to_disk(env,1001,label,.false.)


	!<2>: Enlarge system without truncation
	do while((label+1)<=bg_id)
		!Enlargement for system part
                call dmrg_setup(sys%ham%len,env%ham%len) 
                mshift=0
		call Get_basis(sys%ham,st_sz,sys_bs)
		call Get_hamiltonian_sys(sys,new_sys,sys_bs)


		call deallocate_model(sys)
		call model_transfer(new_sys,sys)
		call deallocate_model(new_sys)

		!Enlargement for environment part
                mshift=sys%len!!!+1
		call Get_basis(env%ham,st_sz,env_bs)
		call Get_hamiltonian_env(env,new_env,env_bs)
		call deallocate_model(env)
		call model_transfer(new_env,env)
		call deallocate_model(new_env)

		!Save to disk
		call model_to_disk(sys,1001,label+1,.true.)
		call model_to_disk(env,1001,label+1,.false.)

		call basis_to_disk(sys_bs,1001,label+1,.true.)
		call basis_to_disk(env_bs,1001,label+1,.false.)

		call deallocate_basis(sys_bs)
		call deallocate_basis(env_bs)

		label=label+1
	end do


	!<3>: Enlarge system with truncation
	do while((label+1)<=end_id)
		call Get_basis(sys%ham,st_sz,sys_bs)
		call Get_basis(env%ham,st_sz,env_bs)
                call dmrg_setup(sys%ham%len,env%ham%len) 

		!tot_pt=sys_bs%len*N_max
		tot_up=int(1.0*(sys_bs%len+env_bs%len)/num_site*tot_num_up+0.001)
		tot_down=int((tot_num_down*1.0d0/Num_site)*(sys_bs%len+env_bs%len)+0.5)
                tot_up=tot_up-mod(tot_up,2)
                !!!!tot_up=sys_bs%len+env_bs%len
                tot_down=0

		if(tot_up>tot_num_up) then
			tot_up=tot_num_up
		endif
		if(tot_down>tot_num_down) then
			tot_down=tot_num_down
		endif


		call Get_super_basis(sys_bs,env_bs,super,tot_up,tot_down)


		!=======================================================================
		!For both ground state and excited states
		do i=1,levels
			call allocate_wave(super,eigvec(i))
			call wave_initialize(eigvec(i))
		end do
		call lanczos_memo_new(sys,env,sys_bs,env_bs,super,eigvec,levels,num_idx)


		call Get_density_matrix_sys_new(sys_den,eigvec,levels)
		call Get_density_matrix_env_new(env_den,eigvec,levels)
		!=======================================================================
		!call Get_sys_den(sys_den,eigvec)

                             if(sys_bs%dim.le.kept_max)then
        call block_transfer_basis(sys_bs, sys_trun)
         do i=1,sys_trun%num
                      sys_trun%sub(i)%mat=0.d0
                        do j=1,sys_trun%sub(i)%sdim
                                sys_trun%sub(i)%mat(j,j)=1.d0
                                enddo
                                   end do

                        else
		call Get_truncation_operator(sys_den,sys_trun,trun_err,num_state)		
                        endif
		call deallocate_block(sys_den)

		!Enlargement for system part
                mshift=0
		call Get_hamiltonian_sys(sys,new_sys,sys_bs, sys_trun)
		!call deallocate_model(sys)
		!call Get_hamiltonian_trun(new_sys,sys,sys_trun)
		call deallocate_model(new_sys)

                             if(env_bs%dim.le.kept_max)then
        call block_transfer_basis(env_bs, env_trun)
         do i=1,env_trun%num
                      env_trun%sub(i)%mat=0.d0
                        do j=1,env_trun%sub(i)%sdim
                                env_trun%sub(i)%mat(j,j)=1.d0
                                enddo
                                   end do
                        else
		call Get_truncation_operator(env_den,env_trun,trun_err,num_state)		
                        endif
		call deallocate_block(env_den)

		!Enlargement for system part
                mshift=sys%len!!+1
		call Get_hamiltonian_env(env,new_env,env_bs, env_trun)
		call deallocate_model(new_env)

		!Save to disk
		call model_to_disk(sys,1001,label+1,.true.)
		call model_to_disk(env,1001,label+1,.false.)

		call truns_to_disk(sys_trun,1001,label+1,.true.)
		call truns_to_disk(env_trun,1001,label+1,.false.)

		call basis_to_disk(sys_bs,1001,label+1,.true.)
		call basis_to_disk(env_bs,1001,label+1,.false.)

		call deallocate_block(sys_trun)
		call deallocate_block(env_trun)

		call deallocate_basis(sys_bs)
		call deallocate_basis(env_bs)
		call deallocate_super_basis(super)

		!For ground state only
		do i=1,levels
			call wave_to_disk(eigvec(i),1001,wave_name(i),7)
			call deallocate_wave(eigvec(i))
		end do

		label=label+1
	end do
	call deallocate_model(sys)
	call deallocate_model(env)

       
end subroutine warm_up_point


!Sweep from left edge to right edge
subroutine left_to_right(bg_id,end_id,trun_err,levels,newflag)
	use pubdata
	implicit none

	real(8),intent(in) :: trun_err
	logical,intent(in) :: newflag
	integer,intent(in) :: bg_id,end_id,levels

	integer :: i,j,x,y,label,num_idx,num_state,tot_up,tot_down
	type(Super_Basis) :: super
	type(Total_Model) :: sys,new_sys,env
	type(Total_Block) :: sys_den,sys_trun,env_trun
	type(Total_Basis) :: sys_bs,env_bs,env_basis
	type(Wavefunction) :: eigvec(levels),tmpvec

	!Start sweep from left to right
	num_idx=2
	label=bg_id
	system_block=.true.
	call model_from_disk(sys,1001,label,.true.)
	call model_from_disk(env,1001,Num_site-label-2,.false.)

        call dmrg_setup(sys%ham%len,env%ham%len) 
        mshift=0
	call Get_basis(sys%ham,st_sz,sys_bs)
	call Get_basis(env%ham,st_sz,env_bs)
	call Get_super_basis(sys_bs,env_bs,super,tot_num_up,tot_num_down)

	!=================================================================================
	!Get initial wavefunction
	do i=1,levels
         if(restart1.eq.0)then
                        call allocate_wave(super,eigvec(i))
                        call wave_initialize(eigvec(i))

                else !Read wavefunction from disk
                        call wave_from_disk(eigvec(i),1001,wave_name(i),7)
                call allocate_wave(super,tmpvec)

                if(eigvec(i)%sys_len.ne.tmpvec%sys_len)then
        call basis_from_disk(env_basis,1001,num_site-label,.false.)
        call truns_from_disk(sys_trun,1001,label,.true.)
    call truns_from_disk(env_trun,1001,num_site-label-1,.false.)
        call wave_trans_LR(sys_bs,env_basis,sys_trun,env_trun,eigvec(i),tmpvec)
                call deallocate_wave(eigvec(i))
                call wave_transfer(tmpvec,eigvec(i))
                call deallocate_block(sys_trun)
                call deallocate_block(env_trun)
                call deallocate_basis(env_basis)
                endif

                call deallocate_wave(tmpvec)

                endif
                enddo

	call lanczos_memo_new(sys,env,sys_bs,env_bs,super,eigvec,levels,num_idx)
	call Get_density_matrix_sys_new(sys_den,eigvec,levels)



	do i=1,levels
		call wave_to_disk(eigvec(i),1001,wave_name(i),7)
	end do
	!=================================================================================

	!call Get_sys_den(sys_den,eigvec)
	call Get_truncation_operator(sys_den,sys_trun,trun_err,num_state)
	call deallocate_block(sys_den)
	call deallocate_model(env)

                mshift=0
	call Get_hamiltonian_sys(sys,new_sys,sys_bs,sys_trun)
	call deallocate_model(new_sys)
	call model_to_disk(sys,1001,label+1,.true.)
	call truns_to_disk(sys_trun,1001,label+1,.true.)
	call basis_to_disk(sys_bs,1001,label+1,.true.)

	call deallocate_basis(sys_bs)
	call deallocate_basis(env_bs)
	call deallocate_super_basis(super)

	label=label+1
	do while((label+1)<=end_id)
		!Get e igen_state
        
		call Get_basis(sys%ham,st_sz,sys_bs)
		call basis_from_disk(env_bs,1001,Num_site-label-1,.false.)
		call basis_from_disk(env_basis,1001,Num_site-label,.false.)
		call model_from_disk(env,1001,Num_site-label-2,.false.)
                call dmrg_setup(sys%ham%len,env%ham%len) 
		call Get_super_basis(sys_bs,env_bs,super,tot_num_up,tot_num_down)

		call truns_from_disk(env_trun,1001,Num_site-label-1,.false.)
		call allocate_wave(super,tmpvec)

		!Wavefunction transformation
		!call wave_initialize(super,tmpvec)
		do i=1,levels
			call wave_trans_LR(sys_bs,env_basis,sys_trun,env_trun,eigvec(i),tmpvec)
			call deallocate_wave(eigvec(i))
			call wave_transfer(tmpvec,eigvec(i))
                        ! call wave_initialize(eigvec(i)) !!!!remove
		end do
		call deallocate_block(sys_trun)
		call deallocate_block(env_trun)
		call deallocate_wave(tmpvec)

		!==========================================================================
		!For both ground state and excited states
		do i=2,levels
			call wave_initialize(eigvec(i))
		end do
	call lanczos_memo_new(sys,env,sys_bs,env_bs,super,eigvec,levels,num_idx)


		call Get_density_matrix_sys_new(sys_den,eigvec,levels)
		call deallocate_model(env)

		do i=1,levels
			call wave_to_disk(eigvec(i),1001,wave_name(i),7)
		end do
		!==========================================================================

		call Get_truncation_operator(sys_den,sys_trun,trun_err,num_state)
		call deallocate_block(sys_den)

                mshift=0
		call Get_hamiltonian_sys(sys,new_sys,sys_bs,sys_trun)
		call deallocate_model(new_sys)

		call model_to_disk(sys,1001,label+1,.true.)
		call truns_to_disk(sys_trun,1001,label+1,.true.)
		call basis_to_disk(sys_bs,1001,label+1,.true.)

		call deallocate_basis(sys_bs)
		call deallocate_basis(env_bs)
		call deallocate_basis(env_basis)
		call deallocate_super_basis(super)

		label=label+1
	end do
	
	call deallocate_model(sys)
	call deallocate_block(sys_trun)

	!For both ground state and excited states
	do i=1,levels
		call deallocate_wave(eigvec(i))
	end do

end subroutine left_to_right


!Sweep from right to left edge with wavefunction transformation
subroutine right_to_left(bg_id,end_id,trun_err,levels,newflag)
	use pubdata
	implicit none

	real(8),intent(in) :: trun_err
	logical,intent(in) :: newflag
	integer,intent(in) :: bg_id,end_id,levels

	integer :: i,j,label,num_idx,num_state,tot_up,tot_down
	type(Super_Basis) :: super
	type(Total_Model) :: new_env,sys,env
	type(Total_Block) :: env_den,sys_trun,env_trun
	type(Total_Basis) :: sys_bs,env_bs,sys_basis
	type(Wavefunction) :: eigvec(levels),tmpvec


	!Enlarge environment until Num_site-2
	num_idx=2
	label=bg_id
	system_block=.false.
	call model_from_disk(sys,1001,label-1,.true.)
	call model_from_disk(env,1001,Num_site-label-1,.false.)

                call dmrg_setup(sys%ham%len,env%ham%len) 
	call Get_basis(sys%ham,st_sz,sys_bs)
	call Get_basis(env%ham,st_sz,env_bs)
	call Get_super_basis(sys_bs,env_bs,super,tot_num_up,tot_num_down)

	!=======================================================================
	!Get initial wavefunction
	do i=1,levels
  if(restart1.eq.0)then
                        call allocate_wave(super,eigvec(i))
                        call wave_initialize(eigvec(i))

                else !Read wavefunction from disk
                        call wave_from_disk(eigvec(i),1001,wave_name(i),7)
                call allocate_wave(super,tmpvec)

                if(eigvec(i)%sys_len.ne.tmpvec%sys_len)then
        write(*,*)'using older wavefunction and wave transformation'
        call basis_from_disk(sys_basis,1001,label+1,.true.)
        call truns_from_disk(sys_trun,1001,label,.true.)
    call truns_from_disk(env_trun,1001,num_site-label-1,.false.)
        call wave_trans_RL(sys_basis,env_bs,sys_trun,env_trun,eigvec(i),tmpvec)

        
                call deallocate_wave(eigvec(i))
                call wave_transfer(tmpvec,eigvec(i))
                call deallocate_block(sys_trun)
                call deallocate_block(env_trun)
                call deallocate_basis(sys_basis)
                endif
                call deallocate_wave(tmpvec)

                endif

	end do

	!For both ground state and excited states
	call lanczos_memo_new(sys,env,sys_bs,env_bs,super,eigvec,levels,num_idx)
	call Get_density_matrix_env_new(env_den,eigvec,levels)

	do i=1,levels
		call wave_to_disk(eigvec(i),1001,wave_name(i),7)
	end do
	!=======================================================================

	!call Get_env_den(env_den,eigvec)
	call Get_truncation_operator(env_den,env_trun,trun_err,num_state)
	call deallocate_block(env_den)
	call deallocate_model(sys)

                mshift=sys%len+1
	call Get_hamiltonian_env(env,new_env,env_bs,env_trun)
	!call deallocate_model(env)
	!call Get_hamiltonian_trun(new_env,env,env_trun)
	call deallocate_model(new_env)

	call model_to_disk(env,1001,Num_site-label,.false.)
	call truns_to_disk(env_trun,1001,Num_site-label,.false.)
	call basis_to_disk(env_bs,1001,Num_site-label,.false.)

	call deallocate_basis(sys_bs)
	call deallocate_basis(env_bs)
	call deallocate_super_basis(super)

	label=label-1
	do while((label-1)>=end_id)
		!Get eigen_state
		call basis_from_disk(sys_bs,1001,label,.true.)
		call basis_from_disk(sys_basis,1001,label+1,.true.)
		call Get_basis(env%ham,st_sz,env_bs)
		call Get_super_basis(sys_bs,env_bs,super,tot_num_up,tot_num_down)
		call model_from_disk(sys,1001,label-1,.true.)

                call dmrg_setup(sys%ham%len,env%ham%len) 

		call truns_from_disk(sys_trun,1001,label,.true.)
		call allocate_wave(super,tmpvec)

		!Wavefunction transformation
		!call wave_initialize(super,tmpvec)
		do i=1,levels
			call wave_trans_RL(sys_basis,env_bs,sys_trun,env_trun,eigvec(i),tmpvec)
			call deallocate_wave(eigvec(i))
			call wave_transfer(tmpvec,eigvec(i))
                      !  call wave_initialize(eigvec(i)) !!!!remove
		end do
		call deallocate_block(sys_trun)
		call deallocate_block(env_trun)
		call deallocate_wave(tmpvec)

		!==========================================================================
		!For both ground state and excited states
		do i=2,levels
			call wave_initialize(eigvec(i))
		end do

		call lanczos_memo_new(sys,env,sys_bs,env_bs,super,eigvec,levels,num_idx)
		call Get_density_matrix_env_new(env_den,eigvec,levels)
		call deallocate_model(sys)

		do i=1,levels
			call wave_to_disk(eigvec(i),1001,wave_name(i),7)
		end do
		!==========================================================================

		!Get env_den and truncate system
		!call Get_env_den(env_den,eigvec)
		call Get_truncation_operator(env_den,env_trun,trun_err,num_state)
		call deallocate_block(env_den)

                mshift=sys%len+1
		call Get_hamiltonian_env(env,new_env,env_bs,env_trun)
		!call deallocate_model(env)
		!call Get_hamiltonian_trun(new_env,env,env_trun)
		call deallocate_model(new_env)

		call model_to_disk(env,1001,Num_site-label,.false.)
		call truns_to_disk(env_trun,1001,Num_site-label,.false.)
		call basis_to_disk(env_bs,1001,Num_site-label,.false.)

		call deallocate_basis(sys_bs)
		call deallocate_basis(env_bs)
		call deallocate_basis(sys_basis)
		call deallocate_super_basis(super)

		label=label-1
	end do

	call deallocate_model(env)
	call deallocate_block(env_trun)

	!For both ground state and excited states
	do i=1,levels
		call deallocate_wave(eigvec(i))
	end do

end subroutine right_to_left


!==========================================================================
!Wavefunction Transformation for left_to_right sweep
!==========================================================================
subroutine wave_trans_LR(sys_bs,env_bs,sys_trun,env_trun,old_wave,new_wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_trun,env_trun
	type(Wavefunction),intent(in) :: old_wave
	type(Wavefunction),intent(inout) :: new_wave

	integer :: i,j,k,x,y,z
	logical :: sys_trun_flag,env_trun_flag,new_wave_flag,sys_bs_flag
	integer :: old_pos,old_sub_dim,new_pos,new_sub_dim
	integer :: old_sys_dim,old_env_dim,new_sys_dim,new_env_dim
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down,st_num_up,st_num_down
	integer :: sys_trun_id,env_trun_id,new_wave_id,sys_bs_id,sys_bs_sub_id

        integer:: st_down_dif, st_num_down1, j1,j2,j12,j3,jtot, j4,j123, j34
        double complex coef1
        real*8  coef11
	double complex,allocatable :: sub(:,:)
        real(8), external :: w6js


	do i=1,old_wave%num
		sys_num_up=old_wave%sub(i)%sys_num_up
		sys_num_down=old_wave%sub(i)%sys_num_down
		env_num_up=old_wave%sub(i)%env_num_up
		env_num_down=old_wave%sub(i)%env_num_down
		old_sys_dim=old_wave%sub(i)%sys_dim
		old_env_dim=old_wave%sub(i)%env_dim

		sys_trun_flag=.false.
		do k=1,sys_trun%num
			if(sys_trun%sub(k)%num_up==sys_num_up) then
			if(sys_trun%sub(k)%num_down==sys_num_down) then
				sys_trun_flag=.true.
				sys_trun_id=k
				goto 101
			endif
			endif
		end do
		101 continue

		if(sys_trun_flag) then
			do k=1,env_bs%num
				if(env_bs%sub(k)%new_num_up==env_num_up) then
				if(env_bs%sub(k)%new_num_down==env_num_down) then

					do x=1,env_bs%sub(k)%num
						st_num_up=env_bs%sub(k)%sub(x)%st_num_up
						st_num_down=env_bs%sub(k)%sub(x)%st_num_down 
						new_env_num_up=env_bs%sub(k)%sub(x)%bl_num_up
						new_env_num_down=env_bs%sub(k)%sub(x)%bl_num_down

                       !!  st_down_dif=down_bias  !!! a site spin quantum number
                      !!  do st_num_down1=-st_down_dif, st_down_dif,  su

						old_pos=env_bs%sub(k)%sub(x)%spos
						old_sub_dim=env_bs%sub(k)%sub(x)%sdim

						env_trun_flag=.false.
						do y=1,env_trun%num
							if(env_trun%sub(y)%num_up==new_env_num_up) then
							if(env_trun%sub(y)%num_down==new_env_num_down) then
								env_trun_flag=.true.
								env_trun_id=y
								goto 102
							endif
							endif
						end do
						102 continue
		do 	new_sys_num_down=abs(sys_num_down-st_num_down), sys_num_down+st_num_down, su
				new_sys_num_up=sys_num_up+st_num_up

						sys_bs_flag=.false.
						do y=1,sys_bs%num
							if(sys_bs%sub(y)%new_num_up==new_sys_num_up) then
							if(sys_bs%sub(y)%new_num_down==new_sys_num_down) then
							do z=1,sys_bs%sub(y)%num
			if((sys_bs%sub(y)%sub(z)%bl_num_up==sys_num_up).and.(sys_bs%sub(y)%sub(z)%bl_num_down==sys_num_down)) then
			if((sys_bs%sub(y)%sub(z)%st_num_up==st_num_up).and.(sys_bs%sub(y)%sub(z)%st_num_down==st_num_down)) then
										sys_bs_flag=.true.
										sys_bs_id=y
										sys_bs_sub_id=z
										goto 103
									endif
									endif
								end do
							endif
							endif
						end do
						103 continue

						new_wave_flag=.false.
				do y=1,new_wave%num
		if((new_wave%sub(y)%sys_num_up==new_sys_num_up).and.(new_wave%sub(y)%sys_num_down==new_sys_num_down)) then
		if((new_wave%sub(y)%env_num_up==new_env_num_up).and.(new_wave%sub(y)%env_num_down==new_env_num_down)) then
								new_wave_flag=.true.
								new_wave_id=y
								goto 104
							endif
							endif
						end do
						104 continue

						if(env_trun_flag.and.sys_bs_flag.and.new_wave_flag) then
							new_sys_dim=new_wave%sub(new_wave_id)%sys_dim
							new_env_dim=new_wave%sub(new_wave_id)%env_dim
							new_pos=sys_bs%sub(sys_bs_id)%sub(sys_bs_sub_id)%spos
							new_sub_dim=sys_bs%sub(sys_bs_id)%sub(sys_bs_sub_id)%sdim

        j12=sys_num_down
        j3=st_num_down
        jtot=0
        j4=new_env_num_down
        j123=new_sys_num_down
        j34=env_num_down
        
                        !! Racah W (j12,j3,jtot,j4,j123,j34)  to W6j( j12, j3,j123, j4, jtot, j34)
        !!coef11=w6js(j12,j3,j123, j4, jtot, j34)
        coef11=(-1)**((j12+j3+j123)/2)/dsqrt((1.0d0+j12)*(1.0d0+j123))
                if(coef11.ne.0.0)then
        coef11=coef11*dsqrt((1.0d0+j123)*(j34+1.0d0))
        if(mod((j12+j3+j4)/2,  2).ne.0)coef11=-coef11
        if(mod((j3+j4-j34)/2,  2).ne.0)coef11=-coef11
        coef1=coef11

					allocate(sub(new_sub_dim,old_sub_dim))

              if(realcode)then
       call DGEMM('T','N',new_sub_dim,old_sub_dim,old_sys_dim,coef11,sys_trun%sub(sys_trun_id)%mat,old_sys_dim,&
                           &old_wave%sub(i)%vec(1:old_sys_dim,old_pos+1:old_pos+old_sub_dim),old_sys_dim,0.0d0,sub,new_sub_dim)
        call DGEMM('N','T',new_sub_dim,new_env_dim,old_sub_dim,1.0d0,sub,new_sub_dim,env_trun%sub(env_trun_id)%mat,&
                                       &new_env_dim,1.0d0,new_wave%sub(new_wave_id)%vec(new_pos+1:new_pos+new_sub_dim,1:new_env_dim),new_sub_dim)

                                else

          call ZGEMM('C','N',new_sub_dim,old_sub_dim,old_sys_dim,coef1,sys_trun%sub(sys_trun_id)%mat,old_sys_dim,&
                        &old_wave%sub(i)%vec(1:old_sys_dim,old_pos+1:old_pos+old_sub_dim),old_sys_dim,czero,sub,new_sub_dim)
                    call ZGEMM('N','C',new_sub_dim,new_env_dim,old_sub_dim,cone,sub,new_sub_dim,env_trun%sub(env_trun_id)%mat,&
                               &new_env_dim,cone,new_wave%sub(new_wave_id)%vec(new_pos+1:new_pos+new_sub_dim,1:new_env_dim),new_sub_dim)


                                endif
        


						deallocate(sub)
						endif
						endif
        enddo
                        
					end do
				endif
				endif
			end do
		endif
	end do

	call wave_normalize(new_wave)

end subroutine wave_trans_LR


!Wavefunction Transformation for right_to_left sweep
subroutine wave_trans_RL(sys_bs,env_bs,sys_trun,env_trun,old_wave,new_wave)
	use pubdata
	implicit none

	type(Total_Basis),intent(in) :: sys_bs,env_bs
	type(Total_Block),intent(in) :: sys_trun,env_trun
	type(Wavefunction),intent(in) :: old_wave
	type(Wavefunction),intent(inout) :: new_wave

	integer :: i,j,k,x,y,z
	logical :: sys_trun_flag,env_trun_flag,new_wave_flag,env_bs_flag
	integer :: old_pos,old_sub_dim,new_pos,new_sub_dim
	integer :: old_sys_dim,old_env_dim,new_sys_dim,new_env_dim
	integer :: sys_num_up,sys_num_down,new_sys_num_up,new_sys_num_down
	integer :: env_num_up,env_num_down,new_env_num_up,new_env_num_down,st_num_up,st_num_down
	integer :: sys_trun_id,env_trun_id,new_wave_id,env_bs_id,env_bs_sub_id
	double complex,allocatable :: sub(:,:)
        integer:: st_down_dif, st_num_down1,j1,j2,j12,j3,jtot, j4,j123, j34
        real*8  coef11
        double complex coef1
        real(8), external :: w6js


	do i=1,old_wave%num
		sys_num_up=old_wave%sub(i)%sys_num_up
		sys_num_down=old_wave%sub(i)%sys_num_down
		env_num_up=old_wave%sub(i)%env_num_up
		env_num_down=old_wave%sub(i)%env_num_down
		old_sys_dim=old_wave%sub(i)%sys_dim
		old_env_dim=old_wave%sub(i)%env_dim

		env_trun_flag=.false.
		do k=1,env_trun%num
			if(env_trun%sub(k)%num_up==env_num_up) then
			if(env_trun%sub(k)%num_down==env_num_down) then
				env_trun_flag=.true.
				env_trun_id=k
				goto 101
			endif
			endif
		end do
		101 continue

		if(env_trun_flag) then
			do k=1,sys_bs%num
				if(sys_bs%sub(k)%new_num_up==sys_num_up) then
				if(sys_bs%sub(k)%new_num_down==sys_num_down) then
					do x=1,sys_bs%sub(k)%num
						st_num_down=sys_bs%sub(k)%sub(x)%st_num_down
						st_num_up=sys_bs%sub(k)%sub(x)%st_num_up
						new_sys_num_up=sys_bs%sub(k)%sub(x)%bl_num_up
						new_sys_num_down=sys_bs%sub(k)%sub(x)%bl_num_down

                        do new_env_num_down=abs(env_num_down-st_num_down), env_num_down+st_num_down, su
						new_env_num_up=env_num_up+st_num_up

						old_pos=sys_bs%sub(k)%sub(x)%spos
						old_sub_dim=sys_bs%sub(k)%sub(x)%sdim

						sys_trun_flag=.false.
						do y=1,sys_trun%num
							if(sys_trun%sub(y)%num_up==new_sys_num_up) then
							if(sys_trun%sub(y)%num_down==new_sys_num_down) then
								sys_trun_flag=.true.
								sys_trun_id=y
								goto 102
							endif
							endif
						end do
						102 continue

						env_bs_flag=.false.
						do y=1,env_bs%num
							if(env_bs%sub(y)%new_num_up==new_env_num_up) then
							if(env_bs%sub(y)%new_num_down==new_env_num_down) then
								do z=1,env_bs%sub(y)%num
									if((env_bs%sub(y)%sub(z)%bl_num_up==env_num_up).and.(env_bs%sub(y)%sub(z)%bl_num_down==env_num_down)) then
									if((env_bs%sub(y)%sub(z)%st_num_up==st_num_up).and.(env_bs%sub(y)%sub(z)%st_num_down==st_num_down)) then
										env_bs_flag=.true.
										env_bs_id=y
										env_bs_sub_id=z
										goto 103
									endif
									endif
								end do
							endif
							endif
						end do
						103 continue

						new_wave_flag=.false.
						do y=1,new_wave%num
							if((new_wave%sub(y)%sys_num_up==new_sys_num_up).and.(new_wave%sub(y)%sys_num_down==new_sys_num_down)) then
							if((new_wave%sub(y)%env_num_up==new_env_num_up).and.(new_wave%sub(y)%env_num_down==new_env_num_down)) then
								new_wave_flag=.true.
								new_wave_id=y
								goto 104
							endif
							endif
						end do
						104 continue

						if(sys_trun_flag.and.env_bs_flag.and.new_wave_flag) then
							new_sys_dim=new_wave%sub(new_wave_id)%sys_dim
							new_env_dim=new_wave%sub(new_wave_id)%env_dim

							new_pos=env_bs%sub(env_bs_id)%sub(env_bs_sub_id)%spos
							new_sub_dim=env_bs%sub(env_bs_id)%sub(env_bs_sub_id)%sdim
        j12=env_num_down
        j3=st_num_down
        jtot=0
        j4=new_sys_num_down
        j123=new_env_num_down
        j34=sys_num_down
        
        coef11=(-1)**((j34+j4+j3)/2)/dsqrt((1.0d0+j12)*(1.0d0+j4))
                if(coef11.ne.0.0)then
        coef11=coef11*dsqrt((1.0d0+j123)*(j34+1.0d0))!!!*(-1)**((j12+j3+j4)/2)
        if(mod((j4+j3+j12)/2,  2).ne.0)coef11=-coef11
        if(mod((j12+j3-j123)/2,  2).ne.0)coef11=-coef11

        coef1=coef11

			
				allocate(sub(new_sys_dim,old_env_dim))

            if(realcode)then
        call DGEMM('N','N',new_sys_dim,old_env_dim,old_sub_dim,coef11,sys_trun%sub(sys_trun_id)%mat,new_sys_dim,&
                &old_wave%sub(i)%vec(old_pos+1:old_pos+old_sub_dim,1:old_env_dim),old_sub_dim,0.0d0,sub,new_sys_dim)
                 call DGEMM('N','N',new_sys_dim,new_sub_dim,old_env_dim,1.0d0,sub,new_sys_dim,env_trun%sub(env_trun_id)%mat,old_env_dim,&
                                                                          &1.0d0,new_wave%sub(new_wave_id)%vec(1:new_sys_dim,new_pos+1:new_pos+new_sub_dim),new_sys_dim)
                                else
         call ZGEMM('N','N',new_sys_dim,old_env_dim,old_sub_dim,coef1,sys_trun%sub(sys_trun_id)%mat,new_sys_dim,&
                         &old_wave%sub(i)%vec(old_pos+1:old_pos+old_sub_dim,1:old_env_dim),old_sub_dim,czero,sub,new_sys_dim)
                call ZGEMM('N','N',new_sys_dim,new_sub_dim,old_env_dim,cone,sub,new_sys_dim,env_trun%sub(env_trun_id)%mat,old_env_dim,&
                                  &cone,new_wave%sub(new_wave_id)%vec(1:new_sys_dim,new_pos+1:new_pos+new_sub_dim),new_sys_dim)
                                        endif


							deallocate(sub)
					endif
                                endif
                                enddo
					end do
				endif
				endif
			end do
		endif
	end do

	call wave_normalize(new_wave)

end subroutine wave_trans_RL

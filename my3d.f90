module tabling
    ! common
    real(kind=8) :: pi = 3.141592653589793, delta = 10.0**(-5)
    complex(kind=8) :: c_one = (1.0,0.0), c_i = (0.0,1.0)
    complex(kind=8), dimension(:), allocatable :: xlst, ylst
    character(1), parameter :: creturn=achar(13), tab=achar(9)          ! ascii character 
    contains
    ! build interpolation table of x
    subroutine buildtable1(nx, xlst)
        integer :: nx, ny, i
        complex(kind=8) :: xlst(1+2*nx*nx)
        do i = 0, 2*nx*nx
            xlst(i+1) = exp(2.0*pi*c_i*i/nx)
        end do
    end subroutine buildtable1
    ! build interpolation table of x,y
    subroutine buildtable2(nx, ny, xlst, ylst)
        integer :: nx, ny, i
        complex(kind=8) :: xlst(1+nx*nx), ylst(1+ny*ny)
        do i = 0, nx*nx
            xlst(i+1) = exp(2.0*pi*c_i*i/nx)
        end do
        do i = 0, ny*ny
            ylst(i+1) = exp(2.0*pi*c_i*i/ny)
        end do
    end subroutine buildtable2
end module tabling
    

! main proram
#include <fintrf.h>
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
	use tabling
	use omp_lib
	!use mkl95_blas
	implicit none
	integer*8 :: nlhs, nrhs
	integer*8 ::  plhs(*), prhs(*)
    ! common variables
    integer :: n, i, j, u, v, w, s, ni, nj, nx, ny, layi, utmp, vtmp, sum1, sum2, tid
    integer :: ncpu, nitr, status, Lx, Ly, nr, nwl, nlay, maskx, masky
	integer, allocatable :: ushift(:,:), vshift(:,:)    ! offset vector
    real(kind=8) :: p, phiVar, dphi, dphis(2)
    real(kind=8) :: rand_str, rand_prd, rand_rng, rnd1, rnd2, buf1
    real(kind=8) :: db1, db2
    real(kind=8), dimension(:,:), allocatable :: target_img
	real(kind=8), dimension(:,:,:), allocatable :: solve_imgs
    real(kind=8), dimension(:,:,:), allocatable :: bufm
    complex(kind=8), dimension(:,:), allocatable :: g0
    complex(kind=8), dimension(:,:,:), allocatable :: gm, gbuf, gbuf2
    real(kind=8) :: dE, Etold, I0, alpha
    real(kind=8), dimension(:), allocatable :: Etnewm, Etoldm
    real(kind=8), dimension(:,:), allocatable :: err_history
    integer*4 :: now(3),iseed, messagen
    integer :: time(0:1), rate
    integer :: time2(0:1), rate2
	character(99) :: msg

    ! rand-seed
	call itime(now)          !now(1)=hour, (2)=minute, (3)=second
    iseed=int(now(2)*now(3)) 
    call random_seed(iseed)

	
	! read vars
    open(3, file = 'tmp1', status = 'old', action = 'read', iostat = status)
    if(status==0) then
        read(3,*) ncpu
		read(3,*) nr
		read(3,*) nwl
		read(3,*) phiVar
		read(3,*) rand_str
		read(3,*) rand_prd
		read(3,*) nitr
		read(3,*) nlay
		read(3,*) nx
		read(3,*) ny
		read(3,*) maskx
		read(3,*) masky
		read(3,*) dphi
		read(3,*) I0
		read(3,*) Etold
		allocate(ushift(nlay,nr),vshift(nlay,nr))
		do i = 1, nlay
            read(3,*) ushift(i,1:nr) 
        end do
		do i = 1, nlay
            read(3,*) vshift(i,1:nr) 
        end do
    end if
	
    close(3)
	! read target_img
	open(4, file = 'tmp2', status = 'old', action = 'read', iostat = status)
    if(status==0) then
        allocate(target_img(nx,ny))
        do i = 1, nx
            read(4,*) target_img(i,1:ny) 
        end do
    end if
    close(4)
	! read g:real part and imaginary part
	allocate(g0(nx,ny),bufm(nx,ny,1))
	open(33, file = 'tmp3', status = 'old', action = 'read', iostat = status)
    if(status==0) then
        do i = 1, nx
            read(33,*) bufm(i,1:ny,1) 
        end do
    end if
    close(33)
    g0 = c_one*bufm(:,:,1)
	open(43, file = 'tmp4', status = 'old', action = 'read', iostat = status)
    if(status==0) then
        do i = 1, nx
            read(43,*) bufm(i,1:ny,1) 
        end do
    end if
    close(43)
    g0 = g0+c_i*bufm(:,:,1)
	deallocate(bufm)
	! read mask
	open(13, file = 'tmp5', status = 'old', action = 'read', iostat = status)
	if(status==0) then
		allocate(solve_imgs(nx,ny,nlay))
		do n=1,nlay
			do i = 1, nx
				read(13,*) solve_imgs(i,1:ny,n) 
			end do
		end do
	end if
    close(13)
	! build table
	if(nx==ny) then
        allocate(xlst(1+2*nx*nx))
        call buildtable1(nx, xlst)
    else
        allocate(xlst(1+nx*nx),ylst(1+ny*ny))
        call buildtable2(nx, ny, xlst, ylst)
    end if
	! allocate other arrays
	allocate(err_history(nitr,2))
	allocate(gbuf(nx,ny,ncpu),gbuf2(nx,ny,ncpu),bufm(nx,ny,ncpu))
    allocate(gm(nx,ny,ncpu))
    allocate(Etoldm(ncpu),Etnewm(ncpu))
	
	
	
	call OMP_SET_NUM_THREADS(ncpu)
    messagen = 12
    ! set initial parameters
    err_history(1,1) = 1.0
    err_history(1,2) = Etold
    do i = 1, ncpu
        gm(:,:,i) = g0
        Etoldm(i) = Etold
    end do
	! solving
	if(nwl==0) then
		do n = 1, nitr
			call system_clock(time(0),rate)
			rand_rng = rand_str*exp(-3.0*n/rand_prd)
			p =  1.0/(1.0+exp(0.1/(rand_rng+0.001)))
			do u = 1, nx, nr
				if(mod(u-1,messagen)==0) call system_clock(time2(0),rate2)
				!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tid,v,buf1,rnd1,rnd2,db1,db2,w,s,alpha,dE,ni,nj,layi,utmp,vtmp,sum1,sum2)
				!$OMP DO
				do v = 1, ny, nr
					tid = OMP_GET_THREAD_NUM()+1
					do layi = 1,nlay
						utmp = ushift(layi,nr)+u
						vtmp = vshift(layi,nr)+v
						buf1 = solve_imgs(utmp,vtmp,layi)
						call random_number(rnd1)
						rnd1=1.0*anint(rnd1)
						forall(w=1:ny) gbuf(:,w,tid) = gm(:,w,tid)
						do ni = 0,nr-1
							do nj = 0,nr-1
								utmp = ushift(layi,nr)+u+ni
								vtmp = vshift(layi,nr)+v+nj
								if (utmp>nx) then
									utmp=utmp-nx
								end if
								if (vtmp>ny) then
									vtmp=vtmp-ny
								end if
								sum1 = 0
								do w=1,nlay
									sum1=sum1+solve_imgs(utmp,vtmp,w)
								end do
								solve_imgs(utmp,vtmp,layi) = rnd1
								sum2 = 0
								do w=1,nlay
									sum2=sum2+solve_imgs(utmp,vtmp,w)
								end do
								! update g0
								if(nx==ny) then
									forall(w=1:nx,s=1:ny) gbuf2(w,s,tid) = xlst((w-1)*(utmp-1)+(s-1)*(vtmp-1)+1)
								else
									forall(w=1:nx,s=1:ny) gbuf2(w,s,tid) = xlst((w-1)*(utmp-1)+1)*ylst((s-1)*(vtmp-1)+1)
								end if
								forall(w=1:ny) gbuf(:,w,tid) = gbuf(:,w,tid)+(1.0/nx/ny)*(exp(c_i*sum2*dphi)-exp(c_i*sum1*dphi))*gbuf2(:,w,tid)
							end do
						end do
						forall(w=1:ny) bufm(:,w,tid) = real(gbuf(:,w,tid))**2+aimag(gbuf(:,w,tid))**2
						alpha = I0/sum(real(bufm(:,:,tid),kind=8))
						Etnewm(tid) = sum( real(target_img-alpha*bufm(:,:,tid),kind=8)**2 )
						dE = Etnewm(tid)-Etoldm(tid)
						call random_number(rnd2)
						if(dE<0.0) then
							if(rnd2>p) then
								Etoldm(tid) = Etnewm(tid)
								forall(w=1:ny) gm(:,w,tid) = gbuf(:,w,tid)
							else
								do ni = 0,nr-1
									do nj = 0,nr-1
										utmp = ushift(layi,nr)+u+ni
										vtmp = vshift(layi,nr)+v+nj
										if (utmp>nx) then
											utmp=utmp-nx
										end if
										if (vtmp>ny) then
											vtmp=vtmp-ny
										end if
										solve_imgs(utmp,vtmp,layi) = buf1
									end do
								end do
								
							end if
						else
							if(rnd2>p) then
								do ni = 0,nr-1
									do nj = 0,nr-1
										utmp = ushift(layi,nr)+u+ni
										vtmp = vshift(layi,nr)+v+nj
										if (utmp>nx) then
											utmp=utmp-nx
										end if
										if (vtmp>ny) then
											vtmp=vtmp-ny
										end if
										solve_imgs(utmp,vtmp,layi) = buf1
									end do
								end do
							else
								Etoldm(tid) = Etnewm(tid)
								forall(w=1:ny) gm(:,w,tid) = gbuf(:,w,tid)
							end if
						end if
					end do
				end do
				!$OMP END DO
				!$OMP BARRIER
				!$OMP END PARALLEL
				
				!!$OMP SINGLE
				gbuf(:,:,1) = g0
				do w = 1, ncpu
					g0 = g0+(gm(:,:,w)-gbuf(:,:,1))
				end do
				bufm(:,:,1) = real(g0(:,:))**2+aimag(g0(:,:))**2
				alpha = I0/sum(real(bufm(:,:,1),kind=8))
				Etold = sum(real(1.0*target_img-alpha*bufm(:,:,1),kind=8)**2)
				forall(w=1:ncpu) gm(:,:,w) = g0
				forall(w=1:ncpu) Etoldm(w) = Etold
				!!$OMP END SINGLE
				! show message
				if(mod(u-1,messagen)==0) then
					call system_clock(time2(1))
					write(msg,'(I4,A15,F7.3,A2,A25,F9.3,A2,F15.0)') n,'-th iteration:',100.0*u/nx,' %',';     t(second)/error: ', (time2(1)-time2(0))/(1.0*rate2), ' /', Etold
					call mexWarnMsgTxt(trim(msg))
				end if
			end do
			write(msg,'(I4,A15,F7.3,A2)') n,'-th iteration:',100.0*nx/nx,' %'
			call mexWarnMsgTxt(trim(msg))
			err_history(n,1) = n
			err_history(n,2) = Etold
			! display message
			call system_clock(time(1))
			write(msg,'(A,A16,I6,A1,I6,A28,F12.5,A2,F15.0)') creturn,'grid-type itr.: ', n, '/', nitr ,';     t(second)/error: ', (time(1)-time(0))/(1.0*rate), ' /', Etold
			call mexWarnMsgTxt(trim(msg))
		end do
	end if
	
	
	
	
	
	
	! export results
    open(unit=11 , file='errhist.txt' , status='replace')
    do i = 1, nitr    
        write(11,'(F7.0,A,E)') err_history(i,1),tab,err_history(i,2)
    end do
    close(11)
    ! export solve_img
    open(unit=21 , file='masko.txt' , status='replace')
	do n=1,nlay
		do i = 1, nx    
			do j = 1, ny-1
				write(21,'(F3.0,A)',ADVANCE='NO') solve_imgs(i,j,n),tab
			end do
			write(21,'(F3.0)') solve_imgs(i,ny,n)
		end do
	end do
    close(21)
	! export g0
	open(unit=31 , file='g0real.txt' , status='replace')
    do i = 1, nx    
		do j = 1, ny-1
			write(31,'(F9.5,A)',ADVANCE='NO') real(g0(i,j)),tab
		end do
		write(31,'(F9.5)') real(g0(i,j))
	end do
    close(31)
	open(unit=51 , file='g0imag.txt' , status='replace')
    do i = 1, nx    
		do j = 1, ny-1
			write(51,'(F9.5,A)',ADVANCE='NO') aimag(g0(i,j)),tab
		end do
		write(51,'(F9.5)') aimag(g0(i,j))
	end do
    close(51)
	
	! free array
	deallocate(target_img,g0)
	deallocate(solve_imgs)
	if(nx==ny) then
		deallocate(xlst)
	else
		deallocate(xlst,ylst)
	end if
	deallocate(err_history)
	deallocate(gbuf,gbuf2,bufm)
    deallocate(gm,ushift,vshift)
    deallocate(Etoldm,Etnewm)
end subroutine mexFunction


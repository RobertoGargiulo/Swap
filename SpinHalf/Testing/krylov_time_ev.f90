Program time_evolution
use exp_sparse

!D dimension of original hilbert space
!


	implicit none
	integer D, i, j, T_max, initial, nz, m_trunc, Trunc, L, nu, distance
	integer, parameter :: nzmax= 500000
	real*8 Utau, Jhop, energy, coeff, bord_prob
	real*8, dimension(:), allocatable :: V
	integer, dimension(:), allocatable:: Hi,Hj,Hi_ex, Hj_ex
	integer, dimension(:,:), allocatable:: n_site
	complex*16, dimension(:), allocatable :: Hop,Hop_ex
	integer argv
	complex*16, dimension(:), allocatable :: state,newstate
	complex, parameter :: im_uni = (0, 1) !imaginary unit
	
	character(len=50):: str
	argv=command_argument_count()
	if (argv.ne.11) then
		write (*,*)  'usage is ./time_evolution initial_state L hilbert_file V_file hopping_file'
		write (*,*)	 'tau J T_max M_trunc Trunc_hilbert output_energy_file'
		write (*,*)  'initial state:'
		write (*,*)	 '1 for larger_trunc_build'
		write (*,*)  'D for all other programs'
		call exit
	end if
	call get_command_argument(2,str)
	read(str,*) L 
	call get_command_argument(3,str)
	open(1,action='read', file=str)
	call get_command_argument(4,str)
	open(2,action='read', file=str)
	call get_command_argument(5,str)
	open(3,action='read', file=str)
	call get_command_argument(6,str)
	read(str,*) Utau
	call get_command_argument(7,str)
	read(str,*) Jhop
	call get_command_argument(8,str)
	read(str,*) T_max
	call get_command_argument(9,str)
	read(str,*) M_trunc
	call get_command_argument(10,str)
	read(str,*) Trunc
	call get_command_argument(11,str)
	open(7,action='write', file=str)

	!read hilbert space dimension
	read(1,*) D
	allocate(n_site(D,L))
	do i=1,D
		read(1,*) n_site(i,:)
	end do
	close(1)
	nu=sum(n_site(1,:))/L
	
	write(*,*) 'hilbert space read'
	
	call get_command_argument(1,str)
	if (str=='1') then
		initial=1
	else if (str=='D') then
		initial =D
	else
		write(*,*) 'non valid initial state'
		return
	end if
	
	!read potential energy matrix
	allocate(V(D))
	do i=1,D
		read(2,*) V(i)
	end do
	close(2)
	write(*,*) 'V read'
	
	!read hopping matrix
	allocate(Hi_ex(nzmax))
	allocate(Hj_ex(nzmax))
	allocate(Hop_ex(nzmax))
	nz=1
	do
		read(3,*, end=10) i,j,coeff
		Hi_ex(nz)=i
		Hj_ex(nz)=j
		Hop_ex(nz)=coeff
		nz=nz+1
		if (Hi_ex(nz-1).ne.Hj_ex(nz-1)) then
			Hj_ex(nz)=Hi_ex(nz-1)
			Hi_ex(nz)=Hj_ex(nz-1)
			Hop_ex(nz)=Hop_ex(nz-1)
			nz=nz+1
		end if
	end do
10	nz=nz-1
	write(*,*) nz
	
	allocate(Hop(nz))
	Hop=-im_uni*Hop_ex(1:nz)
	deallocate(Hop_ex)
	
	allocate(Hi(nz))
	Hi=Hi_ex(1:nz)
	deallocate(Hi_ex)
	
	allocate(Hj(nz))
	Hj=Hj_ex(1:nz)
	deallocate(Hj_ex)	
	
	!calculate time evolution of the state
	allocate(state(D))
	allocate(newstate(D))
	!build initial state
	state=0
	state(initial)=1
	!save the initial potential energy
	energy = V(D)
	
	!evolve and sample potential energy
	write(7,*) energy, 0
	do i=1,T_max
		call evolve(D, nz, M_trunc, Hi, Hj, Hop, state, Jhop, newstate)
		energy=0
		do j=1,D
			energy=energy+V(j)*abs(newstate(j))**2
			state(j)=newstate(j)*exp(-im_uni*Utau*V(j))
		end do
		bord_prob=0
		do j=1,D
			distance=sum(abs(n_site(j,:)-nu))/2
			if (distance.eq.Trunc) then
				bord_prob=bord_prob+abs(state(j))**2
			end if
		end do
		write(7,*) energy/sum(abs(state)**2), bord_prob
	end do
	
	deallocate(V)
	deallocate(Hop)
	deallocate(Hi)
	deallocate(Hj)
	close(7)
	close(3)	
end Program time_evolution

module mod_generate_kpoints
	use mod_f90_kind
	use mod_parameters, only: a0, ncp, runoptions
	use mod_constants
	implicit none
!========================================================================================!
!	wave vector integration variables
	integer :: nkpoints
	real(double),allocatable      :: kbz(:,:),wkbz(:),kbz2d(:,:)

contains
	!-- Generating integration points in the Brillouin Zone --
	subroutine generate_kpoints_bcc110()
		use MPI
		use mod_mpi_pars
		implicit none
		integer         :: i,j,AllocateStatus
		integer					:: m0,m1,imax,ki,jmax,kj,iz
		real(double)    :: beta,beta2,fact
		real(double)    :: akxp(4),akzp(4),akxt,akzt,akx,aky,akz
		!   Generation of the Cunningham points for a centred rectangular
		!   lattice. The units are chosen with the lattice parameter a0.

		allocate( kbz(4**(ncp+1),3),wkbz(4**(ncp+1)),kbz2d(4**(ncp+1),2), STAT=AllocateStatus )
		if (AllocateStatus.ne.0) then
	    write(*,"('[mod_generate_kpoints] Not enough memory for: kbz,wkbz,kbz2d')")
			call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
		end if

!--- k integration parameters ------
		beta  = 0.5d0*sq2
		beta2 = beta*beta
		fact  = sq2/4.d0
		m0 = 2**ncp
		m1 = m0*2
		imax = m0 - 0.5d0
!-----------------------------------
		! Generating k points
		nkpoints   = 0
		do ki = 0,imax
			i      = 2*ki+1
			akxp(1) = tpi*fact*dble(i)/(a0*dble(m0))
			jmax    = nint(beta2*(-2*ki+m0-1))+m0-0.5d0
			do kj=0,jmax
				j = 2*kj+1
				akzp(1) = tpi*fact*dble(j)/(a0*dble(m1)*beta)
	!     Generating k_// inside the full BZ
	!     2nd quadrant
				akxp(2) =-akxp(1)
				akzp(2) = akzp(1)
	!     3rd quadrant
				akxp(3) =-akxp(1)
				akzp(3) =-akzp(1)
	!     4th quadrant
				akxp(4) = akxp(1)
				akzp(4) =-akzp(1)

				do iz=1,4
					nkpoints = nkpoints + 1

					akxt = akxp(iz)
					akzt = akzp(iz)

					! Storing the kpoints in the 2DBZ
					kbz2d(nkpoints,1) = akxt
					kbz2d(nkpoints,2) = akzt

					! Transformation to cartesian axis
					akx = akxt*beta
					aky = akx
					akz = akzt

					! Storing the kpoints in the BZ transformed to the cartesian system
					kbz(nkpoints,1) = akx
					kbz(nkpoints,2) = aky
					kbz(nkpoints,3) = akz
				end do ! iz
			end do ! kj
		end do ! ki

		wkbz(:) = 1.d0/dble(nkpoints)

		if((myrank.eq.0))  write(*,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints

		return
	end subroutine generate_kpoints_bcc110

	subroutine generate_kpoints_fcc100()
		use MPI
		use mod_mpi_pars
		implicit none
		integer         :: iz,AllocateStatus
		integer					:: ki,kj,nkmax,icount
		real(double)		:: auxk,auxw
		real(double)    :: akp(8,2)
!   generation of Cunningham points for a square lattice
!   the units are chosen using the variable for the lattice parameter a0.

		if(ncp.le.0) then
	    write(*,"('[mod_generate_kpoints] ncp must be greater than 0')")
			call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
		end if
		nkpoints = 8*(2**(ncp-1))*(1+2**ncp)
		allocate( kbz(nkpoints,3),wkbz(nkpoints),kbz2d(nkpoints,2), stat=AllocateStatus )
		if (AllocateStatus.ne.0) then
	    write(*,"('[mod_generate_kpoints] Not enough memory for: kbz,wkbz,kbz2d')")
			call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
		end if
		nkmax = 2**ncp
		auxk  = dble(2*nkmax)
		auxw  = 8.d0*dble(2**(ncp+ncp-1))

		icount = 0
		do ki=1,nkmax
			akp(1,1)=pi*(dble(2*ki-1))/(auxk*a0)
			do kj=1,ki
				akp(1,2)=pi*(dble(2*kj-1))/(auxk*a0)

				! generating k_// inside the full bz
				! 1st quadrant
				akp(2,1) = akp(1,2)
				akp(2,2) = akp(1,1)
				! 2nd quadrant
				akp(3,1) =-akp(2,1)
				akp(3,2) = akp(2,2)
				akp(4,1) =-akp(3,2)
				akp(4,2) =-akp(3,1)
				! 3rd quadrant
				akp(5,1) = akp(4,1)
				akp(5,2) =-akp(4,2)
				akp(6,1) = akp(5,2)
				akp(6,2) = akp(5,1)
				! 4th quadrant
				akp(7,1) =-akp(6,1)
				akp(7,2) = akp(6,2)
				akp(8,1) =-akp(7,2)
				akp(8,2) =-akp(7,1)

				do iz=1,8
					icount = icount + 1

					kbz2d(icount,1) = akp(iz,1)
					kbz2d(icount,2) = akp(iz,2)

					kbz(icount,1) = akp(iz,1)+akp(iz,2)
					kbz(icount,2) =-akp(iz,1)+akp(iz,2)
					kbz(icount,3) = 0.d0

					! Weight of k point
					if(ki.eq.kj) then ! Half of the weight if the point is in the diagonal
						wkbz(icount)=0.5d0/auxw
					else
						wkbz(icount)=1.d0/auxw
					end if

				end do
			end do
		end do

		if(icount.ne.nkpoints) then
			 if(myrank.eq.0) write(*,"('[mod_generate_kpoints] Incorrect number of points: icount = ',i0,', nkpoints = ',i0)") icount,nkpoints
			 call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
		end if
		if((myrank.eq.0))  write(*,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints

		return
	end subroutine generate_kpoints_fcc100

  subroutine generate_kpoints_fcc111()
    ! Generation of the k-points points for an hexagonal lattice.
    ! The generated k-points units are chosen using a0.
		use MPI
		use mod_mpi_pars
		use mod_constants, only: pi,sq3
    implicit none
    integer         :: iz,AllocateStatus
    integer         :: ki,kj
    real(kind=8)    :: a1(2),a2(2),b1(2),b2(2),c1(2),c2(2)
    real(kind=8)    :: weight,dnkpoints
    real(kind=8)    :: kv(2),angle,rm(2,2)
    real(double), dimension(3,3) :: rm2,rm3
    real(kind=8), dimension(100000) 	:: wk
    real(kind=8), dimension(100000,2) :: kvm
		real(double), dimension(3) :: qvl,k2anx

!		Lattice vetors in real space
		a1 = [0.5d0, sqrt(3.d0)/2.d0]*a0
		a2 = [0.5d0,-sqrt(3.d0)/2.d0]*a0

!		Reciprocal lattice determination
		b1 =         [a2(2),-a2(1)]
		b1 = 2.d0*pi*[a2(2),-a2(1)]/dot_product(a1,b1)
		b2 =         [a1(2),-a1(1)]
		b2 = 2.d0*pi*[a1(2),-a1(1)]/dot_product(a2,b2)

    angle = -pi/6.d0
    rm(1,1) = cos(angle)
    rm(2,1) = sin(angle)
    rm(1,2) = -rm(2,1)
    rm(2,2) = rm(1,1)

    c1 = (b1+b2)/3.d0
    c2 = b1-c1
    c1 = matmul(rm,c1)
    c2 = matmul(rm,c2)

    kvm = 0.d0

    ! Gamma point
    nkpoints = 1
    kvm(nkpoints,:) = [ 0.d0,0.d0 ]
    wk(nkpoints) = 1.d0
    dnkpoints = wk(nkpoints)

    do ki=1,ncp
      do kj=1,ncp
        weight = 1.d0
        kv = ( ki*c1 + (kj-1.d0)*c2 ) / DBLE(ncp)

        if(kv(1).gt.2.d0*pi/sq3) cycle
        if(kv(1).eq.2.d0*pi/sq3) weight = 0.5d0

        if((ki.eq.ncp).and.(kj.eq.1)) weight = 1.d0/3.d0

        dnkpoints = dnkpoints + weight
        nkpoints = nkpoints + 1
        kvm(nkpoints,:) = kv
        wk(nkpoints) = weight

        do iz=1,5
          angle = dble(iz)*pi/3.d0
          rm(1,1) = cos(angle)
          rm(1,2) = -sin(angle)
          rm(2,1) = sin(angle)
          rm(2,2) = cos(angle)

          dnkpoints = dnkpoints + weight
          nkpoints = nkpoints + 1
          kvm(nkpoints,:) = matmul(rm,kv)
          wk(nkpoints) = weight
        end do
      end do
    end do

		allocate( kbz(nkpoints,3),wkbz(nkpoints),kbz2d(nkpoints,2), stat=AllocateStatus )
		if (AllocateStatus.ne.0) then
	    write(*,"('[mod_generate_kpoints] Not enough memory for: kbz,wkbz,kbz2d')")
			call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
		end if

		! Rotation matrices to transform 2D BZ to 3D BZ
	  rm2 = 0.d0
	  rm2(1,1) = 1.d0/sq3
	  rm2(1,3) = sq2/sq3
	  rm2(3,1) = -sq2/sq3
	  rm2(2,2) = 1.d0
	  rm2(3,3) = 1.d0/sq3

	  rm3 = 0.d0
	  rm3(1,1) = 1.d0/sq2
	  rm3(1,2) = -1.d0/sq2
	  rm3(2,1) = 1.d0/sq2
	  rm3(2,2) = 1.d0/sq2
	  rm3(3,3) = 1.d0

    weight = 0.d0
    do iz = 1,nkpoints
      kbz2d(iz,:) = kvm(iz,:)
      k2anx(1) = kvm(iz,1)
      k2anx(2) = kvm(iz,2)
      k2anx(3) = 0.d0
      qvl = matmul(rm2,k2anx)
      kbz(iz,:) = matmul(rm3,qvl)
      wkbz(iz) = wk(iz)/dnkpoints
      weight = weight + wkbz(iz)
    end do

		if((myrank.eq.0))  write(*,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints

!     write(*,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints
!     write(*,"('[mod_generate_kpoints] ',e9.2,' total weight.')") weight
!     write(*,"('[mod_generate_kpoints] ',f15.2,' dnkpoints.')") dnkpoints

    return
  end subroutine generate_kpoints_fcc111
end module mod_generate_kpoints
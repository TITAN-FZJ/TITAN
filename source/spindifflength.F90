! --------- Sum over wave vectors to calculate spin diffusion length ----------
subroutine spindifflength(e,pos,inn,sdl)
	use mod_f90_kind
	use mod_parameters
	use mod_constants
	use mod_generate_kpoints
	use mod_lattice
	implicit none
	integer         :: i,iz,pos,inn
	real(double)    :: e
	real(double)    :: kp(3)
	complex(double)	:: ikr,expikr,sdl(Npl),sdlgf
	complex(double),dimension(Npl,Npl,18,18)        :: gf
	complex(double),dimension(Npl,Npl,9,9)        	:: gfud,gfdu

	sdl    = zero
	do iz=1,nkpoints
		kp = kbz(iz,:)

		ikr = zi*pos*(kp(1)*r0(inn,1)+kp(2)*r0(inn,2)+kp(3)*r0(inn,3))
		expikr = exp(-ikr)

		! Green function at (k,e+ieta)
		call green(e,eta,kp,gf)

		gfud = gf(:,:, 1: 9,10:18)
		gfdu = gf(:,:,10:18, 1: 9)

		do i=1,Npl
			sdlgf = sum(gfud(i,i,:,:)) - conjg(sum(gfdu(i,i,:,:)))
			sdl(i) = sdl(i) + (expikr*sdlgf*wkbz(iz))
		end do

	end do
	sdl = zi*sdl

	return
end subroutine spindifflength
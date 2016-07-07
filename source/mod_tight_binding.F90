! This module calculates the two center integrals (as given by Slater-Koster)
! and calculates the hamiltonian matrix components in real space
! It also includes all the parameters obtained from first principles calculations
! The system is defined in module mod_define_system, where the parameters
! to be used are chosen in variable mmlayer(Npl+2)

! Material indices:
!------------------------
! 1 - Empty spheres     |
! 2 - Fe (S)            |
! 3 - W (S-1)           |    set 1
! 4 - W (S-2)           |   BCC(110)
! 5 - W (S-3)           |
! 6 - W (blk)           |
!------------------------
! 7  - Empty spheres    |
! 8  - W (S)            |
! 9  - W (S-1)          |    set 2
! 10 - W (S-2)          |   BCC(110)
! 11 - W (S-3)          |
! 12 - W (blk)          |
!------------------------
! 13 - Empty spheres    |
! 14 - Co (S)           |
! 15 - Pt (S-1)         |    set 3
! 16 - Pt (S-2)         |   FCC(100)
! 17 - Pt (S-3)         |
! 18 - Pt (blk)         |
!------------------------
! 19 - Empty spheres    |
! 20 - Pt (S)           |
! 21 - Pt (S-1)         |    set 4
! 22 - Pt (S-2)         |   FCC(100)
! 23 - Pt (S-3)         |
! 24 - Pt (blk)         |
!------------------------

module mod_tight_binding
  use mod_f90_kind
  implicit none
! SOC parameters
  real(double), dimension(:), allocatable  :: lambda
! Total cccupations per layer
  real(double), dimension(:), allocatable  :: npart0
! real-space hoppings
  real(double),allocatable  :: t00(:,:,:,:),t01(:,:,:,:),t02(:,:,:,:)
! SOC L.S matrix
  complex(double), dimension(18,18) :: ls

contains

  subroutine DFT_parameters(cs,cp,cd,ds,dp,dd)
    use mod_f90_kind
    use mod_parameters, only: Npl, Ef, SOC, dfttype, U, Utype, layertype, mmlayer, ry2ev, outputunit
    use mod_mpi_pars
  integer       :: i,j
  integer, parameter :: mmmax=24
  integer, dimension(mmmax) :: layertype_dft
  real(double)              :: shft
  real(double),dimension(mmmax) :: n0s,n0p,n0d,lambdasoc_dft,Ef_dft,U_dft
  real(double),dimension(mmmax) :: cs_ort,cp_ort,cd_ort,ds_ort,dp_ort,dd_ort
  real(double),dimension(mmmax) :: cs_bar,cp_bar,cd_bar,ds_bar,dp_bar,dd_bar
  real(double),dimension(Npl+2),intent(out) :: cs,cp,cd,ds,dp,dd

  n0s    = 999.d0
  n0p    = 999.d0
  n0d    = 999.d0
  cs_bar = 999.d0
  cp_bar = 999.d0
  cd_bar = 999.d0
  ds_bar = 999.d0
  dp_bar = 999.d0
  dd_bar = 999.d0
  cs_ort = 999.d0
  cp_ort = 999.d0
  cd_ort = 999.d0
  ds_ort = 999.d0
  dp_ort = 999.d0
  dd_ort = 999.d0

  layertype_dft = 0
  lambdasoc_dft = 0.d0
  Ef_dft = 0.d0
  U_dft = 0.d0

!----------------------- BEGINNING OF PARAMETERS -----------------------

! Angela's unrelaxed calculations Fe/W(110)
! no SOC, paramagnetic
!
! occupations obtained from spin-polarized (up+down),
! with SOC calculation
!
! Material indices: 1 - Empty spheres
!                   2 - Fe (S)
!                   3 - W (S-1)
!                   4 - W (S-2)
!                   5 - W (S-3)
!                   6 - W (blk)
!

  layertype_dft(1) = 1
  layertype_dft(2) = 2
  layertype_dft(6) = 3

  ! Fermi energy
  Ef_dft(1:6)  = 0.14247d0

  ! Effective intra-site coulomb interaction U (in eV)
  U_dft(2:6) = 1.d0

  ! SOC parameters
  if(SOC) then
    lambdasoc_dft(2)   = 0.004d0  ! Fe
    lambdasoc_dft(3:6) = 0.019d0  ! W
  end if

  ! Empty spheres
  n0s(1)    = 0.047d0
  n0p(1)    = 0.050d0
  n0d(1)    = 0.023d0
  cs_bar(1) = 0.452107d0
  cp_bar(1) = 0.526667d0
  cd_bar(1) = 0.289442d0
  ds_bar(1) = 0.216817d0
  dp_bar(1) = 0.114156d0
  dd_bar(1) = 0.028997d0
!   cs_ort(1) =
!   cp_ort(1) =
!   cd_ort(1) =
!   ds_ort(1) =
!   dp_ort(1) =
!   dd_ort(1) =

  ! Fe (S)
  n0s(2)    = 0.722d0
  n0p(2)    = 0.630d0
  n0d(2)    = 6.569d0
  cs_bar(2) = -0.059955d0
  cp_bar(2) = 0.486489d0
  cd_bar(2) = 0.091421d0
  ds_bar(2) = 0.358579d0
  dp_bar(2) = 0.233078d0
  dd_bar(2) = 0.100813d0
!   cs_ort(2) =
!   cp_ort(2) =
!   cd_ort(2) =
!   ds_ort(2) =
!   dp_ort(2) =
!   dd_ort(2) =

  ! W (S-1)
  n0s(3)    = 0.824d0
  n0p(3)    = 1.005d0
  n0d(3)    = 4.114d0
  cs_bar(3) = -0.231327d0
  cp_bar(3) = 0.504396d0
  cd_bar(3) = 0.118289d0
  ds_bar(3) = 0.390632d0
  dp_bar(3) = 0.234604d0
  dd_bar(3) = 0.181033d0
!   cs_ort(3) =
!   cp_ort(3) =
!   cd_ort(3) =
!   ds_ort(3) =
!   dp_ort(3) =
!   dd_ort(3) =

  ! W (S-2)
  n0s(4)    = 0.828d0
  n0p(4)    = 1.077d0
  n0d(4)    = 4.095d0
  cs_bar(4) =-0.226435d0
  cp_bar(4) = 0.499843d0
  cd_bar(4) = 0.126014d0
  ds_bar(4) = 0.385574d0
  dp_bar(4) = 0.232620d0
  dd_bar(4) = 0.180319d0
!   cs_ort(4) =
!   cp_ort(4) =
!   cd_ort(4) =
!   ds_ort(4) =
!   dp_ort(4) =
!   dd_ort(4) =

  ! W (S-3)
  n0s(5)    = 0.841d0
  n0p(5)    = 1.148d0
  n0d(5)    = 4.067d0
  cs_bar(5) =-0.218636d0
  cp_bar(5) = 0.497383d0
  cd_bar(5) = 0.132964d0
  ds_bar(5) = 0.379463d0
  dp_bar(5) = 0.230094d0
  dd_bar(5) = 0.179174d0
!   cs_ort(5) =
!   cp_ort(5) =
!   cd_ort(5) =
!   ds_ort(5) =
!   dp_ort(5) =
!   dd_ort(5) =

  ! W (bulk)
  n0s(6)    = 0.423d0*2.d0
  n0p(6)    = 0.539d0*2.d0
  n0d(6)    = 2.037d0*2.d0
  cs_bar(6) =-0.255781d0
  cp_bar(6) = 0.558431d0
  cd_bar(6) = 0.145829d0
  ds_bar(6) = 0.422935d0
  dp_bar(6) = 0.250566d0
  dd_bar(6) = 0.191739d0
!   cs_ort(6) =
!   cp_ort(6) =
!   cd_ort(6) =
!   ds_ort(6) =
!   dp_ort(6) =
!   dd_ort(6) =


! Angela's unrelaxed calculations for W (110)
! no SOC, paramagnetic
!
! occupations obtained from paramagnetic,
! with SOC calculation
!
! Material indices:
!                   7  - Empty spheres
!                   8  - W (S)
!                   9  - W (S-1)
!                   10 - W (S-2)
!                   11 - W (S-3)
!                   12 - W (blk)
!
  layertype_dft(7) = 1
  layertype_dft(12) = 3

  ! Fermi energy
  Ef_dft(7:12)  = 0.03313d-3

  ! Effective intra-site coulomb interaction U (in eV)
  U_dft(8:12) = 1.d0

  ! SOC parameters
  if(SOC) then
    lambdasoc_dft(8:12) = 0.019d0  ! W
  end if

  ! Empty spheres
  n0s(7)    = 0.043d0*2.d0
  n0p(7)    = 0.047d0*2.d0
  n0d(7)    = 0.023d0*2.d0
  cs_bar(7) = 0.298842d0
  cp_bar(7) = 0.420951d0
  cd_bar(7) = 0.221018d0
  ds_bar(7) = 0.239901d0
  dp_bar(7) = 0.130442d0
  dd_bar(7) = 0.041150d0
  cs_ort(7) =-0.00649443d0
  cp_ort(7) = 0.92907559d0
  cd_ort(7) = 2.31562212d0
  ds_ort(7) = 0.36953380d0
  dp_ort(7) = 0.39281758d0
  dd_ort(7) = 0.41988460d0

  ! W (S)
  n0s(8)    = 0.401d0*2.d0
  n0p(8)    = 0.413d0*2.d0
  n0d(8)    = 2.037d0*2.d0
  cs_bar(8) =-0.317014d0
  cp_bar(8) = 0.420629d0
  cd_bar(8) = 0.026219d0
  ds_bar(8) = 0.400275d0
  dp_bar(8) = 0.234706d0
  dd_bar(8) = 0.181658d0
  cs_ort(8) =-0.31923450d0
  cp_ort(8) = 0.84254832d0
  cd_ort(8) = 0.03538781d0
  ds_ort(8) = 0.40473197d0
  dp_ort(8) = 0.40770215d0
  dd_ort(8) = 0.19296136d0

  ! W (S-1)
  n0s(9)    = 0.420d0*2.d0
  n0p(9)    = 0.548d0*2.d0
  n0d(9)    = 2.075d0*2.d0
  cs_bar(9) =-0.357662d0
  cp_bar(9) = 0.374054d0
  cd_bar(9) = 0.001349d0
  ds_bar(9) = 0.387063d0
  dp_bar(9) = 0.233989d0
  dd_bar(9) = 0.181245d0
  cs_ort(9) =-0.31943923d0
  cp_ort(9) = 0.84308686d0
  cd_ort(9) = 0.05507815d0
  ds_ort(9) = 0.40704971d0
  dp_ort(9) = 0.40977848d0
  dd_ort(9) = 0.19729412d0

  ! W (S-2)
  n0s(10)    = 0.412d0*2.d0
  n0p(10)    = 0.537d0*2.d0
  n0d(10)    = 2.042d0*2.d0
  cs_bar(10) =-0.337268d0
  cp_bar(10) = 0.384353d0
  cd_bar(10) = 0.011980d0
  ds_bar(10) = 0.384590d0
  dp_bar(10) = 0.231537d0
  dd_bar(10) = 0.179818d0
  cs_ort(10) =-0.32885989d0
  cp_ort(10) = 0.83395538d0
  cd_ort(10) = 0.03552894d0
  ds_ort(10) = 0.40629836d0
  dp_ort(10) = 0.40959958d0
  dd_ort(10) = 0.19608482d0

  ! W (S-3)
  n0s(11)    = 0.413d0*2.d0
  n0p(11)    = 0.542d0*2.d0
  n0d(11)    = 2.045d0*2.d0
  cs_bar(11) =-0.335601d0
  cp_bar(11) = 0.392339d0
  cd_bar(11) = 0.015465d0
  ds_bar(11) = 0.382509d0
  dp_bar(11) = 0.232667d0
  dd_bar(11) = 0.179847d0
  cs_ort(11) =-0.32932672d0
  cp_ort(11) = 0.83282434d0
  cd_ort(11) = 0.03597702d0
  ds_ort(11) = 0.40648561d0
  dp_ort(11) = 0.40947828d0
  dd_ort(11) = 0.19625827d0

  ! W (bulk)
  n0s(12)    = 0.423d0*2.d0
  n0p(12)    = 0.539d0*2.d0
  n0d(12)    = 2.037d0*2.d0
  cs_bar(12) =-0.332984d0
  cp_bar(12) = 0.395000d0
  cd_bar(12) = 0.017591d0
  ds_bar(12) = 0.385482d0
  dp_bar(12) = 0.232944d0
  dd_bar(12) = 0.180038d0
  cs_ort(12) =-0.32766102d0
  cp_ort(12) = 0.83455975d0
  cd_ort(12) = 0.03841926d0
  ds_ort(12) = 0.40642203d0
  dp_ort(12) = 0.40943060d0
  dd_ort(12) = 0.19629882d0


! Angela's unrelaxed calculations for Co/Pt(100)
! no SOC, paramagnetic, with hoh
!
! occupations obtained from spin-polarized,
! with SOC calculation
!
! Material indices:
!                   13 - Empty spheres
!                   14 - Co (S)
!                   15 - Pt (S-1)
!                   16 - Pt (S-2)
!                   17 - Pt (S-3)
!                   18 - Pt (blk)
!
! Experimental lattice parameter of Pt a= 3.92 \AA
! Fermi Level: E_F= -0.05529 Ry

  layertype_dft(13) = 1
  layertype_dft(14) = 2
  layertype_dft(18) = 3

  ! Fermi energy
  Ef_dft(13:18)   = -0.05529d0

  ! Effective intra-site coulomb interaction U (in eV)
  U_dft(14)    = 1.d0
  U_dft(15:18) = 0.6d0

  ! SOC parameters
  if(SOC) then
    lambdasoc_dft(14) = 0.005666d0 ! Co
    lambdasoc_dft(15:18) = 0.044d0  ! Pt
  end if

  ! Empty spheres
  n0s(13)    = 0.068d0
  n0p(13)    = 0.067d0
  n0d(13)    = 0.030d0
  cs_bar(13) = 0.24564797d0
  cp_bar(13) = 0.37056673d0
  cd_bar(13) = 0.15181117d0
  ds_bar(13) = 0.23712713d0
  dp_bar(13) = 0.12927482d0
  dd_bar(13) = 0.03894704d0
  cs_ort(13) = 0.02943149d0
  cp_ort(13) = 0.99865904d0
  cd_ort(13) = 2.43518168d0
  ds_ort(13) = 0.37579750d0
  dp_ort(13) = 0.39933182d0
  dd_ort(13) = 0.42709341d0

  ! Co (S)
  n0s(14)    = 0.632d0
  n0p(14)    = 0.527d0
  n0d(14)    = 7.663d0
  cs_bar(14) =-0.26851717d0
  cp_bar(14) = 0.26004029d0
  cd_bar(14) =-0.10677846d0
  ds_bar(14) = 0.35605769d0
  dp_bar(14) = 0.22665604d0
  dd_bar(14) = 0.09570896d0
  cs_ort(14) =-0.41644908d0
  cp_ort(14) = 0.44261695d0
  cd_ort(14) =-0.26009543d0
  ds_ort(14) = 0.37499902d0
  dp_ort(14) = 0.36824596d0
  dd_ort(14) = 0.09301234d0

  ! Pt (S-1)
  n0s(15)    = 0.819d0
  n0p(15)    = 0.869d0
  n0d(15)    = 8.328d0
  cs_bar(15) =-0.50976104d0
  cp_bar(15) = 0.23571278d0
  cd_bar(15) =-0.29429756d0
  ds_bar(15) = 0.37954151d0
  dp_bar(15) = 0.23937890d0
  dd_bar(15) = 0.15128375d0
  cs_ort(15) =-0.51961616d0
  cp_ort(15) = 0.60434768d0
  cd_ort(15) =-0.30485581d0
  ds_ort(15) = 0.38695479d0
  dp_ort(15) = 0.39646105d0
  dd_ort(15) = 0.15132097d0

  ! Pt (S-2)
  n0s(16)    = 0.798d0
  n0p(16)    = 0.903d0
  n0d(16)    = 8.293d0
  cs_bar(16) =-0.52359992d0
  cp_bar(16) = 0.19479415d0
  cd_bar(16) =-0.30413744d0
  ds_bar(16) = 0.37285715d0
  dp_bar(16) = 0.23361646d0
  dd_bar(16) = 0.15151997d0
  cs_ort(16) =-0.52096271d0
  cp_ort(16) = 0.60572419d0
  cd_ort(16) =-0.30416645d0
  ds_ort(16) = 0.38731090d0
  dp_ort(16) = 0.39770293d0
  dd_ort(16) = 0.15151201d0

  ! Pt (S-3)
  n0s(17)    = 0.786d0
  n0p(17)    = 0.912d0
  n0d(17)    = 8.297d0
  cs_bar(17) =-0.52276751d0
  cp_bar(17) = 0.19382565d0
  cd_bar(17) =-0.30216556d0
  ds_bar(17) = 0.37189534d0
  dp_bar(17) = 0.23320324d0
  dd_bar(17) = 0.15159747d0
  cs_ort(17) =-0.52043991d0
  cp_ort(17) = 0.60646170d0
  cd_ort(17) =-0.30289365d0
  ds_ort(17) = 0.38738172d0
  dp_ort(17) = 0.39781407d0
  dd_ort(17) = 0.15158493d0

  ! Pt (bulk)
  n0s(18)    = 0.399d0*2.d0
  n0p(18)    = 0.441d0*2.d0
  n0d(18)    = 4.160d0*2.d0
  cs_bar(18) =-0.523644d0
  cp_bar(18) = 0.193450d0
  cd_bar(18) =-0.303302d0
  ds_bar(18) = 0.371978d0
  dp_bar(18) = 0.233302d0
  dd_bar(18) = 0.151575d0
  cs_ort(18) =-0.52062844d0
  cp_ort(18) = 0.60622535d0
  cd_ort(18) =-0.30330401d0
  ds_ort(18) = 0.38736493d0
  dp_ort(18) = 0.39778479d0
  dd_ort(18) = 0.15156316d0


! Angela's unrelaxed calculations for Pt(100)
! no SOC, paramagnetic, with hoh
!
! Material indices:
!                   19 - Empty spheres
!                   20 - Pt (S)
!                   21 - Pt (S-1)
!                   22 - Pt (S-2)
!                   23 - Pt (S-3)
!                   24 - Pt (blk)
!
! Experimental lattice parameter of Pt a= 3.92 \AA
! Fermi Level: E_F= -0.05529 Ry

  layertype_dft(19) = 1
  layertype_dft(24) = 3

  ! Fermi energy
  Ef_dft(19:24)   = -0.05529d0

  ! Effective intra-site coulomb interaction U (in eV)
  U_dft(20:24) = 0.6d0

  ! SOC parameters
  if(SOC) then
    lambdasoc_dft(20:24) = 0.044d0  ! Pt
  end if

  ! Empty spheres
  n0s(19)    = 0.0485d0*2.d0
  n0p(19)    = 0.0579d0*2.d0
  n0d(19)    = 0.0329d0*2.d0
  cs_bar(19) = 0.222360d0
  cp_bar(19) = 0.280346d0
  cd_bar(19) = 0.035879d0
  ds_bar(19) = 0.219724d0
  dp_bar(19) = 0.114204d0
  dd_bar(19) = 0.028745d0
  cs_ort(19) = 0.00040176d0
  cp_ort(19) = 0.97478705d0
  cd_ort(19) = 2.41243775d0
  ds_ort(19) = 0.38359822d0
  dp_ort(19) = 0.40527728d0
  dd_ort(19) = 0.43120938d0

  ! Pt (S)
  n0s(20)    = 0.372d0*2.d0
  n0p(20)    = 0.278d0*2.d0
  n0d(20)    = 4.179d0*2.d0
  cs_bar(20) =-0.438911d0
  cp_bar(20) = 0.274259d0
  cd_bar(20) =-0.231783d0
  ds_bar(20) = 0.385410d0
  dp_bar(20) = 0.231370d0
  dd_bar(20) = 0.150482d0
  cs_ort(20) =-0.51471897d0
  cp_ort(20) = 0.61410754d0
  cd_ort(20) =-0.30759283d0
  ds_ort(20) = 0.38519728d0
  dp_ort(20) = 0.39632741d0
  dd_ort(20) = 0.15046392d0

  ! Pt (S-1)
  n0s(21)    = 0.407d0*2.d0
  n0p(21)    = 0.446d0*2.d0
  n0d(21)    = 4.184d0*2.d0
  cs_bar(21) =-0.539709d0
  cp_bar(21) = 0.192392d0
  cd_bar(21) =-0.310160d0
  ds_bar(21) = 0.375983d0
  dp_bar(21) = 0.237037d0
  dd_bar(21) = 0.152617d0
  cs_ort(21) =-0.51109582d0
  cp_ort(21) = 0.61402257d0
  cd_ort(21) =-0.28340503d0
  ds_ort(21) = 0.38801814d0
  dp_ort(21) = 0.39769974d0
  dd_ort(21) = 0.15257825d0

  ! Pt (S-2)
  n0s(22)    = 0.399d0*2.d0
  n0p(22)    = 0.441d0*2.d0
  n0d(22)    = 4.155d0*2.d0
  cs_bar(22) =-0.523213d0
  cp_bar(22) = 0.191546d0
  cd_bar(22) =-0.304713d0
  ds_bar(22) = 0.371251d0
  dp_bar(22) = 0.232710d0
  dd_bar(22) = 0.151367d0
  cs_ort(22) =-0.52265313d0
  cp_ort(22) = 0.60446527d0
  cd_ort(22) =-0.30741084d0
  ds_ort(22) = 0.38724920d0
  dp_ort(22) = 0.39778864d0
  dd_ort(22) = 0.15136209d0

  ! Pt (S-3)
  n0s(23)    = 0.399d0*2.d0
  n0p(23)    = 0.440d0*2.d0
  n0d(23)    = 4.161d0*2.d0
  cs_bar(23) =-0.523688d0
  cp_bar(23) = 0.193381d0
  cd_bar(23) =-0.302989d0
  ds_bar(23) = 0.372072d0
  dp_bar(23) = 0.233317d0
  dd_bar(23) = 0.151614d0
  cs_ort(23) =-0.52025695d0
  cp_ort(23) = 0.60659157d0
  cd_ort(23) =-0.30254739d0
  ds_ort(23) = 0.38738733d0
  dp_ort(23) = 0.39780109d0
  dd_ort(23) = 0.15160067d0

  ! Pt (bulk)
  n0s(24)    = 0.399d0*2.d0
  n0p(24)    = 0.441d0*2.d0
  n0d(24)    = 4.160d0*2.d0
  cs_bar(24) =-0.523644d0
  cp_bar(24) = 0.193450d0
  cd_bar(24) =-0.303302d0
  ds_bar(24) = 0.371978d0
  dp_bar(24) = 0.233302d0
  dd_bar(24) = 0.151575d0
  cs_ort(24) =-0.52062844d0
  cp_ort(24) = 0.60622535d0
  cd_ort(24) =-0.30330401d0
  ds_ort(24) = 0.38736493d0
  dp_ort(24) = 0.39778479d0
  dd_ort(24) = 0.15156316d0

!-------------------------- END OF PARAMETERS --------------------------
  ! Fermi level is obtained for layer 1
  Ef = Ef_dft(mmlayer(1))*ry2ev

  do j=1,Npl+2
    i = mmlayer(j)
    ! Shift in the Fermi level
    shft = Ef - Ef_dft(i)*ry2ev
    ! Layer type: 1 - Empty Sphere; 2 - Magnetic ; 0 - Other
    layertype(j) = layertype_dft(i)
    ! Effective intra-site coulomb interaction U (in Ryd)
    select case (Utype)
    case(0)
      U(j) = 0.d0
    case(1)
      if(layertype(j).eq.2) then
        U(j) = U_dft(i)*ry2ev/13.6d0
      else
        U(j) = 0.d0
      end if
    case(2)
      U(j) = U_dft(i)*ry2ev/13.6d0
    end select
    ! SOC parameters
    lambda(j) = lambdasoc_dft(i)*ry2ev
    ! Occupations
    npart0(j) = n0s(i)+n0p(i)+n0d(i)
    ! C's and Delta's for hopping calculations
    dft_type: select case (dfttype)
    case ("T")
      if((cs_bar(i)+cp_bar(i)+cd_bar(i)+ds_bar(i)+dp_bar(i)+dd_bar(i)).gt.900.d0) then
        if(myrank.eq.0) write(outputunit,"('[DFT_parameters] Missing Tight-binding basis parameters!')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      cs(j) = cs_bar(i) + shft
      cp(j) = cp_bar(i) + shft
      cd(j) = cd_bar(i) + shft
      ds(j) = ds_bar(i)
      dp(j) = dp_bar(i)
      dd(j) = dd_bar(i)
    case ("O")
      if((cs_ort(i)+cp_ort(i)+cd_ort(i)+ds_ort(i)+dp_ort(i)+dd_ort(i)).gt.900.d0) then
        if(myrank.eq.0) write(outputunit,"('[DFT_parameters] Missing Orthogonal basis parameters!')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      cs(j) = cs_ort(i) + shft
      cp(j) = cp_ort(i) + shft
      cd(j) = cd_ort(i) + shft
      ds(j) = ds_ort(i)
      dp(j) = dp_ort(i)
      dd(j) = dd_ort(i)
    case default
      if(myrank.eq.0) write(outputunit,"('[DFT_parameters] Choose between (T)ight-binding or (O)rthogonal DFT parameters!')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end select dft_type
  end do

  cs = cs*ry2ev
  cp = cp*ry2ev
  cd = cd*ry2ev
  ds = ds*sqrt(ry2ev)
  dp = dp*sqrt(ry2ev)
  dd = dd*sqrt(ry2ev)

  return
  end subroutine DFT_parameters

  subroutine rs_hoppings()
    use mod_f90_kind
    use mod_parameters, only: Npl, Utype, mmlayer, nmaglayers, mmlayermag, layertype, outputunit, outputunit_loop
    use mod_lattice
    use mod_mpi_pars
  character(len=30) :: formatvar
  integer       :: i,mu,neighbor
  real(double),dimension(Npl+2) :: cs,cp,cd,ds,dp,dd
  real(double) :: dst,dpt,ddt,w(3),bp(9,9)
  real(double) :: ds2,dsp,dsd,dp2,dpd,dd2
! on-site integrals
  real(double) :: s0,p0,d0t,d0e
! off-site integrals
  real(double) :: sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd

  call DFT_parameters(cs,cp,cd,ds,dp,dd)
  if((layertype(1).ne.1).or.(layertype(Npl+2).ne.1)) then
    if(myrank.eq.0) then
      write(outputunit,"('[mod_tight_binding] Problem defining the system: first or last layer is not Empty Spheres!')")
      write(outputunit,"('[mod_tight_binding] First layer type: ',i0,' , last layer type: ',i0,'.')") layertype(1),layertype(Npl+2)
    end if
    call MPI_Finalize(ierr)
    stop
  end if
  write(formatvar,fmt="(a,i0,a)") '(a,',Npl+2,'(i0,2x))'
  if(myrank_row_hw.eq.0) write(outputunit_loop,fmt=formatvar) '[rs_hoppings] Layer type: ',(mmlayer(i),i=1,Npl+2)

  ! Obtaining the number and list of magnetic layers
  nmaglayers = 0
  do i=1,Npl+2
    if(layertype(i).eq.2) then
      nmaglayers = nmaglayers + 1
      mmlayermag(nmaglayers) = i
    end if
  end do
  if(myrank_row_hw.eq.0) then
    if(nmaglayers.ge.1) then
      write(formatvar,fmt="(a,i0,a)") '(a,i0,a,',nmaglayers,'(i0,2x))'
      write(outputunit_loop,fmt=formatvar) '[rs_hoppings] ',nmaglayers,' magnetic layer(s) in the system: ',(mmlayermag(i),i=1,nmaglayers)
    else
      write(outputunit_loop,fmt="('[rs_hoppings] No magnetic layer in the system. ')")
    end if
  end if
  ! Unifying Utype = 0 and 1 To Utype = 0 if there is no magnetic layer
  if((Utype.le.1).and.(nmaglayers.eq.0)) Utype = 0

  ! Allocating real space hoppings
  inter_plane_hoppings: select case (plnn)
  case(1) ! If inter-plane second nearest neighbors are in n.n. plane
    allocate( t00(Npl+2,0:n0,9,9),t01(Npl+1,n1+n2,9,9) )
  case(2) ! If inter-plane second nearest neighbors are in 2nd. n.n. plane
    allocate( t00(Npl+2,0:n0,9,9),t01(Npl+1,n1,9,9),t02(Npl,n2,9,9) )
  case default
    if(myrank.eq.0) write(outputunit,"('[rs_hoppings] System not defined for more than 2 n.n. planes!')")
    call MPI_Finalize(ierr)
    stop
  end select inter_plane_hoppings

  t00 = 0.d0
! In plane hoppings
  do i=1,Npl+2

    ds2 = ds(i)*ds(i)
    dsp = ds(i)*dp(i)
    dsd = ds(i)*dd(i)
    dp2 = dp(i)*dp(i)
    dpd = dp(i)*dd(i)
    dd2 = dd(i)*dd(i)

    ! on site
    s0  = cs(i) + ds2*3.09d0
    p0  = cp(i) + dp2*2.79d0
    d0t = cd(i) + dd2*2.71d0
    d0e = cd(i) + dd2*1.30d0

    t00(i,0,1,1)  = s0
    do mu=2,4
      t00(i,0,mu,mu) = p0
    end do
    do mu=5,7
      t00(i,0,mu,mu) = d0t
    end do
    do mu=8,9
      t00(i,0,mu,mu) = d0e
    end do

    ! first n.n.
    sss = -ds2*0.593d0
    sps =  dsp*1.18d0
    pps =  dp2*2.36d0
    ppp = -dp2*0.36d0
    sds = -dsd*1.42d0
    pds = -dpd*2.93d0
    pdp =  dpd*0.82d0
    dds = -dd2*3.84d0
    ddp =  dd2*1.85d0
    ddd = -dd2*0.19d0
    do neighbor=1,n01
      w = c0(neighbor,:)
      call intd(sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd,w,bp)
      t00(i,neighbor,:,:) = bp
    end do

    ! second n.n.
    sss = -ds2*0.203d0
    sps =  dsp*0.44d0
    pps =  dp2*0.93d0
    ppp = -dp2*0.05d0
    sds = -dsd*0.60d0
    pds = -dpd*1.29d0
    pdp =  dpd*0.13d0
    dds = -dd2*1.76d0
    ddp =  dd2*0.36d0
    ddd = -dd2*0.02d0
    do neighbor=n01+1,n0
      w = c0(neighbor,:)
      call intd(sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd,w,bp)
      t00(i,neighbor,:,:) = bp
    end do
  end do

  t01 = 0.d0
! Inter plane hoppings
  do i=1,Npl+1

    dst = sqrt(ds(i)*ds(i+1))
    dpt = sqrt(dp(i)*dp(i+1))
    ddt = sqrt(dd(i)*dd(i+1))

    ds2 = dst*dst
    dsp = dst*dpt
    dsd = dst*ddt
    dp2 = dpt*dpt
    dpd = dpt*ddt
    dd2 = ddt*ddt

    ! first n.n.
    sss = -ds2*0.593d0
    sps =  dsp*1.18d0
    pps =  dp2*2.36d0
    ppp = -dp2*0.36d0
    sds = -dsd*1.42d0
    pds = -dpd*2.93d0
    pdp =  dpd*0.82d0
    dds = -dd2*3.84d0
    ddp =  dd2*1.85d0
    ddd = -dd2*0.19d0
    do neighbor=1,n1
      w = c1(neighbor,:)
      call intd(sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd,w,bp)
      t01(i,neighbor,:,:) = bp
    end do
  end do

  second_neighbor_inter_plane_hopping: select case (plnn)
  case(1) ! If inter-plane second nearest neighbors are in n.n. plane
    do i=1,Npl+1

      dst = sqrt(ds(i)*ds(i+1))
      dpt = sqrt(dp(i)*dp(i+1))
      ddt = sqrt(dd(i)*dd(i+1))

      ds2 = dst*dst
      dsp = dst*dpt
      dsd = dst*ddt
      dp2 = dpt*dpt
      dpd = dpt*ddt
      dd2 = ddt*ddt

      ! second n.n.
      sss = -ds2*0.203d0
      sps =  dsp*0.44d0
      pps =  dp2*0.93d0
      ppp = -dp2*0.05d0
      sds = -dsd*0.60d0
      pds = -dpd*1.29d0
      pdp =  dpd*0.13d0
      dds = -dd2*1.76d0
      ddp =  dd2*0.36d0
      ddd = -dd2*0.02d0
      do neighbor=1,n2
        w = c2(neighbor,:)
        call intd(sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd,w,bp)
        t01(i,n1+neighbor,:,:) = bp
      end do
    end do
  case(2) ! If inter-plane second nearest neighbors are in 2nd. n.n. plane
    t02 = 0.d0
    do i=1,Npl

      dst = sqrt(ds(i)*ds(i+2))
      dpt = sqrt(dp(i)*dp(i+2))
      ddt = sqrt(dd(i)*dd(i+2))

      ds2 = dst*dst
      dsp = dst*dpt
      dsd = dst*ddt
      dp2 = dpt*dpt
      dpd = dpt*ddt
      dd2 = ddt*ddt

      ! second n.n.
      sss = -ds2*0.203d0
      sps =  dsp*0.44d0
      pps =  dp2*0.93d0
      ppp = -dp2*0.05d0
      sds = -dsd*0.60d0
      pds = -dpd*1.29d0
      pdp =  dpd*0.13d0
      dds = -dd2*1.76d0
      ddp =  dd2*0.36d0
      ddd = -dd2*0.02d0
      do neighbor=1,n2
        w = c2(neighbor,:)
        call intd(sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd,w,bp)
        t02(i,neighbor,:,:) = bp
      end do
    end do
  end select second_neighbor_inter_plane_hopping

  return
  end subroutine rs_hoppings

  pure subroutine intd(sss,sps,pps,ppp,ss,ps,pp,ds,dp,dd,w,b)
    use mod_f90_kind
    use mod_constants
    implicit none
  real(double), intent(in)  :: sss,sps,pps,ppp,ss,ps,pp,ds,dp,dd,w(3)
  real(double), intent(out) :: b(9,9)

  real(double)  :: x,y,z,xx,xy,yy,yz,zz,zx,xxyy,yyzz,zzxx
  real(double)  :: aux,aux1,aux2,aux3,aux4,r3,f1,f2,f3,f4,f5,f8,g1,g2

  x=w(1)
  y=w(2)
  z=w(3)
  xx=x*x
  xy=x*y
  yy=y*y
  yz=y*z
  zz=z*z
  zx=z*x
  xxyy=xx*yy
  yyzz=yy*zz
  zzxx=zz*xx
  aux=pps-ppp
  r3=sq3
  aux1=r3*ss
  f8=3.d0*zz-1.d0
  f1=xx+yy
  f2=xx-yy
  f3=zz-.5d0*f1
  g1=1.5d0*f2*ds
  g2=r3*f3*ds
  b(1,1)=sss
  b(1,2)=x*sps
  b(1,3)=y*sps
  b(1,4)=z*sps
  b(2,1)=-b(1,2)
  b(2,2)=xx*pps+(1.d0-xx)*ppp
  b(2,3)=xy*aux
  b(2,4)=zx*aux
  b(3,1)=-b(1,3)
  b(3,2)=b(2,3)
  b(3,3)=yy*pps+(1.d0-yy)*ppp
  b(3,4)=yz*aux
  b(4,1)=-b(1,4)
  b(4,2)=b(2,4)
  b(4,3)=b(3,4)
  b(4,4)=zz*pps+(1.d0-zz)*ppp
  b(1,5)=xy*aux1
  b(1,6)=yz*aux1
  b(1,7)=zx*aux1
  b(1,8)=.5d0*f2*aux1
  b(1,9)=.5d0*f8*ss
  b(5,1)=b(1,5)
  b(6,1)=b(1,6)
  b(7,1)=b(1,7)
  b(8,1)=b(1,8)
  b(9,1)=b(1,9)
  f4=.5d0*r3*f2*ps
  f5=.5d0*f8*ps
  aux2=r3*xx*ps+(1.d0-2.d0*xx)*pp
  b(2,5)=aux2*y
  b(2,6)=(r3*ps-2.d0*pp)*xy*z
  b(2,7)=aux2*z
  b(2,8)=(f4+(1.d0-f2)*pp)*x
  b(2,9)=(f5-r3*zz*pp)*x
  aux3=(r3*yy*ps+(1.d0-2.d0*yy)*pp)
  b(3,5)=aux3*x
  b(3,6)=aux3*z
  b(3,7)=b(2,6)
  b(3,8)=(f4-(1.d0+f2)*pp)*y
  b(3,9)=(f5-r3*zz*pp)*y
  aux4=r3*zz*ps+(1.d0-2.d0*zz)*pp
  b(4,5)=b(2,6)
  b(4,6)=aux4*y
  b(4,7)=aux4*x
  b(4,8)=(f4-f2*pp)*z
  b(4,9)=(f5+r3*f1*pp)*z
  b(5,2)=-b(2,5)
  b(6,2)=-b(2,6)
  b(7,2)=-b(2,7)
  b(8,2)=-b(2,8)
  b(9,2)=-b(2,9)
  b(5,3)=-b(3,5)
  b(6,3)=-b(3,6)
  b(7,3)=-b(3,7)
  b(8,3)=-b(3,8)
  b(9,3)=-b(3,9)
  b(5,4)=-b(4,5)
  b(6,4)=-b(4,6)
  b(7,4)=-b(4,7)
  b(8,4)=-b(4,8)
  b(9,4)=-b(4,9)
  b(5,5)=3.d0*xxyy*ds+(f1-4.d0*xxyy)*dp+(zz+xxyy)*dd
  b(5,6)=(3.d0*yy*ds+(1.d0-4.d0*yy)*dp+(yy-1.d0)*dd)*zx
  b(5,7)=(3.d0*xx*ds+(1.d0-4.d0*xx)*dp+(xx-1.d0)*dd)*yz
  b(5,8)=(g1-2.d0*f2*dp+.5d0*f2*dd)*xy
  b(5,9)=(g2-2.d0*r3*zz*dp+.5d0*r3*(1.d0+zz)*dd)*xy
  b(6,5)=b(5,6)
  b(6,6)=3.d0*yyzz*ds+(yy+zz-4.d0*yyzz)*dp+(xx+yyzz)*dd
  b(6,7)=(3.d0*zz*ds+(1.d0-4.d0*zz)*dp+(zz-1.d0)*dd)*xy
  b(6,8)=(g1-(1.d0+2.d0*f2)*dp+(1.d0+.5d0*f2)*dd)*yz
  b(6,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*yz
  b(7,5)=b(5,7)
  b(7,6)=b(6,7)
  b(7,7)=3.d0*zzxx*ds+(zz+xx-4.d0*zzxx)*dp+(yy+zzxx)*dd
  b(7,8)=(g1+(1.d0-2.d0*f2)*dp-(1.d0-.5d0*f2)*dd)*zx
  b(7,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*zx
  b(8,5)=b(5,8)
  b(8,6)=b(6,8)
  b(8,7)=b(7,8)
  b(8,8)=.75d0*f2*f2*ds+(f1-f2*f2)*dp+(zz+.25d0*f2*f2)*dd
  b(8,9)=.5d0*f2*g2-r3*zz*f2*dp+.25d0*r3*(1.d0+zz)*f2*dd
  b(9,5)=b(5,9)
  b(9,6)=b(6,9)
  b(9,7)=b(7,9)
  b(9,8)=b(8,9)
  b(9,9)=f3*f3*ds+3.d0*zz*f1*dp+.75d0*f1*f1*dd

  return
  end subroutine intd

end module mod_tight_binding
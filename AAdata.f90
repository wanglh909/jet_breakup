module data
  use kind
  implicit none



  real(kind=rk):: R, H, Hmin
  real(kind=rk):: Re, Oh, Grav, Kdi = 1.0_rk     !Kdi: k, used before p as the dimensionless group

  integer(kind=ik):: NEX, NEY, NNX, NTN, NTE, NVar, iBW       !NVar: Number of Variables   !NNX: number of nodes in r (abscissa)
  integer(kind=ik):: NNVar = 5    !NNVar: the maximum number of local variables (r,z,u,v,p --> 5)

  integer(kind=ik), parameter:: Nr = 0, Nz = 1, Nu = 2, Nv = 3, Np = 4

  integer(kind=ik):: timestep, step

  real(kind=rk):: error2

  integer(kind=ik), allocatable:: globalNM(:,:), rowNM(:), columnNM(:)
  integer(kind=ik), allocatable:: MDF(:), NOPP(:)

  real(kind=rk):: time
  real(kind=rk), allocatable:: rcoordinate(:), zcoordinate(:), usol(:), vsol(:), psol(:)
  !solp & dtp means previous solution & dt, solpred means the predicted solution
  real(kind=rk), allocatable:: sol(:), solp(:), soldot(:), soldotp(:), soldotpp(:), solpred(:), dsol(:)

  real(kind=rk), allocatable:: fsi_size(:), geta_size(:)

  !for subroutine basis_function
  real(kind=rk):: phi(3,3,9), phisi(3,3,9), phieta(3,3,9), psi(3,3,9), psisi(3,3,9), psieta(3,3,9), &
       phi_1d(3,3), phix_1d(3,3)
  integer(kind=ik),parameter:: convert49(4) = (/1, 3, 7, 9/)


  integer, parameter:: Ng = 3 !number of gausspoints


  !for subroutine prediction
  real(kind=rk):: dt, dtp, CTJ, change, trunerr
  real(kind=rk), parameter:: eps = 1.0e-3_rk

  !for tao_calculation
  real(kind=rk):: Hminlast, Hmin2last, tlast, t2last

  !for frontal solver  
  integer(kind=ik), allocatable:: rNOP(:,:,:)
  integer(kind=ik):: s_mode = 0, ths = 1, bas = 9
  !s_mode: Set to 0 if solving full dynamics and 1 if only mesh (must be 0 if solve mesh by putting BC's on all velocity and pressure vars)
  !bas: number of nodes in element (9)
  !ths: Number of threads, parameter in the kind file (must make clean(wipe)  when changing)

  integer(kind=ik):: RN = 1   !number of regions



end module data

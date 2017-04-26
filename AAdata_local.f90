module Ldata
  use kind
  implicit none


  !for subroutine values_in_an_element
  real(kind=rk):: rlocal(9), zlocal(9), ulocal(9), vlocal(9), plocal(9)
  integer(kind=ik), parameter:: NNr = 1, NNz = 2, NNu = 3, NNv = 4, NNp = 5
  real(kind=rk):: rintfac(3,3), rsi(3,3), reta(3,3), zsi(3,3), zeta(3,3)
  real(kind=rk):: Jp(3,3), s_orth(3,3)
  real(kind=rk):: phir(3,3,9), phiz(3,3,9)
  real(kind=rk):: uintfac(3,3), urintfac(3,3), uzintfac(3,3), &
       vintfac(3,3), vrintfac(3,3), vzintfac(3,3), pintfac(3,3)
  real(kind=rk):: rdotintfac(3,3), zdotintfac(3,3), udotintfac(3,3), vdotintfac(3,3)
  real(kind=rk):: udotlocal(9), vdotlocal(9), rdotlocal(9), zdotlocal(9)
  real(kind=rk):: rsi_down(3),  zsi_down(3), rsi_up(3), zsi_up(3), reta_left(3), zeta_left(3), &
       reta_right(3), zeta_right(3)
  real(kind=rk):: rintfac_right(3), uintfac_right(3), vintfac_right(3), &
       rdotintfac_right(3), zdotintfac_right(3)


  !for subroutine define_sf
  real(kind=rk):: Aterm(3,3), Bterm(3,3)
  ! real(kind=rk), allocatable:: sf(:,:)
  real(kind=rk):: M1 = 0.0_rk, M2 = 0.0_rk
  real(kind=rk), parameter:: eps1 = 1.0_rk, eps2 = 1.0_rk, epss = 0.1_rk
  


  !for subroutine values_in_sj
  real(kind=rk):: s_orth_r(3,3), s_orth_z(3,3)
  real(kind=rk):: Aterm_r(3,3),Aterm_z(3,3), Bterm_r(3,3), Bterm_z(3,3)
  real(kind=rk):: Jp_r(3,3), Jp_z(3,3), rJp_r(3,3)
  real(kind=rk)::  phir_r(3,3,9), phir_z(3,3,9), phiz_r(3,3,9), phiz_z(3,3,9)
  real(kind=rk):: urintfac_r(3,3), urintfac_z(3,3), uzintfac_r(3,3), uzintfac_z(3,3), &
       vrintfac_r(3,3), vrintfac_z(3,3), vzintfac_r(3,3), vzintfac_z(3,3)




  !for subroutine VI_in_sj
  ! real(kind=rk), allocatable:: sj(:,:,:,:)


  !for subroutine SVI_in_sj



end module Ldata

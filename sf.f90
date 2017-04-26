
subroutine define_sf(m,i, sf, LNVar, LNOPP)
  use kind
  use data
  use Ldata
  use NOP_mod, only: gaussian_quadrature, gaussian_quadrature_1d

  implicit none

  integer(kind=ik), intent(in):: m, i, LNVar, LNOPP(9)
  real(kind=rk), intent(out):: sf(LNVar)

  integer(kind=ik):: k, l, ipp   !no i (i is the i in main program)
  real(kind=rk):: integrandRsi_V(Ng,Ng), integrandReta_V(Ng,Ng), integrandRu_V(Ng,Ng), integrandRv_V(Ng,Ng), integrandRp(Ng,Ng)
  real(kind=rk):: integrandRsi_S(Ng), integrandReta_S(Ng), integrandRu_S(Ng), integrandRv_S(Ng)


  do k = 1, Ng, 1
     do l = 1, Ng, 1

        !Aterm(k,l,i)
        Aterm(k,l) = ( zeta(k,l)**2 + reta(k,l)**2 )*phisi(k,l,i) - &
             ( zeta(k,l)*zsi(k,l) + reta(k,l)*rsi(k,l) )*phieta(k,l,i)
        Bterm(k,l) = -( zsi(k,l)*zeta(k,l) + rsi(k,l)*reta(k,l) )*phisi(k,l,i) &
             + ( zsi(k,l)**2 + rsi(k,l)**2 )*phieta(k,l,i)

        integrandRsi_V(k,l) = ( s_orth(k,l) + epss )*Aterm(k,l)/Jp(k,l) - &
             eps1*phisi(k,l,i)*fsi_size(m)*log( rsi(k,l)**2 + zsi(k,l)**2 )

        integrandReta_V(k,l) = ( 1.0_rk/s_orth(k,l) + epss )*Bterm(k,l)/Jp(k,l) - &
             eps2*phieta(k,l,i)*geta_size(m)*log( reta(k,l)**2 + zeta(k,l)**2 )


        integrandRu_V(k,l) = ( Re*phi(k,l,i)*( udotintfac(k,l) + ( uintfac(k,l) - rdotintfac(k,l) )*urintfac(k,l) +  &
             ( vintfac(k,l) - zdotintfac(k,l) )*uzintfac(k,l) )  - &
             Kdi*pintfac(k,l)*Oh*( phir(k,l,i) + phi(k,l,i)/rintfac(k,l) ) + &
             Oh*( 2.0_rk*urintfac(k,l)*phir(k,l,i) + ( uzintfac(k,l) + vrintfac(k,l) )*phiz(k,l,i) + &
             phi(k,l,i)*2.0_rk/rintfac(k,l)**2 *uintfac(k,l) ) )  *rintfac(k,l)*abs(Jp(k,l))

        integrandRv_V(k,l) = ( Re*phi(k,l,i)*( vdotintfac(k,l) + ( uintfac(k,l) - rdotintfac(k,l) )*vrintfac(k,l) + &
             ( vintfac(k,l) - zdotintfac(k,l) )*vzintfac(k,l) + Grav ) - Kdi*pintfac(k,l)*Oh*phiz(k,l,i) + &
             Oh*( ( uzintfac(k,l) + vrintfac(k,l) ) *phir(k,l,i) + 2.0_rk*vzintfac(k,l)*phiz(k,l,i) ) ) &
             *rintfac(k,l)*abs(Jp(k,l))

        integrandRp(k,l) = psi(k,l,i)* ( urintfac(k,l) + uintfac(k,l)/rintfac(k,l) + vzintfac(k,l) ) &
             *rintfac(k,l)*abs(Jp(k,l))
     end do
  end do

  sf(LNOPP(i)+Nr) = gaussian_quadrature(integrandRsi_V)
  sf(LNOPP(i)+Nz) = gaussian_quadrature(integrandReta_V)
  sf(LNOPP(i)+Nu) = gaussian_quadrature(integrandRu_V)
  sf(LNOPP(i)+Nv) = gaussian_quadrature(integrandRv_V)
  if( MDF( globalNM(m,i) ).eq.NNp ) then
     sf(LNOPP(i)+Np) = gaussian_quadrature(integrandRp)
  end if

  ! sf(i,NNr) = gaussian_quadrature(integrandRsi_V)
  ! sf(i,NNz) = gaussian_quadrature(integrandReta_V)
  ! sf(i,NNu) = gaussian_quadrature(integrandRu_V)
  ! sf(i,NNv) = gaussian_quadrature(integrandRv_V)
  ! sf(i,NNp) = gaussian_quadrature(integrandRp)    !sf( (2,4,5,6,8) ,5) = 0.0_rk


  !*******************************************adding SI to sf********************************************************

  if( (rowNM(m).eq.1) .and. ( (i.eq.1) .or. (i.eq.2) .or. (i.eq.3) ) )  then
     do k = 1, Ng, 1    !three gausspoints
        ipp = i  !phix_1d(k,ipp)
        integrandRsi_S(k) = phix_1d(k,ipp) * fsi_size(m) * log( rsi_down(k)**2 + zsi_down(k)**2 )
     end do
     sf(LNOPP(i)+Nr) = sf(LNOPP(i)+Nr)  -  M1*gaussian_quadrature_1d(integrandRsi_S)
  end if

  if( (rowNM(m).eq.NEY) .and. ( (i.eq.7) .or. (i.eq.8) .or. (i.eq.9) ) )  then
     do k = 1, Ng, 1    !three gausspoints
        ipp = i - 6  !phix_1d(k,ipp)
        integrandRsi_S(k) = phix_1d(k,ipp) * fsi_size(m) * log( rsi_up(k)**2 + zsi_up(k)**2 )
     end do
     sf(LNOPP(i)+Nr) = sf(LNOPP(i)+Nr)  -  M1*gaussian_quadrature_1d(integrandRsi_S)
  end if

  if( (columnNM(m).eq.1) .and. ( (i.eq.1) .or. (i.eq.4) .or. (i.eq.7) ) )  then
     do k = 1, Ng, 1    !three gausspoints
        ipp = i/3 + 1  !phix_1d(k,ipp)
        integrandReta_S(k) = phix_1d(k,ipp) * geta_size(m) * log( reta_left(k)**2 + zeta_left(k)**2 )
     end do
     sf(LNOPP(i)+Nz) = sf(LNOPP(i)+Nz)  -  M2*gaussian_quadrature_1d(integrandReta_S)
  end if

  if( (columnNM(m).eq.NEX) .and. ( (i.eq.3) .or. (i.eq.6) .or. (i.eq.9) ) )  then
     do k = 1, Ng, 1    !three gausspoints
        ipp = i/3  !phix_1d(k,ipp)
        integrandReta_S(k) = phix_1d(k,ipp) * geta_size(m) * log( reta_right(k)**2 + zeta_right(k)**2 )


        integrandRsi_S(k) = phi_1d(k,ipp)*( -zeta_right(k)*( uintfac_right(k) - rdotintfac_right(k) ) + &
             reta_right(k)*( vintfac_right(k) - zdotintfac_right(k) ) )*rintfac_right(k)

        integrandRu_S(k) = ( reta_right(k)*phix_1d(k,ipp) / ( reta_right(k)**2 + zeta_right(k)**2 ) + &
             phi_1d(k,ipp)/rintfac_right(k) )*rintfac_right(k)*( reta_right(k)**2 + zeta_right(k)**2 )**0.5_rk

        integrandRv_S(k) = zeta_right(k)*phix_1d(k,ipp) / ( reta_right(k)**2 + zeta_right(k)**2 )**0.5_rk &
             *rintfac_right(k)

     end do

     sf(LNOPP(i)+Nz) = sf(LNOPP(i)+Nz)  -  M2*gaussian_quadrature_1d(integrandReta_S)

     sf(LNOPP(i)+Nr) = gaussian_quadrature_1d(integrandRsi_S)          !directly replace Rsi_V with Rsi_S, viz KBC
     sf(LNOPP(i)+Nu) = sf(LNOPP(i)+Nu)  + gaussian_quadrature_1d(integrandRu_S)
     sf(LNOPP(i)+Nv) = sf(LNOPP(i)+Nv)  + gaussian_quadrature_1d(integrandRv_S)
  end if

  !    sf(i,NNz) = sf(i,NNz)  -  M2*gaussian_quadrature_1d(integrandReta_S)

  !    sf(i,NNr) = gaussian_quadrature_1d(integrandRsi_S)          !directly replace Rsi_V with Rsi_S, viz KBC
  !    sf(i,NNu) = sf(i,NNu)  + gaussian_quadrature_1d(integrandRu_S)
  !    sf(i,NNv) = sf(i,NNv)  + gaussian_quadrature_1d(integrandRv_S)
  ! end if

  !*****************************************************************************************************************



end subroutine define_sf


subroutine values_in_an_element(m)
  use kind
  use data
  use Ldata



  implicit none
  integer(kind=ik):: m, j, k, l, n, npp

  !give value to r1-r9 & z1-z9
  !define ulocal(9), vlocal(9), plocal(4)
  do j = 1, 9, 1

     rlocal(j) = sol( NOPP( globalNM(m,j) ) + Nr )
     zlocal(j) = sol( NOPP( globalNM(m,j) ) + Nz )
     ulocal(j) = sol( NOPP( globalNM(m,j) ) + Nu )
     vlocal(j) = sol( NOPP( globalNM(m,j) ) + Nv )
     if( MDF( globalNM(m,j) ) .gt. Np )  plocal(j) = sol( NOPP( globalNM(m,j) ) + Np )   ! plocal(2, 4, 5, 6, 8) remain unchanged

     ! rlocal(j) = rcoordinate( globalNM(m,j) )
     ! zlocal(j) = zcoordinate( globalNM(m,j) )
     ! ulocal(j) = usol( globalNM(m,j) )
     ! vlocal(j) = vsol( globalNM(m,j) )
     ! !if( MDF( globalNM(m,j) ) .gt. Np ) 
     ! plocal(j) = psol( globalNM(m,j) )   ! plocal(2, 4, 5, 6, 8) = 0.0_rk

     rdotlocal(j) = soldot( NOPP( globalNM(m,j) ) + Nr )
     zdotlocal(j) = soldot( NOPP( globalNM(m,j) ) + Nz )
     udotlocal(j) = soldot( NOPP( globalNM(m,j) ) + Nu )
     vdotlocal(j) = soldot( NOPP( globalNM(m,j) ) + Nv )

  end do


  rintfac(:,:) = 0.0_rk
  rsi(:,:) = 0.0_rk
  reta(:,:) = 0.0_rk
  zsi(:,:) = 0.0_rk
  zeta(:,:) = 0.0_rk
  do k = 1, Ng, 1          !relate to a(Ng) (the value of si for gaussian_quadrature)
     do l = 1, Ng, 1           !relate to a(Ng) (the value of eta for gaussian_quadrature)
        !define rsi, reta, zsi, zeta
        do n = 1, 9, 1
           rintfac(k,l) = rintfac(k,l) + rlocal(n)*phi(k,l,n)
           rsi(k,l) = rsi(k,l) + rlocal(n)*phisi(k,l,n)
           reta(k,l) = reta(k,l) + rlocal(n)*phieta(k,l,n)
           zsi(k,l) = zsi(k,l) + zlocal(n)*phisi(k,l,n)
           zeta(k,l) = zeta(k,l) + zlocal(n)*phieta(k,l,n)
        end do
        !define Jp(Ng,Ng)
        Jp(k,l) = rsi(k,l)*zeta(k,l) - reta(k,l)*zsi(k,l)
        s_orth(k,l) = ( ( rsi(k,l)**2 + zsi(k,l)**2 )/( reta(k,l)**2 + zeta(k,l)**2 ) )**0.5_rk

        !define phir, phiz
        do n = 1, 9, 1           !n is the notation of phi (phi1 - phi9)
           phir(k,l,n) = 1.0_rk/Jp(k,l) * &
                ( zeta(k,l)*phisi(k,l,n) - zsi(k,l)*phieta(k,l,n) )
           phiz(k,l,n) = 1.0_rk/Jp(k,l) * &
                ( -reta(k,l)*phisi(k,l,n) + rsi(k,l)*phieta(k,l,n) )
        end do

     end do
  end do



  !define uintfac, urintfac, uzintfac, vintfac, vrintfac, vzintfac
  uintfac(:,:) = 0.0_rk
  urintfac(:,:) = 0.0_rk
  uzintfac(:,:) = 0.0_rk
  vintfac(:,:) = 0.0_rk
  vrintfac(:,:) = 0.0_rk
  vzintfac(:,:) = 0.0_rk
  pintfac(:,:) = 0.0_rk

  rdotintfac(:,:) = 0.0_rk
  zdotintfac(:,:) = 0.0_rk
  udotintfac(:,:) = 0.0_rk
  vdotintfac(:,:) = 0.0_rk

  do k = 1, Ng
     do l = 1, Ng

        do n = 1, 9, 1
           uintfac(k,l) = uintfac(k,l) +  ulocal(n)*phi(k,l,n)
           urintfac(k,l) = urintfac(k,l) +  ulocal(n)*phir(k,l,n)
           uzintfac(k,l) = uzintfac(k,l) +  ulocal(n)*phiz(k,l,n)
           vintfac(k,l) = vintfac(k,l) +  vlocal(n)*phi(k,l,n)
           vrintfac(k,l) = vrintfac(k,l) +  vlocal(n)*phir(k,l,n)
           vzintfac(k,l) = vzintfac(k,l) +  vlocal(n)*phiz(k,l,n)
           pintfac(k,l) =  pintfac(k,l) +  plocal(n)*psi(k,l,n)

           rdotintfac(k,l) = rdotintfac(k,l) + rdotlocal(n)*phi(k,l,n)
           zdotintfac(k,l) = zdotintfac(k,l) + zdotlocal(n)*phi(k,l,n)
           udotintfac(k,l) = udotintfac(k,l) + udotlocal(n)*phi(k,l,n)
           vdotintfac(k,l) = vdotintfac(k,l) + vdotlocal(n)*phi(k,l,n)
        end do

     end do
  end do   !end loop for k&l, preparation for gaussian_quadrature




!!!!!!!!!!!!!!!!!!!!!!!!!!preparation for the SIs, define the temrs in SI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !other than free surface, SIs only exit in elliptic mesh equations
  rsi_down(:) = 0.0_rk
  zsi_down(:) = 0.0_rk
  rsi_up(:) = 0.0_rk
  zsi_up(:) = 0.0_rk
  reta_left(:) = 0.0_rk
  zeta_left(:) = 0.0_rk
  reta_right(:) = 0.0_rk
  zeta_right(:) = 0.0_rk

  rintfac_right(:) = 0.0_rk
  uintfac_right(:) = 0.0_rk
  vintfac_right(:) = 0.0_rk
  rdotintfac_right(:) = 0.0_rk
  zdotintfac_right(:) = 0.0_rk

  do k = 1, Ng, 1  !the three gausspoints for SI

     do n = 1, 3, 1  !the summation of three terms, eg: rsi = sum( rlocal(7,8,9)*phix_1d(1,2,3) )

        if ( rowNM(m).eq.1 ) then     
           !reta, zeta not appear in the SI of Rsi, and the SI of Reta is 0 because d(eta)=0
           !so reta, zeta are not needed
           npp = n  !only need r1, r2, r3 -->r(npp)
           rsi_down(k) = rsi_down(k) + rlocal(npp)*phix_1d(k,n)
           zsi_down(k) = zsi_down(k) + zlocal(npp)*phix_1d(k,n)
        end if

        if ( rowNM(m).eq.NEY ) then
           npp = n + 6  !only need r7, r8, r9 -->r(npp)
           rsi_up(k) = rsi_up(k) + rlocal(npp) * phix_1d(k,n)
           zsi_up(k) = zsi_up(k) + zlocal(npp) * phix_1d(k,n)
        end if

        if ( columnNM(m).eq.1 ) then
           npp = 3*n - 2  !only need r1, r4, r7 -->r(npp)
           reta_left(k) = reta_left(k) + rlocal(npp) * phix_1d(k,n)
           zeta_left(k) = zeta_left(k) + zlocal(npp) * phix_1d(k,n)
        end if

        if ( columnNM(m).eq.NEX ) then
           npp = 3*n  !only need r3, r6, r9 -->r(npp)
           reta_right(k) = reta_right(k) + rlocal(npp) * phix_1d(k,n)
           zeta_right(k) = zeta_right(k) + zlocal(npp) * phix_1d(k,n)

           rintfac_right(k) = rintfac_right(k) + rlocal(npp) * phi_1d(k,n)
           uintfac_right(k) = uintfac_right(k) + ulocal(npp) * phi_1d(k,n)
           vintfac_right(k) = vintfac_right(k) + vlocal(npp) * phi_1d(k,n)
           rdotintfac_right(k) = rdotintfac_right(k) + rdotlocal(npp) * phi_1d(k,n)
           zdotintfac_right(k) = zdotintfac_right(k) + zdotlocal(npp) * phi_1d(k,n)

        end if

     end do  !end for n
  end do    !end for k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine values_in_an_element

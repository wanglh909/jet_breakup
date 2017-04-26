

subroutine values_in_sj(m,i,j)
  use kind
  use data
  use Ldata

  implicit none
  integer(kind=ik), intent(in):: m,i,j
  integer(kind=ik):: k, l, n

  do k = 1, Ng, 1          !relate to a(Ng) (the value of si for gaussian_quadrature)
     do l = 1, Ng, 1           !relate to a(Ng) (the value of eta for gaussian_quadrature)

        !Actually, these terms are (k,l,i,j), but can be refurbished in each cycle of i&j, so just need (k,l)
        !s_orth_r(k,l,i,j)
        s_orth_r(k,l) = 1.0_rk / ( reta(k,l)**2 + zeta(k,l)**2 ) * ( rsi(k,l)*phisi(k,l,j)/s_orth(k,l) &
             - reta(k,l)*phieta(k,l,j)*s_orth(k,l) )

        !s_orth_z(k,l,i,j)
        s_orth_z(k,l) = 1.0_rk / ( reta(k,l)**2 + zeta(k,l)**2 ) * ( zsi(k,l)*phisi(k,l,j)/s_orth(k,l) &
             - zeta(k,l)*phieta(k,l,j)*s_orth(k,l) )

        !Aterm_r(k,l,i,j)
        Aterm_r(k,l) = 2.0_rk*reta(k,l)*phieta(k,l,j)*phisi(k,l,i) &
             - ( rsi(k,l)*phieta(k,l,j) + reta(k,l)*phisi(k,l,j) )*phieta(k,l,i)

        !Aterm_z(k,l,i,j)
        Aterm_z(k,l) = 2.0_rk*zeta(k,l)*phieta(k,l,j)*phisi(k,l,i) &
             - ( zsi(k,l)*phieta(k,l,j) + zeta(k,l)*phisi(k,l,j) )*phieta(k,l,i)

        !Bterm_r(k,l,i,j)
        Bterm_r(k,l) = - ( reta(k,l)*phisi(k,l,j) + rsi(k,l)*phieta(k,l,j) )*phisi(k,l,i) &
             + 2.0_rk*rsi(k,l)*phisi(k,l,j)*phieta(k,l,i) 

        !Bterm_z(k,l,i,j)
        Bterm_z(k,l) = - ( zeta(k,l)*phisi(k,l,j) + zsi(k,l)*phieta(k,l,j) )*phisi(k,l,i) &
             + 2.0_rk*zsi(k,l)*phisi(k,l,j)*phieta(k,l,i) 

        !Jp_r((k,l,i,j)
        Jp_r(k,l) = phisi(k,l,j)*zeta(k,l) - phieta(k,l,j)*zsi(k,l)

        !Jp_r((k,l,i,j)
        Jp_z(k,l) = rsi(k,l)*phieta(k,l,j) - reta(k,l)*phisi(k,l,j)

        !rJp_r(k,l,i,j)
        rJp_r(k,l) = phi(k,l,j)*Jp(k,l) + rintfac(k,l)*Jp_r(k,l)


        do n = 1, 9            !phir_r(k,l,n,j)
           phir_r(k,l,n) = -1.0_rk/Jp(k,l)**2 *Jp_r(k,l)*( phisi(k,l,n)*zeta(k,l) - phieta(k,l,n)*zsi(k,l) )     !d phir(n) / d rj

           phir_z(k,l,n) = -1.0_rk/Jp(k,l)**2 *Jp_z(k,l)*( phisi(k,l,n)*zeta(k,l) - phieta(k,l,n)*zsi(k,l) ) + &
                1.0_rk/Jp(k,l)*( phisi(k,l,n)*phieta(k,l,j) - phieta(k,l,n)*phisi(k,l,j) )

           phiz_r(k,l,n) = -1.0_rk/Jp(k,l)**2 *Jp_r(k,l)*( -phisi(k,l,n)*reta(k,l) + phieta(k,l,n)*rsi(k,l) ) + &
                1.0_rk/Jp(k,l)*( -phisi(k,l,n)*phieta(k,l,j) + phieta(k,l,n)*phisi(k,l,j) )

           phiz_z(k,l,n) = -1.0_rk/Jp(k,l)**2 *Jp_z(k,l)*( -phisi(k,l,n)*reta(k,l) + phieta(k,l,n)*rsi(k,l) )
        end do


        !urintfac_r(k,l,j)
        urintfac_r(k,l) = 0.0_rk
        urintfac_z(k,l) = 0.0_rk
        uzintfac_r(k,l) = 0.0_rk
        uzintfac_z(k,l) = 0.0_rk
        vrintfac_r(k,l) = 0.0_rk
        vrintfac_z(k,l) = 0.0_rk
        vzintfac_r(k,l) = 0.0_rk
        vzintfac_z(k,l) = 0.0_rk
        do n = 1, 9
           urintfac_r(k,l) = urintfac_r(k,l) + ulocal(n)*phir_r(k,l,n)
           urintfac_z(k,l) = urintfac_z(k,l) + ulocal(n)*phir_z(k,l,n)
           uzintfac_r(k,l) = uzintfac_r(k,l) + ulocal(n)*phiz_r(k,l,n)
           uzintfac_z(k,l) = uzintfac_z(k,l) + ulocal(n)*phiz_z(k,l,n)
           vrintfac_r(k,l) = vrintfac_r(k,l) + vlocal(n)*phir_r(k,l,n)
           vrintfac_z(k,l) = vrintfac_z(k,l) + vlocal(n)*phir_z(k,l,n)
           vzintfac_r(k,l) = vzintfac_r(k,l) + vlocal(n)*phiz_r(k,l,n)
           vzintfac_z(k,l) = vzintfac_z(k,l) + vlocal(n)*phiz_z(k,l,n)
        end do



     end do
  end do            !end loop for k&l


end subroutine values_in_sj




subroutine VI_in_sj(m,i,j, sj, LNVar, LNOPP)
  use kind
  use data
  use Ldata
  use NOP_mod, only: gaussian_quadrature

  implicit none

  integer(kind=ik), intent(in):: m,i,j, LNVar, LNOPP(9)
  real(kind=rk), intent(out):: sj(LNVar, LNVar)

  integer(kind=ik):: k, l
  real(kind=rk):: integrandRsi_r_V(Ng,Ng), integrandRsi_z_V(Ng,Ng), integrandReta_r_V(Ng,Ng), integrandReta_z_V(Ng,Ng)
  real(kind=rk):: integrandRu_r_V(Ng,Ng), integrandRu_z_V(Ng,Ng), integrandRu_u(Ng,Ng), integrandRuv(Ng,Ng), integrandRu_p(Ng,Ng)
  real(kind=rk):: integrandRv_r_V(Ng,Ng), integrandRv_z_V(Ng,Ng), integrandRv_u(Ng,Ng), integrandRvv(Ng,Ng),integrandRv_p(Ng,Ng) 
  real(kind=rk):: integrandRp_r(Ng,Ng), integrandRp_z(Ng,Ng), integrandRp_u(Ng,Ng), integrandRp_v(Ng,Ng)


  !define integrand(k,l,p,q):
  do k = 1, Ng, 1          !relate to a(Ng) (the value of si for gaussian_quadrature)
     do l = 1, Ng, 1           !relate to a(Ng) (the value of eta for gaussian_quadrature)

        !integrandRsi_u,v,p = integrandReta_u,v,p = 0.0_rk

        integrandRsi_r_V(k,l) = s_orth_r(k,l) * Aterm(k,l) / Jp(k,l)  +  &
             ( s_orth(k,l) + epss ) * Aterm_r(k,l) / Jp(k,l)  &
             -  ( s_orth(k,l) + epss ) * Aterm(k,l) / Jp(k,l)**2 * Jp_r(k,l)  &
             -  eps1 * phisi(k,l,i) * fsi_size(m) / ( rsi(k,l)**2 + zsi(k,l)**2 ) &
             * 2.0_rk*rsi(k,l) * phisi(k,l,j)

        integrandRsi_z_V(k,l) = s_orth_z(k,l) * Aterm(k,l) / Jp(k,l)  +  &
             ( s_orth(k,l) + epss ) * Aterm_z(k,l) / Jp(k,l)  &
             -  ( s_orth(k,l) + epss ) * Aterm(k,l) / (Jp(k,l)**2) * Jp_z(k,l)  &
             -  eps1 * phisi(k,l,i) * fsi_size(m) / ( rsi(k,l)**2 + zsi(k,l)**2 ) &
             * 2.0_rk*zsi(k,l) * phisi(k,l,j) 

        integrandReta_r_V(k,l) = -1.0_rk/(s_orth(k,l)**2) * s_orth_r(k,l) * Bterm(k,l) / Jp(k,l)  &
             +  ( 1.0_rk/s_orth(k,l) + epss ) * Bterm_r(k,l) / Jp(k,l)  &
             -  ( 1.0_rk/s_orth(k,l) + epss ) * Bterm(k,l) / (Jp(k,l)**2) * Jp_r(k,l)  &
             -  eps2 * phieta(k,l,i) * geta_size(m) / ( reta(k,l)**2 + zeta(k,l)**2 ) &
             * 2.0_rk*reta(k,l) * phieta(k,l,j)

        integrandReta_z_V(k,l) =  -1.0_rk/s_orth(k,l)**2 * s_orth_z(k,l) * Bterm(k,l) / Jp(k,l)  &
             +  ( 1.0_rk/s_orth(k,l) + epss ) * Bterm_z(k,l) / Jp(k,l)  &
             -  ( 1.0_rk/s_orth(k,l) + epss ) * Bterm(k,l) / (Jp(k,l)**2) * Jp_z(k,l)  &
             -  eps2 * phieta(k,l,i) * geta_size(m) / ( reta(k,l)**2 + zeta(k,l)**2 ) &
             * 2.0_rk*zeta(k,l) * phieta(k,l,j)


        ! inertial = 
        ! viscous = 
        ! integrandRu_r(k,l) = inertial + viscous
        integrandRu_r_V(k,l) = ( Re*phi(k,l,i)* ( -CTJ*phi(k,l,j)/dt*urintfac(k,l) + &
             ( uintfac(k,l) - rdotintfac(k,l) )*urintfac_r(k,l) + ( vintfac(k,l)- zdotintfac(k,l) )*uzintfac_r(k,l) ) - &
             Kdi*pintfac(k,l)*Oh*( phir_r(k,l,i) - phi(k,l,i)/rintfac(k,l)**2 *phi(k,l,j) ) + &
             Oh*( 2.0_rk*urintfac_r(k,l)*phir(k,l,i) + 2.0_rk*urintfac(k,l)*phir_r(k,l,i) + &
             ( uzintfac_r(k,l) + vrintfac_r(k,l) )*phiz(k,l,i) + ( uzintfac(k,l) + vrintfac(k,l) )*phiz_r(k,l,i) - &
             4.0_rk*phi(k,l,i)*phi(k,l,j)/rintfac(k,l)**3 *uintfac(k,l) ) ) *rintfac(k,l)*abs(Jp(k,l)) + &
             
             ( Re*phi(k,l,i)*( udotintfac(k,l) + ( uintfac(k,l) - rdotintfac(k,l) )*urintfac(k,l) +  &
             ( vintfac(k,l) - zdotintfac(k,l) )*uzintfac(k,l) )  - &
             Kdi*pintfac(k,l)*Oh*( phir(k,l,i) + phi(k,l,i)/rintfac(k,l) ) + &
             Oh*( 2.0_rk*urintfac(k,l)*phir(k,l,i) + ( uzintfac(k,l) + vrintfac(k,l) )*phiz(k,l,i) + &
             phi(k,l,i)*2.0_rk/rintfac(k,l)**2 *uintfac(k,l) ) ) *rJp_r(k,l) 


        integrandRu_z_V(k,l) = ( Re*phi(k,l,i)* ( ( uintfac(k,l) - rdotintfac(k,l) )*urintfac_z(k,l) - &
             CTJ*phi(k,l,j)/dt*uzintfac(k,l) + ( vintfac(k,l)- zdotintfac(k,l) )*uzintfac_z(k,l) ) - &
             Kdi*pintfac(k,l)*Oh*phir_z(k,l,i) + &
             Oh*( 2.0_rk*urintfac_z(k,l)*phir(k,l,i) + 2.0_rk*urintfac(k,l)*phir_z(k,l,i) + &
             ( uzintfac_z(k,l) + vrintfac_z(k,l) )*phiz(k,l,i) + ( uzintfac(k,l) + vrintfac(k,l) )*phiz_z(k,l,i) ) ) &
             *rintfac(k,l)*abs(Jp(k,l)) + &
             
             ( Re*phi(k,l,i)*( udotintfac(k,l) + ( uintfac(k,l) - rdotintfac(k,l) )*urintfac(k,l) +  &
             ( vintfac(k,l) - zdotintfac(k,l) )*uzintfac(k,l) )  - &
             Kdi*pintfac(k,l)*Oh*( phir(k,l,i) + phi(k,l,i)/rintfac(k,l) ) + &
             Oh*( 2.0_rk*urintfac(k,l)*phir(k,l,i) + ( uzintfac(k,l) + vrintfac(k,l) )*phiz(k,l,i) + &
             phi(k,l,i)*2.0_rk/rintfac(k,l)**2 *uintfac(k,l) ) )  *rintfac(k,l)* Jp_z(k,l)



        integrandRu_u(k,l) = ( Re*phi(k,l,i)* ( CTJ*phi(k,l,j)/dt + phi(k,l,j)*urintfac(k,l) +  &
             ( uintfac(k,l) - rdotintfac(k,l) )*phir(k,l,j) + ( vintfac(k,l)- zdotintfac(k,l) )*phiz(k,l,j) ) + &
             Oh*( 2.0_rk*phir(k,l,i)*phir(k,l,j) + phiz(k,l,i)*phiz(k,l,j) + 2.0_rk/rintfac(k,l)**2 *phi(k,l,i)*phi(k,l,j) ) ) &
             *rintfac(k,l)*abs(Jp(k,l))  

        integrandRuv(k,l) = ( Re*phi(k,l,i)*phi(k,l,j)*uzintfac(k,l) + Oh*phir(k,l,j)*phiz(k,l,i) ) *rintfac(k,l)*abs(Jp(k,l)) 

        integrandRu_p(k,l) = -Kdi*psi(k,l,j)*Oh*( phir(k,l,i) + phi(k,l,i)/rintfac(k,l) ) *rintfac(k,l)*abs(Jp(k,l))
        !j=2,4,5,6,8,  integrandRu_p(k,l) = 0.0_rk




        integrandRv_r_V(k,l) = ( Re*phi(k,l,i)* ( -CTJ*phi(k,l,j)/dt*vrintfac(k,l) + &
             ( uintfac(k,l) - rdotintfac(k,l) )*vrintfac_r(k,l) + ( vintfac(k,l)- zdotintfac(k,l) )*vzintfac_r(k,l) ) - &
             Kdi*pintfac(k,l)*Oh*phiz_r(k,l,i) + &
             Oh*( ( uzintfac_r(k,l) + vrintfac_r(k,l) )*phir(k,l,i) + ( uzintfac(k,l) + vrintfac(k,l) )*phir_r(k,l,i) + &
             2.0_rk*vzintfac_r(k,l)*phiz(k,l,i) + 2.0_rk*vzintfac(k,l)*phiz_r(k,l,i) ) ) *rintfac(k,l)*abs(Jp(k,l)) + &
             
             ( Re*phi(k,l,i)*( vdotintfac(k,l) + ( uintfac(k,l) - rdotintfac(k,l) )*vrintfac(k,l) + &
             ( vintfac(k,l) - zdotintfac(k,l) )*vzintfac(k,l) + Grav ) - Kdi*pintfac(k,l)*Oh*phiz(k,l,i) + &
             Oh*( ( uzintfac(k,l) + vrintfac(k,l) ) *phir(k,l,i) + 2.0_rk*vzintfac(k,l)*phiz(k,l,i) ) ) *rJp_r(k,l) 


        integrandRv_z_V(k,l) = ( Re*phi(k,l,i)* ( ( uintfac(k,l) - rdotintfac(k,l) )*vrintfac_z(k,l) - &
             CTJ*phi(k,l,j)/dt*vzintfac(k,l) + ( vintfac(k,l)- zdotintfac(k,l) )*vzintfac_z(k,l) ) - &
             Kdi*pintfac(k,l)*Oh*phiz_z(k,l,i) + &
             Oh*( ( uzintfac_z(k,l) + vrintfac_z(k,l) )*phir(k,l,i) + ( uzintfac(k,l) + vrintfac(k,l) )*phir_z(k,l,i) + &
             2.0_rk*vzintfac_z(k,l)*phiz(k,l,i) + 2.0_rk*vzintfac(k,l)*phiz_z(k,l,i) ) ) *rintfac(k,l)*abs(Jp(k,l)) + &
             
             ( Re*phi(k,l,i)*( vdotintfac(k,l) + ( uintfac(k,l) - rdotintfac(k,l) )*vrintfac(k,l) + &
             ( vintfac(k,l) - zdotintfac(k,l) )*vzintfac(k,l) + Grav ) - Kdi*pintfac(k,l)*Oh*phiz(k,l,i) + &
             Oh*( ( uzintfac(k,l) + vrintfac(k,l) ) *phir(k,l,i) + 2.0_rk*vzintfac(k,l)*phiz(k,l,i) ) ) *rintfac(k,l)* Jp_z(k,l)


        integrandRv_u(k,l) = ( Re*phi(k,l,i)*phi(k,l,j)*vrintfac(k,l) + Oh*phiz(k,l,j)*phir(k,l,i) ) *rintfac(k,l)*abs(Jp(k,l)) 

        integrandRvv(k,l) = ( Re*phi(k,l,i)*( CTJ*phi(k,l,j)/dt + ( uintfac(k,l) - rdotintfac(k,l) )*phir(k,l,j) + &
             phi(k,l,j)*vzintfac(k,l) + (vintfac(k,l) - zdotintfac(k,l) )*phiz(k,l,j) ) + &
             Oh*( phir(k,l,i)*phir(k,l,j) + 2.0_rk*phiz(k,l,i)*phiz(k,l,j) ) ) *rintfac(k,l)*abs(Jp(k,l)) 


        integrandRv_p(k,l) = -Kdi*psi(k,l,j)*Oh*phiz(k,l,i) *rintfac(k,l)*abs(Jp(k,l))  




        integrandRp_r(k,l) = psi(k,l,i)* ( urintfac_r(k,l) - uintfac(k,l)/rintfac(k,l)**2 *phi(k,l,j) + vzintfac_r(k,l) ) &
             *rintfac(k,l)*abs(Jp(k,l)) + &
             
             psi(k,l,i)* ( urintfac(k,l) + uintfac(k,l)/rintfac(k,l) + vzintfac(k,l) ) *rJp_r(k,l)


        integrandRp_z(k,l) = psi(k,l,i)* ( urintfac_z(k,l) + vzintfac_z(k,l) ) *rintfac(k,l)*abs(Jp(k,l)) + &
             
             psi(k,l,i)* ( urintfac(k,l) + uintfac(k,l)/rintfac(k,l) + vzintfac(k,l) ) *rintfac(k,l)* Jp_z(k,l)


        integrandRp_u(k,l) = ( phi(k,l,j)/rintfac(k,l) + phir(k,l,j) )*psi(k,l,i) *rintfac(k,l)*abs(Jp(k,l))

        integrandRp_v(k,l) = phiz(k,l,j)*psi(k,l,i) *rintfac(k,l)*abs(Jp(k,l)) 




     end do
  end do   !end loop for k&l, preparation for gaussian_quadrature


  sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = gaussian_quadrature(integrandRsi_r_V)            ! sjRsi_r(i,j)   
  sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = gaussian_quadrature(integrandRsi_z_V)            ! sjRsi_z(i,j)  
  sj(LNOPP(i)+Nr,LNOPP(j)+Nu) = 0.0_rk                                            ! sjRsi_u(i,j)   
  sj(LNOPP(i)+Nr,LNOPP(j)+Nv) = 0.0_rk                                            ! sjRsi_v(i,j)  
  if( MDF( globalNM(m,j) ).eq.NNp ) then
     sj(LNOPP(i)+Nr,LNOPP(j)+Np) = 0.0_rk                                          ! sjRsi_p(i,j)  
  end if

  sj(LNOPP(i)+Nz,LNOPP(j)+Nr) = gaussian_quadrature(integrandReta_r_V)           ! sjReta_r(i,j)   
  sj(LNOPP(i)+Nz,LNOPP(j)+Nz) = gaussian_quadrature(integrandReta_z_V)           ! sjReta_z(i,j)  
  sj(LNOPP(i)+Nz,LNOPP(j)+Nu) = 0.0_rk                                          ! sjReta_u(i,j)   
  sj(LNOPP(i)+Nz,LNOPP(j)+Nv) = 0.0_rk                                          ! sjReta_v(i,j)  
  if( MDF( globalNM(m,j) ).eq.NNp ) then
     sj(LNOPP(i)+Nz,LNOPP(j)+Np) = 0.0_rk                                          ! sjReta_p(i,j)  
  end if

  sj(LNOPP(i)+Nu,LNOPP(j)+Nr) = gaussian_quadrature(integrandRu_r_V)           ! sjRur(i,j)   
  sj(LNOPP(i)+Nu,LNOPP(j)+Nz) = gaussian_quadrature(integrandRu_z_V)           ! sjRuz(i,j) 
  sj(LNOPP(i)+Nu,LNOPP(j)+Nu) = gaussian_quadrature(integrandRu_u)           ! sjRuu(i,j)   
  sj(LNOPP(i)+Nu,LNOPP(j)+Nv) = gaussian_quadrature(integrandRuv)           ! sjRuv(i,j)  
  if( MDF( globalNM(m,j) ).eq.NNp ) then
     sj(LNOPP(i)+Nu,LNOPP(j)+Np) = gaussian_quadrature(integrandRu_p)           ! sjRup(i,j) 
  end if

  sj(LNOPP(i)+Nv,LNOPP(j)+Nr) = gaussian_quadrature(integrandRv_r_V)           ! sjRvr(i,j)   
  sj(LNOPP(i)+Nv,LNOPP(j)+Nz) = gaussian_quadrature(integrandRv_z_V)           ! sjRvz(i,j) 
  sj(LNOPP(i)+Nv,LNOPP(j)+Nu) = gaussian_quadrature(integrandRv_u)           ! sjRvu(i,j)   
  sj(LNOPP(i)+Nv,LNOPP(j)+Nv) = gaussian_quadrature(integrandRvv)           ! sjRvv(i,j)  
  if( MDF( globalNM(m,j) ).eq.NNp ) then
     sj(LNOPP(i)+Nv,LNOPP(j)+Np) = gaussian_quadrature(integrandRv_p)           ! sjRvp(i,j) 
  end if


  if( MDF( globalNM(m,i) ).eq.NNp ) then
     sj(LNOPP(i)+Np,LNOPP(j)+Nr) = gaussian_quadrature(integrandRp_r)           ! sjRpr(i,j)   
     sj(LNOPP(i)+Np,LNOPP(j)+Nz) = gaussian_quadrature(integrandRp_z)           ! sjRpz(i,j) 
     sj(LNOPP(i)+Np,LNOPP(j)+Nu) = gaussian_quadrature(integrandRp_u)           ! sjRpu(i,j)   
     sj(LNOPP(i)+Np,LNOPP(j)+Nv) = gaussian_quadrature(integrandRp_v)           ! sjRpv(i,j)  
     if( MDF( globalNM(m,j) ).eq.NNp ) then
        sj(LNOPP(i)+Np,LNOPP(j)+Np) = 0.0_rk                                        ! sjRpp(i,j) 
     end if
  end if
                       

end subroutine VI_in_sj




subroutine SI_in_sj(m,i,j, sj, LNVar, LNOPP)                     !adding SI to sj
  use kind
  use data
  use Ldata
  use NOP_mod, only: gaussian_quadrature_1d

  implicit none

  integer(kind=ik), intent(in):: m,i,j, LNVar, LNOPP(9)
  real(kind=rk), intent(out):: sj(LNVar, LNVar)

  integer(kind=ik):: k, ipp, jpp
  real(kind=rk):: integrandRsi_r_S(Ng), integrandRsi_z_S(Ng), integrandReta_r_S(Ng), integrandReta_z_S(Ng)
  real(kind=rk):: integrandRsi_u_S(Ng), integrandRsi_v_S(Ng), &
       integrandRu_r_S(Ng), integrandRu_z_S(Ng), integrandRv_r_S(Ng), integrandRv_z_S(Ng)



  if( (rowNM(m).eq.1) .and. ( (i.eq.1) .or. (i.eq.2) .or. (i.eq.3) ) .and. ( (j.eq.1) .or. (j.eq.2) .or. (j.eq.3) ) ) then
     do k = 1, Ng, 1    !three gausspoints
        ipp = i  !phix_1d(k,ipp)
        jpp = j  !phix_1d(k,jpp)
        integrandRsi_r_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( rsi_down(k)**2 + zsi_down(k)**2 ) &
             * 2.0_rk*rsi_down(k) * phix_1d(k,jpp)
        integrandRsi_z_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( rsi_down(k)**2 + zsi_down(k)**2 ) &
             * 2.0_rk*zsi_down(k) * phix_1d(k,jpp)
     end do
     sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = sj(LNOPP(i)+Nr,LNOPP(j)+Nr) - M1*gaussian_quadrature_1d(integrandRsi_r_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = sj(LNOPP(i)+Nr,LNOPP(j)+Nz) - M1*gaussian_quadrature_1d(integrandRsi_z_S)
  end if

  if( (rowNM(m).eq.NEY) .and. ( (i.eq.7) .or. (i.eq.8) .or. (i.eq.9) ) .and. ( (j.eq.7) .or. (j.eq.8) .or. (j.eq.9) ) ) then
     do k = 1, Ng, 1    !three gausspoints
        ipp = i - 6  !phix_1d(k,l)
        jpp = j - 6  !phix_1d(k,jpp)
        integrandRsi_r_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( rsi_up(k)**2 + zsi_up(k)**2 ) &
             * 2.0_rk*rsi_up(k) * phix_1d(k,jpp)
        integrandRsi_z_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( rsi_up(k)**2 + zsi_up(k)**2 ) &
             * 2.0_rk*zsi_up(k) * phix_1d(k,jpp)
     end do
     sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = sj(LNOPP(i)+Nr,LNOPP(j)+Nr) - M1*gaussian_quadrature_1d(integrandRsi_r_S)
     sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = sj(LNOPP(i)+Nr,LNOPP(j)+Nz) - M1*gaussian_quadrature_1d(integrandRsi_z_S)
  end if

  if( (columnNM(m).eq.1) .and. ( (i.eq.1) .or. (i.eq.4) .or. (i.eq.7) ) .and. ( (j.eq.1) .or. (j.eq.4) .or. (j.eq.7) ) ) then
     do k = 1, Ng, 1    !three gausspoints
        ipp = i/3 + 1  !phix_1d(k,ipp)
        jpp = j/3 + 1  !phix_1d(k,jpp)
        integrandReta_r_S(k) = phix_1d(k,ipp) * geta_size(m) / ( reta_left(k)**2 + zeta_left(k)**2 ) &
             * 2.0_rk*reta_left(k) * phix_1d(k,jpp)
        integrandReta_z_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( reta_left(k)**2 + zeta_left(k)**2 ) &
             * 2.0_rk*zeta_left(k) * phix_1d(k,jpp)
     end do
     sj(LNOPP(i)+Nz,LNOPP(j)+Nr) = sj(LNOPP(i)+Nz,LNOPP(j)+Nr) - M1*gaussian_quadrature_1d(integrandReta_r_S)
     sj(LNOPP(i)+Nz,LNOPP(j)+Nz) = sj(LNOPP(i)+Nz,LNOPP(j)+Nz) - M1*gaussian_quadrature_1d(integrandReta_z_S)
  end if



  if( (columnNM(m).eq.NEX) .and. ( (i.eq.3) .or. (i.eq.6) .or. (i.eq.9) ) )  then      

     if( (j.eq.3) .or. (j.eq.6) .or. (j.eq.9) ) then
        do k = 1, Ng, 1    !three gausspoints
           ipp = i/3  !phix_1d(k,l)
           jpp = j/3  !phix_1d(k,l)
           integrandReta_r_S(k) = phix_1d(k,ipp) * geta_size(m) / ( reta_right(k)**2 + zeta_right(k)**2 ) &
                * 2.0_rk*reta_right(k) * phix_1d(k,jpp)
           integrandReta_z_S(k) = phix_1d(k,ipp) * fsi_size(m) / ( reta_right(k)**2 + zeta_right(k)**2 ) &
                * 2.0_rk*zeta_right(k) * phix_1d(k,jpp)


           integrandRsi_r_S(k) = phi_1d(k,ipp)*( zeta_right(k)*CTJ/dt*phi_1d(k,jpp) + &
                phix_1d(k,jpp)*( vintfac_right(k) - zdotintfac_right(k) ) )*rintfac_right(k) + &
                
                phi_1d(k,ipp)*( -zeta_right(k)*( uintfac_right(k) - rdotintfac_right(k) ) + &
                reta_right(k)*( vintfac_right(k) - zdotintfac_right(k) ) )*phi_1d(k,jpp)

           integrandRsi_z_S(k) = phi_1d(k,ipp)*( -phix_1d(k,jpp)*( uintfac_right(k) - rdotintfac_right(k) ) - &
                reta_right(k)*CTJ/dt*phi_1d(k,jpp) )*rintfac_right(k)

           integrandRsi_u_S(k) = -phi_1d(k,jpp)*phi_1d(k,ipp)*zeta_right(k)*rintfac_right(k)

           integrandRsi_v_S(k) = phi_1d(k,jpp)*phi_1d(k,ipp)*reta_right(k)*rintfac_right(k)


           integrandRu_r_S(k) = ( phix_1d(k,ipp)*phix_1d(k,jpp) / ( reta_right(k)**2 + zeta_right(k)**2 ) - &
                reta_right(k)*phix_1d(k,ipp) / ( reta_right(k)**2 + zeta_right(k)**2 )**2 * 2.0_rk*reta_right(k)*phix_1d(k,jpp) - &
                1.0_rk/rintfac_right(k)**2 *phi_1d(k,jpp)*phi_1d(k,ipp) ) * &
                rintfac_right(k)*( reta_right(k)**2 + zeta_right(k)**2 )**0.5_rk + &
                
                ( reta_right(k)*phix_1d(k,ipp) / ( reta_right(k)**2 + zeta_right(k)**2 ) + phi_1d(k,ipp)/rintfac_right(k) )&
                *( phi_1d(k,jpp)*( reta_right(k)**2 + zeta_right(k)**2 )**0.5_rk + &
                rintfac_right(k) / ( reta_right(k)**2 + zeta_right(k)**2 )**0.5_rk *reta_right(k)*phix_1d(k,jpp) )

           integrandRu_z_S(k) = -reta_right(k)*phix_1d(k,ipp)*( reta_right(k)**2 + zeta_right(k)**2 )**(-1.5_rk) *&
                2.0_rk*zeta_right(k)*phix_1d(k,jpp)*rintfac_right(k) + &
                
                ( reta_right(k)*phix_1d(k,ipp) / ( reta_right(k)**2 + zeta_right(k)**2 ) + phi_1d(k,ipp)/rintfac_right(k) ) *&
                rintfac_right(k) / ( reta_right(k)**2 + zeta_right(k)**2 )**0.5_rk *zeta_right(k)*phix_1d(k,jpp)


           integrandRv_r_S(k) = zeta_right(k)*phix_1d(k,ipp)*&
                ( -( reta_right(k)**2 + zeta_right(k)**2 )**(-1.5_rk) *reta_right(k)*phix_1d(k,jpp)*rintfac_right(k) + &
                ( reta_right(k)**2 + zeta_right(k)**2 )**(-0.5_rk) *phi_1d(k,jpp) )

           integrandRv_z_S(k) = phix_1d(k,ipp)*phix_1d(k,jpp) / sqrt( reta_right(k)**2 + zeta_right(k)**2 ) *rintfac_right(k) - &
                zeta_right(k)**2 *phix_1d(k,ipp)*phix_1d(k,jpp) / ( reta_right(k)**2 + zeta_right(k)**2 )**1.5_rk *rintfac_right(k)


        end do
        sj(LNOPP(i)+Nz,LNOPP(j)+Nr) = sj(LNOPP(i)+Nz,LNOPP(j)+Nr) - M1*gaussian_quadrature_1d(integrandReta_r_S)
        sj(LNOPP(i)+Nz,LNOPP(j)+Nz) = sj(LNOPP(i)+Nz,LNOPP(j)+Nz) - M1*gaussian_quadrature_1d(integrandReta_z_S)!?should be M2, doesn't matter because M1=M2=0

        sj(LNOPP(i)+Nr,LNOPP(j)+Nr) = gaussian_quadrature_1d(integrandRsi_r_S)     !directly replace Rsi_V with Rsi_S, viz KBC
        sj(LNOPP(i)+Nr,LNOPP(j)+Nz) = gaussian_quadrature_1d(integrandRsi_z_S)
        sj(LNOPP(i)+Nr,LNOPP(j)+Nu) = gaussian_quadrature_1d(integrandRsi_u_S)
        sj(LNOPP(i)+Nr,LNOPP(j)+Nv) = gaussian_quadrature_1d(integrandRsi_v_S)

        sj(LNOPP(i)+Nu,LNOPP(j)+Nr) = sj(LNOPP(i)+Nu,LNOPP(j)+Nr) + gaussian_quadrature_1d(integrandRu_r_S)
        sj(LNOPP(i)+Nu,LNOPP(j)+Nz) = sj(LNOPP(i)+Nu,LNOPP(j)+Nz) + gaussian_quadrature_1d(integrandRu_z_S)

        sj(LNOPP(i)+Nv,LNOPP(j)+Nr) = sj(LNOPP(i)+Nv,LNOPP(j)+Nr) + gaussian_quadrature_1d(integrandRv_r_S)
        sj(LNOPP(i)+Nv,LNOPP(j)+Nz) = sj(LNOPP(i)+Nv,LNOPP(j)+Nz) + gaussian_quadrature_1d(integrandRv_z_S)

     else
        do k = Nr, Np
           if( MDF( globalNM(m,j) ).lt.NNp .and. k.eq.Np ) cycle
           sj(LNOPP(i)+Nr,LNOPP(j)+k) = 0.0_rk        !d Rsi/ d xj = 0.0_rk other than j = 3,6,9
        end do

     end if
  end if



end subroutine SI_in_sj





subroutine jacobian_check(m, i, sj, sf, LNOPP, LNVar, id)
  use kind
  use data

  implicit none

  integer(kind=ik), intent(in):: m, i, LNVar, LNOPP(LNVar), id
  real(kind=rk), intent(in):: sj(LNVar, LNVar), sf(LNVar)

  integer(kind=ik):: k,p,q
  real(kind=rk):: sj_check, sf_check(LNVar), diffsj, rel, delta, TOL

  delta = 10e-8_rk
  TOL = 10e-6_rk

  !change sol( NOPP( globalNM(m,k) ) + p )
  do k = 1, bas
     do p = 0, MDF( globalNM(m,k) )-1    !p = 0 --> si, p = MDF( globalNM(m,k) ) --> v or p

        !add delta to xj
        sol( NOPP( globalNM(m,k) ) + p ) = sol( NOPP( globalNM(m,k) ) + p ) + delta       
        !define soldot
        if (timestep.le.5) then
           soldot = ( sol - solp )/dt
        else
           soldot = 2.0_rk*( sol - solp )/dt - soldotp
        end if


        call values_in_an_element(m,id)
! if(m.eq.1 .and. i.eq.3) then
! write(*,*) 'check jacobian'
! pause
! end if
        call define_sf(m,i,sf_check, LNVar, LNOPP, id)

        
        !"do i = 1, bas" is in assemble subroutine
        do q = 0, MDF( globalNM(m,i) )-1    !q = 0 --> Rsi, q = MDF( globalNM(m,i) ) --> Rv or Rp
           
           if(q.eq.0) then
              if(i.ne.3 .and. i.ne.6 .and. i.ne.9) cycle
           end if
           if(q.eq.1) cycle

           sj_check = ( sf_check( LNOPP(i)+q ) - sf( LNOPP(i)+q ) )/delta
           diffsj = abs( sj_check - sj( LNOPP(i)+q, LNOPP(k)+p ) ) 
           rel = abs( diffsj / sj( LNOPP(i)+q, LNOPP(k)+p ) )

           if( rel.gt.TOL ) then
              if( diffsj.gt.TOL) then
                 write(*,*) 'ele =', m, ', i =', i,q+1, ', j =', k,p+1
                 write(*,*) 'sf =', sf( LNOPP(i)+q ), ', sf_check =', sf_check( LNOPP(i)+q )
                 write(*,*) 'sj =', sj( LNOPP(i)+q , LNOPP(k)+p ), ', sj_check =', sj_check
                !  pause
              ! else
              !    write(*,*) 'relative error large, but diffsj small'
              end if
           ! else
           !    write(*,*) 'ele =', m, ', i =', i,q+1, ', j =', k,p+1
           end if
        end do

        !substract deltax to uj
        sol( NOPP( globalNM(m,k) ) + p ) = sol( NOPP( globalNM(m,k) ) + p ) - delta
        !define soldot
        if (timestep.le.5) then
           soldot = ( sol - solp )/dt
        else
           soldot = 2.0_rk*( sol - solp )/dt - soldotp
        end if
        call values_in_an_element(m,id)

     end do
  end do

  return

end subroutine jacobian_check

subroutine break_check
use kind
use data, only: R, NTE, NEX, NNX, globalNM, columnNM, rcoordinate, timestep, time, Hmin, zcoordinate, vsol, Oh,&
     Hminlast, Hmin2last, tlast, t2last


implicit none

integer(kind=ik):: i, ipp, m, Ehmin, Ehmin1     
!Ehmin: the element number where Hmin appears
!Ehmin1: the element number where Hmin1 appears
real(kind=rk):: Hmin1, zhmin, vhmin, zhmin1, vhmin1, dz, dv      
!Hmin1 is 1.2*Hmin; zhmin, vhmin are at Hmin; zhmin1, vhmin1 are at Hmin1
real(kind=rk):: rr(3), z(3), v(3), phi(3)
real(kind=rk):: A, B, C, x1, x2, x, ReLocal


     Hmin = R
     zhmin = 0.0_rk
     Ehmin = NEX
     zhmin1 = 0.0_rk
     vhmin1 = 0.0_rk
     do m = 1, NTE
        if( columnNM(m) .eq. NEX ) then
           if( rcoordinate( globalNM(m,3) ) .lt. Hmin ) then
              Hmin = rcoordinate( globalNM(m,3) )
              zhmin = zcoordinate( globalNM(m,3) )
              Ehmin = m
           end if
           if( rcoordinate( globalNM(m,6) ) .lt. Hmin ) then
              Hmin = rcoordinate( globalNM(m,6) )
              zhmin = zcoordinate( globalNM(m,6) )
              Ehmin = m
           end if
           if( rcoordinate( globalNM(m,9) ) .lt. Hmin ) then
              Hmin = rcoordinate( globalNM(m,9) )
              zhmin = zcoordinate( globalNM(m,9) )
              Ehmin = m
           end if
        end if
     end do
     write(*,*) 'Hmin:', Hmin
     
     Hmin2last = Hminlast
     Hminlast = Hmin
     t2last = tlast
     tlast = time
     
     if (timestep.eq.0) then
        open(unit = 10, file = 'Hmin.dat', status = 'replace')
        write(10,'(A)') 'variables = "t", "Hmin" '
        write(10,'(A)') 'zone T  = "solution" ' 
     else
        open(unit = 10, file = 'Hmin.dat', status = 'old', access = 'append')
     end if

     write(10,'(2es14.7)') time, Hmin
     close (10)


     Hmin1 = 1.2_rk*Hmin
     do m = NEX, Ehmin, NEX       !only consider the situation when Hmin1 is under Hmin
        if( ( rcoordinate( globalNM(m,9) ) - Hmin1 )*( rcoordinate( globalNM(m,3) ) - Hmin1 ) .le. 0.0_rk ) then
           Ehmin1 = m
           exit
        end if
     end do
        
     do i = 1, 3
        ipp = i*3
        rr(i) = rcoordinate( globalNM(Ehmin1,ipp) ) 
        z(i) = zcoordinate( globalNM(Ehmin1,ipp) ) 
        v(i) = vsol( globalNM(Ehmin1,ipp) ) 
     end do

     A = 2.0_rk*rr(1) - 4.0_rk*rr(2) + 2.0_rk*rr(3)
     B = -3.0_rk*rr(1) + 4.0_rk*rr(2) - rr(3)
     C = rr(1) - Hmin1
     x1 = 1.0_rk/(2.0_rk*A) * ( -B - sqrt( B**2 - 4.0_rk*A*C ) )
     x2 = 1.0_rk/(2.0_rk*A) * ( -B + sqrt( B**2 - 4.0_rk*A*C ) )
     x = 0.0_rk
     if( x1.gt.0.0_rk .and. x1.lt.1.0_rk )  x=x1
     if( x2.gt.0.0_rk .and. x2.lt.1.0_rk )  x=x2
     if(x.ne.x1 .and. x.ne.x2) write(*,*) "error in finding si of Hmin1"
     
     phi(1) = 1.0_rk - 3.0_rk*x + 2.0_rk*x**2
     phi(2) = 4.0_rk*( x - x**2 )
     phi(3) = -x + 2.0_rk*x**2

     zhmin1 = 0.0_rk
     vhmin1 = 0.0_rk
     do i = 1, 3
        zhmin1 = zhmin1 + z(i)*phi(i)
        vhmin1 = vhmin1 + v(i)*phi(i)
     end do

     dz = zhmin - zhmin1
     ! dv = vhmin - vhmin1
     
     ReLocal =  abs( dz*vhmin1/Oh )

     if (timestep.eq.0) then
        open(unit = 10, file = 'ReLocal.dat', status = 'replace')
        write(10,'(A)') 'variables = "Hmin", "ReLocal", "x", "dz", "z", "v" '
        write(10,'(A)') 'zone T  = "solution" ' 
     else
        open(unit = 10, file = 'ReLocal.dat', status = 'old', access = 'append')
     end if

     write(10,'(6es15.7)') Hmin, ReLocal, x, dz, zhmin1, vhmin1
     close (10)

   end subroutine break_check

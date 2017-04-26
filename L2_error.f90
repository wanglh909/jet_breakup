subroutine L2_error
use kind
use data, only: NVar, NTN, timestep, step, time, dt, error2, dsol, Nr, Nz, Nu, Nv, Np
use NOP_mod, only: NOPP, MDF

implicit none
integer(kind=ik):: i
real(kind=rk):: dr(NTN), dz(NTN), du(NTN), dv(NTN), dp(NTN), drmax, dzmax, dumax, dvmax, dpmax

        !define error2
        error2 = 0.0_rk
        do i = 1, NVar, 1      
           error2 = error2 + dsol(i)**2
        end do
        error2 = sqrt(error2)             !absolute value of dsol


        dp(:) = 0.0_rk
        do i = 1, NTN, 1
           dr(i) = abs( dsol( NOPP(i) + Nr ) )
           dz(i) = abs( dsol( NOPP(i) + Nz ) )
           du(i) = abs( dsol( NOPP(i) + Nu ) )
           dv(i) = abs( dsol( NOPP(i) + Nv ) )
           if ( MDF(i) .gt. Np )     dp(i) = abs( dsol( NOPP(i) + Np ) )
        end do
        drmax = 0.0d0
        dzmax = 0.0d0
        dumax = 0.0d0
        dvmax = 0.0d0
        dpmax = 0.0d0
        do i = 1, NTN
           if( dr(i) .gt. drmax )   drmax = dr(i)
           if( dz(i) .gt. dzmax )   dzmax = dz(i)
           if( du(i) .gt. dumax )   dumax = du(i)
           if( dv(i) .gt. dvmax )   dvmax = dv(i)
           if( dp(i) .gt. dpmax )   dpmax = dp(i)
        end do
        
     if (timestep.eq.0 .and. step.eq.1) then
        open(unit = 10, file = 'L2_error.dat', status = 'replace')
        write(10, '(A)') 'variables = "dr", "dz", "du", "dv", "dp", "total" '
     else
        open(unit = 10, file = 'L2_error.dat', status = 'old', access = 'append')
     end if

     if(step.eq.1) then
        write(10, '(A)') ' '
        write(10, '(A, 2es15.7, i8, A)') '----------------', time, dt, timestep, '---------------------'
        write(10, '(A)') ' '
     end if

     write(10,'(6es15.7)') drmax, dzmax, dumax, dvmax, dpmax, error2
     close(10)

   end subroutine L2_error

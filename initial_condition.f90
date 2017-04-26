
subroutine initial_condition
  use kind
  use data
  use NOP_mod, only: NOP



  implicit none

  integer(kind=ik):: i, j, rowN, columnN  !the rowN, columnN for a node, different from rowNM

  rcoordinate(:) = 0.0_rk
  zcoordinate(:) = 0.0_rk


  do i = 1, NTN, 1            !globalN

     if( mod( real(i,rk),real(NNX,rk) ) .eq. 0.0_rk ) then
        rowN = i/NNX           !no dble, rowN of nodes not elements (different from rowN in NOP)
        columnN = NNX
     else
        rowN = i/NNX + 1         !no dble
        columnN = mod(i,NNX)
     end if

     zcoordinate( i ) = H/real(2*NEY,rk)*real(rowN-1,rk)
     rcoordinate( i ) = ( R - 0.1_rk*cos( zcoordinate(i)*3.14_rk/H ) )/real(NNX-1,rk)*real(columnN-1,rk)

  end do

  usol(:) = 0.0_rk
  vsol(:) = 0.0_rk
  psol(:) = 0.0_rk    !include p at the node and p needed to be interpolated



  !**********************************************define sol(:)********************************************
  do i = 1, NTN, 1
     sol( NOPP(i) + Nr ) = rcoordinate(i)
     sol( NOPP(i) + Nz ) = zcoordinate(i)
     sol( NOPP(i) + Nu ) = usol(i)
     sol( NOPP(i) + Nv ) = vsol(i)
     if ( MDF(i) .gt. Np )    sol( NOPP(i) + Np ) = psol(i)
  end do
  !*******************************************************************************************************

Hminlast = 0.0_rk
Hmin2last = 0.0_rk
tlast = 0.0_rk
t2last = 0.0_rk

end subroutine initial_condition

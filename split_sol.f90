
subroutine split_sol
  use kind
  use data, only: NTN, NNX, NEY, rcoordinate, zcoordinate, usol, vsol, psol, sol, Nr, Nz, Nu, Nv, Np, NOPP, MDF

  implicit none
  integer(kind=ik):: i, j

  do i = 1, NTN, 1
     rcoordinate(i) = sol( NOPP(i) + Nr ) 
     zcoordinate(i) = sol( NOPP(i) + Nz ) 
     usol(i) = sol( NOPP(i) + Nu ) 
     vsol(i) = sol( NOPP(i) + Nv ) 
     if ( MDF(i) .gt. Np )     psol(i) = sol( NOPP(i) + Np ) 
  end do

  !interpolate p(:)
  do i = 1, 2*NEY + 1, 2
     do j = 2, NNX-1, 2
        psol( j + NNX*(i-1) ) = 0.5_rk*( psol( j + NNX*(i-1) - 1 ) + psol( j + NNX*(i-1) + 1 ) )
     end do
  end do
  do i = 2, 2*NEY , 2
     do j = 1, NNX, 1
        psol( j + NNX*(i-1) ) = 0.5_rk*( psol( j + NNX*(i-1) - NNX ) + psol( j + NNX*(i-1) + NNX ) )
     end do
  end do




end subroutine split_sol


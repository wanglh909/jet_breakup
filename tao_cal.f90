subroutine tao_calculation
  use kind
  use data, only: Hminlast, Hmin2last, tlast, t2last

  implicit none
  real(kind=rk):: t0


  t0 = ( t2last - Hmin2last/Hminlast*tlast ) / ( 1.0_rk - Hmin2last/Hminlast )
  write(*,*) 't0 =', t0

end subroutine tao_calculation

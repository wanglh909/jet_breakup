!-*- mode: f90;-*-  

!p = 15 is double, p = 17 is extended, p = 19 is quadruple
!if you change p, you have to make clean 

!it has to be the first thing compiled
!use this module in every subroutine and function and module (if you use it in a module it applies to all subs and functions contained in the module automatically).

  !integer:: -> integer(kind=ik)::
  !double precision:: -> real(kind=rk)::
  !Technically integers like 5 should be 5_ik but don't worry about that since you will always use the default kind of 4 for integers probably, I do and didn't change them.

  ! k = 1.0d0 -> 1.0_rk
  ! k = 1.0d-9 -> 1.0e-9_rk
  ! k = DBLE(m) -> REAL(m,rk) ; inverse: INT(m,ik)

Module kind
  !Specifies kind values for integer and real types
  integer*4, parameter :: rk = selected_real_kind(p=15), rk2 = selected_real_kind(p=15), ik = 4

end Module kind

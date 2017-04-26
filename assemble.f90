
  !call assemble_local(ele,local,loc,NOPPl,NB,id,dum)
  !ele: element to be assembled (DONT CHANGE)
  !local(NB,NB): local matrix (Warning assemble as (j,i), we normally assemble i,j as i = row, j = col)
  !NB: Integer, defines sizes (# number of variables in element (DONT CHANGE)
  !loc(NB): local rhs
  !NOPPl(bas): Local NOPPl, mut be filled and returned
  !id: Thread id, integer (DONT CHANGE)
  !dum: integer no purpose (DONT CHANGE)

subroutine assemble_local(m, locJacTr, locRHS, LNOPP, LNVar, id, dum)

  use kind
  use data
  use Ldata, only: NNp

  implicit none

  integer(kind=ik), intent(in):: m, LNVar, id, dum
  real(kind=rk), intent(out):: locJacTr(LNVar, LNVar), locRHS(LNVar)

  integer(kind=ik):: LNOPP(bas)
  real(kind=rk):: locJac(LNVar, LNVar), locRes(LNVar)

  integer(kind=ik):: i,j, k,l, p,q
  integer(kind=ik):: ele_c

  ele_c = 80
  
  LNOPP = 0
  do i = 1, bas
     do j = 1, i-1
        LNOPP(i) = LNOPP(i) + MDF( globalNM(m,j) )
     end do
     LNOPP(i) = LNOPP(i) + 1
  end do




  
  
  locJac = 0.0_rk
  locJacTr = 0.0_rk
  locRes = 0.0_rk !locRes = 0.0_rk
  locRHS = 0.0_rk  
  
  call values_in_an_element(m)


  do i = 1, bas, 1            !i is the i in the residual Res(i) & Jac(ij)

     !define locRes: R(r,z,c)i --> Res
     call define_sf(m,i, locRes, LNVar, LNOPP)

     !define locJac: dR(r,z)i/d(u,v,p)j --> Jac
     do j = 1, bas, 1            !stand for (u,v,p)j, j is the j in the notation of uj, vj, pj
        call values_in_sj(m,i,j)
        call VI_in_sj(m,i,j, locJac, LNVar, LNOPP)
        call SI_in_sj(m,i,j, locJac, LNVar, LNOPP)
     end do         !end loop for j

  if( timestep.eq.1 .and. m.eq.ele_c)  call jacobian_check(m, i, locJac, locRes, LNOPP, LNVar)

  end do             !end loop for i, Res(i) & Jac(i,j)

  call Dirichlet_BC(m, locJac, locRes, LNVar, LNOPP)
  
  
  do i = 1, LNVar
     do j = 1, LNVar
        locJacTr(i,j) = locJac(j,i)
     end do
  end do

  locRHS = -locRes  !locRes


  if( timestep.eq.1 .and. m.eq.ele_c ) pause

end subroutine assemble_local





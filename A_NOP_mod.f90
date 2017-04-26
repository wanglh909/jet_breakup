module NOP_mod
  use kind
  use data

  implicit none




contains

  subroutine variableN      !defines MDF & NOPP & rNOP
    implicit none

    integer(kind=ik):: i,j,k,l

    allocate( MDF(NTN), NOPP(NTN) )

    MDF(:) = NNVar-1
    do i = 1, NNX, 2
       do j = 0, NEY ,1
          MDF( i + j*2*NNX ) = NNVar
       end do
    end do

    NOPP(:) = 0
    do i = 1, NTN, 1
       do j = 1, i-1, 1
          NOPP(i) = NOPP(i) + MDF(j)
       end do
       NOPP(i) = NOPP(i) + 1
    end do

    NVar = NOPP(NTN) + MDF(NTN) - 1


    allocate( rNOP(NTN,4,2) )
    rNOP(:,:,:) = 0
    Do i = 1, NTE, 1
       Do j = 1, 9, 1
          l = 0
          k = 1
152       if (l.eq.0.and.k.le.4) then
             if (rNOP(globalNM(i,j),k,1).eq.0) then
                rNOP(globalNM(i,j),k,1) = i
                rNOP(globalNM(i,j),k,2) = j
                l = 1
             end if
             k = k + 1
             goto 152

          end if
       end do
    end do

    ! open(unit = 10, file = 'rNOP.dat', status = 'replace')
    ! do i = 1, NTN, 1
    !    write(10,'(A, i8, A)') '----------------------globalN:', i, '------------------------'
    !    do j = 1, 4
    !       write(10,'(2i8)') rNOP(i,j,1), rNOP(i,j,2)
    !    end do
    ! end do
    ! close(10)
    ! pause


  end subroutine variableN
  



subroutine NOP
    implicit none
    integer(kind=ik) :: eleN, locN
    integer(kind=ik):: globalNb  !globalN basis for each element

    NNX = 2*NEX + 1
    NTN = NNX*( 2*NEY + 1 )
    NTE = NEX*NEY

    allocate(  globalNM(NTE,9), rowNM(NTE), columnNM(NTE) )

    do eleN = 1, NTE
       !calculate rowN & columnN
       if( mod( dble(eleN),dble(NEX) ) .eq. 0.0_rk ) then
          rowNM(eleN) = eleN/NEX              !no dble, rowN of elements
          columnNM(eleN) = NEX
       else
          rowNM(eleN) = eleN/NEX + 1          !no dble
          columnNM(eleN) = mod(eleN,NEX)
       end if


       globalNb = 2*NNX*(rowNM(eleN)-1) + 2*(columnNM(eleN)-1)

       do locN = 1, bas

          !calculate globalN
          select case (locN)
          case(1)
             globalNM(eleN, locN) = globalNb + 1
          case(2)
             globalNM(eleN, locN) = globalNb + 2
          case(3)
             globalNM(eleN, locN) = globalNb + 3
          case(4)
             globalNM(eleN, locN) = globalNb + NNX + 1
          case(5)
             globalNM(eleN, locN) = globalNb + NNX + 2
          case(6)
             globalNM(eleN, locN) = globalNb + NNX + 3
          case(7)
             globalNM(eleN, locN) = globalNb + 2*NNX + 1
          case(8)
             globalNM(eleN, locN) = globalNb + 2*NNX + 2
          case(9)
             globalNM(eleN, locN) = globalNb + 2*NNX + 3
          end select

       end do

    end do

  end subroutine NOP



  

  
  subroutine size_function

    implicit none
    integer(kind=ik)::i,j
    real(kind=rk):: x   !the value of fsi_size & geta_size at the most dense area

    allocate( fsi_size(NTE), geta_size(NTE) )

    x = 0.1_rk
    x = 1.0_rk - x
    
    do i = 1, NTE

       ! fsi_size(i) = x/real(1-NEX,rk)*real(columnNM(i),rk) + 1.0_rk - x/real(1-NEX,rk)

       ! if( rowNM(i).le.NEY/2 ) then
       !    geta_size(i) = x/real(1-NEY/2,rk)*real(rowNM(i),rk) + 1.0_rk - x/real(1-NEY/2,rk)
       ! else 
       !    geta_size(i) = x/real(NEY/2-1,rk)*real(rowNM(i),rk) + 1.0_rk - x/(0.5_rk - 1/real(NEY,rk) )
       ! end if

       if( rowNM(i).le.3*NEY/4 ) then
          geta_size(i) = x/real(1-3*NEY/4,rk)*real(rowNM(i),rk) + 1.0_rk - x/real(1-3*NEY/4,rk)
       else 
          geta_size(i) = x/real(NEY/4,rk)*real(rowNM(i),rk) + 1.0_rk - x/0.25_rk
       end if

    end do

    fsi_size(:) = 1.0_rk
    ! geta_size(:) = 1.0_rk
   
  end subroutine size_function


  double precision function damfac(x)
    implicit none

    double precision:: x


    if(x.gt.20.0_rk) then
       damfac = 0.1_rk

    else if(x.gt.10.0_rk) then
       damfac = 0.3_rk

    else if(x.gt.1.0_rk) then
       damfac = 0.6_rk

    else
       damfac = 1.0_rk
    end if

  end function damfac







double precision function gaussian_quadrature(f) 
  implicit none

  !gaussian quadrature integration in 0-1 with 3 mesh points.
  ! f(Ng,Ng)(integrand) needs to be defined:
  !use "do" loop to give value to f(Ng): eg, f(i) = a(i)**2 (This stands for f(x) = x**2 )

  integer(kind=ik):: i, j
  real(kind=rk):: f(Ng,Ng)
  real(kind=rk), parameter, dimension(Ng):: w = (/0.2777777777777777777777777777778_rk,&
       0.44444444444444444444444444444444_rk, 0.2777777777777777777777777777778_rk/)

  gaussian_quadrature = 0.0_rk
  do i = 1,Ng,1
     do j = 1,Ng,1
     gaussian_quadrature = gaussian_quadrature + w(i)*w(j)*f(i,j)
     end do
  end do

end function gaussian_quadrature

double precision function gaussian_quadrature_1d(f)
  implicit none

  !gaussian quadrature integration in 0-1 with Ng mesh points.
  ! f(Ng)(integrand) needs to be defined
  !use "do" loop to give value to f(Ng): eg, f(i) = a(i)**2 (This stands for f(x) = x**2 )

  integer(kind=ik):: i
  real(kind=rk):: f(Ng)
  real(kind=rk), parameter, dimension(Ng):: w = (/0.2777777777777777777777777777778_rk,&
       0.44444444444444444444444444444444_rk, 0.2777777777777777777777777777778_rk/)

  gaussian_quadrature_1d = 0.0_rk
  do i = 1,Ng,1
     gaussian_quadrature_1d = gaussian_quadrature_1d + w(i)*f(i)
  end do

  return
end function gaussian_quadrature_1d






end module NOP_mod


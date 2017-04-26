


subroutine graph
  use kind
  use data, only: timestep, time, dt, NTN, NTE, rcoordinate, zcoordinate, usol, vsol, psol, globalNM, step


  implicit none

  integer(kind=ik):: i


  if (timestep.eq.0) then
     ! if (step.eq.0) then
     open(unit = 10, file = 'mesh.dat', status = 'replace')
  else
     open(unit = 10, file = 'mesh.dat', status = 'old', access = 'append')
  end if

  
  ! write(10, '(A)') ' '
  ! write(10, '(A, 2es15.7, i8, A)') '----------------next time step:', time, dt, timestep, '---------------------'
  ! write(10, '(A)') ' '
  write(10, '(A)') 'variables = "r", "z", "u", "v", "p"'
  write(10,'(A,es14.7,A,i8,A,i8,A)') 'Zone T = "mesh", STRANDID = 1, SOLUTIONTIME =', time, &
       ', Datapacking = Point, Zonetype = FEQuadrilateral, N =', NTN, ', E =', NTE, &
       ', DT = (double,double,double,double,double)'



  ! write(10,'(A,i8,A,i8,A,i8,A)') 'Zone T = "mesh", STRANDID = 1, SOLUTIONTIME =', step, &
  !      ', Datapacking = Point, Zonetype = FEQuadrilateral, N =', NTN, ', E =', NTE, &
  !      ', DT = (double,double,double,double,double)'


  do i = 1, NTN, 1
     write(10,'(5es15.7)') rcoordinate(i), zcoordinate(i), usol(i), vsol(i), psol(i)
  end do

  !form the quadrilateral with nodes
  do i = 1, NTE, 1
     write(10,'(4i8)') globalNM(i,1), globalNM(i,3), globalNM(i,9), globalNM(i,7)
  end do

  close(10)

end subroutine graph

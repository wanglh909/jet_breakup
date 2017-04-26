!-*- mode: f90;-*-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module front_mod
  use kind, only: rk, ik
  use omp_lib
  ! use data, only: NV, NN, iBW, rNOP, NE, NOP, MDF, NOPP, s_mode, BN, EOI, load, FNOP, DNOP, NEr, NEX, NEY
  use data, only: NV=>NVar, NN=>NTN, iBW, rNOP, NE=>NTE, NOP=>globalNM, MDF, NOPP, s_mode, load=>dsol, NEX, NEY, ths, bas


  !Search for "USER_SPECIFIED" to find parts in need of modification
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!USER_SPECIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !ths: Number of threads, parameter in the kind file (must make clean(wipe)  when changing)

  !iBW: Bandwidth (warning will be changed by init_front

  !s_mode: Set to 0 if solving full dynamics and 1 if only mesh 
  !(must be zero if solve mesh by putting BC's on all velocity and pressure vars)
  !bas: number of nodes in element (9)




  !LIST OF THINGS TO DO
  ! make sure rNOP exist and is filled
  ! make the assemble_local subroutine
  ! make sure to initialize the frontal solver before using


  
  !rNOP(NN,4,2): Reverse NOP, create if you don't have it using below code
  


  !REQUIRED SUBROUTINES FROM YOU
  !Must have a callable assembly subroutine:
  !call assemble_local(ele,local,loc,NOPPl,NB,id,dum)
  !ele: element to be assembled (DONT CHANGE)
  !local(NB,NB): local matrix (Warning assemble as (j,i), we normally assemble i,j as i = row, j = col)
  !NB: Integer, defines sizes (# number of variables in element (DONT CHANGE)
  !loc(NB): local rhs
  !NOPPl(bas): Local NOPPl, mut be filled and returned
  !id: Thread id, integer (DONT CHANGE)
  !dum: integer no purpose (DONT CHANGE)

  !!!!!!!!!!!!!!!!!!!!USER_SPECIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=rk), allocatable :: loadc(:), LT(:,:), load_dum(:)
  integer(kind=ik), allocatable ::  IT(:,:), lt_i(:), NODf(:,:), NODf_n(:,:,:), NOA(:), NCN(:),&
       ele_list(:,:), lt_i2(:), dm_list(:,:), dm_type(:), Etest(:), eprs(:,:,:), dprs(:,:,:)
  integer(kind=ik) :: NB, LT_l, m2, m1, mDM, rad_ele, z_ele, iDM, ths2, ths3, ptn_g, ptn2_g, NSi,&
       prm, prm2, max_prs, piv(1), s_frac, f_cut   !?piv(ths)
  real(kind=rk) :: tr

  !LT: Stores the pivotal rows
  !IT: Stores the information on how to use LT to back substitute
  !loadc: b vector used during solution, but not where the answer goes
  !iPg: global pointers i
  !jPg: global pointers j
  !rhead: Where the local matrix goes rowwise
  !chead: Where the local goes column wise
  !NB: Number of variables in a local block
  !front: the front
  !Lt_l: Max size of any given fully summed pivotal row to be held
  
  !Layout of sub domains
  !Arrows indicate the direction that elimination  occurs in
  !Middle domains may also be subdivided into smaller domains and 
  !eliminated via nested dissection
  !----------------------------------------------------------------------------------------------|
  !               X|--->      |           |            |            |--->        |       * <-----|
  !                |    7...  |           |            |            |   iDM-1    |               |
  !                |          |           |            |            |            |               |
  !                |----------|-----------|------------|------------|------------|               |
  !                |          |           |            |            |            |               |
  !                |    6     |           |            |            |            |               |
  !                |          |           |            |            |            |               |
  !                |----------|---------- |------------|------------|------------| -             |
  !                |          |           |            |            |            | |             |
  !                |    5     |           |            |            |            | |=ptn2        |
  !                |          |           |            |            |            | |             |
  !       d1       |----------|-----------|------------|------------|------------| -             |
  !                |          |           |            |            |            |       iDM     |
  !                |    4     |           |            |            |            |               |
  !                |          |           |            |            |            |               |
  !                |----------|-----------|------------|------------|------------|               |
  !                |          |           |            |            |            |               |
  !                |    3     |           |            |            |            |               |
  !                |          |           |            |            |            |               |
  !                |----------|-----------|------------|------------|------------|               |
  !                |          |           |            |            |            |               |
  !                |    2     |           |            |            |            |               |
  !                |          |           |            |            |            |               |
  !------>*________|Z--->_____|--->______ |____________|____________|____________|Y______________|
  !                                       |------------|
  !                                             =ptn
  !X = end_ele_d1 
  !Y = end_ele_dn
  !Z = first_ele_di
  !*Default setup is lazy domain 1 and iDM numbering (Number from bottom up)
  !Must manully create these ordering for max efficiency






contains
  subroutine init_front(init,ptn,ptn2)
    !init = 1 for first entry, 2 for changing setup after first, 0 for exit
    !ptn = size of sub rows in middle (-1 for overide single, -2 for overide double, -3 optimal nat/nest), -4 drops
    !ptn2 = size of sub columns in middle (if ptn<0 make ptn2 1)

    !Recommended ptn2 = rad_ele/3ish (want three regions close to equal, 
    !so for 28 pick 10 for 2 regions of 10 and 1 of 8)
    !ptn = 4-7, 5 in general but test it out
    !If you intend to use nested ptn = ptn2 is better
    use data, only: RN
    ! use locals, only: RN
    implicit none
    integer(kind=ik) :: init, h, i, j, k, l, m, n, o, p, q, r, g, NNf(NN), last_e, last_n, flag1, l_i,&
         ele, ptn, rp, zp, ptn2, end_ele_d1, end_ele_dn, first_ele_di
    logical :: inc

    !Deallocate arrays
    if (init.ne.1) then
       deallocate(loadc,LT,IT,lt_i,load_dum,NODf,NOA,ele_list,dm_list,lt_i2,dm_type,Etest,NCN)
       call setup_nested(.FALSE.)
    end if
    
    if (init.ge.1) then
       f_cut = 1
       call omp_set_num_threads(ths)
       tr = 0.0_rk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!USER_SPECIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Last element in domain 1 and iDM (remember domain iDM is ordered from the lart element back)
       !Only need for multifront_nat or multifront2
       end_ele_d1 = 10!EOP(2,1)-1
       end_ele_dn = 10!BNOP(2,3,BN(2,3)-2)

       !Should be +1 to end_ele_d1
       first_ele_di = end_ele_d1 + 1

       !Number of radial elements
       rad_ele = NEX!

       !Number of axial elements (Not counting weird regions)
       ! z_ele = 0
       ! do i = 2, 2, 2
       !    z_ele = z_ele + BN(i,3)
       ! end do
       ! z_ele = z_ele - 3
       z_ele = NEY
!
       write(*,*) 'rad, axl:', rad_ele, z_ele
       write(*,*) 'eled1, eledn:', end_ele_d1, end_ele_dn
!!!!!!!!!!!!!!!!!!!!!!!!!USER_SPECIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       ptn_g = ptn
       ptn2_g = ptn2
       if (ptn_g.eq.-3) then
          ptn2_g = rad_ele/2
          !if (mod(rad_ele,2).ne.0) then
             do
                ptn2_g = ptn2_g - 1
                write(*,*) ptn2_g, rad_ele/ptn2_g, mod(rad_ele,ptn2_g)
                if (rad_ele/ptn2_g.ge.3) exit
             end do
             if (mod(rad_ele,3).ne.0) ptn2_g = ptn2_g + 1
          !end if
          ptn_g = (rad_ele - 2* ptn2_g)/2 + mod((rad_ele - 2* ptn2_g),2)
          if (mod((rad_ele - 2* ptn2_g),2).eq.0) ptn_g = ptn_g + 1
          write(*,*) "SELECTED PTN'S FOR MODE 1:", ptn_g, ptn2_g
       end if

!
      



       !Calculate frontwidth needed
       if (ptn_g.ne.-2.and.ptn_g.ne.-1.and.ptn_g.ne.-4) then
          
          if (s_mode.eq.1) then
             iBW = (2*(rad_ele+2*ptn_g+2*ptn2_g)+1)*2
          else      !r,z,u,v                      !pa                      !pb
             iBW = ((2*(rad_ele+2*ptn_g+2*ptn2_g)+1)*4 + (rad_ele+2*ptn_g+ptn2_g+1) + 1 + ptn_g + 1)
          end if
       end if
       !write(*,*) 'Max Front width:', iBW

       !Number of vars in an element
       allocate(NCN(NE))
       NCN(:) = 0
       do j = 1, NE, 1
          do i = 1, bas, 1
             NCN(j) = NCN(j) + MDF(NOP(j,i))

          end do
       end do

       !Save partition numbers and max subfront
       if (ptn_g.ge.ptn2_g) then
          if (mod(ptn_g,2).ne.0) then
             NSi = (ptn_g+1)/2
          else
             NSi = ptn_g/2
          end if
       else
          if (mod(ptn2_g,2).ne.0) then
             NSi = (ptn2_g+1)/2
          else
             NSi = ptn2_g/2
          end if
       end if
       NSi = NSi*NSi

       !Find max NCN => NB
       NB = 0
       do i = 1, NE, 1
          if (NCN(i).gt.NB) NB = INT(NCN(i),ik)
       end do

       !Size of radial and axial partitions
       rp = rad_ele/ptn2_g
       zp = z_ele/ptn_g
       
       LT_l = iBW
       !Determine number of domains
       if (mod(rad_ele,ptn2_g).ne.0) rp = rp + 1
       iDM = rp*zp+2
       if (mod(z_ele,ptn_g).ne.0) then
          iDM = iDM + rp
       end if
       if (mod(rad_ele,ptn2_g).ne.0) rp = rp - 1
      
       !Override (2 for dual domain, 1 for single domain)
       if (ptn_g.eq.-1) then
          iDM = 1
       else if (ptn_g.eq.-2) then
          iDM = 2
       else if (ptn_g.eq.-4) then
          iDM = 2*(RN/4+1)
          if (iDM.eq.2) then
             ptn_g = -2
          end if
       end if
       
       if (ptn_g.lt.0) then
          NSi = 1
       end if

       allocate(loadc(NV),load_dum(NV),LT(NV,LT_l),IT(NV,3+LT_l),lt_i(iDM+1),dm_type(iDM),Etest(NE),&
            NODf(iDM+1,NV),NOA(NV),ele_list(NE+1,iDM),dm_list(iDM,2))
       if (ptn.ne.-4) then
          allocate(lt_i2(3))
       else
          allocate(lt_i2(iDM/2+1))
       end if
       
       !Only use two threads if 2 available
       if (ths.gt.1) then
          ths2 = 2
          ths3 = (ths - 2)/2
          if (ths3.lt.1) ths3 = 1
       else
          ths2 = 1
          ths3 = 1
       end if

       if (ths3*2+2.gt.ths) write(*,*) 'WARNING NESTED THREADS EXCEEDS TOTAL (multifront2)'
     
       
       ele_list = 0
       if (iDM.eq.1) then !Single domain standard order
          do i = 1, NE, 1
             ele_list(i,1) = i
          end do
       else if (iDM.eq.2) then !Dual domain standard order, and reversed standard
          m = 0
          do i = 1, NE, 1
             call exclude_check(i,inc)
             if (inc) m = m + 1
          end do

          
          k = 0
          i = 1
          do
             
             call exclude_check(i,inc)
             if (inc) then
                k = k + 1
                ele_list(k,1) = i
             end if
             i = i + 1
             if (k.eq.m/2) exit
          end do
          
          k = 0
          i = NE
          do 
             call exclude_check(i,inc)
             if (inc) then
                k = k + 1
                ele_list(k,2) = i
             end if
             i = i - 1
             if (k+m/2.eq.m) exit
          end do
       else if (ptn.eq.-4) then
       write(*,*) 'error using drop version' 
          ! j = 1
          ! do i = 1, iDM-1, 2
          
          !    m = 1!EOP(j,BN(j,1))    !start
          !    if (j+3.le.RN) then
          !       n = EOP(j+3,NEr(j+3)) !end
          !    else
          !       n = EOP(j+2,NEr(j+2))
          !    end if
             
          !    o = (n-m)/2 + m         !middle
          !    k = 0
          !    do p = m, o, 1
          !       k = k + 1
          !       call exclude_check(p,inc)
          !       if (inc) ele_list(k,i) = p
          !    end do

          !    k = 0
          !    do p = n, o+1, -1
          !       k = k + 1
          !       call exclude_check(p,inc)
          !       if (inc) ele_list(k,i+1) = p
          !    end do

          !    !write(*,*) m, n
          !    !write(*,*) j, j+3, i, i+1
          !    !pause
          !    j = j + 4
          ! end do
       else
             
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!USER_SPECIFIED!!!!!!!!!!!!!PROBABLY SKIP
          
          ! Lazy Way, leave as is or create your own order
          k = 0
          do i = 1, end_ele_d1, 1
             k = k + 1 
             ele_list(k,1) = i
          end do
          
          k = 0
          do i = NE, end_ele_dn, -1
             k = k + 1
             ele_list(k,iDM) = i
          end do
          
          !Non-lazy example
          !First region order natural, or switch out for above for lazy way
          ! k = 0
          ! o = 15
          ! do i = 1, BN(o,3), 1
          !    do j = BN(o,1), 1, -1
          !       k = k + 1
          !       ele_list(k,1) = EOI(o,i,j)
          !    end do
          ! end do
          ! o = 14
          ! do i = 1, BN(o,3), 1
          !    do j = BN(o,1), 1, -1
          !       k = k + 1
          !       ele_list(k,1) = EOI(o,i,j)
          !    end do
          ! end do
          ! o = 1
          ! do i = 1, BN(o,3), 1
          !    do j = BN(o,1), 1, -1
          !       k = k + 1
          !       ele_list(k,1) = EOI(o,i,j)
          !    end do
          ! end do
          
         
          ! !First region order natural, or switch out for above for lazy way
          ! k = 0
          ! o = 12
          ! do i = BN(o,3), 1, -1
          !    do j = 1, BN(o,1), 1
          !       k = k + 1
          !       ele_list(k,iDM) = EOI(o,i,j)
          !    end do
          ! end do
          ! o = 13
          ! do i = BN(o,3), 1, -1
          !    do j = 1, BN(o,1), 1
          !       k = k + 1
          !       ele_list(k,iDM) = EOI(o,i,j)
          !    end do
          ! end do
          ! o = 26
          ! do i = BN(o,3), 1, -1
          !    do j = 1, BN(o,1), 1
          !       k = k + 1
          !       ele_list(k,iDM) = EOI(o,i,j)
          !    end do
          ! end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!USER_SPECIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          
          !If even division into 2 not possible, make sure it is 2's leftover 1
          !Not 1's leftover huge number
          m = mod(z_ele,ptn_g)
          r = mod(rad_ele,ptn2_g)
          if (r.ne.0) then
             l = 1
          else
             l = 0
          end if
         
          !Size of remainder partitions 
          prm = m
          prm2 = r
          

          !k local count indice in domain
          !j starting element
          !i domain number
          !ptn_g partition number
          !o loop upwards axial
          !n loop sideways radially
          !iDM total # of domains
          !if remainder in rp l = 1 else 0
          !if remainder in rp r = remainder else 0
          !if remainder in zp m = remainder else 0
          !rp total number of cuts upwards (not including additionally remainder l)
          !zp total number of cuts sideways (not including remainder)
          !g keeps track of jumps up in current domain
          !q keeps track of jumps right in current domain
          k = 0
          j = first_ele_di-1
          q = 0
          i = 2
          
          do 


             k = 0
             g = 0
             do p = 1, rp+l, 1
                !write(*,*) '1 Filling domain:', i
                if (p.lt.rp) then
                   do o = 1, ptn2_g, 1
                      do n = 1, ptn_g, 1
                         k = k + 1
                         ele_list(k,i) = grid(q+n,o+g,j,rad_ele)
                      end do

                   end do
                   dm_type(i) = 1
                   i = i + 1
                   k = 0
                   g = g + ptn2_g
                else if (p.eq.rp.and.r.eq.0) then

                   do o = ptn2_g, 1, -1
                      do n = 1, ptn_g, 1
                         k = k + 1
                         ele_list(k,i) = grid(q+n,o+g,j,rad_ele)
                      end do

                   end do
                   dm_type(i) = 1
                   i = i + 1
                   k = 0
                   g = g + ptn2_g
                else if (p.eq.rp.and.r.ne.0) then
                   do o = 1, ptn2_g, 1
                      do n = 1, ptn_g, 1
                         k = k + 1
                         ele_list(k,i) = grid(q+n,o+g,j,rad_ele)
                      end do

                   end do
                   dm_type(i) = 1
                   i = i + 1
                   k = 0
                   g = g + ptn2_g
                else
                   do o = r, 1, -1
                      do n = 1, ptn_g, 1
                         k = k + 1
                         ele_list(k,i) = grid(q+n,o+g,j,rad_ele)
                      end do

                   end do
                   dm_type(i) = 2 !2 is short ptn2_g
                   i = i + 1
                   k = 0
                   g = g + r
                end if
                !write(*,*) 'g:', g


             end do

             q = q + ptn_g
             !write(*,*) 'q:', q
             !pause
             !o = o + 1
             if (((i-2)/(rp+l)).eq.zp) exit
          end do
          
          if (m.ne.0) then


             k = 0
             g = 0
             do p = 1, rp+l, 1
                !write(*,*) '2 Filling domain:', i
                if (p.lt.rp) then
                   do o = 1, ptn2_g, 1
                      do n = 1, m, 1
                         k = k + 1
                         ele_list(k,i) = grid(q+n,o+g,j,rad_ele)
                      end do

                   end do
                   dm_type(i) = 3 !3 is short ptn_g
                   i = i + 1
                   k = 0
                   g = g + ptn2_g
                else if (p.eq.rp.and.r.eq.0) then

                   do o = ptn2_g, 1, -1
                      do n = 1, m, 1
                         k = k + 1
                         ele_list(k,i) = grid(q+n,o+g,j,rad_ele)
                      end do

                   end do
                   dm_type(i) = 3
                   i = i + 1
                   k = 0
                   g = g + ptn2_g
                else if (p.eq.rp.and.r.ne.0) then
                   do o = 1, ptn2_g, 1
                      do n = 1, m, 1
                         k = k + 1
                         ele_list(k,i) = grid(q+n,o+g,j,rad_ele)
                      end do

                   end do
                   dm_type(i) = 3
                   i = i + 1
                   k = 0
                   g = g + ptn2_g
                else
                   do o = r, 1, -1
                      do n = 1, m, 1
                         k = k + 1
                         ele_list(k,i) = grid(q+n,o+g,j,rad_ele)
                      end do

                   end do
                   dm_type(i) = 4 !4 is short both
                   i = i + 1
                   k = 0
                   g = g + r
                end if

                !write(*,*) 'g:', g

                !write(*,*) 'q:', q

             end do
             q = q + m
          end if
          !pause
          do i = 1, NE, 1
             o = 0
             do l = 1, iDM, 1
                k = 0
                do
                   k = k + 1
                   if (ele_list(k,l).eq.0) exit
                   if (ele_list(k,l).eq.i) then
                      o = o + 1
                   end if
                end do
             end do

             if (o.ne.1) then

                write(*,*) 'Error in list:', o, i
             end if
          end do


          !Domain list (Used for level 2 elimination two fronts)
          dm_list = 0
          do i = 1, iDM/2, 1
             dm_list(i,1) = i
          end do

          k = 0
          do i = iDM, iDM/2+1, -1
             k = k + 1
             dm_list(k,2) = i
          end do


          ! if (m.ne.0) then
          !    i = iDM-2
          !    k = 0
          !    h = 0
          !    do n = 1, m, 1
          !       do o = 1, rad_ele/2+mod(rad_ele,2), 1
          !          j = j + 1
          !          k = k + 1
          !          ele_list(k,i) = j
          !       end do
          !       do o = 1, rad_ele/2, 1
          !          j = j + 1
          !          h = h + 1
          !          ele_list(h,i+1) = j
          !       end do
          !    end do
          ! end if

          ! k = 0
          ! j = BNOP(2,1,1)-1
          ! do i = 2, iDM - 2 - 2*m, 2
          !    k = 0
          !    h = 0
          !    do n = 1, l, 1
          !       do o = 1, rad_ele/2+mod(rad_ele,2), 1
          !          j = j + 1
          !          k = k + 1
          !          ele_list(k,i) = j
          !       end do
          !       do o = 1, rad_ele/2, 1
          !          j = j + 1
          !          h = h + 1
          !          ele_list(h,i+1) = j
          !       end do
          !    end do
          ! end do

          ! if (m.ne.0) then
          !    i = iDM-2
          !    k = 0
          !    h = 0
          !    do n = 1, m, 1
          !       do o = 1, rad_ele/2+mod(rad_ele,2), 1
          !          j = j + 1
          !          k = k + 1
          !          ele_list(k,i) = j
          !       end do
          !       do o = 1, rad_ele/2, 1
          !          j = j + 1
          !          h = h + 1
          !          ele_list(h,i+1) = j
          !       end do
          !    end do
          ! end if
          !write(*,*) j, BNOP(12,1,1)
          !pause
       end if

       !NOA stores number of appearance per variable
       NOA(:) = 0
       NODf = 0
       do i = 1, NN, 1
          do j = 1, 4, 1
             if (rNOP(i,j,1).eq.0) cycle
             call exclude_check(rNOP(i,j,1),inc)
             if (.not.inc) cycle
             do k = 0, MDF(i)-1, 1
                NOA(NOPP(i)+k) = NOA(NOPP(i)+k) + 1
             end do
          end do
       end do

      

       !Use NODf to find the amount of variables that will be
       !deleted in each level and sublevels
       if (iDM.gt.1) then
          !Determine l_i offsets LEVEL 1
          l_i = 0
          lt_i(1) = 0
          do l = 1, iDM, 1
             j = 0
             do 
                j = j + 1
                ele = ele_list(j,l)
                if (ele.eq.0) exit
                do i = 1, 9, 1

                   o = NOP(ele,i)
                   p = NOPP(o)
                   !Add new vars
                   NODf(l,p) = NODf(l,p) + 1
                   do k = 1, MDF(o)-1, 1
                      NODf(l,p+k) = NODf(l,p+k) + 1

                   end do

                   !Mark and count fully summed
                   if (NODf(l,p).eq.NOA(p)) then
                      l_i = l_i + MDF(o)
                      do k = 0, MDF(o)-1, 1
                         NODf(l,p+k) = -NODf(l,p+k)
                      end do
                   end if

                end do
             end do
             lt_i(l+1) = l_i
             !write(*,*) 'l_i1:', l_i
          end do

          !Determine l_i offsets LEVEL 2
          if (iDM.gt.2.and.ptn.ne.-4) then
             lt_i2(1) = lt_i(iDM+1)




             l = 1
             j = 1

             do
                j = j + 1
                ele = dm_list(j,l)
                if (ele.eq.0) exit
                NODf(ele,:) = NODf(ele,:) + NODF(ele-1,:)
                do i = 1, NV, 1
                   if (NODf(ele,i).eq.NOA(i)) then
                      l_i = l_i + 1
                      NODf(ele,i) = -NODf(ele,i)
                   end if
                end do
                !write(*,*) 'l_i:', l_i
             end do
             lt_i2(2) = l_i
             !write(*,*) 'l_i2:', l_i
             j = 1
             l = 2
             do
                j = j + 1
                ele = dm_list(j,l)
                if (ele.eq.0) exit
                NODf(ele,:) = NODf(ele,:) + NODF(ele+1,:)
                do i = 1, NV, 1
                   if (NODf(ele,i).eq.NOA(i)) then
                      l_i = l_i + 1
                      NODf(ele,i) = -NODf(ele,i)
                   end if
                end do
             end do
             lt_i2(3) = l_i
             !write(*,*) 'l_i2:', l_i
          end if
       end if

       if (ptn_g.eq.-4) then
          lt_i2(1) = lt_i(iDM+1)

          j = 1
          do k = 1, iDM-1, 2
             j = j + 1
             NODf(k+1,:) = NODf(k+1,:) + NODF(k,:)
             do i = 1, NV, 1
                if (NODf(k+1,i).eq.NOA(i)) then
                   l_i = l_i + 1
                   NODf(k+1,i) = -NODf(k+1,i)
                end if
             end do
             lt_i2(j) = l_i
          end do
          !write(*,*) l_i, NV
          !pause
       end if

      

       !Setup arrays for nested ordering
       !if (iDM.gt.2) then
       call setup_nested(.TRUE.)
       ! else
       !    allocate(eprs(1,1,1),dprs(1,1,1))
       ! end if
       !write(*,*) 'SETUP DONE'

       call ETESTer()
    end if


    
  end subroutine init_front

  subroutine exclude_check(ele,inc)
    ! use locals, only: exc_list
    implicit none
    integer(kind=ik) :: ele
    logical inc

  
    
    ! if (exc_list(DNOP(ele)).eq.1) then
       inc = .TRUE.
    ! else
    !    inc = .FALSE.
    ! end if      
    !

  end subroutine exclude_check

  subroutine ETESTer()
    implicit none
    integer(kind=ik) :: E(NE), i,j, ele
    logical :: inc
    E = 0
    do ele = 1, NE, 1
       do i = 1, iDM, 1
          do j = 1, NE, 1
             if (ele_list(j,i).eq.0) exit
             if (ele.eq.ele_list(j,i)) E(ele) = E(ele)+1
          end do
       end do
    end do

    do ele = 1, NE, 1
       if (E(ele).ne.1) then
          call exclude_check(ele,inc)
          if (ele.eq.0.and.inc) then
             write(*,*) 'failed', ele,'with',E(ele)
          end if
          if (E(ele).gt.1) then
             write(*,*) 'failed', ele,'with',E(ele)
          end if
          !pause
       end if
    end do
    write(*,*) 'ETEST DONE'

  end subroutine ETESTER

  integer(kind=ik) function grid(i,j,s_ele,r_ele)
    implicit none
    integer(kind=ik) :: i, j, s_ele, r_ele

    grid = s_ele + (i-1)*(r_ele) + j
    return
    
  end function grid

  subroutine multifront_single(L2res,t)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose

    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    implicit none
    real(kind=rk) :: front_m(iBW,iBW), loadf_m(iBW), loadf_dumm(iBW), L2res, t1, t
    integer(kind=ik) :: e_n(4), iPg_m(iBW), jPg_m(iBW), fs_m, i, j, l_i, id

    if (iDM.ne.1) then
       write(*,*) 'Error iDM != 1'
       return
    end if
    !mode = 1 is nested natural ordering
    !mode = 2 is nested dissection (ptn should equal ptn2)
    t1 = REAL(omp_get_wtime(),rk)
    piv = 0

    
    !Blank global arrays
    loadc = 0.0_rk      !working array
    load_dum = 0.0_rk   !dummy right hand side
    load = 0.0_rk       !Solution vector

    front_m(:,:) = 0.0_rk  !Single front matrix
    loadf_m(:) = 0.0_rk    !single right hand
    loadf_dumm(:) = 0.0_rk !single dummy
    iPg_m = 0              !local to global rows
    jPg_m = 0              !local to glbal col
    fs_m = 0               !Size of front
    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances



    !Single front
    call Front_solve(fs_m,front_m,loadf_m,loadf_dumm,iPg_m,jPg_m,l_i,1,1,1,e_n)

    !Perform back substitution
    call back_sub(1)

    !Calculate L2res
    L2res = 0.0_rk
    do i = 1, NV, 1
       L2res = L2res + (load_dum(i))**2
    end do
    L2res = sqrt(L2res)

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - t1
    j = 0
    do i = 1, ths, 1
       j = piv(i) + j
    end do
    write(*,*) 'Small Pivots:', j

     
    
  end subroutine multifront_single

  subroutine multifront_dual(L2res,t)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose

    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    implicit none
    real(kind=rk) :: front(iBW,iBW,iDM),  loadf(iBW,iDM), loadf_dum(iBW,iDM), L2res, t1, t
    integer(kind=ik) :: iPg(iBW,iDM), jPg(iBW,iDM), fs(iDM), e_n(4), i, j, k, l_i, id, s

    if (iDM.ne.2) then
       write(*,*) 'Error iDM != 2'
       return
    end if
    !mode = 1 is nested natural ordering
    !mode = 2 is nested dissection (ptn should equal ptn2)
    t1 = REAL(omp_get_wtime(),rk)
    piv = 0

    
    !Blank global arrays
    IT = 0
    loadc = 0.0_rk      !working array
    load_dum = 0.0_rk   !dummy right hand side
    load = 0.0_rk       !Solution vector
    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances

    front = 0.0_rk
    loadf = 0.0_rk
    loadf_dum = 0.0_rk
    iPg = 0
    jPg = 0
    fs = 0
      
    !PHASE 1: Solve each half
    !$omp parallel private(id,j,l_i) num_threads(ths2)
    !$omp do
    do j = 1, iDM, 1
       id = omp_get_thread_num()+1  !Thread id used to determine assembly arrays

       l_i = lt_i(j)                !initilaize var counter to predetermined starting point


       !Solve domain j and save front
       call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n)

    end do
    !$omp end do
    !$omp end parallel
  

    !PHASE 3: Stitch halves
    l_i = lt_i(3)
    call add_front(front(:,:,iDM),front(:,:,1),loadf(:,iDM),loadf(:,1),loadf_dum(:,iDM),&
         loadf_dum(:,1),iPg(:,iDM),iPg(:,1),jPg(:,iDM),jPg(:,1),fs(iDM),fs(1),iDM,1,1,1)
    call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
         iPg(:,1),jPg(:,1),l_i,fs(1),id)

    write(*,*) 'fs:', fs(1), l_i

    !Perform back substitution
    call back_sub(1)

    !Calculate L2res
    L2res = 0.0_rk
    !$omp parallel private(i) reduction(+:L2res)
    !$omp do
    do i = 1, NV, 1
       L2res = L2res + (load_dum(i))**2
    end do
    !$omp end do
    !$omp end parallel
    L2res = sqrt(L2res)

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - t1
    j = 0
    do i = 1, ths, 1
       j = piv(i) + j
    end do
    write(*,*) 'Small Pivots:', j

     
    
  end subroutine multifront_dual


 subroutine multifront_drops(L2res,t)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose

    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    implicit none
    real(kind=rk) :: front(iBW,iBW,iDM),  loadf(iBW,iDM), loadf_dum(iBW,iDM), L2res, t1, t
    integer(kind=ik) :: iPg(iBW,iDM), jPg(iBW,iDM), fs(iDM), e_n(4), i, j, k, l_i, id, s

    ! if (iDM.ne.2) then
    !    write(*,*) 'Error iDM != 2'
    !    return
    ! end if
    !mode = 1 is nested natural ordering
    !mode = 2 is nested dissection (ptn should equal ptn2)
    t1 = REAL(omp_get_wtime(),rk)
    piv = 0

    
    !Blank global arrays
    IT = 0
    loadc = 0.0_rk      !working array
    load_dum = 0.0_rk   !dummy right hand side
    load = 0.0_rk       !Solution vector
    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances

    front = 0.0_rk
    loadf = 0.0_rk
    loadf_dum = 0.0_rk
    iPg = 0
    jPg = 0
    fs = 0
      
    !PHASE 1: Solve each half
    !$omp parallel private(id,j,l_i) num_threads(ths)
    !$omp do
    do j = 1, iDM, 1
       id = omp_get_thread_num()+1  !Thread id used to determine assembly arrays

       l_i = lt_i(j)                !initilaize var counter to predetermined starting point


       !Solve domain j and save front
       call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n)

    end do
    !$omp end do
    !$omp end parallel
  

    !PHASE 3: Stitch halves
    !$omp parallel private(id,i,l_i) num_threads(ths)
    !$omp do
    do i = 1, iDM-1, 2
       l_i = lt_i2((i+1)/2)
       id =  omp_get_thread_num()+1
       call add_front(front(:,:,i+1),front(:,:,i),loadf(:,i+1),loadf(:,i),loadf_dum(:,i+1),&
            loadf_dum(:,i),iPg(:,i+1),iPg(:,i),jPg(:,i+1),jPg(:,i),fs(i+1),fs(i),i+1,i,1,1)
       call elim_stitch(front(:,:,i),loadf(:,i),loadf_dum(:,i),&
            iPg(:,i),jPg(:,i),l_i,fs(i),id)
    end do
    !$omp end do
    !$omp end parallel
       
    !write(*,*) 'fs:', fs(:), l_i

    !Perform back substitution
    call back_sub(1)

    !Calculate L2res
    L2res = 0.0_rk
    !$omp parallel private(i) reduction(+:L2res)
    !$omp do
    do i = 1, NV, 1
       L2res = L2res + (load_dum(i))**2
    end do
    !$omp end do
    !$omp end parallel
    L2res = sqrt(L2res)

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - t1
    j = 0
    do i = 1, ths, 1
       j = piv(i) + j
    end do
    write(*,*) 'Small Pivots:', j

     
    
  end subroutine multifront_drops

  
  subroutine multifront_nat(L2res,t)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose

    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    implicit none
    real(kind=rk) :: front(iBW,iBW,iDM), loadf(iBW,iDM), loadf_dum(iBW,iDM), L2res, t1, t
    integer(kind=ik) :: iPg(iBW,iDM), jPg(iBW,iDM), fs(iDM), mode, e_n(4), i, j, k, l_i, id

    if (iDM.lt.3) then
       write(*,*) 'Error iDM lt 3'
       return
    end if
    
    !mode = 1 is nested natural ordering
    !mode = 2 is nested dissection (ptn should equal ptn2)
    t1 = REAL(omp_get_wtime(),rk)
    piv = 0

    
    !Blank global arrays
    loadc = 0.0_rk      !working array
    load_dum = 0.0_rk   !dummy right hand side
    load = 0.0_rk       !Solution vector

                           !Size of front
    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances

    !Blank arrays for front that are domain wise
    front = 0.0_rk
    loadf = 0.0_rk
    loadf_dum = 0.0_rk
    iPg = 0
    jPg = 0
    fs = 0

    !PHASE 1: All smallest level Natural done using arbitray threads!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    !$omp parallel private(id,j,l_i)
    !$omp do
    do j = 1, iDM, 1
       id = omp_get_thread_num()+1  !Thread id used to determine assembly arrays

       l_i = lt_i(j)                !initilaize var counter to predetermined starting point


       !Solve domain j and save front
       call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n)

    end do
    !$omp end do
    !$omp end parallel


    !PHASE 2: Sweep through domains from two fronts in natural order of domains (2 threads at most)!!!!!!!!!!!!!
    !Not required for single domain front or dual domain front

    !$omp parallel sections private(j,k,l_i,i,id) num_threads(ths2)
    !$omp section
    id = omp_get_thread_num() + 1
    k = 1
    j = 1
    l_i = lt_i2(j)


    do 
       k = k + 1
       i = dm_list(k,j)
       if (i.eq.0) exit

       call add_front(front(:,:,i),front(:,:,1),loadf(:,i),loadf(:,1),loadf_dum(:,i),&
            loadf_dum(:,1),iPg(:,i),iPg(:,1),jPg(:,i),jPg(:,1),fs(i),fs(1),i,1,1,1)
       call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
            iPg(:,1),jPg(:,1),l_i,fs(1),id)

    end do
    !$omp section
    id = omp_get_thread_num() + 1
    k = 1
    j = 2
    l_i = lt_i2(j)

    do 
       k = k + 1
       i = dm_list(k,j)
       if (i.eq.0) exit

       call add_front(front(:,:,i),front(:,:,iDM),loadf(:,i),loadf(:,iDM),loadf_dum(:,i),&
            loadf_dum(:,iDM),iPg(:,i),iPg(:,iDM),jPg(:,i),jPg(:,iDM),fs(i),fs(iDM),i,iDM,1,1)
       call elim_stitch(front(:,:,iDM),loadf(:,iDM),loadf_dum(:,iDM),&
            iPg(:,iDM),jPg(:,iDM),l_i,fs(iDM),id)
    end do
    !$omp end parallel sections

    
   
    !PHASE 3: Stitch and eliminate final fronts!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    l_i = lt_i2(3)
    id = 1
    call add_front(front(:,:,iDM),front(:,:,1),loadf(:,iDM),loadf(:,1),loadf_dum(:,iDM),&
         loadf_dum(:,1),iPg(:,iDM),iPg(:,1),jPg(:,iDM),jPg(:,1),fs(iDM),fs(1),iDM,1,1,1)
    call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
         iPg(:,1),jPg(:,1),l_i,fs(1),id)

    !Perform back substitution
    call back_sub(2)

    !Calculate L2res
    L2res = 0.0_rk
    !$omp parallel private(i) reduction(+:L2res)
    !$omp do
    do i = 1, NV, 1
       L2res = L2res + (load_dum(i))**2
    end do
    !$omp end do
    !$omp end parallel
    L2res = sqrt(L2res)

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - t1
    j = 0
    do i = 1, ths, 1
       j = piv(i) + j
    end do
    write(*,*) 'Small Pivots:', j

     
    
  end subroutine multifront_nat

  subroutine multifront_nest(L2res,t)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose
    
    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    implicit none
    real(kind=rk) :: front(iBW,iBW,iDM), loadf(iBW,iDM),&
         tp, ts, ti, loadf_dum(iBW,iDM), L2res, t1, t, front_n(iBW,iBW,NSi,ths),&
         loadf_n(iBW,NSi,ths), loadf_dumn(iBW,NSi,ths)
    integer(kind=ik) :: iPg(iBW,iDM), jPg(iBW,iDM), fs(iDM), mode, &
         e_n(4), i, j, k, l, m, n, l_i, id,&
    iPg_n(iBW,NSi,ths), jPg_n(iBW,NSi,ths), fs_n(NSi,ths), o, p
    write(*,*) 'Nested Broken'
    pause
    return
    
    if (iDM.lt.3) then
       write(*,*) 'Error iDM lt 3'
       return
    end if

    mode = 2
    t1 = REAL(omp_get_wtime(),rk)
    piv = 0

    
    !Blank global arrays
    loadc = 0.0_rk      !working array
    load_dum = 0.0_rk   !dummy right hand side
    load = 0.0_rk       !Solution vector

    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances

    !Blank arrays for front that are domain wise
    front = 0.0_rk
    loadf = 0.0_rk
    loadf_dum = 0.0_rk
    iPg = 0
    jPg = 0
    fs = 0

    !Phase 1.1: Outer two domains not nested order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    l_i = 0
    !$omp parallel private(id,j,l_i)
    !$omp do
    do j = 1, iDM, iDM-1
       id = omp_get_thread_num()+1

       fs(j) = 0
       l_i = lt_i(j)
       !write(*,*) 'Solving Domain:', j, fs(j)
       !write(*,*) 'Hi from thread 2:', omp_get_thread_num()
       call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n)
       !write(*,*) 'Solved Domain:', j, fs(j)
    end do
    !$omp end do
    !$omp end parallel


    !PHASE 1.2: All smallest level Nested!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    !$omp parallel private(id,j,l,l_i,m,o,p)
    !$omp do
    do j = 2, iDM-1, 1
       id = omp_get_thread_num()+1
       l_i = lt_i(j)
       !Initialize nested fronts within domains
       front_n(:,:,:,id) = 0.0_rk
       loadf_n(:,:,id) = 0.0_rk
       loadf_dumn(:,:,id) = 0.0_rk
       iPg_n(:,:,id) = 0
       jPg_n(:,:,id) = 0
       fs_n(:,id) = 0
       NODf_n(:,:,id) = 0

       !Determine domain type (shape, most are type 1)
       l = dm_type(j)

       !Assemble all smallest level nested regions (consist of 2x2 elements)
       do m = 1, NSi, 1
          if (eprs(1,m,j).eq.0) exit
          !write(*,*) eprs(:,m,j)
          if (m.eq.1) then
             !write(*,*) 'hi', fs(j)
             call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),&
                  jPg(:,j),l_i,id,j,3,eprs(:,m,j))
             !write(*,*) 'bye'
          else
             call Front_solve(fs_n(m,id),front_n(:,:,m,id),loadf_n(:,m,id),loadf_dumn(:,m,id),iPg_n(:,m,id),&
                  jPg_n(:,m,id),l_i,id,m,mode,eprs(:,m,j))
          end if
       end do

       !Use nested dissection to combine sub domains into 1 final level
       !Once done the remaining front resides in the domain front
       do m = 1, max_prs, 1
          o = dprs(1,m,l)
          p = dprs(2,m,l)
          if (o.eq.0) exit
          if (o.eq.1) then
             call add_front(front_n(:,:,p,id),front(:,:,j),&
                  loadf_n(:,p,id),loadf(:,j),&
                  loadf_dumn(:,p,id),loadf_dum(:,j),&
                  iPg_n(:,p,id),iPg(:,j),&
                  jPg_n(:,p,id),jPg(:,j),fs_n(p,id),fs(j),&
                  p,j,3,id)
             call elim_stitch(front(:,:,j),loadf(:,j),&
                  loadf_dum(:,j),&
                  iPg(:,j),jPg(:,j),l_i,fs(j),id)
          else
             call add_front(front_n(:,:,p,id),front_n(:,:,o,id),&
                  loadf_n(:,p,id),loadf_n(:,o,id),&
                  loadf_dumn(:,p,id),loadf_dumn(:,o,id),&
                  iPg_n(:,p,id),iPg_n(:,o,id),&
                  jPg_n(:,p,id),jPg_n(:,o,id),fs_n(p,id),fs_n(o,id),&
                  p,o,mode,id)
             call elim_stitch(front_n(:,:,o,id),loadf_n(:,o,id),&
                  loadf_dumn(:,o,id),&
                  iPg_n(:,o,id),jPg_n(:,o,id),l_i,fs_n(o,id),id)
          end if

       end do



    end do
    !$omp end do
    !$omp end parallel





    !PHASE 2: Sweep through domains from two fronts in natural order of domains (2 threads at most)
    !Not required for single domain front or dual domain front

    !$omp parallel sections private(j,k,l_i,i,l,id) num_threads(ths2)

    !$omp section
    id = omp_get_thread_num() + 1
    k = 1
    j = 1
    l_i = lt_i2(j)


    do 
       k = k + 1
       i = dm_list(k,j)
       if (i.eq.0) exit

       call add_front(front(:,:,i),front(:,:,1),loadf(:,i),loadf(:,1),loadf_dum(:,i),&
            loadf_dum(:,1),iPg(:,i),iPg(:,1),jPg(:,i),jPg(:,1),fs(i),fs(1),i,1,1,1)
       call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
            iPg(:,1),jPg(:,1),l_i,fs(1),id)

    end do
    !$omp section
    id = omp_get_thread_num() + 1
    k = 1
    j = 2
    l_i = lt_i2(j)



    do 
       k = k + 1
       i = dm_list(k,j)
       if (i.eq.0) exit

       call add_front(front(:,:,i),front(:,:,iDM),loadf(:,i),loadf(:,iDM),loadf_dum(:,i),&
            loadf_dum(:,iDM),iPg(:,i),iPg(:,iDM),jPg(:,i),jPg(:,iDM),fs(i),fs(iDM),i,iDM,1,1)
       call elim_stitch(front(:,:,iDM),loadf(:,iDM),loadf_dum(:,iDM),&
            iPg(:,iDM),jPg(:,iDM),l_i,fs(iDM),id)
    end do



    !$omp end parallel sections


    l_i = lt_i2(3)
    id = 1
    call add_front(front(:,:,iDM),front(:,:,1),loadf(:,iDM),loadf(:,1),loadf_dum(:,iDM),&
         loadf_dum(:,1),iPg(:,iDM),iPg(:,1),jPg(:,iDM),jPg(:,1),fs(iDM),fs(1),iDM,1,1,1)
    call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
         iPg(:,1),jPg(:,1),l_i,fs(1),id)




    !Perform back substitution
    call back_sub(1)

    !Calculate L2res
    L2res = 0.0_rk
    !$omp parallel private(i) reduction(+:L2res)
    !$omp do
    do i = 1, NV, 1
       L2res = L2res + (load_dum(i))**2
    end do
    !$omp end do
    !$omp end parallel
    L2res = sqrt(L2res)

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - t1
    j = 0
    do i = 1, ths, 1
       j = piv(i) + j
    end do
    write(*,*) 'Small Pivots:', j

     
    
  end subroutine multifront_nest

  subroutine multifront2(L2res,t)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose

    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    implicit none
    real(kind=rk) :: front(iBW,iBW,iDM), loadf(iBW,iDM), loadf_dum(iBW,iDM),&
         L2res, t1, t2, t3, t, ti
    integer(kind=ik) :: iPg(iBW,iDM), jPg(iBW,iDM), fs(iDM), e_n(4,ths), &
         i, j, k, l_i, id, o, p, f1, f2, f3, l_i2, l_i3, s_stop, thsb
    
    if (ths3*2+2.gt.ths) then
       write(*,*) 'WARNING NESTED THREADS EXCEEDS TOTAL'
       return
    end if
    if (iDM.lt.3) then
       write(*,*) 'Error iDM lt 3'
       return
    end if
   
    ti = REAL(omp_get_wtime(),rk)
    piv = 0
    !s_frac = iDM/10
    s_stop = 4
    thsb = 2 + ths3
    
    !Blank global arrays
    loadc = 0.0_rk          !working array
    load_dum = 0.0_rk       !dummy right hand side
    load = 0.0_rk           !Solution vector
           
    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances



    !Single front
    !Dual front and nested fronts
    !Blank arrays for front that are domain wise
    front = 0.0_rk
    loadf = 0.0_rk
    loadf_dum = 0.0_rk
    iPg = 0
    jPg = 0
    fs = 0
    if (tr.eq.0.0_rk) then
       !do o = 1, 2, 1
       t1 = REAL(omp_get_wtime(),rk)
       id = 1
          do j = 1, iDM/4, 1
             !omp_get_thread_num()+1  !Thread id used to determine assembly arrays

             l_i = lt_i(j)                !initilaize var counter to predetermined starting point


             !Solve domain j and save front
             call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n(:,id))

          end do

          t2 = REAL(omp_get_wtime(),rk)



          l_i = lt_i2(1)
          do k = 2, iDM/4, 1

             i = dm_list(k,1)

             call add_front(front(:,:,i),front(:,:,1),loadf(:,i),loadf(:,1),loadf_dum(:,i),&
                  loadf_dum(:,1),iPg(:,i),iPg(:,1),jPg(:,i),jPg(:,1),fs(i),fs(1),i,1,1,1)
             call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
                  iPg(:,1),jPg(:,1),l_i,fs(1),id)

          end do
          t3 = REAL(omp_get_wtime(),rk)
          ! write(*,*) "sec1:", t3 - t2
          ! write(*,*) "sec11:", t2-t1
          tr = (t2-t1)/(t3-t2)*1.0_rk/REAL(ths3,rk)
          write(*,*) 'tr:', tr
          !Blank global arrays
          loadc = 0.0_rk      !working array
          load_dum = 0.0_rk   !dummy right hand side
          load = 0.0_rk       !Solution vector

          l_i = 0                !Variables removed from front
          NODf = 0               !Records variable appearances



          !Single front
          !Dual front and nested fronts
          !Blank arrays for front that are domain wise
          front = 0.0_rk
          loadf = 0.0_rk
          loadf_dum = 0.0_rk
          iPg = 0
          jPg = 0
          fs = 0

          p = 3
          o = 1000
          do i = 1, iDM/4, 1
             k = 1
             f1 = i
             f2 = INT((1.0_rk+1.0_rk/tr)*REAL(f1,rk),ik) + f1
 
             f3 = INT(REAL(f2 - f1,rk)/tr,ik) + f2
                
             do
               
                if (f3.eq.iDM/2) exit
                k = k + 1
                j = f3 - f2
                f1 = f2
                f2 = f3
                f3 = INT(REAL((f2-f1),rk)/tr,ik) + f2
                if (f3.gt.iDM/2) f3 = iDM/2
                if (f2-f1.lt.s_stop) f3 = iDM/2
                !write(*,*) f3
             end do
             write(*,*) i, k, f3-f2, j
             if (j.lt.ths3) cycle
             if (f3-f2.lt.1) cycle
             if (f3-f2-j.lt.-2) cycle
             if (f3-f2+j.lt.o) then
                o = f3-f2+j
                p = i
             end if
                
            
            
          end do
          s_frac = p
          write(*,*) 's_frac selected:', s_frac
    end if
 
    
    
    !Mode 1 is nested natural order
    !$omp parallel sections private(f1,f2,f3,i,j,k,id,l_i,l_i2) num_threads(2)

    !$omp section
    f1 = s_frac
    f2 = INT((1.0_rk+1.0_rk/tr)*REAL(f1,rk),ik) + f1
    f3 = 0
    l_i2 = 0
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PHASE 1.0:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$omp parallel sections private(i,j,k,id,l_i) shared(l_i2,f1,f2,f3) num_threads(2)

    !$omp section

    id = 1
    do j = 1, f1, 1
         !Thread id used to determine assembly arrays

       l_i = lt_i(j)                !initilaize var counter to predetermined starting point


       !Solve domain j and save front
       call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n(:,id))

    end do
   
    l_i2 = lt_i2(1)
    do k = 2, f1, 1

       i = dm_list(k,1)


       call add_front(front(:,:,i),front(:,:,1),loadf(:,i),loadf(:,1),loadf_dum(:,i),&
            loadf_dum(:,1),iPg(:,i),iPg(:,1),jPg(:,i),jPg(:,1),fs(i),fs(1),i,1,1,1)
       call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
            iPg(:,1),jPg(:,1),l_i2,fs(1),id)

    end do

    !$omp section

    !$omp parallel do private(i,j,id,l_i) shared(f1,f2) num_threads(ths3)
    do j = f1+1, f2, 1
       !Thread id used to determine assembly arrays
       id = omp_get_thread_num()+2
       !write(*,*) 'id1:', id
       l_i = lt_i(j)                !initilaize var counter to predetermined starting point


       !Solve domain j and save front
       call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n(:,id))

    end do
    !$omp end parallel do
    !$omp end parallel sections

    f3 = INT(REAL(f2 - f1,rk)/tr,ik) + f2
    if (f3.gt.iDM/2) f3 = iDM/2
    if (f3-f2.lt.s_stop) f3 = iDM/2
    do
       !$omp parallel sections private(i,j,k,id,l_i) shared(l_i2,f1,f2,f3) num_threads(2)
       !$omp section
       id = 1
      
       do k = f1+1, f2, 1

          i = dm_list(k,1)


          call add_front(front(:,:,i),front(:,:,1),loadf(:,i),loadf(:,1),loadf_dum(:,i),&
               loadf_dum(:,1),iPg(:,i),iPg(:,1),jPg(:,i),jPg(:,1),fs(i),fs(1),i,1,1,1)
          call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
               iPg(:,1),jPg(:,1),l_i2,fs(1),id)

       end do
       
       !$omp section
       
       !$omp parallel do private(i,j,id,l_i) shared(f3,f2) num_threads(ths3)
       do j = f2+1, f3, 1
            !Thread id used to determine assembly arrays
          id = omp_get_thread_num() + 2
          l_i = lt_i(j)                !initilaize var counter to predetermined starting point
          !write(*,*) 'id2:', id

          !Solve domain j and save front
          call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n(:,id))

       end do
       !$omp end parallel do
       
       !$omp end parallel sections
       if (f3.eq.iDM/2) exit
       !j = f3 - f2
       f1 = f2
       f2 = f3
       f3 = INT(REAL(f2-f1,rk)/tr,ik) + f2
       if (f3.gt.iDM/2) f3 = iDM/2
       if (f2-f1.lt.s_stop) f3 = iDM/2

    end do

    id = 1
    do k = f2+1, f3, 1

       i = dm_list(k,1)
       !if (i.eq.0) write(*,*) 'Error 1'

       call add_front(front(:,:,i),front(:,:,1),loadf(:,i),loadf(:,1),loadf_dum(:,i),&
            loadf_dum(:,1),iPg(:,i),iPg(:,1),jPg(:,i),jPg(:,1),fs(i),fs(1),i,1,1,1)
       call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
            iPg(:,1),jPg(:,1),l_i2,fs(1),id)

    end do
   
    !$omp section
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PHASE 2.0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    f1 = s_frac
    f2 = INT((1.0_rk+1.0_rk/tr)*REAL(f1,rk),ik) + f1
    f3 = 0
    l_i2 = 0
    !$omp parallel sections private(i,j,k,id,l_i) shared(l_i2,f1,f2,f3) num_threads(2)

    !$omp section
    id = thsb
    do i = 1, f1, 1
        !Thread id used to determine assembly arrays
       j = dm_list(i,2)
       l_i = lt_i(j)                !initilaize var counter to predetermined starting point


       !Solve domain j and save front
       call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n(:,id))

    end do

   
    l_i2 = lt_i2(2)
    do k = 2, f1, 1

       i = dm_list(k,2)


       call add_front(front(:,:,i),front(:,:,iDM),loadf(:,i),loadf(:,iDM),loadf_dum(:,i),&
            loadf_dum(:,iDM),iPg(:,i),iPg(:,iDM),jPg(:,i),jPg(:,iDM),fs(i),fs(iDM),i,iDM,1,1)
       call elim_stitch(front(:,:,iDM),loadf(:,iDM),loadf_dum(:,iDM),&
            iPg(:,iDM),jPg(:,iDM),l_i2,fs(iDM),id)

    end do
   
    !$omp section
    
    !$omp parallel do private(i,j,id,l_i) shared(f1,f2) num_threads(ths3)
    do i = f1+1, f2, 1
        !Thread id used to determine assembly arrays
       j = dm_list(i,2)
       l_i = lt_i(j)                !initilaize var counter to predetermined starting point
       id = thsb+1+omp_get_thread_num()
       
       !Solve domain j and save front
       call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n(:,id))

    end do
    !$omp end parallel do
    !$omp end parallel sections

    f3 = INT(REAL(f2 - f1,rk)/tr,ik) + f2
    if (f3.gt.iDM/2) f3 = iDM/2
    if (f3-f2.lt.s_stop) f3 = iDM/2
    do
       !$omp parallel sections private(i,j,k,id,l_i) shared(l_i2,f1,f2,f3) num_threads(2)
       !$omp section
       id = thsb
       
       do k = f1+1, f2, 1

          i = dm_list(k,2)


          call add_front(front(:,:,i),front(:,:,iDM),loadf(:,i),loadf(:,iDM),loadf_dum(:,i),&
            loadf_dum(:,iDM),iPg(:,i),iPg(:,iDM),jPg(:,i),jPg(:,iDM),fs(i),fs(iDM),i,iDM,1,1)
          call elim_stitch(front(:,:,iDM),loadf(:,iDM),loadf_dum(:,iDM),&
            iPg(:,iDM),jPg(:,iDM),l_i2,fs(iDM),id)

       end do
      
       !$omp section
       
       !$omp parallel do private(i,j,id,l_i) shared(f3,f2) num_threads(ths3)
       do i = f2+1, f3, 1
            !Thread id used to determine assembly arrays
          j = dm_list(i,2)
          l_i = lt_i(j)                !initilaize var counter to predetermined starting point
          id = thsb+1+omp_get_thread_num()

          !Solve domain j and save front
          call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n(:,id))

       end do
       !$omp end parallel do

       !$omp end parallel sections
       if (f3.eq.iDM-iDM/2) exit
       !j = f3 - f2
       f1 = f2
       f2 = f3
       f3 = INT(REAL(f2-f1,rk)/tr,ik) + f2
       if (f3.gt.iDM-iDM/2) f3 = iDM - iDM/2
       if (f2-f1.lt.s_stop) f3 = iDM-iDM/2

    end do

    id = thsb
    do k = f2+1, f3, 1

       i = dm_list(k,2)

       call add_front(front(:,:,i),front(:,:,iDM),loadf(:,i),loadf(:,iDM),loadf_dum(:,i),&
            loadf_dum(:,iDM),iPg(:,i),iPg(:,iDM),jPg(:,i),jPg(:,iDM),fs(i),fs(iDM),i,iDM,1,1)
       call elim_stitch(front(:,:,iDM),loadf(:,iDM),loadf_dum(:,iDM),&
            iPg(:,iDM),jPg(:,iDM),l_i2,fs(iDM),id)

    end do
     !write(*,*) l_i2, 'END PHASE 2'
    !$omp end parallel sections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PHASE 3.0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !PHASE 3: Sweep through domains from two fronts in natural order of domains (2 threads at most)
    !Not required for single domain front or dual domain front

   
 
    id = 1
    l_i = lt_i2(3)
    !Combine and eliminate final fronts
    call add_front(front(:,:,iDM),front(:,:,1),loadf(:,iDM),loadf(:,1),loadf_dum(:,iDM),&
         loadf_dum(:,1),iPg(:,iDM),iPg(:,1),jPg(:,iDM),jPg(:,1),fs(iDM),fs(1),iDM,1,1,1)
    call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
         iPg(:,1),jPg(:,1),l_i,fs(1),id)


    !Perform back substitution
    call back_sub(2)

    !Calculate L2res
    L2res = 0.0_rk
    !$omp parallel private(i) reduction(+:L2res)
    !$omp do
    do i = 1, NV, 1
       L2res = L2res + (load_dum(i))**2
    end do
    !$omp end do
    !$omp end parallel
    L2res = sqrt(L2res)
 

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - ti
    j = 0
    do i = 1, ths, 1
       j = piv(i) + j
    end do
    write(*,*) 'Small Pivots:', j

 
    
  end subroutine multifront2

   subroutine back_sub(mode)
    implicit none
    integer(kind=ik) :: i, j, k, l, mode, last

    do i = NV, 1, -1
       if (IT(i,2).ne.0) then
          last = i
          exit
       end if
    end do

    !For single, dual or nested single thread front
    if (mode.eq.1) then
       load(IT(last,3)) = loadc(IT(last,2))
       do i = last-1, 1, -1
          do j = 1, IT(i,1), 1
             loadc(IT(i,2)) = loadc(IT(i,2)) - LT(i,j)*load(IT(i,j+3))
          end do

          load(IT(i,3)) = loadc(IT(i,2))

       end do
       
    else !For mode 1 parallel implementation
       
       load(IT(NV,3)) = loadc(IT(NV,2))


       do i = NV-1, lt_i2(3)+1, -1
          do j = 1, IT(i,1), 1
             loadc(IT(i,2)) = loadc(IT(i,2)) - LT(i,j)*load(IT(i,j+3))
          end do
          load(IT(i,3)) = loadc(IT(i,2))
       end do

       !$omp parallel private(k,i,j) num_threads(ths2)
       !$omp do
       do k = 2, 1, -1
          do i = lt_i2(k+1), lt_i2(k)+1, -1
             do j = 1, IT(i,1), 1
                loadc(IT(i,2)) = loadc(IT(i,2)) - LT(i,j)*load(IT(i,j+3))
             end do
             load(IT(i,3)) = loadc(IT(i,2))
          end do
       end do
       !$omp end do
       !$omp end parallel


       !$omp parallel private(k,i,j)
       !$omp do
       do k = iDM, 1, -1
          do i = lt_i(k+1), lt_i(k)+1, -1
             do j = 1, IT(i,1), 1
                loadc(IT(i,2)) = loadc(IT(i,2)) - LT(i,j)*load(IT(i,j+3))
             end do
             load(IT(i,3)) = loadc(IT(i,2))
          end do
       end do
       !$omp end do
       !$omp end parallel
    end if

  end subroutine back_sub

  !Determines the ordering for nested dissection mode 2
  subroutine setup_nested(opt)
    implicit none
    integer(kind=ik) :: i, j, k, l, o, p, m, n, enn(NSi), en(NSi), e_n(4), type, pairs
    logical :: opt

    write(*,*) 'Skipping Nested allocations'
    return
    
    if (opt) then
       !o = (ptn_g + mod(ptn_g,2) + ptn2_g + mod(ptn2_g,2))
       allocate(eprs(4,NSi,2:iDM-1),NODf_n(NSi,NV,ths)) !element ordering for smallist level
       eprs = 0

       !PHASE 1: All smallest level Nested
       do j = 2, iDM-1, 1
          !Type 1: 2x2
          !Type 2: 1x2
          !Type 3: 2x1
          !Type 4: 1x1
          select case(dm_type(j))
          case(1)
             o = ptn_g
             p = ptn2_g
          case(2)
             o = ptn_g
             p = prm2
          case(3)
             o = prm
             p = ptn2_g
          case(4)
             o = prm
             p = prm2
          case default
             write(*,*) 'Missing case:', dm_type(j), j
          end select

          l = -1 !Keeps count on element
          m = 1  !Keeps count of sub-domain
          !Find groups of 4
          do k = 1, p/2, 1

             do i = 1, o/2, 1
                l = l + 2
                eprs(1,m,j) = ele_list(l,j)
                eprs(2,m,j) = ele_list(l+1,j)
                eprs(3,m,j) = ele_list(l+o,j)
                eprs(4,m,j) = ele_list(l+o+1,j)

                m = m + 1
             end do
             if (mod(o,2).ne.0) then !if type 3
                l = l + 1
                eprs(1,m,j) = ele_list(l+1,j)
                eprs(2,m,j) = ele_list(l+o+1,j)
                !e_n(3) = 0
                !e_n(4) = 0
                m = m + 1
             end if
             l = l + o

          end do
          if (mod(p,2).ne.0) then !if type 2 
             do i = 1, o/2, 1
                l = l + 2
                eprs(1,m,j) = ele_list(l,j)
                eprs(2,m,j) = ele_list(l+1,j)
                !e_n(3) = 0
                !e_n(4) = 0
                m = m + 1
             end do
             if (mod(o,2).ne.0) then !if type 4
                l = l + 1
                eprs(1,m,j) = ele_list(l+1,j)
                !e_n(2) = 0
                !e_n(3) = 0
                !e_n(4) = 0
                m = m + 1
             end if
          end if

       end do


       !Run through dry to see number of pairs required
       !Reducing type 1 is 3 pairs, type 2,3 1 pair, type 4 0 pairs
       n = 0
       do type = 1, 4, 1 
          select case(type)
          case(1)
             o = ptn_g
             p = ptn2_g
          case(2)
             o = ptn_g
             p = prm2
          case(3)
             o = prm
             p = ptn2_g
          case(4)
             o = prm
             p = prm2
          case default
             write(*,*) 'Missing case:', dm_type(j), j
          end select
          if (o.eq.0.or.p.eq.0) cycle
          o = o + mod(o,2)
          p = p + mod(p,2)
          o = o/2
          p = p/2

          pairs = 0
          do 
             l = -1
             m = 1
             do k = 1, p/2, 1
                do i = 1, o/2, 1
                   l = l + 2
                   m = m + 1
                   pairs = pairs + 3
                end do
                if (mod(o,2).ne.0) then
                   l = l + 1
                   pairs = pairs + 1
                   m = m + 1
                end if
                l = l + o
             end do
             if (mod(p,2).ne.0) then
                do i = 1, o/2, 1
                   l = l + 2
                   m = m + 1
                   pairs = pairs + 1
                end do
                if (mod(o,2).ne.0) then
                   l = l + 1
                   m = m + 1
                end if
             end if
             o = o + mod(o,2)
             p = p + mod(p,2)
             o = o/2
             p = p/2

             if (o.eq.1.and.p.eq.1) exit
          end do

          if (pairs.gt.n) then
             n = pairs
          end if

       end do
       max_prs = n
       !Allocate array for nested ordering
       !Only needs to store pairing for the 4 types not every sub domain
       allocate(dprs(2,max_prs,4))
       dprs = 0



       !Run through and determine pairs
       do type = 1, 4, 1 
          select case(type)
          case(1)
             o = ptn_g
             p = ptn2_g
          case(2)
             o = ptn_g
             p = prm2
          case(3)
             o = prm
             p = ptn_g
          case(4)
             o = prm
             p = prm2
          case default
             write(*,*) 'Missing case:', dm_type(j), j
          end select
          if (o.eq.0.or.p.eq.0) cycle
          o = o + mod(o,2)
          p = p + mod(p,2)
          o = o/2
          p = p/2
          do i = 1, o*p, 1
             en(i) = i
          end do


          !pause
          pairs = 0
          do 



             !write(*,*) 'o,p:', o, p
             l = -1
             m = 1
             do k = 1, p/2, 1

                do i = 1, o/2, 1
                   l = l + 2
                   e_n(1) = l
                   e_n(2) = l+1
                   e_n(3) = l+o
                   e_n(4) = l+o+1

                   pairs = pairs + 1  
                   dprs(1,pairs,type) = en(e_n(1))
                   dprs(2,pairs,type) = en(e_n(2))
                   ! write(*,*) '1e_n:', e_n(1:2)
                   ! write(*,*) '1dprs:', dprs(1,pairs,type), dprs(2,pairs,type)
                   pairs = pairs + 1
                   dprs(1,pairs,type) = en(e_n(1))
                   dprs(2,pairs,type) = en(e_n(3))
                   pairs = pairs + 1
                   dprs(1,pairs,type) = en(e_n(1))
                   dprs(2,pairs,type) = en(e_n(4))
                   ! write(*,*) '2e_n:', e_n(3:4)
                   ! write(*,*) '2dprs:', dprs(1,pairs,type), dprs(2,pairs,type)
                   enn(m) = en(e_n(1))
                   m = m + 1

                end do
                if (mod(o,2).ne.0) then
                   l = l + 1
                   e_n(1) = l+1
                   e_n(2) = l+o+1
                   e_n(3) = 0
                   e_n(4) = 0

                   pairs = pairs + 1
                   dprs(1,pairs,type) = en(e_n(1))
                   dprs(2,pairs,type) = en(e_n(2))

                   enn(m) = en(e_n(1))
                   m = m + 1
                end if
                l = l + o

             end do
             if (mod(p,2).ne.0) then
                do i = 1, o/2, 1
                   l = l + 2
                   e_n(1) = l
                   e_n(2) = l+1
                   e_n(3) = 0
                   e_n(4) = 0

                   pairs = pairs + 1
                   dprs(1,pairs,type) = en(e_n(1))
                   dprs(2,pairs,type) = en(e_n(2))

                   enn(m) = en(e_n(1))
                   m = m + 1

                end do
                if (mod(o,2).ne.0) then
                   l = l + 1
                   e_n(1) = l+1
                   e_n(2) = 0
                   e_n(3) = 0
                   e_n(4) = 0
                   !write(*,*) 'e_n(:):', e_n(:)
                   enn(m) = en(e_n(1))
                   m = m + 1
                   !MOVE ONLY LAST FRONT
                end if
             end if



             o = o + mod(o,2)
             p = p + mod(p,2)
             o = o/2
             p = p/2
             en = enn
             if (o.eq.1.and.p.eq.1) exit
          end do

       end do
    else
       deallocate(eprs,dprs,NODf_n)
    end if
       
  end subroutine setup_nested
    
  !Eliminates a front which was the result of combining two fronts
  subroutine elim_stitch(front,loadf,loadf_dum,iPg,jPg,l_i,fs,id)
    implicit none
    integer(kind=ik) :: fs, iPg(iBW),jPg(iBW), NOPf(NE,bas)
    real(kind=rk) :: front(iBW,iBW), loadf(iBW), loadf_dum(iBW)
    integer(kind=ik) ::  i, j, k, flag1, last_e, last_n,&
         RSUM(iBW), CSUM(iBW), rvs, cvs, CPIV, RPIV, l_i, elim, j_max, k_max, id
    real(kind=rk) ::  PIVOT, FAC


    rvs = 0
    cvs = 0

    !Extract fully summed from add_front
    do i = 1, fs, 1
       if (iPg(i).lt.0) then
          rvs = rvs + 1
          RSUM(rvs) = i
          iPg(i) = -iPg(i)
       end if
       if (jPg(i).lt.0) then
          cvs = cvs + 1
          CSUM(cvs) = i
          jPg(i) = -jPg(i)
       end if
    end do

    ! write(*,*) 'Bout to call elim'
    ! write(*,*) 'l_i', l_i
    ! write(*,*) 'Calling elim'
    !Call elimination subroutine
    call elim_summed(front,loadf,loadf_dum,iPg,jPg,rvs,cvs,RSUM,CSUM,l_i,fs,id)

    
    
  end subroutine elim_stitch

  !Adds two partial fronts together
  subroutine add_front(front,front_m,loadf,loadf_m,loadf_dum,loadf_dumm,iPg,iPg_m,jPg,jPg_m,fs,fs_m,id2,id1,mode,th)
    implicit none
    real(kind=rk) :: front(iBW,iBW), front_m(iBW,iBW), loadf(iBW), loadf_m(iBW),&
         loadf_dum(iBW), loadf_dumm(iBW)
    integer(kind=ik) :: iPg(iBW), jPg(iBW), fs,&
         iPg_m(iBW), jPg_m(iBW), fs_m, i, j, k, chead, rhead, flag1, id2, id1, mode, th
    if (iDM.eq.1) return
    !Should never be called for single domain

    
    
    !Modes determine which NDF to use
    !Combine domain id2 and id1
    if (mode.eq.1) then
       NODf(id1,:) = NODf(id2,:) + NODf(id1,:)

       !Add nonexistant rows and columns
       do i = 1, fs, 1
          !Find matching row
          flag1 = 0
          do k = 1, fs_m, 1
             if (iPg(i).eq.iPg_m(k)) then
                !rhead = k
                flag1 = 1

             end if
          end do
          if (flag1.eq.0) then
             fs_m = fs_m + 1
            
             !rhead = fs_m
             iPg_m(fs_m) = iPg(i)
             jPg_m(fs_m) = iPg(i)

             !write(*,*) 'added row'
          end if
       end do
       !write(*,*) 'fs:', fs_m, id1, id2


       !Transfer front to master front and flag fully summed in iPg/jPg
       do i = 1, fs, 1
          !Find matching row
          flag1 = 0
          do k = 1, fs_m, 1
             if (iPg(i).eq.ABS(iPg_m(k))) then
                rhead = k
                if (NODf(id1,iPg_m(k)).eq.NOA(iPg_m(k))) iPg_m(k) = -iPg_m(k)
                flag1 = 1
             end if
          end do
          if (flag1.eq.0) then
             fs_m = fs_m + 1
             rhead = fs_m
             iPg_m(fs_m) = iPg(i)
             jPg_m(fs_m) = iPg(i)

             write(*,*) 'error'
          end if


          do j = 1, fs, 1
             !Find Matching column
             flag1 = 0
             do k = 1, fs_m, 1
                if (jPg(j).eq.ABS(jPg_m(k))) then
                   chead = k
                   if (jPg_m(k).gt.0) then
                      if (NODf(id1,jPg_m(k)).eq.NOA(jPg_m(k))) then
                         jPg_m(k) = -jPg_m(k)
                         !write(*,*) 'Found fully summed:', jPg_m(k)
                      end if
                   end if
                   flag1 = 1
                end if
             end do
             if (flag1.eq.0) then
                fs_m = fs_m + 1
                chead = fs_m
                jPg_m(fs_m) = iPg(i)
                write(*,*) 'error'
             end if
             front_m(chead,rhead) = front(j,i) + front_m(chead,rhead)
          end do
          loadf_m(rhead) = loadf_m(rhead) + loadf(i)
          loadf_dumm(rhead) = loadf_dumm(rhead) + loadf_dum(i)
       end do
       
       !combine nested
    else if (mode.eq.2) then
       
       NODf_n(id1,:,th) = NODf_n(id2,:,th) + NODf_n(id1,:,th)
       
       !Add nonexistant rows and columns
       do i = 1, fs, 1
          !Find matching row
          flag1 = 0
          do k = 1, fs_m, 1
             if (iPg(i).eq.iPg_m(k)) then
                !rhead = k
                flag1 = 1

             end if
          end do
          if (flag1.eq.0) then
             fs_m = fs_m + 1
             !rhead = fs_m
             iPg_m(fs_m) = iPg(i)
             jPg_m(fs_m) = iPg(i)

             !write(*,*) 'added row'
          end if
       end do



       !Transfer front to master front and flag fully summed in iPg/jPg
       do i = 1, fs, 1
          !Find matching row
          flag1 = 0
          do k = 1, fs_m, 1
             if (iPg(i).eq.ABS(iPg_m(k))) then
                rhead = k
                if (NODf_n(id1,iPg_m(k),th).eq.NOA(iPg_m(k))) iPg_m(k) = -iPg_m(k)
                flag1 = 1
             end if
          end do
          if (flag1.eq.0) then
             fs_m = fs_m + 1
             rhead = fs_m
             iPg_m(fs_m) = iPg(i)
             jPg_m(fs_m) = iPg(i)

             write(*,*) 'error'
          end if


          do j = 1, fs, 1
             !Find Matching column
             flag1 = 0
             do k = 1, fs_m, 1
                if (jPg(j).eq.ABS(jPg_m(k))) then
                   chead = k
                   if (jPg_m(k).gt.0) then
                      if (NODf_n(id1,jPg_m(k),th).eq.NOA(jPg_m(k))) then
                         jPg_m(k) = -jPg_m(k)
                         !write(*,*) 'Found fully summed:', jPg_m(k)
                      end if
                   end if
                   flag1 = 1
                end if
             end do
             if (flag1.eq.0) then
                fs_m = fs_m + 1
                chead = fs_m
                jPg_m(fs_m) = iPg(i)
                write(*,*) 'error'
             end if
             front_m(chead,rhead) = front(j,i) + front_m(chead,rhead)
          end do
          loadf_m(rhead) = loadf_m(rhead) + loadf(i)
          loadf_dumm(rhead) = loadf_dumm(rhead) + loadf_dum(i)
       end do

       !Combine domain and nested
    else 
       NODf(id1,:) = NODf_n(id2,:,th) + NODf(id1,:)

        !Add nonexistant rows and columns
       do i = 1, fs, 1
          !Find matching row
          flag1 = 0
          do k = 1, fs_m, 1
             if (iPg(i).eq.iPg_m(k)) then
                !rhead = k
                flag1 = 1

             end if
          end do
          if (flag1.eq.0) then
             fs_m = fs_m + 1
             !rhead = fs_m
             iPg_m(fs_m) = iPg(i)
             jPg_m(fs_m) = iPg(i)

             !write(*,*) 'added row'
          end if
       end do



       !Transfer front to master front and flag fully summed in iPg/jPg
       do i = 1, fs, 1
          !Find matching row
          flag1 = 0
          do k = 1, fs_m, 1
             if (iPg(i).eq.ABS(iPg_m(k))) then
                rhead = k
                if (NODf(id1,iPg_m(k)).eq.NOA(iPg_m(k))) iPg_m(k) = -iPg_m(k)
                flag1 = 1
             end if
          end do
          if (flag1.eq.0) then
             fs_m = fs_m + 1
             rhead = fs_m
             iPg_m(fs_m) = iPg(i)
             jPg_m(fs_m) = iPg(i)

             write(*,*) 'error'
          end if


          do j = 1, fs, 1
             !Find Matching column
             flag1 = 0
             do k = 1, fs_m, 1
                if (jPg(j).eq.ABS(jPg_m(k))) then
                   chead = k
                   if (jPg_m(k).gt.0) then
                      if (NODf(id1,jPg_m(k)).eq.NOA(jPg_m(k))) then
                         jPg_m(k) = -jPg_m(k)
                         !write(*,*) 'Found fully summed:', jPg_m(k)
                      end if
                   end if
                   flag1 = 1
                end if
             end do
             if (flag1.eq.0) then
                fs_m = fs_m + 1
                chead = fs_m
                jPg_m(fs_m) = iPg(i)
                write(*,*) 'error'
             end if
             front_m(chead,rhead) = front(j,i) + front_m(chead,rhead)
          end do
          loadf_m(rhead) = loadf_m(rhead) + loadf(i)
          loadf_dumm(rhead) = loadf_dumm(rhead) + loadf_dum(i)
       end do
    end if



  end subroutine add_front

 
  !Solves front from element list
  Subroutine Front_solve(fs,front,loadf,loadf_dum,iPg,jPg,l_i,id,dm,mode,e_n)
    implicit none
    integer(kind=ik) :: ele1, ele2, fs, NOPf(NE,bas), iPg(iBW),jPg(iBW),e1, e2, es
    real(kind=rk) :: front(iBW,iBW), loadf(iBW), loadf_dum(iBW)
    integer(kind=ik) :: ele, i, j, k, rhead(NB), chead(NB), NOPPl(bas), flag1,&
         RSUM(iBW), CSUM(iBW), rvs, cvs, CPIV, RPIV, l_i, elim, j_max, k_max, &
         last_e, last_n, o, p, l, NK(NB), id, dm, mode, e_n(4)
    real(kind=rk) ::  local(NB,NB), loc(NB), PIVOT, FAC, t2, t1
    
    !Initialize local to global
    CPIV = 0 !Pivotal column
    RPIV = 0 !Pivotal row
    cvs = 0 !#summed col
    rvs = 0 !#summed rows (same as cvs)
    !piv = 0 !#small pivots (global now)

    !Find last appearance of node using NODf
    NOPf = NOP
    flag1 = 0
   

    !if (iDM.ne.dm) then
       
    if (mode.eq.1) then !Use domain ele_list to find summed vars and change NOPf
      
       
       e1 = 0
       do
          e1 = e1 + 1
          ele = ele_list(e1,dm)
          if (ele.eq.0) exit
         


          do i = 1, 9, 1
             o = NOP(ele,i)
             p = NOPP(o)
             NODf(dm,p) = NODf(dm,p) + 1
             if (NODf(dm,p).eq.NOA(p)) then
                NOPf(ele,i) = -NOP(ele,i)

             end if
             do k = 1, MDF(o)-1, 1
                NODf(dm,p+k) = NODf(dm, p+k) + 1
             end do
          end do
       end do
    else if (mode.eq.2) then !Use passed in list to find summed vars and change NOPf and nested NDF
       e1 = 0
       do
          e1 = e1 + 1
          if (e1.eq.5) exit
          ele = e_n(e1)
          if (ele.eq.0) exit
          

          do i = 1, 9, 1
             o = NOP(ele,i)
             p = NOPP(o)
             NODf_n(dm,p,id) = NODf_n(dm,p,id) + 1
             if (NODf_n(dm,p,id).eq.NOA(p)) then
                NOPf(ele,i) = -NOP(ele,i)

             end if
             do k = 1, MDF(o)-1, 1
                NODf_n(dm,p+k,id) = NODf_n(dm, p+k,id) + 1
             end do
          end do
       end do
    else if (mode.eq.3) then !Use passed in list to find summed vars and change NOPf and domain NDF
       e1 = 0
       do
          e1 = e1 + 1
          if (e1.eq.5) exit
          ele = e_n(e1)
          if (ele.eq.0) exit
          

          do i = 1, 9, 1
             o = NOP(ele,i)
             p = NOPP(o)
             NODf(dm,p) = NODf(dm,p) + 1
             if (NODf(dm,p).eq.NOA(p)) then
                NOPf(ele,i) = -NOP(ele,i)

             end if
             do k = 1, MDF(o)-1, 1
                NODf(dm,p+k) = NODf(dm, p+k) + 1
             end do
          end do
       end do
    end if
       
    e1 = 0
    do
       !Extract next element from proper list
       e1 = e1 + 1
       if (mode.eq.1) then
          ele = ele_list(e1,dm)
          if (ele.eq.0) exit
       else
          if (e1.eq.5) exit
          ele = e_n(e1)
          if (ele.eq.0) exit
       end if
          
       !USER_SPECIFIED
       !Assemble local using proper thread (id)
       call assemble_local(ele,local,loc,NOPPl,NB,id,0)
       !USER_SPECIFIED

       !Create NK vector
       do i = 1, bas, 1
          j = NOPP(NOP(ele,i))
          o = NOPPl(i)
          do k = 0, MDF(NOP(ele,i))-1, 1
             NK(k+o) = k+j
          end do
       end do
       
       
       !Create heading victors
       do i = 1, NCN(ele), 1
          flag1 = 0
          do j = 1, fs, 1
             if (NK(i).eq.iPg(j)) then
                rhead(i) = j
                flag1 = 1
                exit
             end if
          end do
          !Create new row and column
          if (flag1.eq.0) then
             fs = fs + 1
             iPg(fs) = NK(i)
             jPg(fs) = NK(i)
             rhead(i) = fs
             chead(i) = fs
          else !find column
             do j = 1, fs, 1
                if (NK(i).eq.jPg(j)) then                
                   chead(i) = j
                   exit
                end if

             end do
          end if

       end do

       !Warning front width exceeded, will probably seg fault
       if (fs.gt.iBW) then
          write(*,*) 'Error front exceeds max size, fs:', fs
          pause
          return
       end if
  
       !Add local to front
       do i = 1, NCN(ele), 1
          do j = 1, NCN(ele), 1
             front(chead(j),rhead(i)) = local(j,i) + front(chead(j),rhead(i))
          end do
          loadf(rhead(i)) = loc(i) + loadf(rhead(i))
          loadf_dum(rhead(i)) = loc(i) + loadf_dum(rhead(i))
       end do
       ! do i = 1, NCN(ele), 1
       !    loadf(rhead(i)) = loc(i) + loadf(rhead(i))
       ! end do

       
       !Find fully summed using NOPf
       do i = 1, bas, 1
          if (NOPf(ele,i).lt.0) then
             do k = 0, MDF(NOP(ele,i))-1, 1

                rvs = rvs + 1
                RSUM(rvs) = rhead(k+NOPPl(i))
                cvs = cvs + 1
                CSUM(cvs) = chead(k+NOPPl(i))

             end do
          end if
       end do

       !Call elimination phase
       call elim_summed(front,loadf,loadf_dum,iPg,jPg,rvs,cvs,RSUM,CSUM,l_i,fs,id)

    

    end do
   
    !if (piv.ne.0) write(*,*) piv, 'small pivots'
    !write(*,*) 'rem. fs:', fs

  End Subroutine Front_solve

 

  subroutine elim_summed(front,loadf,loadf_dum,iPg,jPg,rvs,cvs,RSUM,CSUM,l_i,fs,id)
    implicit none
    integer(kind=ik) :: fs, iPg(iBW),jPg(iBW)
    real(kind=rk) :: front(iBW,iBW), loadf(iBW), loadf_dum(iBW)
    integer(kind=ik) ::  i, j, k, id,&
         RSUM(iBW), CSUM(iBW), rvs, cvs, CPIV, RPIV, l_i, elim, j_max, k_max, o, p, l,n,m, q, r
    real(kind=rk) ::  PIVOT, FAC, QQ(iBW), c_piv, GG(iBW), QG(iBW,iBW)
    !write(*,*) 'Called'
    !eliminate fully summed full pivotal choice
    elim = rvs
    if (l_i.lt.0.or.l_i.gt.NV) write(*,*) 'Error l_i', l_i
    do i = 1, elim, 1

       !Find sufficiently large pivot
       j_max = 1
       k_max = 1
       PIVOT = ABS(front(CSUM(k_max),RSUM(j_max)))
       if (PIVOT.lt.0.0001_rk) then
          do k = 1, rvs, 1
             
             do j = 1, cvs, 1
                if (PIVOT.lt.ABS(front(CSUM(k),RSUM(j)))) then
                   j_max = j
                   k_max = k
                   PIVOT = ABS(front(CSUM(k),RSUM(j)))
                end if
             end do
          end do
       end if
       RPIV = RSUM(j_max)
       CPIV = CSUM(k_max)
       
       !Reorder list
       !if (j_max.ne.rvs) then
          do j = j_max, rvs-1
             RSUM(j) = RSUM(j+1)
          end do
       !end if
       do j = 1, rvs-1, 1
          if (RSUM(j).gt.RPIV) RSUM(j) = RSUM(j) - 1
       end do

       !Reorder list
       !if (k_max.ne.cvs) then
          do j = k_max, cvs-1
             CSUM(j) = CSUM(j+1)
          end do
       !end if
       do j = 1, cvs-1, 1
          if (CSUM(j).gt.CPIV) CSUM(j) = CSUM(j) - 1
       end do

       !Decrement rows and columns to delete
       rvs = rvs - 1
       cvs = cvs - 1
       
       !Set PIVOT
       PIVOT = front(CPIV,RPIV)
       !Check if pivot too small
       if (ABS(PIVOT).le.1.0e-18_rk) then
          !write(*,*) 'Small Pivot:', PIVOT
          piv(id) = piv(id)+1
       end if
         
       !Normalize Pivotal Row
       do j = 1, fs, 1
          front(j,RPIV) = front(j,RPIV)/PIVOT
       end do
       loadf(RPIV) = loadf(RPIV)/PIVOT
       

       !Save pivotal row
       l_i = l_i+1           !increment eliminated var count
       IT(l_i,1) = fs-1      !Save number in row
       IT(l_i,2) = iPg(RPIV) !Save global headers
       IT(l_i,3) = jPg(CPIV) !Save column header for pivot
       loadc(iPg(RPIV)) = loadf(RPIV) 
       load_dum(iPg(RPIV)) = loadf_dum(RPIV)
       ! if (fs.eq.1) then
       !    write(*,*) iPg(RPIV)
       !    fs = 0
       !    exit
       ! end if
       !
       
       o = CPIV-1
       p = CPIV+1
       do j = 1, o, 1
          LT(l_i,j) = front(j,RPIV)
          IT(l_i,3+j) = jPg(j)
       end do
       k = 1
       do j = p, fs, 1
          LT(l_i,j-k) = front(j,RPIV)
          IT(l_i,3+j-k) = jPg(j)
       end do
       !Cogy out pivotal row and column
       do j = 1, fs, 1
          QQ(j) = front(j,RPIV)
          GG(j) = front(CPIV,j)
       end do
       c_piv = loadf(RPIV) !Save rhs
       
       !Perform actual elimination
       ![                            |           ]
       ![                            |           ]
       ![              1             |     2     ]
       ![                            |           ]
       ![                            |           ]
       ![                            |           ]
       ![                            |           ]
       ![----------------------------------------]
       ![                            |           ]
       ![              3             |     4     ]
       ![                            |           ]
       

       
       !Precalculate QQ*GG, slower
       
       ! do k = 1, RPIV-1, 1
       !    do j = 1, o, 1
       !       QG(j,k) = QQ(j)*GG(k)
       !    end do
       !    do j = p, fs, 1
       !       QG(j,k) = QQ(j)*GG(k)
       !    end do
       ! end do

       ! do k = RPIV+1, fs, 1
       !    do j = 1, o, 1
       !       QG(j,k) = QQ(j)*GG(k)
       !    end do
       !    do j = p, fs, 1
       !       QG(j,k) = QQ(j)*GG(k)
       !    end do
       ! end do

       ! if (ths.eq.1) then
          do k = 1, RPIV-1, 1
             !FAC = front(CPIV,k)
             if (GG(k).ne.0.0_rk) then
                do j = 1, o, 1 !Doing Sector 1 update
                   front(j,k) = front(j,k) - QQ(j)*GG(k)
                end do
                !front(1:o,k) = front(1:o,k) - QQ(1:o)*GG(k)

                loadf(k) = loadf(k) - c_piv*GG(k)
             end if

          end do

          do k = 1, RPIV-1, 1
             !FAC = front(CPIV,k)
             do j = p, fs, 1 !Doing Sector 2 update
                front(j-1,k) = front(j,k) - QQ(j)*GG(k)
             end do
             !front(p-1:fs-1,k) = front(p:fs,k) - QQ(p:fs)*GG(k)
          end do

          do k = RPIV+1, fs, 1
             !FAC = front(CPIV,k)
             do j = 1, o, 1 !Doing Sector 3 update
                front(j,k-1) = front(j,k) - QQ(j)*GG(k)
             end do
             !front(1:o,k-1) = front(1:o,k) - QQ(1:o)*GG(k)
             loadf_dum(k-1) = loadf_dum(k)

          end do

          do k = RPIV+1, fs, 1
             !FAC = front(CPIV,k)
             do j = p, fs, 1 !Doing Sector 4 update
                front(j-1,k-1) = front(j,k) - QQ(j)*GG(k)
             end do
             !front(p-1:fs-1,k-1) = front(p:fs,k) - QQ(p:fs)*GG(k)
             loadf(k-1) = loadf(k) - c_piv*GG(k)
          end do
       ! else
          ! !$omp parallel sections private(j,k) num_threads(ths3)
          
          ! !$omp section
          
          
          
          
          ! do k = 1, RPIV-1, 1
          !    !FAC = front(CPIV,k)
          !    if (GG(k).ne.0.0_rk) then
          !       do j = 1, o, 1 !Doing Sector 1 update
          !          front(j,k) = front(j,k) - QQ(j)*GG(k)
          !       end do
          !       !front(1:o,k) = front(1:o,k) - QQ(1:o)*GG(k)

          !       loadf(k) = loadf(k) - c_piv*GG(k)
          !    end if

          ! end do

          ! !$omp section

          ! do k = 1, RPIV-1, 1
          !    !FAC = front(CPIV,k)
          !    do j = p, fs, 1 !Doing Sector 2 update
          !       front(j-1,k) = front(j,k) - QQ(j)*GG(k)
          !    end do
          !    !front(p-1:fs-1,k) = front(p:fs,k) - QQ(p:fs)*GG(k)
          ! end do

          ! !$omp section
          ! do k = RPIV+1, fs, 1
          !    !FAC = front(CPIV,k)
          !    do j = 1, o, 1 !Doing Sector 3 update
          !       front(j,k-1) = front(j,k) - QQ(j)*GG(k)
          !    end do
          !    !front(1:o,k-1) = front(1:o,k) - QQ(1:o)*GG(k)
          !    loadf_dum(k-1) = loadf_dum(k)

          ! end do

          ! !$omp section
          ! do k = RPIV+1, fs, 1
          !    !FAC = front(CPIV,k)
          !    do j = p, fs, 1 !Doing Sector 4 update
          !       front(j-1,k-1) = front(j,k) - QQ(j)*GG(k)
          !    end do
          !    !front(p-1:fs-1,k-1) = front(p:fs,k) - QQ(p:fs)*GG(k)
          !    loadf(k-1) = loadf(k) - c_piv*GG(k)
          ! end do

          ! !$omp end parallel sections
       ! end if
       
       !Old serial
       ! do k = 1, RPIV-1, 1
       !    !FAC = front(CPIV,k)
       !    if (GG(k).ne.0.0_rk) then
       !       do j = 1, CPIV-1, 1 !Doing Sector 1 update
       !          front(j,k) = front(j,k) - QQ(j)*GG(k)!LT(l_i,j)*FAC
       !       end do
       !       loadf(k) = loadf(k) - c_piv*GG(k)!load((iPg(RPIV))*FAC
       !    end if
       !    do j = CPIV+1, fs, 1 !Doing Sector 2 update
       !       front(j-1,k) = front(j,k) - QQ(j)*GG(k)!LT(l_i,j)*FAC
       !    end do
       ! end do

       

       ! do k = RPIV+1, fs, 1
       !    !FAC = front(CPIV,k)
       !    do j = 1, CPIV-1, 1 !Doing Sector 3 update
       !       front(j,k-1) = front(j,k) - QQ(j)*GG(k)!LT(l_i,j)*FAC
       !    end do
          
          
       !    do j = CPIV+1, fs, 1 !Doing Sector 4 update
       !       front(j-1,k-1) = front(j,k) - QQ(j)*GG(k)!LT(l_i,j)*FAC
       !    end do
       !    loadf(k-1) = loadf(k) - c_piv*GG(k)!loadc(iPg(RPIV))*FAC
       !    loadf_dum(k-1) = loadf_dum(k)
       ! end do
       

       !Blank out the area outside the front
       do j = 1, fs, 1
          front(fs,j) = 0.0_rk
       end do
       do j = 1, fs-1, 1
          front(j,fs) = 0.0_rk
       end do
       loadf(fs) = 0.0_rk
       loadf_dum(fs) = 0.0_rk
       fs  = fs - 1

       !Shift global indices
       do j = RPIV, fs, 1
          iPg(j) = iPg(j+1)
       end do
       do j = CPIV, fs, 1
          jPg(j) = jPg(j+1)
       end do

    


    end do
  
  end subroutine elim_summed

  subroutine elim_summedp(front,loadf,loadf_dum,iPg,jPg,rvs,cvs,RSUM,CSUM,l_i,fs,piv)
    implicit none
    integer(kind=ik) :: fs, iPg(iBW),jPg(iBW)
    real(kind=rk) :: front(iBW,iBW), loadf(iBW), loadf_dum(iBW)
    integer(kind=ik) ::  i, j, k,&
         RSUM(iBW), CSUM(iBW), rvs, cvs, CPIV, RPIV, l_i, elim, j_max, k_max, o, p, l,n,m, q, r, piv
    real(kind=rk) ::  PIVOT, FAC, QQ(iBW), c_piv, GG(iBW)!, QG(iBW,iBW)
    
    !eliminate fully summed full pivotal choice
    elim = rvs
    do i = 1, elim, 1

       !Find sufficiently large pivot
       j_max = 1
       k_max = 1
       PIVOT = ABS(front(CSUM(k_max),RSUM(j_max)))
       if (PIVOT.lt.0.0001_rk) then
          do k = 1, rvs, 1
             
             do j = 1, cvs, 1
                if (PIVOT.lt.ABS(front(CSUM(k),RSUM(j)))) then
                   j_max = j
                   k_max = k
                   PIVOT = ABS(front(CSUM(k),RSUM(j)))
                end if
             end do
          end do
       end if
       RPIV = RSUM(j_max)
       CPIV = CSUM(k_max)
       
       !Reorder list
       if (j_max.ne.rvs) then
          do j = j_max, rvs-1
             RSUM(j) = RSUM(j+1)
          end do
       end if
       do j = 1, rvs-1, 1
          if (RSUM(j).gt.RPIV) RSUM(j) = RSUM(j) - 1
       end do

       !Reorder list
       if (k_max.ne.cvs) then
          do j = k_max, cvs-1
             CSUM(j) = CSUM(j+1)
          end do
       end if
       do j = 1, cvs-1, 1
          if (CSUM(j).gt.CPIV) CSUM(j) = CSUM(j) - 1
       end do

       !Decrement rows and columns to delete
       rvs = rvs - 1
       cvs = cvs - 1
       
       !Set PIVOT
       PIVOT = front(CPIV,RPIV)
       !Check if pivot too small
       if (ABS(PIVOT).le.1.0e-18_rk) then
          !write(*,*) 'Small Pivot:', PIVOT
          piv = piv+1
       end if
         
       !Normalize Pivotal Row
       do j = 1, fs, 1
          front(j,RPIV) = front(j,RPIV)/PIVOT
       end do
       loadf(RPIV) = loadf(RPIV)/PIVOT
       

       !Save pivotal row
       l_i = l_i+1           !increment eliminated var count
       IT(l_i,1) = fs-1      !Save number in row
       IT(l_i,2) = iPg(RPIV) !Save global headers
       IT(l_i,3) = jPg(CPIV) !Save column header for pivot
       loadc(iPg(RPIV)) = loadf(RPIV) 
       load_dum(iPg(RPIV)) = loadf_dum(RPIV)
       if (fs.eq.1) then
          fs = 0
          exit
       end if
       !
       do j = 1, CPIV-1, 1
          LT(l_i,j) = front(j,RPIV)
          IT(l_i,3+j) = jPg(j)
       end do
       k = 1
       do j = CPIV+1, fs, 1
          LT(l_i,j-k) = front(j,RPIV)
          IT(l_i,3+j-k) = jPg(j)
       end do
       !Cogy out pivotal row and column
       do j = 1, fs, 1
          QQ(j) = front(j,RPIV)
          GG(j) = front(CPIV,j)
       end do
       c_piv = loadf(RPIV) !Save rhs
       
       !Perform actual elimination
       ![                            |           ]
       ![                            |           ]
       ![              1             |     2     ]
       ![                            |           ]
       ![                            |           ]
       ![                            |           ]
       ![                            |           ]
       ![----------------------------------------]
       ![                            |           ]
       ![              3             |     4     ]
       ![                            |           ]
       

       o = CPIV-1
       p = CPIV+1

     
       !$omp parallel sections private(j,k) num_threads(ths3)
       
       !$omp section
       do k = 1, RPIV-1, 1
          !FAC = front(CPIV,k)
          if (GG(k).ne.0.0_rk) then
             do j = 1, o, 1 !Doing Sector 1 update
                front(j,k) = front(j,k) - QQ(j)*GG(k)
             end do
             !front(1:o,k) = front(1:o,k) - QQ(1:o)*GG(k)

             loadf(k) = loadf(k) - c_piv*GG(k)
          end if

       end do

       !$omp section

       do k = 1, RPIV-1, 1
          !FAC = front(CPIV,k)
          do j = p, fs, 1 !Doing Sector 2 update
             front(j-1,k) = front(j,k) - QQ(j)*GG(k)
          end do
          !front(p-1:fs-1,k) = front(p:fs,k) - QQ(p:fs)*GG(k)
       end do

       !$omp section
       do k = RPIV+1, fs, 1
          !FAC = front(CPIV,k)
          do j = 1, o, 1 !Doing Sector 3 update
             front(j,k-1) = front(j,k) - QQ(j)*GG(k)
          end do
          !front(1:o,k-1) = front(1:o,k) - QQ(1:o)*GG(k)
          loadf_dum(k-1) = loadf_dum(k)

       end do

       !$omp section
       do k = RPIV+1, fs, 1
          !FAC = front(CPIV,k)
          do j = p, fs, 1 !Doing Sector 4 update
             front(j-1,k-1) = front(j,k) - QQ(j)*GG(k)
          end do
          !front(p-1:fs-1,k-1) = front(p:fs,k) - QQ(p:fs)*GG(k)
          loadf(k-1) = loadf(k) - c_piv*GG(k)
       end do

       !$omp end parallel sections

       
       

       !Blank out the area outside the front
       do j = 1, fs, 1
          front(fs,j) = 0.0_rk
       end do
       do j = 1, fs-1, 1
          front(j,fs) = 0.0_rk
       end do
       loadf(fs) = 0.0_rk
       loadf_dum(fs) = 0.0_rk
       fs  = fs - 1

       !Shift global indices
       do j = RPIV, fs, 1
          iPg(j) = iPg(j+1)
       end do
       do j = CPIV, fs, 1
          jPg(j) = jPg(j+1)
       end do

    


    end do
  
  end subroutine elim_summedp

   subroutine multifront_old(L2res,t,mode)
    !L2res: returns this value
    !t: returns solve time, useful for comparing different params
    !mode: determines natural(1) or nested(2) sub ordering
    !If you overrode init_front with -1 or -2 ptn then no purpose

    !DO NOT CALCULATE L2RES outside of this subroutine
    !Just feed in L2res and it will be returned
    !If you must use the variable load_dum to calculate yourself
    implicit none
    real(kind=rk) :: front(iBW,iBW,iDM), front_m(iBW,iBW), loadf(iBW,iDM), loadf_m(iBW),&
         tp, ts, ti, loadf_dum(iBW,iDM), loadf_dumm(iBW), L2res, t1, t2, t3, t4, t,&
         front_n(iBW,iBW,NSi,ths), loadf_n(iBW,NSi,ths), loadf_dumn(iBW,NSi,ths)
    integer(kind=ik) :: NOPf(NE,bas), iPg(iBW,iDM), jPg(iBW,iDM), fs(iDM), mode, e_n(4),&
         iPg_m(iBW), jPg_m(iBW), fs_m, i, j, k, l, m, n, chead, rhead, l_i, i_ms, j_ms, id, step,&
         iPg_n(iBW,NSi,ths), jPg_n(iBW,NSi,ths), fs_n(NSi,ths), o, p

    !mode = 1 is nested natural ordering
    !mode = 2 is nested dissection (ptn should equal ptn2)
    t1 = REAL(omp_get_wtime(),rk)
    piv = 0

    
    !Blank global arrays
    loadc = 0.0_rk      !working array
    load_dum = 0.0_rk   !dummy right hand side
    load = 0.0_rk       !Solution vector

    front_m(:,:) = 0.0_rk  !Single front matrix
    loadf_m(:) = 0.0_rk    !single right hand
    loadf_dumm(:) = 0.0_rk !single dummy
    iPg_m = 0              !local to global rows
    jPg_m = 0              !local to glbal col
    fs_m = 0               !Size of front
    l_i = 0                !Variables removed from front
    NODf = 0               !Records variable appearances



    !Single front
    if (iDM.eq.1) then
       call Front_solve(fs_m,front_m,loadf_m,loadf_dumm,iPg_m,jPg_m,l_i,1,1,1,e_n)

    else !Dual front and nested fronts
       !Blank arrays for front that are domain wise
       front = 0.0_rk
       loadf = 0.0_rk
       loadf_dum = 0.0_rk
       iPg = 0
       jPg = 0
       fs = 0
      

       !Mode 1 is nested natural order
       if (mode.eq.1) then

          !PHASE 1: All smallest level Natural done using arbitray threads
          !$omp parallel private(id,j,l_i)
          !$omp do
          do j = 1, iDM, 1
             id = omp_get_thread_num()+1  !Thread id used to determine assembly arrays

             l_i = lt_i(j)                !initilaize var counter to predetermined starting point
             
             
             !Solve domain j and save front
             call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n)

          end do
          !$omp end do
          !$omp end parallel
          
       else

          !Phase 1.1: Outer two domains not nested order
          l_i = 0
          !$omp parallel private(id,j,l_i)
          !$omp do
          do j = 1, iDM, iDM-1
             id = omp_get_thread_num()+1

             fs(j) = 0
             l_i = lt_i(j)
             !write(*,*) 'Solving Domain:', j, fs(j)
             !write(*,*) 'Hi from thread 2:', omp_get_thread_num()
             call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),jPg(:,j),l_i,id,j,1,e_n)
             !write(*,*) 'Solved Domain:', j, fs(j)
          end do
          !$omp end do
          !$omp end parallel
         

          !PHASE 1.2: All smallest level Nested
          id = 1

          !$omp parallel private(id,j,l,l_i,m,o,p)
          !$omp do
          do j = 2, iDM-1, 1
             id = omp_get_thread_num()+1
             l_i = lt_i(j)
             !Initialize nested fronts within domains
             front_n(:,:,:,id) = 0.0_rk
             loadf_n(:,:,id) = 0.0_rk
             loadf_dumn(:,:,id) = 0.0_rk
             iPg_n(:,:,id) = 0
             jPg_n(:,:,id) = 0
             fs_n(:,id) = 0
             NODf_n(:,:,id) = 0

             !Determine domain type (shape, most are type 1)
             l = dm_type(j)

             !Assemble all smallest level nested regions (consist of 2x2 elements)
             do m = 1, NSi, 1
                if (eprs(1,m,j).eq.0) exit
                !write(*,*) eprs(:,m,j)
                if (m.eq.1) then
                   !write(*,*) 'hi', fs(j)
                   call Front_solve(fs(j),front(:,:,j),loadf(:,j),loadf_dum(:,j),iPg(:,j),&
                        jPg(:,j),l_i,id,j,3,eprs(:,m,j))
                   !write(*,*) 'bye'
                else
                   call Front_solve(fs_n(m,id),front_n(:,:,m,id),loadf_n(:,m,id),loadf_dumn(:,m,id),iPg_n(:,m,id),&
                        jPg_n(:,m,id),l_i,id,m,mode,eprs(:,m,j))
                end if
             end do

             !Use nested dissection to combine sub domains into 1 final level
             !Once done the remaining front resides in the domain front
             do m = 1, max_prs, 1
                o = dprs(1,m,l)
                p = dprs(2,m,l)
                if (o.eq.0) exit
                if (o.eq.1) then
                   call add_front(front_n(:,:,p,id),front(:,:,j),&
                        loadf_n(:,p,id),loadf(:,j),&
                        loadf_dumn(:,p,id),loadf_dum(:,j),&
                        iPg_n(:,p,id),iPg(:,j),&
                        jPg_n(:,p,id),jPg(:,j),fs_n(p,id),fs(j),&
                        p,j,3,id)
                   call elim_stitch(front(:,:,j),loadf(:,j),&
                        loadf_dum(:,j),&
                        iPg(:,j),jPg(:,j),l_i,fs(j),id)
                else
                   call add_front(front_n(:,:,p,id),front_n(:,:,o,id),&
                        loadf_n(:,p,id),loadf_n(:,o,id),&
                        loadf_dumn(:,p,id),loadf_dumn(:,o,id),&
                        iPg_n(:,p,id),iPg_n(:,o,id),&
                        jPg_n(:,p,id),jPg_n(:,o,id),fs_n(p,id),fs_n(o,id),&
                        p,o,mode,id)
                   call elim_stitch(front_n(:,:,o,id),loadf_n(:,o,id),&
                        loadf_dumn(:,o,id),&
                        iPg_n(:,o,id),jPg_n(:,o,id),l_i,fs_n(o,id),id)
                end if

             end do



          end do
          !$omp end do
          !$omp end parallel

       end if



       !PHASE 2: Sweep through domains from two fronts in natural order of domains (2 threads at most)
       !Not required for single domain front or dual domain front
     if (iDM.gt.2) then
        !$omp parallel sections private(j,k,l_i,i,l,id) num_threads(ths2)
       
        !$omp section
        id = omp_get_thread_num() + 1
        k = 1
        j = 1
        l_i = lt_i2(j)
        
        write(*,*) 'hi'
        do 
           k = k + 1
           i = dm_list(k,j)
           if (i.eq.0) exit

           call add_front(front(:,:,i),front(:,:,1),loadf(:,i),loadf(:,1),loadf_dum(:,i),&
                loadf_dum(:,1),iPg(:,i),iPg(:,1),jPg(:,i),jPg(:,1),fs(i),fs(1),i,1,1,1)
           call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
                iPg(:,1),jPg(:,1),l_i,fs(1),id)

        end do
        !$omp section
        id = omp_get_thread_num() + 1
        k = 1
        j = 2
        l_i = lt_i2(j)
       

        write(*,*) 'hi1'
        do 
           k = k + 1
           i = dm_list(k,j)
           if (i.eq.0) exit

           call add_front(front(:,:,i),front(:,:,iDM),loadf(:,i),loadf(:,iDM),loadf_dum(:,i),&
                loadf_dum(:,iDM),iPg(:,i),iPg(:,iDM),jPg(:,i),jPg(:,iDM),fs(i),fs(iDM),i,iDM,1,1)
           call elim_stitch(front(:,:,iDM),loadf(:,iDM),loadf_dum(:,iDM),&
                iPg(:,iDM),jPg(:,iDM),l_i,fs(iDM),id)
        end do

              write(*,*) 'hi2'
           
          !$omp end parallel sections
       end if
    
       !Determine which domains contain the remaining fronts
       ! k = 0
       ! do i = 1, iDM, 1
       !    k = k + 1
       !    if (dm_list(k,1).eq.0) exit
       ! end do
       if (iDM.eq.2) then
          !k = 2
          l_i = lt_i(3)
       else
          l_i = lt_i2(3)
       end if
       ! i = k - 1

       !Combine and eliminate final fronts
       ! call add_front(front(:,:,i+1),front(:,:,i),loadf(:,i+1),loadf(:,i),loadf_dum(:,i+1),&
       !      loadf_dum(:,i),iPg(:,i+1),iPg(:,i),jPg(:,i+1),jPg(:,i),fs(i+1),fs(i),i+1,i,1,1)
       ! call elim_stitch(front(:,:,i),loadf(:,i),loadf_dum(:,i),&
       !      iPg(:,i),jPg(:,i),l_i,fs(i),id)
         write(*,*) 'hi3'
       call add_front(front(:,:,iDM),front(:,:,1),loadf(:,iDM),loadf(:,1),loadf_dum(:,iDM),&
            loadf_dum(:,1),iPg(:,iDM),iPg(:,1),jPg(:,iDM),jPg(:,1),fs(iDM),fs(1),iDM,1,1,1)
       call elim_stitch(front(:,:,1),loadf(:,1),loadf_dum(:,1),&
            iPg(:,1),jPg(:,1),l_i,fs(1),id)

       
    end if
     
    !Perform back substitution
    call back_sub(mode)

    !Calculate L2res
    L2res = 0.0_rk
    !$omp parallel private(i) reduction(+:L2res)
    !$omp do
    do i = 1, NV, 1
       L2res = L2res + (load_dum(i))**2
    end do
    !$omp end do
    !$omp end parallel
    L2res = sqrt(L2res)

    !Save total runtime
    t = REAL(omp_get_wtime(),rk) - t1
    j = 0
    do i = 1, ths, 1
       j = piv(i) + j
    end do
    write(*,*) 'Small Pivots:', j

     
    
  end subroutine multifront_old



end module front_mod

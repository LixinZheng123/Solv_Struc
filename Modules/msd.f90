!
!  This module is to find first solvation shell (FSS) of a defect from a MD simulation trajectory.
!  Lixin Zheng, lixin.zheng@temple.edu
!  May 2014.
!  Modified 20140901, add error message.
!  Modified 20141014, refactor with rdf codes.
!  Refactored 20141118
!----------------------------------------
! 
MODULE find_fss

  !
  USE main_variables,    ONLY : DP
  !
  !
  SUBROUTINE find_hydroxide(spc1,spc2,ion1,ion2,cnt_ion,i_ion1,i_ion2)
    !**************************************************************
    !**************************************************************
    !
    !  ONLY numION, INDEX iIO(1), iiH(1:numION,1), AND COORDINATES rIO(1:3,1) NEED TO BE PASSED TO NEXT MODULE.
    !
    !**************************************************************
    !**************************************************************
    IMPLICIT NONE
    !
    real, intent(inout) :: spc1(3,nsp(1))
    real, intent(inout) :: spc2(3,nsp(2))
    INTEGER, intent(out):: cnt_ion=0
    INTEGER             :: numW=0
    INTEGER             :: e_step=0
    INTEGER             :: i1,i2
    INTEGER             :: p=0
    REAL                :: d

    do i1=1,nsp(1)
      q=0
      do i2=1, nsp(2)
        CALL get_distance(spc1(1:3,i1),spc2(1:3,i2),d)
        if (d*convertBA .lt. r_cov) then
          q=q+1
          iiH(q,p+1)=iH
          bond(q,p+1)=iH
        endif
      enddo
      !*************************
      if (q .eq. 2) then
        numW=numW+1
      else if (q .eq. numH) then
        p=p+1
        iIO(p)=iO
      end if
      !**************************
    enddo
    numION=p
    do i=1, numION
      iIO_record(i)=iIO(i)
    enddo
    !**************************
    !1. The most straight forward situation
    !**************************
    if ((numION+numW) .eq. nsp(1) .and. numION .eq. 1) then
      do j=1,3
        rIO(j,1) = rO(j,iIO(1))
      enddo
    endif
    !**************************
    !2. The different situations that would probably go wrong.
    !Just skip the judgement part, and adopt the previous ion index.
    !**************************
    if ((numION+numW) .ne. nsp(1) .or. numION .gt. 2 .or. numION .lt. 1) then 
      !if ((numION+numW) .ne. nsp(1) .or. numION .gt. 2) then
      !  write(*,*) time, numION, numION+numW
      !endif
      iIO(1)=iIOp(1)
      do i=1,numH
        iiH(i,1)=iiHp(i)
      enddo
      do j=1,3
        rIO(j,1)=rO(j,iIO(1))
      enddo
    endif
    !**************************
    !3. When numION=2, try to find the real ion
    !**************************
    if ((numION+numW) .eq. nsp(1) .and. numION .eq. 2) then
      do p=1,2
        do j=1,3
          rIO(j,p) = rO(j,iIO(p))
        enddo
      enddo
      !**************************
      !compare and change order
      !so that the total distances from O to H of iIO(1) is the shortest
      !If cs=2, OH-
      if (cs .eq. 2) then
        do p=1,numION
          CALL get_distance(rIO(1:3,p),rH(1:3,iiH(1,p)),dr(p))
        enddo
        do p=2,numION
          if(dr(p) .lt. dr(1)) then
           iIO(1)=iIO(p)
            dr(1)=dr(p)
            do j=1,3
              rIO(j,1)=rIO(j,p)
              iiH(1,1)=iiH(1,p)
            enddo
          endif
        enddo
      !IF cs=3, Find the mutual iiH that ion1 and ion2 has
      else if (cs .eq. 3) then
        do q=1,3
          do q1=1,3
            if (iiH(q,1) .eq. iiH(q1,2)) then
              iiH_real=iiH(q,1)
              exit
            endif
          enddo
        enddo
        do p=1,2
          CALL get_distance(rIO(1:3,p),rH(1:3,iiH_real),dr(p))
        enddo
        do p=2,numION
          if(dr(p) .lt. dr(1)) then
           iIO(1)=iIO(p)
            dr(1)=dr(p)
            do j=1,3
              rIO(j,1)=rIO(j,p)
              do q=1,numH
                iiH(q,1)=iiH(q,p)
              enddo
            enddo
          endif
        enddo
      endif
    endif
    !
    !
    !***
    !Modified 20140901
    !To determine if the ion has transfered too fast
    error=1
    !
    if (ncount .gt. 1) then
      do i=1,new_nsp(1)
        if (iIO(1) .eq. iO_FSS(i) ) then
          error=0
          exit
        endif
      enddo
    endif
    if (error .eq. 1 .and. ncount .gt. 1) then
      write(5,*) "Error1: The new ion has gone out of the FSS of previous ion."
      e_step=step
      write(5,*) step,time,iIO(1) 
    endif
    !***
    !
    !
    !**************************
    !Write in .ion
    !**************************
    write(2,*) step,time,numION
    write(2,*) iIO(1),iiH(1:numH,1)
    write(2,*) rIO(1:3,1)
    !
    !**************************
    !**************************
    if (iIO(1) .ne. iIOp(1)) then
      !Write in .trans
      write(4,*) step, iIO(1), iiH(1:numH,1)
      !Modified 20141014
      !We also want to write out the non-rattling-transfer index in the same time.
      if (iIO(1) .ne. iIOp(2)) then
        write(10,*) step, iIO(1)
      endif
      iIOp(2)=iIOp(1)
      iIOp(1)=iIO(1)
      !End of modification
    endif
    do i=1,numH
      iiHp(i)=iiH(i,1)
    enddo
    !
    !
    !
  end SUBROUTINE find_hydroxide

  SUBROUTINE find_hydronium()
  
  end SUBROUTINE find_hydronium
  
  SUBROUTINE FSS()
    !**************************************************************
    !**************************************************************
    !
    !MODULE 2: FIND THE 1ST SOLVATION SHELL 
    !ONLY numION, INDEX iIO(1), iiH(1:numION,1), AND COORDINATES rIO(1:3,1) NEED TO BE PASSED FROM LAST MODULE.
    !
    !**************************************************************
    !**************************************************************
    !
    !**************************
    !This if statement is to make sure that the 1st ion has been found;
    !If not, we will write nothing to the .fss file.
    !**************************
    if (start .eq. 0 .and. numION .eq. 0) then
      cycle
    else if (start .eq. 0 .and. numION .ne. 0) then
      start=1
    endif
    !**************************
    !
    !**************************
    !Find FSS oxygens
    !**************************
    p=0
    q=0
    cov_bond=0
    do i=1,nsp(1)
      CALL get_distance(rO(1:3,i),rIO(1:3,1),d)
      if (d*convertBA .lt. r_FSS) then
        p=p+1  
        iO_FSS(p)=i
        do j=1,3
          rO_FSS(j,p)=rO(j,i)
        enddo
      endif
    enddo
    new_nsp(1)=p
    !**************************
    !Find the hydrogens relates to FSS oxygens
    !**************************
    do i=1,nsp(2)
      CALL get_distance(rH(1:3,i),rIO(1:3,1),dr(1))
      do p=1,new_nsp(1)
        CALL get_distance(rH(1:3,i),rO_FSS(1:3,p),dr(2))
        if (dr(2)*convertBA .lt. r_cov .or. dr(1)*convertBA .lt. (r_FSS-r_cov)) then
          q=q+1
          iH_FSS(q)=i
          do j=1,3
            rH_FSS(j,q)=rH(j,i)
          enddo
          exit
        endif     
      enddo
    enddo
    new_nsp(2)=q
    !**************************
    !Record the link between O and H
    !**************************
    do p=1,new_nsp(1)
      k=0
      do q=1, new_nsp(2) 
        CALL get_distance(rH_FSS(1:3,q),rO_FSS(1:3,p),d)
        if (d*convertBA .lt. r_cov) then
          k=k+1
          cov_bond(k,p)=iH_FSS(q)
        endif
      enddo
    enddo
    !**************************
    !Write index into file.fss 
    !**************************
    if (numION .eq. 1) write(3,*) step, time, numION, new_nsp(1),new_nsp(2)
    if (numION .ne. 1) write(3,*) step, time, numION, new_nsp(1),new_nsp(2), iIO_record(1:numION)
    write(3,*) iIO(1),iiH(1:numH,1)
    do p=1,new_nsp(1)
      write(3,*) iO_FSS(p),rO_FSS(1:3,p),cov_bond(1:4,p)
    enddo
    do q=1,new_nsp(2)
      write(3,*) iH_FSS(q),rH_FSS(1:3,q)
    enddo
    !
    !
    !Modified 20140901
    !if (e_step .eq. step) then
      if (cs .eq. 2) then
        if ((new_nsp(1)*2-1) .ne. new_nsp(2) .and. numION .ne. 0) then
          write(5,*) "Error2: (nsp(1)*2-1) .ne. nsp(2)!"
          write(5,*) step, time, "numION:", numION, new_nsp(1),new_nsp(2)
          write(5,*) "iIO:",iIO_record(1:numION)
        endif
      endif
      if (cs .eq. 3) then
        if ((new_nsp(1)*2+1) .ne. new_nsp(2) .and. numION .ne. 0) then
          write(5,*) "Error2: (nsp(1)*2+1) .ne. nsp(2)!"
          write(5,*) step, time, "numION:", numION, new_nsp(1),new_nsp(2) 
          write(5,*) "iIO:",iIO_record(1:numION)
        endif
      endif
    !endif
    ncount=ncount+1
  enddo
  
  !
  END SUBROUTINE FSS
  !
END MODULE find_fss

!
!  This module is to find first solvation shell (FSS) of a defect from a MD simulation trajectory.
!  Lixin Zheng, lixin.zheng@temple.edu
!  May 2014.
!  Modified 20140901, add error message.
!  Modified 20141014, refactor with rdf codes.
!  Refactored 20141118
!----------------------------------------
! 
MODULE defect
  !
  !
  USE main_variables,  ONLY : DP
  USE gen_lib
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !---------------------------------------------------
    SUBROUTINE find_hydroxide(spc1,spc2,cov_bond,Ion1,Ion2,num_I,i_Ion1,i_Ion2,error)
      !
      USE main_variables,  ONLY : nsp,r_FSS
      !
      IMPLICIT NONE
      !
      REAL(DP), intent(in)           :: spc1(3,nsp(1))
      REAL(DP), intent(in)           :: spc2(3,nsp(2))
      !REAL(DP), intent(out),  save   :: Ion1(3)         !Coordinate of 1st species ion
      REAL(DP), intent(out)          :: Ion1(3)         !Coordinate of 1st species ion
      REAL(DP), intent(out)          :: Ion2(3)         !Coordinate of 2nd species ion
      REAL(DP)                       :: Ion1_prvs(3)    !Coordinate of 1st species ion in the previous step
      REAL(DP)                       :: Ion1_record(3,4)!Coordinate of 1st species ion candidates
      REAL(DP)                       :: Ion2_record(3,4)!Coordinate of 2nd species ion candidates
      INTEGER,  intent(out)          :: num_I           !Count of 1st species ion
      !INTEGER,  intent(out),  save   :: i_Ion1          !Index of 1st species ion
      INTEGER,  intent(out)          :: i_Ion1          !Index of 1st species ion
      !INTEGER,  intent(out),  save   :: i_Ion2          !Index of 2nd species ion
      INTEGER,  intent(out)          :: i_Ion2          !Index of 2nd species ion
      INTEGER                        :: i_Ion1_record(4)!Index of 1st species ion candidates
      INTEGER                        :: i_Ion2_record(4)!Index of 2nd species ion candidates
      INTEGER,  intent(out)          :: cov_bond(4,nsp(1))
      INTEGER                        :: num_cov
      INTEGER                        :: num_W=0         !Count of normal waters
      INTEGER                        :: i,j,ii
      INTEGER                        :: error
      REAL(DP)                       :: d, dr(3)
      !
      error=0
      num_I=0
      !
      do j=1,3
        Ion1_prvs(j)=Ion1(j)
      enddo
      !
      !
      do i=1,nsp(1)
        CALL find_cov(spc1(1:3,i),spc2(1:3,1:nsp(2)),cov_bond(1:4,i),num_cov)
        if (num_cov .eq. 2) num_W=num_W+1
        if (num_cov .eq. 1) then
          num_I=num_I+1
          i_Ion1_record(num_I)=i
          i_Ion2_record(num_I)=cov_bond(num_I,i)
        end if
      enddo
      !
      !
      if ((num_I+num_W) .eq. nsp(1)) then
        !
        do i=1,num_I
          do j=1,3
            Ion1_record(j,i) = spc1(j,i_Ion1_record(i))
          enddo
        enddo
        !
        !
        !1. The most straight forward situation
        if (num_I .eq. 1) ii=1
        !
        !2. When number of ion is 2, try to find the real ion
        if (num_I .eq. 2) then
          ii=1
          do i=1,num_I
            CALL get_distance(Ion1_record(1:3,i),Ion2_record(1:3,i),dr(i))
          enddo
          if(dr(2) .lt. dr(1)) ii=2
          !
        endif
        !
        !
        i_Ion1=i_Ion1_record(ii)
        i_Ion2=i_Ion2_record(ii)
        do j=1,3
          Ion1(j)=Ion1_record(j,ii)
        enddo
        !
      endif
      !
      !3. The different situations that would probably go wrong.
      !Just skip the judgement part, and adopt the previous ion index.
      if ((num_I+num_W) .ne. nsp(1) .or. num_I .gt. 2 .or. num_I .lt. 1) then
        do j=1,3
          Ion1(j)=spc1(j,i_Ion1)
          Ion2(j)=spc2(j,i_Ion2)
        enddo
      endif
      !
      !
      !Determine if the ion has transfered too fast
      CALL get_distance(Ion1(1:3),Ion1_prvs(1:3),d)
      if (d .gt. r_FSS) then
        error=1
      endif
      !
      !
      return
      !
    end SUBROUTINE find_hydroxide
    !---------------------------------------------------
    !
    !---------------------------------------------------
    SUBROUTINE find_hydronium(spc1,spc2,cov_bond,Ion1,Ion2,num_I,i_Ion1,i_Ion2,error)
      !
      USE main_variables,  ONLY : nsp,r_FSS
      !
      IMPLICIT NONE
      !
      REAL(DP), intent(in)           :: spc1(3,nsp(1))
      REAL(DP), intent(in)           :: spc2(3,nsp(2))
      !REAL(DP), intent(out),  save   :: Ion1(3)         !Coordinate of 1st species ion
      REAL(DP), intent(out)          :: Ion1(3)         !Coordinate of 1st species ion
      REAL(DP), intent(out)          :: Ion2(3,3)       !Coordinate of 2nd species ion
                                                        ! (num_H, dimension)
      REAL(DP)                       :: Ion1_prvs(3)    !Coordinate of 1st species ion in the previous step
      REAL(DP)                       :: Ion1_record(3,4)!Coordinate of 1st species ion candidates
                                                        ! (dimension,num_I)
      REAL(DP)                       :: Ion2_record(3,3,4)!Coordinate of 2nd species ion candidates
                                                        ! (num_H, dimension,num_I)
      INTEGER,  intent(out)          :: num_I           !Count of 1st species ion
      !INTEGER,  intent(out),  save   :: i_Ion1          !Index of 1st species ion
      INTEGER,  intent(out)          :: i_Ion1          !Index of 1st species ion
      !INTEGER,  intent(out),  save   :: i_Ion2(3)       !Index of 2nd species ion
      INTEGER,  intent(out)          :: i_Ion2(3)       !Index of 2nd species ion
      INTEGER                        :: i_Ion2_trans    !Index of the transfering 2nd species ion
      INTEGER                        :: i_Ion1_record(4)!Index of 1st species ion candidates
      INTEGER                        :: i_Ion2_record(3,4)!Index of 2nd species ion candidates
                                                        ! (num_H, num_I)
      INTEGER,  intent(out)          :: cov_bond(4,nsp(1))
      INTEGER                        :: num_cov
      INTEGER                        :: num_W=0         !Count of normal waters
      INTEGER                        :: i,j,k,ii,k1
      INTEGER                        :: error
      REAL(DP)                       :: d, dr(3)
      !
      error=0
      num_I=0
      !
      do j=1,3
        Ion1_prvs(j)=Ion1(j)
      enddo
      !
      !
      do i=1,nsp(1)
        CALL find_cov(spc1(1:3,i),spc2(1:3,1:nsp(2)),cov_bond(1:4,i),num_cov)
        if (num_cov .eq. 2) num_W=num_W+1
        if (num_cov .eq. 3) then
          num_I=num_I+1
          i_Ion1_record(num_I)=i
          do k=1,3
            i_Ion2_record(k,num_I)=cov_bond(num_I,i)
          enddo
        end if
      enddo
      !
      !
      if ((num_I+num_W) .eq. nsp(1)) then
        !
        do i=1,num_I
          do j=1,3
            Ion1_record(j,i) = spc1(j,i_Ion1_record(i))
          enddo
        enddo
        !
        !
        !1. The most straight forward situation
        if (num_I .eq. 1) ii=1
        !
        !2. When number of ion is 2
        if (num_I .eq. 2) then
          ii=1
          !Find the mutual iiH that ion1 and ion2 has
          do k=1,3
            do k1=1,3
              if (i_Ion2_record(k,1) .eq. i_Ion2_record(k1,2)) then
                i_Ion2_trans=i_Ion2_record(k,1)
                exit
              endif
            enddo
          enddo
          !
          do i=1,2
            CALL get_distance(Ion1_record(1:3,i),spc2(1:3,i_Ion2_trans),dr(i))
          enddo
          if(dr(2) .lt. dr(1)) ii=i
          !
        endif
        !
        !
        i_Ion1=i_Ion1_record(ii)
        do k=1,3
          i_Ion2(k)=i_Ion2_record(k,ii)
        enddo
        !
        do j=1,3
          Ion1(j)=Ion1_record(j,ii)
          do k=1,3
            Ion2(k,j)=Ion2_record(k,j,ii)
          enddo
        enddo
        !
      endif
      !
      !3. The different situations that would probably go wrong.
      !Just skip the judgement part, and adopt the previous ion index.
      if ((num_I+num_W) .ne. nsp(1) .or. num_I .gt. 2 .or. num_I .lt. 1) then
        do j=1,3
          Ion1(j)=spc1(j,i_Ion1)
          do k=1,3
            Ion2(k,j)=spc2(j,i_Ion2(k))
          enddo
        enddo
      endif
      !
      !
      !Modified 20140901
      !To determine if the ion has transfered too fast
      CALL get_distance(Ion1(1:3),Ion1_prvs(1:3),d)
      if (d .gt. r_FSS) then
        error=1
      endif
      !
      !
      return
      !
    end SUBROUTINE find_hydronium
    !---------------------------------------------------
    
    !---------------------------------------------------
    SUBROUTINE FSS(Ion1,spc1,spc2,cov_bond,indexs,num_fss)
      !
      !
      USE main_variables,  ONLY : nsp,r_FSS,convertBA
      !
      IMPLICIT NONE
      !
      REAL(DP), intent(in)           :: spc1(3,nsp(1))
      REAL(DP), intent(in)           :: spc2(3,nsp(2))
      REAL(DP), intent(in)           :: Ion1(3)         !Coordinate of 1st species ion
                                                        ! (num_H, dimension)
      INTEGER,  intent(in)           :: cov_bond(1:4,nsp(1))
      INTEGER,  intent(out)          :: indexs(2,30)    !index of 2 fss species
      INTEGER,  intent(out)          :: num_fss(2)      !number of 2 fss species
      INTEGER                        :: i,i1,k
      INTEGER                        :: error=0
      INTEGER                        :: switch=0
      REAL(DP)                       :: d, dr(3)
      !
      num_fss=0
      !
      !
      do i=1,nsp(1)
        CALL get_distance(spc1(1:3,i),Ion1(1:3),d)
        if (d*convertBA .lt. r_FSS) then
          num_fss(1)=num_fss(1)+1
          indexs(1,num_fss(1))=i
        endif
      enddo
      !
      !
      !Find the hydrogens relates to FSS oxygens
      do i=1,nsp(2)
        CALL get_distance(spc2(1:3,i),Ion1(1:3),d)
        !2 cases a hydrogen can be included in fss
        !1st, distance is within cutoff
        if (d*convertBA .lt. r_FSS) then
          switch=1
        !2nd, hydrogen is covelent-bonded to a fss oxygen
        else
          do i1=1,num_fss(1)
            do k=1,4
              if (cov_bond(k,indexs(1,i1)) .eq. 0) exit
              if (i .eq. cov_bond(k,indexs(1,i1))) switch=1
            enddo
          enddo
        endif
        if (switch .eq. 1) then
          num_fss(2)=num_fss(2)+1
          indexs(2,num_fss(2))=i
        endif
      enddo
      !
      !
      return
      !
    END SUBROUTINE FSS
    !
END MODULE defect

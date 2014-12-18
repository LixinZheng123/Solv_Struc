!
!  Refactored 20141118
!  Lixin Zheng, lixin.zheng@temple.edu
!----------------------------------------
! 
MODULE gen_lib
  !
  !
  USE main_variables,    ONLY : DP
  IMPLICIT NONE
  !
  CONTAINS
  !
  !---------------------------------------------------
  SUBROUTINE find_cov(spc1,spc2,cov,num)
    !
    USE main_variables,  ONLY : nsp,convertBA,r_cov
    !
    IMPLICIT NONE
    !
    REAL(DP), intent(in)    :: spc1(3)            !Coordinate of the species that want to find covalent bond
    REAL(DP), intent(in)    :: spc2(3,nsp(2))     !Coordinate of the 2nd species
    INTEGER,  intent(out)   :: cov(4)             !Index of spc2 that are covalent-bonded to spc1
    INTEGER,  intent(out)   :: num                !Number of spc2 that are covalent-bonded to spc1
    REAL(DP)                :: dis                !Distance
    INTEGER                 :: i,n
    !
    cov=0
    n=0
    !
    do i=1, nsp(2)
      CALL get_distance(spc1(1:3),spc2(1:3,i),dis)
      if (dis*convertBA .lt. r_cov) then
          n=n+1
          cov(n)=i
        endif
    enddo
    !
    return
    !
  end SUBROUTINE find_cov
  !---------------------------------------------------
  !
  !
  !
  !---------------------------------------------------
  SUBROUTINE get_distance(pos1,pos2,dist)
    
    USE main_variables,    ONLY : celldm
    IMPLICIT NONE
    !
    REAL(DP), intent(in)    :: pos1(3),pos2(3)
    REAL(DP), intent(inout) :: dist
    REAL(DP)                :: delta(3)
    INTEGER                 :: i
    !
    do i=1,3
      delta(i)=pos1(i)-pos2(i)
      delta(i)=delta(i)-nint(delta(i)/celldm)*celldm
    enddo
    dist=sqrt(delta(1)**2+delta(2)**2+delta(3)**2)
    !
    return
    !
  END SUBROUTINE get_distance
  !---------------------------------------------------
  !
  !
  !
  !---------------------------------------------------
  SUBROUTINE EOF_PrintOut(err,k)
    !
    IMPLICIT NONE
    !
    INTEGER, intent(in)    :: err,k
    !
    write(*,*) '  End of File Reached'
    write(*,fmt='(1X, "  Total number of Samples: ", I7)' )(k)
    ! 
    return
    !
  END SUBROUTINE EOF_PrintOut
  !---------------------------------------------------
  !
  !
  !
  !---------------------------------------------------
END MODULE gen_lib

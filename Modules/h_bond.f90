!
!  This module is to find hydrogen bonds (HB) of a given defect, or of the whone system.
!  Lixin Zheng, lixin.zheng@temple.edu
!  May 2014.
!  Refactored 20141118
!----------------------------------------
! 
MODULE hydrogen_bond
  !
  !
  USE main_variables,  ONLY : DP
  !
  !
  !---------------------------------------------------
  SUBROUTINE find_hb((spc1,spc2,hb,num)
    !
    USE main_variables,  ONLY : convertBA,r_hb,angle_hb
    USE fss,             ONLY : new_nsp
    !
    IMPLICIT NONE
    !
    REAL(DP), intent(in)    :: spc1(3)            !Coordinate of the species that want to find covalent bond
    REAL(DP), intent(in)    :: spc2(3,new_nsp(2)) !Coordinate of the 2nd species
    INTEGER,  intent(out)   :: hb(4)              !Index of spc2 that are covalent-bonded to spc1
    INTEGER,  intent(out)   :: num                !Number of spc2 that are covalent-bonded to spc1
    REAL(DP)                :: dis                !Distance
    INTEGER                 :: i
    !
    !cov=0
    !n=0
    !!
    !do i=1, nsp(2)
    !  CALL get_distance(spc1(1:3),spc2(1:3,i),dis)
    !  if (dis*convertBA .lt. r_cov) then
    !      n=n+1
    !      cov(n)=i
    !    endif
    !enddo
    !
  end SUBROUTINE find_hb
  !---------------------------------------------------
  !
  !
end MODULE hydrogen_bond

!
!  Lixin Zheng, lixin.zheng@temple.edu
!  Refactored 20141118
!----------------------------------------
!
PROGRAM msd_main
  !
  USE main_variables
  USE gen_lib
  USE find_fss,        only : find_hydroxide, find_hydronium, FSS
  !
  IMPLICIT NONE
  !
  !======Global variables======
  INTEGER              :: ncount=0
  INTEGER              :: step 
  INTEGER              :: ierror
  REAL                 :: time
  !======Variables concerning .pos======
  REAL,ALLOCATABLE     :: rO(:,:),rH(:,:)
  !======Variables concerning .fss======
  !REAL                 :: rO_FSS(3,15)         !atomic positions (dimension,atomic-species)
  !REAL                 :: rH_FSS(3,30)         !atomic positions (dimension,atomic-species)
  !REAL                 :: rIO_record(3,5)      !Ion position (dimension,numION)
  !REAL                 :: rIH_record(3,3)      !ion positions (dimension,numH)
  !REAL                 :: rIO(3)               !Real ion position (dimension)
  !REAL                 :: rIH(3,3)             !Real ion position (dimension,numH)
  !INTEGER              :: new_nsp(2)           !number of each atomic species
  !INTEGER              :: iO_FSS(15)
  !INTEGER              :: iH_FSS(30)
  !INTEGER              :: iIH_record(5,4)      !index of Ion hydrogen (numH,numION)
  !INTEGER              :: iIO_record(5)        !index of O* read from Ion
  !INTEGER              :: iIO                  !index of Ion oxygen
  !INTEGER              :: iIH(3)               !index of Ion oxygen
  !INTEGER              :: numION               !Number of Ions in each timestep. 
  !INTEGER              :: numH                 !Number of Hydrogen in each ion. cs=2, numH=1; cs=3, numH=3
  !INTEGER              :: iIOp(2),iiHp(4)
  !INTEGER              :: iiH_trans            !Find the hydrogen, through which PT occurs.
  !INTEGER              :: cov_bond(4,15)       !The covalent bonded hydrogen to each oxygen
  !======Local variables======
  REAL                 :: d,dr(5)
  INTEGER              :: e_step
  INTEGER              :: j,k,p,q,q1
  INTEGER              :: start=0
  INTEGER              :: numW                 !number of water
  INTEGER              :: error
  !
  namelist /input/ cs,filename, nsp, stepstart, stepstop
  !
  !initialization
  CALL init
  !
  !Open files
  CALL files(1)
  !
  !************************************
  !main loop
  main_loop: do 
    !end of file
    read(1,*,iostat=ierror) step, time 
    if (ierror .lt. 0) then
      CALL EOF_PrintOut(ierror,ncount)
      exit main_loop
    endif
    !
    !READ FILE .POS
    CALL read_file(rO(1:3,nsp(1)),rH(1:3,nsp(2)))
    !
    if (cs .eq. 2) then
      CALL find_hydroxide(rO(1:3,nsp(1)),rH(1:3,nsp(2)),rIO(1:3),rIH(1:3,1),numION,iIO,iIH(1))
    else if (cs .eq. 3) then
      CALL find_hydronium(rO(1:3,nsp(1)),rH(1:3,nsp(2)),rIO(1:3),rIH(1:3,1:3),numION,iIO,iIH(1:3))
    endif
    !
    CALL find_FSS(rIO(1:3),rIH(1:3,1:numH),rO(1:3,nsp(1)),rH(1:3,nsp(2)))
    !
  enddo main_loop
  !************************************
  !
  !close files
  CALL files(2)
  !
  !
  ! 
  contains
  !
  !=========================================================================================
  !=========================================================================================
  !
  SUBROUTINE init
    !
    IMPLICIT NONE
    !read Namelist input for stdin
    read(*,input)
    !
    !print Intro
    write(*,*) '--------------------------------------------------'
    write(*,*) '|     Find the Ion and First Solvation Shell      |'
    write(*,*) '--------------------------------------------------'
    write(*,*) 'Position file      : ', trim(filename)//'.pos'
    write(*,*) 'Defect file        : ', trim(filename)//'.ion'
    write(*,*) '1st solvation file : ', trim(filename)//'.fss'
    write(*,*) 'Ion transfer file  : ', trim(filename)//'.trans ',trim(filename)//'.trans_no_r'
    write(*,*) 'Ion error file     : ', trim(filename)//'.error'
    !
    !allocation
    allocate(rO(3,nsp(1)))
    allocate(rH(3,nsp(2)))
    !Initiation
    iIO=0
    iIO_record=0
    iIH_record=0
    iIOp=0
    iiHp=0
    if (cs .eq. 2) numH=1
    if (cs .eq. 3) numH=3
    !
  end SUBROUTINE init
  !
  !=========================================================================================
  !=========================================================================================
  !
  SUBROUTINE files(io)
    !
    IMPLICIT NONE
    !
    INTEGER, intent(in)  :: io
    !
    select case(io)
      case(1)
        open(unit=1, file=(trim(filename)//'.pos'), status='old')
        open(unit=2, file=(trim(filename)//'.ion'), status='unknown')
        open(unit=3, file=(trim(filename)//'.fss'), status='unknown')
        open(unit=4, file=(trim(filename)//'.trans'), status='unknown')
        open(unit=5, file=(trim(filename)//'.error'), status='unknown')
        open(unit=10, file=(trim(filename)//'.trans_no_r'), status='unknown')
        !
      case(2)
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(10)
      case default
    end select
    !
  end SUBROUTINE files
  !
END PROGRAM find_fss

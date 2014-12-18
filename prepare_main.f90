!
!  Lixin Zheng, lixin.zheng@temple.edu
!  Refactored 20141118
!----------------------------------------
!
PROGRAM prepare
  !
  USE main_variables
  USE gen_lib
  USE defect
  !USE h_bond,          only :
  !
  IMPLICIT NONE
  !
  !======Global variables======
  INTEGER              :: ncount=0
  INTEGER              :: step 
  INTEGER              :: ierror
  REAL                 :: time
  !======Local variables======
  REAL                 :: d
  !INTEGER              :: e_step
  INTEGER              :: iIOp(2)              !Previous ion oxygen index
  INTEGER              :: i,j,k
  INTEGER              :: start=0              !To determine if should write fss.
  INTEGER              :: numW                 !number of water
  INTEGER              :: numH                 !number of hydrogen in a ion
                                               !For cs=2, numH=1
                                               !For cs=3, numH=3
  INTEGER              :: errmsg               !Error message
  !INTEGER              :: i
  !
  namelist /input/ cs,filename, nsp, stepstart, stepstop
  !
  !initialization
  CALL init
  !
  !Open files
  CALL files(1)
  !
  !========================================
  !main loop
  main_loop: do 
    !end of file
    read(pos_id,*,iostat=ierror) step, time 
    if (ierror .lt. 0) then
      CALL EOF_PrintOut(ierror,ncount)
      exit main_loop
    endif
    !
    !-------------------------------------
    !read file .pos
    do i=1,nsp(1)
      read(pos_id,*) rO(1:3,i)
    enddo
    do i=1,nsp(2)
      read(pos_id,*) rH(1:3,i)
    enddo
    !-------------------------------------
    !
    !-------------------------------------
    !find defect
    if (cs .eq. 2) then
      CALL find_hydroxide(rO(1:3,1:nsp(1)),rH(1:3,1:nsp(2)),covalent(1:4,1:nsp(1)), &
                          rIO(1:3),rIH(1:3,1),numION,iIO,iIH(1),errmsg)
    else if (cs .eq. 3) then
      CALL find_hydronium(rO(1:3,1:nsp(1)),rH(1:3,1:nsp(2)),covalent(1:4,1:nsp(1)), &
                          rIO(1:3),rIH(1:3,1:3),numION,iIO,iIH(1:3),errmsg)
    endif
    !
    if (errmsg .eq. 1) then
      write(err_id,*) "Error1: Distance of new and old ion too big!"
      !e_step=step
      write(err_id,*) step,time,iIO
    endif
    !-------------------------------------
    !
    !
    !
    !-------------------------------------
    !write file .ion
    write(ion_id,*) step,time,numION
    write(ion_id,*) iIO,iIH(1:numH)
    write(ion_id,*) rIO(1:3)
    !write file .trans
    if (iIO .ne. iIOp(2)) then
      write(trans_id,*)
    endif
    iIOp(2)=iIOp(1)
    iIOp(1)=iIO
    !-------------------------------------
    !
    !
    !
    !-------------------------------------
    !find first solvation shell atoms
    CALL find_FSS(rIO(1:3),rO(1:3,1:nsp(1)),rH(1:3,1:nsp(2)),covalent(1:4,1:nsp(1)), &
                  i_FSS(1:2,1:30),new_nsp(1:2))
    !
    if (errmsg .eq. 1) then
      if (cs .eq. 2) then
        write(err_id,*) "Error2: (nsp(1)*2-1) .ne. nsp(2)!"
        write(err_id,*) step, time, "numION:", numION, new_nsp(1),new_nsp(2)
        write(err_id,*) "iIO:",iIO
      else if (cs .eq. 3) then
        write(err_id,*) "Error2: (nsp(1)*2+1) .ne. nsp(2)!"
        write(err_id,*) step, time, "numION:", numION, new_nsp(1),new_nsp(2)
        write(err_id,*) "iIO:",iIO
      endif
    endif
    !-------------------------------------
    !
    !
    !
    !-------------------------------------
    !write file .fss
    if (start .eq. 0 .and. numION .eq. 0) then
      !This is to make sure that the 1st ion has been found;
      !If not, we will write numION=0, new_nsp=0 in .fss file.
      write(fss_id,*) step,time,"0","0","0"
      cycle
    else if (start .eq. 0 .and. numION .ne. 0) then
      start=1
    endif
    !
    write(fss_id,*) step, time, numION, new_nsp(1:2)
    if (cs .eq. 2) write(fss_id,*) iIO,iIH(1)
    if (cs .eq. 3) write(fss_id,*) iIO,iIH(1:3)
    do i=1,new_nsp(1)
      write(fss_id,*) i_FSS(1,i),rO(1:3,i_FSS(1,i),covalent(1:4,i_FSS(1,i))
    enddo
    do i=1,new_nsp(2)
      write(fss_id,*) i_FSS(2,i),rH(1:3,i_FSS(2,i)
    enddo
    !-------------------------------------
    !
    !
    !
  enddo main_loop
  !========================================
  !
  !close files
  CALL files(2)
  !
  !
  !
  !
  !
  !
  !
  !
  ! 
  CONTAINS
  !
  !
  !
  SUBROUTINE initialization
    !
    IMPLICIT NONE
    !read Namelist input for stdin
    read(*,input)
    !
    !print Intro
    write(*,*) '-------------------------------------------'
    write(*,*) '|     Preparation for furture analysis     |'
    write(*,*) '-------------------------------------------'
    write(*,*) 'Position file      : ', trim(filename)//'.pos'
    write(*,*) 'Defect file        : ', trim(filename)//'.ion'
    write(*,*) '1st solvation file : ', trim(filename)//'.fss'
    write(*,*) 'Ion transfer file  : ', trim(filename)//'.trans ',trim(filename)//'.trans_no_r'
    write(*,*) 'HB file            : ', trim(filename)//'.hbcase'
    write(*,*) 'Error file         : ', trim(filename)//'.error'
    !
    if (cs .eq. 2) numH=1
    if (cs .eq. 3) numH=3
    
  end SUBROUTINE initialization
  !
  !
  !
  SUBROUTINE files(io)
    !
    IMPLICIT NONE
    !
    INTEGER, intent(in)  :: io
    !
    select case(io)
      case(1)
        pos_id=1
        open(unit=pos_id, file=(trim(filename)//'.pos'), status='old')
        ion_id=2
        open(unit=ion_id, file=(trim(filename)//'.ion'), status='unknown')
        fss_id=3
        open(unit=fss_id, file=(trim(filename)//'.fss'), status='unknown')
        trans_id=4
        open(unit=trans_id, file=(trim(filename)//'.trans'), status='unknown')
        err_id=5
        open(unit=err_id, file=(trim(filename)//'.error'), status='unknown')
        trans2_id=10
        open(unit=trans2_id, file=(trim(filename)//'.trans_no_r'), status='unknown')
        hb_id=11
        open(unit=hb_id, file=(trim(filename)//'.hbcase'), status='unknown')
        !
      case(2)
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(10)
        close(11)
      case default
    end select
    !
  end SUBROUTINE files
  !
  !
  !
  !
  !
  !
  !
END PROGRAM prepare

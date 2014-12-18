!
!  Refactored 20141118
!  Lixin Zheng, lixin.zheng@temple.edu
!----------------------------------------
!
MODULE main_variables
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER   :: DP = selected_real_kind(14,200) !double-precision kind
  !
  !======Parameters======
  REAL(DP), PARAMETER  :: pi=4.d0*atan(1.d0)
  REAL, PARAMETER      :: celldm = 23.5170     !dimension of the supercell 
  REAL, PARAMETER      :: convertBA= 0.529     !Transfer Bohr to Angestrom9
  REAL, PARAMETER      :: r_FSS=3.5            !distance (Angestrom) between O and O for first solvation
  REAL, PARAMETER      :: r_HB=3.5
  REAL, PARAMETER      :: angle_HB=4.d0*atan(1.d0)/6
  REAL, PARAMETER      :: r_cov=1.24           !covalent bond length of OH in a regular H2O 
  !======Global variables======
  CHARACTER(LEN=40)    :: filename
  INTEGER              :: nsp(2)               !number of each atomic species
  INTEGER              :: cs
  INTEGER              :: stepstart
  INTEGER              :: stepstop
  INTEGER              :: pos_id, wfc_id,  &    !ID for each files  
                          ion_id, fss_id, trans_id, trans2_id, hb_id, err_id, &
                          msd_id, rdf_id
  REAL,ALLOCATABLE     :: covalent(:,:)
  !======Variables concerning .pos======
  REAL,ALLOCATABLE     :: rO(:,:),rH(:,:)
  !======Variables concerning .fss======
  REAL                 :: rO_FSS(3,15)         !atomic positions (dimension,atomic-species)
  REAL                 :: rH_FSS(3,30)         !atomic positions (dimension,atomic-species)
  REAL                 :: rIO(3)               !Real ion position (dimension)
  REAL                 :: rIH(3,3)             !Real ion position (dimension,numH)
  INTEGER              :: new_nsp(2)           !number of each atomic species
  INTEGER              :: i_FSS(2:30)
  INTEGER              :: iIO                  !index of Ion oxygen
  INTEGER              :: iIH(3)               !index of Ion oxygen
  INTEGER              :: numION               !Number of Ions in each timestep. 
  INTEGER              :: cov_bond(4,15)       !The covalent bonded hydrogen to each oxygen
  !
  !
  namelist /global/ cs,filename,nsp
  namelist /control/ stepstart,stepstop
  !namelist /rdf/
  !
  !
  contains
  !
  !
  SUBROUTINE init
    !
    IMPLICIT NONE
    !read Namelist input for stdin
    read(*,global)
    read(*,control)
    !
    !allocation
    allocate(rO(3,nsp(1)))
    allocate(rH(3,nsp(2)))
    allocate(covalent(4,nsp(1)))
    !
  end SUBROUTINE init
  !
  !
  !
  !
END MODULE main_variables

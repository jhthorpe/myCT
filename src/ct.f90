PROGRAM ct
  USE myUtils

  IMPLICIT NONE

  !interface to subroutine calls
  INTERFACE
    SUBROUTINE input(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13)
      INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a9
      REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a4 
      INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: a1
      CHARACTER(LEN=2), DIMENSION(0:), INTENT(IN) :: a5
      REAL(KIND=4), DIMENSION(0:), INTENT(INOUT) :: a2 
      INTEGER, DIMENSION(0:), INTENT(INOUT) :: a8
      REAL(KIND=4), INTENT(INOUT) :: a3,a6,a10,a11,a12
      INTEGER, INTENT(INOUT) :: a7,a13
    END SUBROUTINE input

    SUBROUTINE levels_HO(a1,a2,a3,a4,a5,a6,a7)
      INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: a7
      CHARACTER(LEN=2), DIMENSION(0:), INTENT(IN) :: a5 
      INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: a1
      REAL(KIND=4), DIMENSION(0:,0:), INTENT(IN) :: a4
      REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: a2
      REAL(KIND=4), INTENT(IN) :: a3,a6
    END SUBROUTINE levels_HO

    SUBROUTINE stat_model(a1,a2,a3,a4,a5,a6,a7,a8,a9)
      INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: a5,a6
      REAL(KIND=4), DIMENSION(0:,0:), INTENT(IN) :: a7
      REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: a8
      INTEGER(KIND=4), INTENT(IN) :: a3
      REAL(KIND=4), INTENT(IN) :: a1,a2,a4,a9
    END SUBROUTINE stat_model
  END INTERFACE 
 
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ids,Ngues
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: Elvl
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: Wi
  CHARACTER(LEN=2), DIMENSION(0:3) :: names
  INTEGER(KIND=4), DIMENSION(0:3) :: nvib
  REAL(KIND=4), DIMENSION(0:3) :: Eelc
  INTEGER, DIMENSION(0:4) :: options
  INTEGER(KIND=4) :: mqm,nX
  REAL(KIND=4) :: Etol,Eint,t1,t2,T,E_rel,dE
  INTEGER :: i,j,k

  CALL CPU_TIME(t1)

  WRITE(*,*) "=============================================================="
  WRITE(*,*) "			myCT, v0.0			    "
  WRITE(*,*) ""
  WRITE(*,*) "		James H. Thorpe				    "
  WRITE(*,*) "		University of Florida			    "
  WRITE(*,*) "		https://github.com/jhthorpe/myCT	    "
  WRITE(*,*) 
  WRITE(*,*) "Enumerating potential overlaps for the reaction"
  WRITE(*,*) "A(+) + B -> A + B(+)"
  WRITE(*,*)
  WRITE(*,*) "Please do not include ZPE in your reported energies"
  WRITE(*,*) 

  names = ['A+', 'B ', 'A ', 'B+']

  CALL input(nvib,Eelc,Etol,Wi,names,Eint,mqm,options,Ngues,T,E_rel,dE,nX)
  WRITE(*,*) "=============================================================="

  IF (options(2) .NE. 0) CALL stat_model(E_rel,dE,nX,T,options,nvib,Wi,Eelc,Eint)
  IF (options(1) .NE. 0) THEN
    IF (options(0) .EQ. 0) THEN
      CALL levels_HO(nvib,Eelc,Etol,Wi,names,Eint,Ngues)
    END IF
  END IF

  DEALLOCATE(Wi)
  IF (ALLOCATED(Elvl)) DEALLOCATE(Elvl)
  IF (ALLOCATED(ids)) DEALLOCATE(ids)

  CALL CPU_TIME(t2)
  WRITE(*,*) "ct finished in (s)", (t2 - t1)

END PROGRAM

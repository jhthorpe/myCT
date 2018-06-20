PROGRAM ct
  USE myUtils

  IMPLICIT NONE

  !interface to subroutine calls
  INTERFACE
    SUBROUTINE input(a1,a2,a3,a4,a5,a6,a7,a8)
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a4 
      INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: a1
      CHARACTER(LEN=2), DIMENSION(0:), INTENT(IN) :: a5
      REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: a2 
      INTEGER, DIMENSION(0:), INTENT(INOUT) :: a8
      REAL(KIND=8), INTENT(INOUT) :: a3,a6
      INTEGER, INTENT(INOUT) :: a7
    END SUBROUTINE input
  END INTERFACE 
 
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ids
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: Elvl
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Wi
  CHARACTER(LEN=2), DIMENSION(0:3) :: names
  INTEGER(KIND=4), DIMENSION(0:3) :: nvib
  REAL(KIND=8), DIMENSION(0:3) :: Eelc
  INTEGER, DIMENSION(0:1) :: options
  REAL(KIND=8) :: Etol,Eint,t1,t2
  INTEGER :: mqm
  INTEGER :: i,j,k

  CALL CPU_TIME(t1)

  WRITE(*,*)
  WRITE(*,*) "Enumerating potential overlaps for the reaction"
  WRITE(*,*) "A(+) + B -> A + B(+)"
  WRITE(*,*)
  WRITE(*,*) "Please do not include ZPE in your reported energies"

  names = ['A+', 'B ', 'A ', 'B+']

  CALL input(nvib,Eelc,Etol,Wi,names,Eint,mqm,options)
  IF (options(0) .EQ. 0) THEN
!    CALL levels_HO(nvib,Eelc,Etol,Wi,names,Eint)
  END IF

  DEALLOCATE(Wi)
  DEALLOCATE(Elvl)
  DEALLOCATE(ids)

  CALL CPU_TIME(t2)
  WRITE(*,*) "ct finished in (s)", (t2 - t1)

END PROGRAM

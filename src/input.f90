!---------------------------------------------------------------------
!               Get input and print it 
!---------------------------------------------------------------------
    !Values
    ! Wi        :       2D real4, array of vibrational frequencies of
    ! molecules
    ! nvib      :       1D int, array of number of vibrational modes
    ! Eelc      :       1D real4, array of elc energies in cm-1
    ! names     :       1D chr*2, array of molecule names for
    ! convenience
    ! Etol      :       real4, tolerance energy in cm-1
    ! Eint      :       real4, intitial internal energy
    ! mqm       :       int, max quantum number
    ! options   :       1D int, list of options
    ! Ngues     :       2D int4, initial guess array [rcts/prds,levels]

  SUBROUTINE input(nvib,Eelc,Etol,Wi,names,Eint,mqm,options,Ngues)

    IMPLICIT NONE

    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Ngues
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Wi
    INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: nvib
    CHARACTER(LEN=2), DIMENSION(0:), INTENT(IN) :: names
    REAL(KIND=4), DIMENSION(0:), INTENT(INOUT) :: Eelc
    INTEGER, DIMENSION(0:), INTENT(INOUT) :: options
    CHARACTER(LEN=6), DIMENSION(0:3) :: fnames
    REAL(KIND=4), INTENT(INOUT) :: Etol,Eint
    INTEGER, INTENT(INOUT) :: mqm

    LOGICAL :: ex
    INTEGER :: i,j,k

    fnames = ['A+.dat', 'B.dat ', 'A.dat ', 'B+.dat']

    OPEN(unit=1,file='input.dat',access='sequential',status='old')
      READ(1,*) nvib(0)
      READ(1,*) nvib(1)
      READ(1,*) nvib(2)
      READ(1,*) nvib(3)
      READ(1,*)
      READ(1,*) Eelc(0)
      READ(1,*) Eelc(1)
      READ(1,*) Eelc(2)
      READ(1,*) Eelc(3)
      READ(1,*)
      READ(1,*) Etol
      READ(1,*) Eint
      READ(1,*)
      READ(1,*) options(0)
      READ(1,*) mqm
    CLOSE(unit=1)

    !setup arrays for fast indexing
    ALLOCATE(Wi(0:3,0:MAXVAL(nvib)-1))

    !allocate Ngues and get it if exists
    INQUIRE(file="Ngues", EXIST=ex)
    IF (ex) THEN
      WRITE(*,*) "Reading inital guess from Ngues"
      ALLOCATE(Ngues(0:1,0:2*MAXVAL(nvib)-1))
      OPEN(unit=1,file="Ngues",status="old",access="sequential")
      READ(1,*) Ngues(0,0:nvib(0)-1)
      READ(1,*) Ngues(0,nvib(0):nvib(0)+nvib(1)-1)
      READ(1,*) Ngues(1,0:nvib(2)-1)
      READ(1,*) Ngues(1,nvib(2):nvib(2)+nvib(3)-1)
      CLOSE(unit=1)
    END IF

    !get harmonic frequencies
    IF (options(0) .EQ. 0) THEN
      DO i=0,3
        IF (nvib(i) .NE. 0) THEN
          OPEN(unit=1,file=fnames(i),status='old',access='sequential')
            DO j=0,nvib(i)-1
              READ(1,*) Wi(i,j)
            END DO
          CLOSE(unit=1)
        END IF
      END DO
    ELSE
      WRITE(*,*) "Sorry, I can only do harmonic frequencies for now"
      WRITE(*,*) "VPT2 will be implemented soon"
      STOP
    END IF


    !print output
    WRITE(*,*)
    WRITE(*,*) "Initial Internal Energy (cm-1) : ", Eint
    WRITE(*,*) "Energy Tolerance (cm-1) :", Etol
    WRITE(*,*)
    DO i=0,3
      WRITE(*,*) names(i)
      WRITE(*,*) "Elc: ", Eelc(i)
      DO j=0,nvib(i)-1
        WRITE(*,*) Wi(i,j)
      END DO
      WRITE(*,*)
    END DO
    IF (ALLOCATED(Ngues)) THEN
      WRITE(*,*) "Initial Guess vectors"
      WRITE(*,*) "A+ levels"
      WRITE(*,*) Ngues(0,0:nvib(0)-1)
      WRITE(*,*) "B levels"
      WRITE(*,*) Ngues(0,nvib(0):nvib(0)+nvib(1)-1)
      WRITE(*,*) "A levels"
      WRITE(*,*) Ngues(1,0:nvib(2)-1)
      WRITE(*,*) "B+ levels"
      WRITE(*,*) Ngues(1,nvib(2):nvib(2)+nvib(3)-1)
    END IF

  END SUBROUTINE input


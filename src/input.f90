!---------------------------------------------------------------------
!               Get input and print it 
!---------------------------------------------------------------------
    !Values
    ! Wi        :       2D dp, array of vibrational frequencies of
    ! molecules
    ! nvib      :       1D int, array of number of vibrational modes
    ! Eelc      :       1D dp, array of elc energies in cm-1
    ! names     :       1D chr*2, array of molecule names for
    ! convenience
    ! Etol      :       dp, tolerance energy in cm-1
    ! Eint      :       dp, intitial internal energy
    ! mqm       :       int, max quantum number
    ! options   :       1D int, list of options

  SUBROUTINE input(nvib,Eelc,Etol,Wi,names,Eint,mqm,options)

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Wi
    INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: nvib
    CHARACTER(LEN=2), DIMENSION(0:), INTENT(IN) :: names
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Eelc
    INTEGER, DIMENSION(0:), INTENT(INOUT) :: options
    CHARACTER(LEN=6), DIMENSION(0:3) :: fnames
    REAL(KIND=8), INTENT(INOUT) :: Etol,Eint
    INTEGER, INTENT(INOUT) :: mqm

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

  END SUBROUTINE input


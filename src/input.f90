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
    ! mqn       :       int, max quantum number
    ! options   :       1D int, list of options
    ! Ngues     :       2D int4, initial guess array [rcts/prds,levels]
    ! T		:	real4, temperature in K

  SUBROUTINE input(nvib,Eelc,Etol,Wi,names,Eint,mqn,options,Ngues,T,Emax,dE,nX)
    IMPLICIT NONE

    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Ngues
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Wi
    INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: nvib
    CHARACTER(LEN=2), DIMENSION(0:), INTENT(IN) :: names
    REAL(KIND=4), DIMENSION(0:), INTENT(INOUT) :: Eelc
    INTEGER, DIMENSION(0:), INTENT(INOUT) :: options
    CHARACTER(LEN=6), DIMENSION(0:3) :: fnames
    INTEGER(KIND=4), INTENT(INOUT) :: mqn,nX
    REAL(KIND=4), INTENT(INOUT) :: Etol,Eint,T,Emax,dE

    LOGICAL :: ex
    INTEGER :: i,j,k

    fnames = ['A+.dat', 'B.dat ', 'A.dat ', 'B+.dat']

    OPEN(unit=1,file='input.dat',access='sequential',status='old')
      READ(1,*)
      READ(1,*) nvib(0)
      READ(1,*) nvib(1)
      READ(1,*) nvib(2)
      READ(1,*) nvib(3)
      READ(1,*)
      READ(1,*) options(0)
      READ(1,*)
      READ(1,*)
      READ(1,*) options(1)
      READ(1,*)
      READ(1,*) Eelc(0)
      READ(1,*) Eelc(1)
      READ(1,*) Eelc(2)
      READ(1,*) Eelc(3)
      READ(1,*)
      READ(1,*) Etol
      READ(1,*) Eint
      READ(1,*)
      READ(1,*) mqn
      READ(1,*) 
      READ(1,*) 
      READ(1,*) options(2)
      READ(1,*) T
      READ(1,*) Emax
      READ(1,*) dE
      READ(1,*) nX
      READ(1,*) options(3)
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

    !get frequencies
    IF (options(0) .EQ. 0) THEN !harmonic only
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

    !quantum print 
    WRITE(*,*) "---------- Quantum Treatment ----------"
    IF (options(1) .NE. 0) THEN
      
      WRITE(*,*)
      WRITE(*,*) "Initial Internal Energy (cm-1) : ", Eint
      WRITE(*,*) "Energy Tolerance (cm-1) :", Etol
      WRITE(*,*)
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
    ELSE
      WRITE(*,*) "No quantum treatment requested"
    END IF
    WRITE(*,*) 

    !statistical print
    WRITE(*,*) "--------- Statistical Treatment ---------" 
    WRITE(*,*) 
    IF (options(2) .NE. 0) THEN
      IF (options(2) .EQ. 1) THEN 
        WRITE(*,*) "Statistical treatment	: Average Energy" 
        WRITE(*,*) "Temperature (K)	:", T
      END IF
      IF (options(2) .EQ. 2) THEN
        WRITE(*,*) "Statistical treatment	: Energy Distribution"
        WRITE(*,*) "Temperature (K)	:", T
        WRITE(*,*) "Max energy (cm-1)	:", Emax
        WRITE(*,*) "Energy binsize (cm-1)	:", dE
        WRITE(*,*) "Fraction bins		:", nX 
        IF (options(3) .EQ. 0) WRITE(*,*) "Density treatment	: External"
        IF (options(3) .EQ. 1) WRITE(*,*) "Density treatment	: HO"
      END IF
    ELSE
      WRITE(*,*) "No statistical treatment requested"
    END IF

    !energy level print
    WRITE(*,*) 
    WRITE(*,*) "------------- Energy Levels -------------" 
    WRITE(*,*) 
    IF (options(0) .EQ.0 ) WRITE(*,*) "Vibrational Frequencies : Harmonic"  
    WRITE(*,*) 
    DO i=0,3
      WRITE(*,*) names(i)
      WRITE(*,*) "Electronic Energy (cm-1) :", Eelc(i)
      WRITE(*,*) 
      IF (options(0) .EQ. 0 .OR. options(0) .EQ. 1) THEN
        WRITE(*,*) "Harmonic Frequencies"   
        DO j=0,nvib(i)-1
          WRITE(*,*) j, Wi(i,j)
        END DO
      END IF
      WRITE(*,*) 
    END DO

  END SUBROUTINE input


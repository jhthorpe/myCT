!---------------------------------------------------------------------
!       !Find all levels that have overlapping energys 
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
    ! numcombo  :       int*8, number of valid combinations to consider
    ! Etot      :       real4, total energy of A(+) + B in cm-1
    ! Esys      :       real4, total energy of current combination of A +
    ! B(+)
    ! tol       :       real4, tolerance for equivlence of two reals
    ! a,b       :       int, number considered, number accepted
    ! mqn       :       int, max quantum number
    ! Ngues     :       1D int4, guess array of quantum numbers
    ! l1,l2     :       int4, l1 - length of hash table, l2- length of key
    ! A         :       1D bool, hash table truth table
    ! B         :       2D int4, hash table keys
    ! C         :       1D bool, hash table vals

  SUBROUTINE levels_HO(nvib,Eelc,Etol,Wi,names,Eint,Ngues)
    USE myUtils

    IMPLICIT NONE

    INTERFACE
      SUBROUTINE init_HO(Wi,m,Et,Es,N0,Ez)
        INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: N0
        REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: Wi 
        INTEGER(KIND=4), INTENT(IN) :: m
        REAL(KIND=4), INTENT(IN) :: Et,Es,Ez
      END SUBROUTINE init_HO
      
      SUBROUTINE recall_HO(N,Wi,ids,A,B,C,Eu,El,Esys,idx,nall,ngood,l1,l2)
         INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B
         LOGICAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: A,C
         INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: N
         INTEGER(KIND=4),DIMENSION(0:), INTENT(IN) :: ids
         REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: Wi
         INTEGER(KIND=4), INTENT(INOUT) :: nall,ngood,l1
         INTEGER(KIND=4), INTENT(IN) :: l2,idx
         REAL(KIND=4), INTENT(IN) :: Eu,El,Esys
      END SUBROUTINE recall_HO

      REAL(KIND=4) FUNCTION energy_HO(N,W,m)
        INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: N
        REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: W
        INTEGER(KIND=4), INTENT(IN) ::m
      END FUNCTION energy_HO

    END INTERFACE

    REAL(KIND=4), PARAMETER :: tol=1.0E-8

    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Ngues
    CHARACTER(LEN=2), DIMENSION(0:), INTENT(IN) :: names
    INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: nvib
    REAL(KIND=4), DIMENSION(0:,0:), INTENT(IN) :: Wi
    REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: Eelc
    REAL(KIND=4), INTENT(IN) :: Etol,Eint

    INTEGER(KIND=4), DIMENSION(0:nvib(0)+nvib(1)-1) :: N0
    INTEGER(KIND=4), DIMENSION(0:nvib(2)+nvib(3)-1) :: N1
    REAL(KIND=4), DIMENSION(0:nvib(0)+nvib(1)-1) :: W0
    REAL(KIND=4), DIMENSION(0:nvib(2)+nvib(3)-1) :: W1
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ids,B
    LOGICAL, DIMENSION(:), ALLOCATABLE :: A,C
    INTEGER(KIND=4) :: l1,l2
    INTEGER(KIND=4) :: nall,ngood
    REAL(KIND=4) :: Etot, Esys,Ezpe,t1,t2,t3
    INTEGER :: i,j,k

    !-------------------------
    ! Let's start with  A(+) + B

    CALL CPU_TIME(t1)

    WRITE(*,*) 

    nall = 0
    ngood = 0

    !assign IDs and move frequencies
    ALLOCATE(ids(0:1,0:2*MAXVAL(nvib)-1))
    W0(0:nvib(0)-1) = Wi(0,0:nvib(0)-1) 
    W0(nvib(0):nvib(0)+nvib(1)-1) = Wi(1,0:nvib(1)-1)
    ids(0,0:nvib(0)-1) = (/ (0, i=0, nvib(0)-1) /)
    ids(0,nvib(0):nvib(0)+nvib(1)-1) = (/ (1, i=0, nvib(1)-1) /)

    !Find starting point for the calculation of A+ & B levels
    WRITE(*,*) "-----------------------------------------"
    WRITE(*,*) "Starting search for levels of A+ and B"

    Esys = Eelc(0) + Eelc(1)

    ! If a user guess was given...
    N0 = 0
    Ezpe = energy_HO(N0,W0,nvib(0)+nvib(1))
    IF (ALLOCATED(Ngues) ) THEN
      WRITE(*,*) "Using Ngues as intial guess vectors"
      N0(0:nvib(0)+nvib(1)-1) = Ngues(0,0:nvib(0)+nvib(1)-1)
    ELSE   
      CALL init_HO(W0,nvib(0)+nvib(1),Eint,Esys,N0,Ezpe)
    END IF

    !total energy of reactants
    Etot = Esys + Ezpe + Eint

    !setup hash tables
    l1 = 10*(nvib(0)+nvib(1))
    l2 = nvib(0)+nvib(1)
    CALL hash_qinit_1Dint4_bool(A,B,C,l1,l2)

    !open unformatted binary file for writing, and do search
    WRITE(*,*) 
    WRITE(*,*) "-----------------------------------------"
    WRITE(*,*) "Starting recursive search" 
    OPEN(unit=2,file="A+B.bin",status='replace',access='sequential',form='unformatted')
    CALL recall_HO(N0,W0,ids(0,:),A,B,C,Esys+Ezpe+Eint+Etol,Esys+Ezpe+Eint-Etol,Esys,0,nall,ngood,l1,l2)
    CLOSE(unit=2)

    !write results to non binary file
    OPEN(unit=3,file="A+Blevels",status='replace',access='sequential')
    OPEN(unit=2,file="A+B.bin",status='old',access='sequential',form='unformatted')
    WRITE(3,*)  "A+ #vibs :", nvib(0), "B #vibs :", nvib(1) 
    DO i=0,ngood-1
      READ(2) N0
      WRITE(3,*) N0, Esys+energy_HO(N0,W0,nvib(0)+nvib(1))
    END DO
    CLOSE(unit=3)
    CLOSE(unit=2)

    DEALLOCATE(A)
    DEALLOCATE(B)
    DEALLOCATE(C)
    CALL CPU_TIME(t2)

    WRITE(*,*)
    WRITE(*,*) "A(+) + B finished in (s)", t2-t1
    WRITE(*,*) "At a rate of (levels/second)", nall/(t2-t1)
    WRITE(*,*) "Total levels considered :", nall
    WRITE(*,*) "Total levels in range   :", ngood
    WRITE(*,*) "-----------------------------------------"

    !-------------------------
    ! Now for A + B(+)
    nall = 0
    ngood = 0

    !assign IDs and move frequencies
    ids(1,0:nvib(2)-1) = (/ (2, i=0, nvib(2)-1) /)
    ids(1,nvib(2):nvib(2)+nvib(3)-1) = (/ (3, i=0, nvib(3)-1) /)
    W1(0:nvib(0)-1) = Wi(2,0:nvib(2)-1) 
    W1(nvib(2):nvib(2)+nvib(3)-1) = Wi(3,0:nvib(3)-1)

    !Find starting point for the calculation of A+ & B levels
    WRITE(*,*) "-----------------------------------------"
    WRITE(*,*) "Starting search for levels of A and B+"

    Esys = Eelc(2) + Eelc(3)

    ! If a user guess was given...
    N1 = 0
    Ezpe = energy_HO(N1,W1,nvib(2)+nvib(3))
    IF (ALLOCATED(Ngues) ) THEN
      WRITE(*,*) "Using Ngues as intial guess vectors"
      N1(0:nvib(2)+nvib(3)-1) = Ngues(1,0:nvib(2)+nvib(3)-1)
    ELSE   
      CALL init_HO(W1,nvib(2)+nvib(3),Etot-Esys-Ezpe,Esys,N1,Ezpe)
    END IF

    !setup hash tables
    l1 = 10*(nvib(2)+nvib(3))
    l2 = nvib(2)+nvib(3)
    CALL hash_qinit_1Dint4_bool(A,B,C,l1,l2)

    !open unformatted binary file for writing, and do search
    WRITE(*,*) 
    WRITE(*,*) "-----------------------------------------"
    WRITE(*,*) "Starting recursive search" 
    OPEN(unit=2,file="AB+.bin",status='replace',access='sequential',form='unformatted')
    CALL recall_HO(N1,W1,ids(1,:),A,B,C,Etot+Etol,Etot-Etol,Esys,0,nall,ngood,l1,l2)
    CLOSE(unit=2)

    !write results to non binary file
    OPEN(unit=3,file="AB+levels",status='replace',access='sequential')
    OPEN(unit=2,file="AB+.bin",status='old',access='sequential',form='unformatted')
    WRITE(3,*)  "A #vibs :", nvib(2), "B+ #vibs :", nvib(3) 
    DO i=0,ngood-1
      READ(2) N1
      WRITE(3,*) N1, Esys+energy_HO(N1,W1,nvib(2)+nvib(3))
    END DO
    CLOSE(unit=3)
    CLOSE(unit=2)

    DEALLOCATE(A)
    DEALLOCATE(B)
    DEALLOCATE(C)
    CALL CPU_TIME(t3)

    WRITE(*,*)
    WRITE(*,*) "A + B(+) finished in (s)", t3-t2
    WRITE(*,*) "At a rate of (levels/second)", nall/(t3-t2)
    WRITE(*,*) "Total levels considered :", nall
    WRITE(*,*) "Total levels in range   :", ngood
    WRITE(*,*) "-----------------------------------------"


    DEALLOCATE(ids)

  END SUBROUTINE levels_HO


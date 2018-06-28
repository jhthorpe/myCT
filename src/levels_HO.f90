!---------------------------------------------------------------------
!       !Find all levels that have overlapping energys 
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
    ! numcombo  :       int*8, number of valid combinations to consider
    ! Etot      :       dp, total energy of A(+) + B in cm-1
    ! Esys      :       dp, total energy of current combination of A +
    ! B(+)
    ! tol       :       dp, tolerance for equivlence of two reals
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
      SUBROUTINE init_HO(Wi,nvib,Et,N0)
        INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: N0
        REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Wi 
        INTEGER(KIND=4), INTENT(IN) :: nvib
        REAL(KIND=8), INTENT(IN) :: Et
      END SUBROUTINE init_HO
      
      SUBROUTINE enumerate_HO(N,Wi,ids,A,B,C,Eu,El,Esys,idx,nall,ngood,l1,l2)
         INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B
         LOGICAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: A,C
         INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: N
         INTEGER(KIND=4),DIMENSION(0:), INTENT(IN) :: ids
         REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Wi
         INTEGER(KIND=4), INTENT(INOUT) :: nall,ngood,l1
         INTEGER(KIND=4), INTENT(IN) :: l2,idx
         REAL(KIND=8), INTENT(IN) :: Eu,El,Esys
      END SUBROUTINE enumerate_HO
    END INTERFACE

    REAL(KIND=8), PARAMETER :: tol=1.0D-16

    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) :: Ngues
    CHARACTER(LEN=2), DIMENSION(0:), INTENT(IN) :: names
    INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: nvib
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Wi
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Eelc
    REAL(KIND=8), INTENT(IN) :: Etol,Eint

    INTEGER(KIND=4), DIMENSION(0:nvib(0)+nvib(1)-1) :: N0
    INTEGER(KIND=4), DIMENSION(0:nvib(2)+nvib(3)-1) :: N1
    REAL(KIND=8), DIMENSION(0:nvib(0)+nvib(1)-1) :: W0
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ids,B
    LOGICAL, DIMENSION(:), ALLOCATABLE :: A,C
    INTEGER(KIND=4) :: l1,l2
    INTEGER(KIND=4) :: nall,ngood
    REAL(KIND=8) :: Etot, Esys,t1,t2
    INTEGER :: i,j,k

    CALL CPU_TIME(t1)

    WRITE(*,*) 

    nall = 0
    ngood = 0

    !assign IDs and move frequencies
    W0(0:nvib(0)-1) = Wi(0,0:nvib(0)-1) 
    W0(nvib(0):nvib(0)+nvib(1)-1) = Wi(1,0:nvib(1)-1)

    ALLOCATE(ids(0:1,0:2*MAXVAL(nvib)-1))
    ids(0,0:nvib(0)-1) = (/ (0, i=0, nvib(0)-1) /)
    ids(0,nvib(0):nvib(0)+nvib(1)-1) = (/ (1, i=0, nvib(1)-1) /)
    ids(1,0:nvib(2)-1) = (/ (2, i=0, nvib(2)-1) /)
    ids(1,nvib(2):nvib(2)+nvib(3)-1) = (/ (3, i=0, nvib(3)-1) /)

    !Find starting point for the calculation of A+ & B levels
    WRITE(*,*) "-----------------------------------------"
    WRITE(*,*) "Starting search for levels of A+ and B"

    Esys = Eelc(0) + Eelc(1)

    ! If a user guess was given...
    IF (ALLOCATED(Ngues) ) THEN
      WRITE(*,*) "Using Ngues as intial guess vectors"
      N0(0:nvib(0)+nvib(1)-1) = Ngues(0,0:nvib(0)+nvib(1)-1)
    ELSE   
      CALL init_HO(W0,nvib(0)+nvib(1),Eelc(0)+Eelc(1)+Etol-Eint,N0)
    END IF

    !setup hash tables
    l1 = 1000
    l2 = nvib(0)+nvib(1)
    CALL hash_qinit_1Dint4_bool(A,B,C,l1,l2)

    !open unformatted binary file for writing
    OPEN(unit=2,file='A+Blevels',status='replace',access='sequential',form='unformatted')
    CALL enumerate_HO(N0,W0,ids(0,:),A,B,C,Esys+Eint+Etol,Etot+Eint-Etol,Esys,0,nall,ngood,l1,l2)
    CLOSE(unit=2)

    CALL CPU_TIME(t2)

    WRITE(*,*)
    WRITE(*,*) "A(+) + B finished in (s)", t2-t1
    WRITE(*,*) "Total levels considered :", nall
    WRITE(*,*) "Total levels in range   :", ngood
    WRITE(*,*) "-----------------------------------------"

    DEALLOCATE(ids)

  END SUBROUTINE levels_HO


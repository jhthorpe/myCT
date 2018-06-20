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
    ! Y0        :       2D bool, tracks if we've found point yet or not

  SUBROUTINE levels_HO(nvib,Eelc,Etol,Wi,names,Eint)

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: tol=1.0D-16

    CHARACTER(LEN=2), DIMENSION(0:), INTENT(IN) :: names
    INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: nvib
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Wi
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Eelc
    REAL(KIND=8), INTENT(IN) :: Etol,Eint

    INTEGER(KIND=4), DIMENSION(0:nvib(0)+nvib(1)-1) :: N0,N1
    INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ids
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: Y0
    INTEGER(KIND=4) :: a,b
    REAL(KIND=8) :: Etot, Esys
    INTEGER :: i,j,k

    a = 0
    b = 0

    !assign IDs
    ALLOCATE(ids(0:1,0:2*MAXVAL(nvib)-1))
    ids(0,0:nvib(0)-1) = (/ (0, i=0, nvib(0)-1) /)
    ids(0,nvib(0):nvib(0)+nvib(1)-1) = (/ (1, i=0, nvib(1)-1) /)
    ids(1,0:nvib(2)-1) = (/ (2, i=0, nvib(2)-1) /)
    ids(1,nvib(2):nvib(2)+nvib(3)-1) = (/ (3, i=0, nvib(3)-1) /)

    !Find combinations of A(+) & B with the input internal energy and
    !around tol
    N0(0:nvib(0)+nvib(1)-1) = (/ (0, i=0, nvib(0)+nvib(1)-1) /)
    !Etot = Eelc(0) + Eelc(1) + SUM(Elvl(0,0:nvib(0)+nvib(1)-1,0))

    !Find starting point for the calculation
    !CALL start_HO(Wi,nvib,)
    !k = MAXLOC(Elvl(0,0,:),3, MASK = (Elvl(0,0,:) - Etot+Eint+Etol)
    !.LT. &
    !      tol)-1 

    STOP

    !setup initial arrays
    ALLOCATE(Y0(k:k,0:nvib(0)+nvib(1)-1))
    Y0(k,0:nvib(0)+nvib(1)-1) = (/ (.FALSE., i=0,nvib(0)+nvib(1)-1) /)
    N0(0) = k

    OPEN(unit=2,file='A+Blevels',status='replace',access='sequential',form='unformatted')
    CALL enumerate_HO(N0,Wi(0,:),ids(0,:),Y0,Etot+Eint+Etol,Etot+Eint-Etol,Eelc(0)+Eelc(1),0,a,b,nvib(0)+nvib(1))
    CLOSE(unit=2)

    WRITE(*,*)
    WRITE(*,*) "A(+) + B finished..."
    WRITE(*,*) "Total levels considered :", a
    WRITE(*,*) "Total levels in range   :", b

    DEALLOCATE(ids)

  END SUBROUTINE levels_HO


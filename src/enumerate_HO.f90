!---------------------------------------------------------------------
!       !recursively aquire all the possible energy values
!---------------------------------------------------------------------
    !Values
    ! N         :       1D int, list of quantum numbers
    ! Wi        :       1D dp, list of HO fundementals 
    ! X         :       1D int, list of ids
    ! Y         :       2D bool, truth table for if we've considered a
    ! point
    ! Eu        :       dp, upper energy limit
    ! El        :       dp, lower energy limit
    ! Ep        :       dp, energy passed down to us
    ! idx       :       int, which index in W we are working on now
    ! a         :       int, number of values considered
    ! b         :       int, number of energies accepted
    ! m         :       int, number of vibrational modes 

  RECURSIVE SUBROUTINE enumerate_HO(N,Wi,X,Y,Eu,El,Ep,idx,a,b,m)
    IMPLICIT NONE
    EXTERNAL energy_HO
  
    
    REAL(KIND=8), PARAMETER :: tol=1.0D-16

    LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Y
    INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: N
    INTEGER(KIND=4),DIMENSION(0:), INTENT(IN) :: X
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Wi
    INTEGER(KIND=4), INTENT(INOUT) :: a,b
    INTEGER(KIND=4), INTENT(IN) :: m,idx
    REAL(KIND=8), INTENT(IN) :: Eu,El,Ep

    INTEGER(KIND=4), DIMENSION(0:m-1) :: newN
    INTEGER(KIND=4) :: nmax,ulim,llim
    INTEGER(KIND=4) :: i,j,k
    REAL(KIND=8) :: Er,energy_HO

    ulim = UBOUND(Y,2)
    llim = LBOUND(Y,2)
    Er = energy_HO(N,Wi,m)

    !base cases


    !1) we are below zero anywhere
    IF (MINVAL(N) .LT. 0) THEN
      WRITE(*,*) "lt 0 at N=",N
      RETURN

    !2) we have seen this index before 

    !1) we are out of bounds of Y0


    !3) we are outside the upper energy bounds


    !4) we are outside the lower energy bounds


    ELSE

    !go through all index
      DO i=0,m-1
        newN(:) = N(:)
        newN(i) = N(i)+1
        CALL enumerate_HO(newN,Wi,X,Y,Eu,El,Ep,i,a,b,m)
        newN(i) = N(i)-1
        CALL enumerate_HO(newN,Wi,X,Y,Eu,El,Ep,i,a,b,m)
      END DO


    END IF


    RETURN

  END SUBROUTINE enumerate_HO

!---------------------------------------------------------------------
!       !Find a good initial guess for the quantum number vector
!---------------------------------------------------------------------
    !Values
    ! Wi        :       2D dp, array of vibrational frequencies of
    ! nvib      :       int4, number of vibrational modes
    ! Et        :       dp, total energy of A(+) + B in cm-1
    ! Ngues     :       1D int4, guess array of quantum numbers
    ! N0        :       1D int4, output array of quantum numbers 

SUBROUTINE init_HO(Wi,nvib,Et,N0)
  IMPLICIT NONE

  REAL(KIND=8), PARAMETER :: tol=1.0D-16

  INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: N0
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Wi
  INTEGER(KIND=4), INTENT(IN) :: nvib
  REAL(KIND=8), INTENT(IN) :: Et

  WRITE(*,*) 
  WRITE(*,*) "Starting Recursive search for an initial solution guess"

  WRITE(*,*) "sorry, this is not implimented yet, please suply guess vector"
  STOP

END SUBROUTINE init_HO

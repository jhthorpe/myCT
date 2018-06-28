!---------------------------------------------------------------------
!       !Find a good initial guess for the quantum number vector
!---------------------------------------------------------------------
  ! As this is a harmonic oscillator, we know exactly how the energy
  !   levels will change from one to another. So we start with the 
  !   smallest vibration, and get to within +/- 1 number of the 
  !   range of answers, then let the recursive program do the rest 
  !   as we need to explore all that solution space anyways

  !Values
  ! W           :       2D dp, array of vibrational frequencies of
  ! m	        :       int4, number of vibrational modes
  ! Ei          :       dp,  internal energy of molecule (above ZPE)
  ! Ngues       :       1D int4, guess array of quantum numbers
  ! N0          :       1D int4, output array of quantum numbers 
  ! Ez		:	real8, zero point energy		

SUBROUTINE init_HO(W,m,Ei,Es,N0,Ez)
  IMPLICIT NONE

  INTERFACE
    REAL(KIND=8) FUNCTION energy_HO(N,W,m)
      INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: N
      REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: W
      INTEGER(KIND=4), INTENT(IN) ::m
    END FUNCTION energy_HO
  END INTERFACE

  REAL(KIND=8), PARAMETER :: tol=1.0D-16

  INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: N0
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: W
  REAL(KIND=8), INTENT(INOUT) :: Ez
  INTEGER(KIND=4), INTENT(IN) :: m
  REAL(KIND=8), INTENT(IN) :: Ei,Es

  INTEGER(KIND=4) :: i,j,n
  REAL(KIND=8) :: Ev,Et

  WRITE(*,*) 
  WRITE(*,*) "Generating initial guess for HO solutions" 

  !First, find the zero point energy
  N0 = 0
  Ez = energy_HO(N0,W,m) 
  Et = Es + Ez + Ei

  !Second, find the smallest vib mode
  i = MINLOC(W,1)-1

  !Remove the zero point energy of the smallest mode
  Ev = Ez - W(i)*0.5D0

  WRITE(*,*) "Target energy is :", Et
  WRITE(*,*) "Internal energy is : ", Ei
  WRITE(*,*) "Single Point energy is:", Es
  WRITE(*,*) "ZPE is :", Ez

  !Calculate the remaining energy, and which QN of 
  !  the smallest mode gets us there
  n = FLOOR((Et-Es-Ev)/W(i) - 0.5D0)
  IF (n .LT. 0) THEN !catch negative case
    n = 0
  END IF
  N0(i) = n

  WRITE(*,*)
  WRITE(*,*) "Initial Guess vector is..."
  WRITE(*,*) N0
  WRITE(*,*) "energy is:", Es+energy_HO(N0,W,m)


END SUBROUTINE init_HO

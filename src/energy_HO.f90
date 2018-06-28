!---------------------------------------------------------------------
!       get energy of combination for harmonic oscillator
!---------------------------------------------------------------------
    !values
    ! N         :       1D list of quanutm numbers
    ! W         :       1D list of fundemental frequencies
    ! m         :       int, number of modes

    !vector vector operations?

REAL(KIND=8) FUNCTION energy_HO(N,W,m)

  IMPLICIT NONE

  INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: N
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: W
  INTEGER(KIND=4), INTENT(IN) ::m
  REAL(KIND=8) :: temp
  INTEGER :: i

  temp = 0.0D0
  DO i=0,m-1
    temp = temp + (N(i)+0.5D0)*W(i)
  END DO

  energy_HO = temp

 END FUNCTION energy_HO



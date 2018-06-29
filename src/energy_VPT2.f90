!---------------------------------------------------------------------
!       get energy of combination for VPT2
!---------------------------------------------------------------------
    !values
    ! N         :       1D list of quanutm numbers
    ! G         :       real4, G0 term
    ! W         :       1D list of Wi frequencies 
    ! X         :       2D list of Xij coupling constants
    ! m         :       int, number of modes

    ! Can I use matrix vector multiplication here?
    ! check terms

  REAL(KIND=4) FUNCTION energy_VPT2(N,G,W,X,m)

    IMPLICIT NONE

    INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: N
    REAL(KIND=4), DIMENSION(0:,0:), INTENT(IN) :: X
    REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: W
    INTEGER(KIND=4), INTENT(IN) ::m,G
    REAL(KIND=4) :: temp
    INTEGER :: i,j

    temp = G

    DO j=0,m-1
      DO i=0,j
        temp = temp + 2*(N(j)+0.5E0)*(N(i)+0.5E0)*X(i,j)
      END DO
      temp = temp + (N(j) + 0.5E0)*W(j)
    END DO

    energy_VPT2 = temp

  END FUNCTION energy_VPT2


!---------------------------------------------------------------------
!      	stat_model 
!			James H. Thorpe
!			Sep 10, 2018
!	-  Uses general statistical theory to calculate ... 
!		... average energ of products, ...
!		...  or energy distribution of products
!---------------------------------------------------------------------
  !Values
  ! Emax	:	int4, max energy to calculate (cm-1)
  ! dE		:	int4, energy bin stepsize (cm-1)	
  ! nX		:	int4, number of fraction bins
  ! T		:	real4, temperature (K)
  ! options	:	1d int4, options array
  ! Wi		:	2d real4, harmonic frequencies
  ! dens	:	2D real4, density of states vectors

SUBROUTINE stat_model(Emax,dE,nX,T,options,nvib,Wi)
  IMPLICIT NONE

  !Inout
  INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: options,nvib
  REAL(KIND=4), DIMENSION(0:,0:), INTENT(IN) :: Wi
  INTEGER(KIND=4), INTENT(IN) :: nX
  REAL(KIND=8), INTENT(IN) :: Emax,dE,T

  !Internal
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dens
  INTEGER(KIND=4) :: Ebins

  !average energy
    

  !energy distribution


END SUBROUTINE stat_model 

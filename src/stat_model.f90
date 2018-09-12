!---------------------------------------------------------------------
!      	stat_model 
!			James H. Thorpe
!			T. Lam Nguyen
!			Sep 10, 2018
!	-  Uses general statistical theory to calculate ... 
!		... average energ of products, ...
!		...  or energy distribution of products
!---------------------------------------------------------------------
  !Values
  ! E_rel	:	int4, max energy relative to rcts to calculate (cm-1)
  ! dE		:	int4, energy bin stepsize (cm-1)	
  ! nX		:	int4, number of fraction bins
  ! T		:	real4, temperature (K)
  ! options	:	1d int4, options array
  ! Wi		:	2d real4, harmonic frequencies
  ! dens	:	2D real4, density of states vectors
  ! E_elc	:	1D real4, electronic energy
  ! E_int	:	real4, internal energy of products (cm-1)
  ! E_tot	:	real4, total energy of system
  ! Eint	:	real4, internal energy of reactants

SUBROUTINE stat_model(E_rel,dE,nX,T,options,nvib,Wi,E_elc,Eint)
  IMPLICIT NONE

  INTERFACE
    REAL(KIND=4) FUNCTION energy_HO(N,W,m)
      INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: N
      REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: W
      INTEGER(KIND=4), INTENT(IN) ::m
    END FUNCTION energy_HO
  END INTERFACE
  
  !Inout
  INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: options,nvib
  REAL(KIND=4), DIMENSION(0:,0:), INTENT(IN) :: Wi
  REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: E_elc
  INTEGER(KIND=4), INTENT(IN) :: nX
  REAL(KIND=4), INTENT(IN) :: E_rel,dE,T,Eint

  !Internal
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: N
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: dens
  INTEGER(KIND=4) :: max_idx,dummy
  REAL(KIND=4) :: t1,t2,A_conv,B_conv,A_enr,B_enr,E_int,E_tot,dum
  INTEGER :: i,j,k

  CALL CPU_TIME(t1)

  WRITE(*,*) "Statistical Energy Average"
  WRITE(*,*) 

  IF (options(3) .NE. 0) THEN
    WRITE(*,*) "Sorry, can only use external density of states currently"
    RETURN 
  END IF

  !get internal energies
  E_tot = E_elc(0) + E_elc(1)

  !reactants
  IF (options(0) .EQ. 0) THEN
    ALLOCATE(N(0:nvib(0)-1)) 
    N=0
    E_tot = E_tot + energy_HO(N,Wi(0,:),nvib(0))
    DEALLOCATE(N)
    ALLOCATE(N(0:nvib(1)-1))
    N=0
    E_tot = E_tot + energy_HO(N,Wi(1,:),nvib(1))
    DEALLOCATE(N)
  ELSE
    WRITE(*,*) "Cannot do VPT2 yet, exit stat_model"
    RETURN    
  END IF
  E_tot = E_tot + Eint

  !Products
  E_int = E_tot
  IF (options(0) .EQ. 0) THEN
    ALLOCATE(N(0:nvib(2)-1))  
    N=0
    E_int = E_int - (E_elc(2) + energy_HO(N,Wi(2,:),nvib(2)))
    DEALLOCATE(N)
    ALLOCATE(N(0:nvib(3)-1))
    N=0
    E_int = E_int - (E_elc(3) + energy_HO(N,Wi(3,:),nvib(3)))
    DEALLOCATE(N)
  END IF
  WRITE(*,*) "Total internal (cm-1)	:", E_int

  !The following statistical models were developed by T. Lam Nguyen

  !get density
  max_idx = FLOOR(E_int/dE)-1
  ALLOCATE(dens(0:1,0:max_idx))

  OPEN(unit=50,file='A.dens',access='sequential',status='old')
  OPEN(unit=51,file='B+.dens',access='sequential',status='old')
  DO i=0,max_idx
    READ(50,*) dummy,dum,dens(0,i),dum 
    READ(51,*) dummy,dum,dens(1,i),dum
  END DO
  CLOSE(unit=51,status='keep')
  CLOSE(unit=50,status='keep')

  !average energy
  A_conv = 0.0D0
  B_conv = 0.0D0
  A_enr = 0.0D0 
  B_enr = 0.0D0

  !convelute density of states
  DO i=0,max_idx
    A_conv = A_conv + dens(0,i)*dens(1,max_idx-i)*dE
    B_conv = B_conv + dens(1,i)*dens(0,max_idx-i)*dE
  END DO
  !sum energies
  DO i=0,max_idx
    A_enr = A_enr + i*dE*dens(0,i)*dens(1,max_idx-i)*dE
    B_enr = B_enr + i*dE*dens(1,i)*dens(0,max_idx-i)*dE
  END DO
  
  !Average energy analysis
  WRITE(*,*) "< E_A > (cm-1) 		:", A_enr/A_conv
  WRITE(*,*) "< E_B > (cm-1) 		:", B_enr/B_conv
  WRITE(*,*) " E_A %			:", A_enr/A_conv/E_int*100
  WRITE(*,*) " E_B %			:",B_enr/B_conv/E_int*100
!  WRITE(*,*) "A conv, B conv", A_conv,B_conv
!  WRITE(*,*) "from differences :", E_int-(A_enr/A_conv)

  !energy distribution


  DEALLOCATE(dens)
  CALL CPU_TIME(t2)
  WRITE(*,*) 
  WRITE(*,*) "stat_model finished in (s)",t2-t1

END SUBROUTINE stat_model 

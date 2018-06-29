!---------------------------------------------------------------------
!       !recursively aquire all the possible energy values
!---------------------------------------------------------------------
    !Values
    ! N         :       1D int, list of current quantum numbers
    ! Wi        :       1D real4, list of HO fundementals 
    ! ids       :       1D int, list of ids
    ! A         :       1D bool, hash table truth table
    ! B         :       2D int4, hash table keys
    ! C         :       1D bool, hash table vals
    ! Eu        :       real4, upper energy limit
    ! El        :       real4, lower energy limit
    ! Es      :       real4, sum of electronic energy 
    ! idx       :       int, which index in W we are working on now
    ! nall      :       int, number of values considered
    ! ngoo      :       int, number of energies accepted
    ! m         :       int, number of vibrational modes 

RECURSIVE SUBROUTINE recall_HO(N,Wi,ids,A,B,C,Eu,El,Es,idx,nall,ngood,l1,l2)
  USE myUtils

  IMPLICIT NONE

  INTERFACE
    REAL(KIND=4) FUNCTION energy_HO(N,W,m)
      INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: N
      REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: W
      INTEGER(KIND=4), INTENT(IN) ::m
    END FUNCTION energy_HO
    
    INTEGER(KIND=4) FUNCTION hash_1Dint4(A)
      INTEGER(KIND=4),DIMENSION(0:),INTENT(IN) :: A
    END FUNCTION hash_1Dint4
  END INTERFACE
  
  REAL(KIND=4), PARAMETER :: tol=1.0E-8

  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: B
  LOGICAL(KIND=4), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: A,C 
    
  INTEGER(KIND=4), DIMENSION(0:), INTENT(INOUT) :: N
  INTEGER(KIND=4),DIMENSION(0:), INTENT(IN) :: ids
  REAL(KIND=4), DIMENSION(0:), INTENT(IN) :: Wi
  INTEGER(KIND=4), INTENT(INOUT) :: nall,ngood,l1
  INTEGER(KIND=4), INTENT(IN) :: l2,idx
  REAL(KIND=4), INTENT(IN) :: Eu,El,Es

  INTEGER(KIND=4), DIMENSION(0:l2-1) :: newN
  INTEGER(KIND=4) :: i,j, hash
  REAL(KIND=4) :: Ev
  LOGICAL :: val

  !base cases
  hash =  hash_1Dint4(N)

  !1) we have seen this index before 
  CALL hash_qsearch_1Dint4_bool(A,B,C,N,val,hash,l1,val,hash_1Dint4)
  IF (val) THEN
    RETURN
  END IF
    
  !2) we are below zero anywhere
  IF (MINVAL(N) .LT. 0) THEN
    RETURN
  END IF

  !3) check if we lost energy from vpt breaking

  Ev = energy_HO(N,Wi,l2)

  !4) we are outside the upper energy bounds
  IF ( Eu .LT. (Es + Ev) ) THEN
    val = .FALSE.
    CALL hash_qinsert_1Dint4_bool(A,B,C,N,val,hash,l1,nall,hash_1Dint4)
    newN = N
    DO i=0,l2-1
      newN(i) = N(i)-1 !search down
      CALL recall_HO(newN,Wi,ids,A,B,C,Eu,El,Es,i,nall,ngood,l1,l2)
      newN(i) = newN(i) + 1 
    END DO
    RETURN

  !5) we are outside the lower energy bounds
  ELSE IF ( El .GT. (Es + Ev)) THEN
    val = .FALSE.
    CALL hash_qinsert_1Dint4_bool(A,B,C,N,val,hash,l1,nall,hash_1Dint4)
    newN = N
    DO i=0,l2-1
      newN(i) = N(i) + 1  !search up
      CALL recall_HO(newN,Wi,ids,A,B,C,Eu,El,Es,i,nall,ngood,l1,l2)
      newN(i) = N(i) - 1
    END DO
    RETURN

  ELSE

  !6) we're in range, explore whole area
    ngood = ngood + 1
    val = .TRUE.
    CALL hash_qinsert_1Dint4_bool(A,B,C,N,val,hash,l1,nall,hash_1Dint4)
    newN = N
    WRITE(2) newN
    DO i=0,l2-1
      newN(i) = N(i)+1  !search up
      CALL recall_HO(newN,Wi,ids,A,B,C,Eu,El,Es,i,nall,ngood,l1,l2)
      newN(i) = newN(i)-2 !search down
      CALL recall_HO(newN,Wi,ids,A,B,C,Eu,El,Es,i,nall,ngood,l1,l2)
      newN(i) = newN(i) + 1
    END DO
    RETURN

  END IF


END SUBROUTINE recall_HO

! Hash function, based off of FNV-1a algorithm created by
!   Fowler-Noll-Vo. Find information here:
!    http://www.isthe.com/chongo/tech/comp/fnv/
!
!    Note that I use a different prime offset, to comply
!      with fortran 90 SIGNED 32 bit ints

INTEGER(KIND=4) FUNCTION hash_1Dint4(A)
  IMPLICIT NONE

  INTEGER(KIND=4), DIMENSION(0:), INTENT(IN) :: A
  INTEGER(KIND=4) :: hash,i,val

  val = 2147483647

  DO i=0,SIZE(A)-1 
    val = IEOR(val,IBITS(A(i),1,7))
    val = val * 16777619
    val = IEOR(val,IBITS(A(i),9,7))
    val = val * 16777619
    val = IEOR(val,IBITS(A(i),17,7))
    val = val * 16777619
    val = IEOR(val,IBITS(A(i),25,7))
    val = val * 16777619
  END DO

  hash_1Dint4 = val

END FUNCTION hash_1Dint4 

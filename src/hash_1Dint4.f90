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
  BYTE :: ish 

  val = 2147483647

  ! I think that either IEOR or IBITS is not acting as intended.
  DO i=0,SIZE(A)-1 
    !val = IEOR(val,IBITS(A(i),0,7))*16777619
    !val = IEOR(val,IBITS(A(i),8,7))*16777619
    !val = IEOR(val,IBITS(A(i),16,7))*16777619
    !val = IEOR(val,IBITS(A(i),24,7))*16777619
!hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24)
    val = IEOR(val,IBITS(A(i),0,7))
    val = val+ISHFT(val,1)+ISHFT(val,4)+ISHFT(val,7)+ISHFT(val,8)+ISHFT(val,24)
    val = IEOR(val,IBITS(A(i),8,7))
    val = val+ISHFT(val,1)+ISHFT(val,4)+ISHFT(val,7)+ISHFT(val,8)+ISHFT(val,24)
    val = IEOR(val,IBITS(A(i),16,7))
    val = val+ISHFT(val,1)+ISHFT(val,4)+ISHFT(val,7)+ISHFT(val,8)+ISHFT(val,24)
    val = IEOR(val,IBITS(A(i),24,7))
    val = val+ISHFT(val,1)+ISHFT(val,4)+ISHFT(val,7)+ISHFT(val,8)+ISHFT(val,24)
  END DO

  val = ABS(val)
 
  hash_1Dint4 = val

END FUNCTION hash_1Dint4 

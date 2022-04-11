subroutine cramer(n,a,b,x)

  use determinant

  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(INOUT) :: x(n)
  REAL, ALLOCATABLE, DIMENSION(:,:) :: ak
  INTEGER :: k
  REAL :: d

  d = det(n, a)

  ALLOCATE(ak(n,n))
  do k = 1, n
    ak = a
    ak(:,k) = b
    x(k) = det(n, ak)/d
  end do

end subroutine cramer

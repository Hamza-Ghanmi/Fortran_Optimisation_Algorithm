subroutine Jacobi(n, a, b, x)

  use alg_tools
  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(INOUT) :: x(n)
  REAL :: AM(n,n), AN(n,n), IAM(n,n)
  REAL :: epsilon
  INTEGER :: i, j

! Jacobi Iterative Method algorithm
  do i = 1, n
    AM(i,i) = a(i,i)
  end do
  AN = AM - a
  IAM = inv(n,AM)
  epsilon = 1
  do while ( epsilon .gt. (10.0**(-3)) )
    x = MATMUL(IAM, MATMUL(AN, x) + b)
    epsilon = SQRT(SUM((MATMUL(a, x) - b)**2, DIM=1))
  end do

end subroutine Jacobi

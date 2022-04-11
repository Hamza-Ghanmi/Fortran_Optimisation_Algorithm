subroutine Relaxation(n, a, b, x)

  use alg_tools
  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(INOUT) :: x(n)
  REAL :: AM(n,n), AN(n,n), E(n,n), D(n,n), IAM(n,n)
  REAL :: epsilon, omega = 0.7
  INTEGER :: i, j

! Relaxation Iterative Method algorithm
  do i = 1, n
    do j = 1, i-1
      E(i,j) = -a(i,j)
    end do
    D(i,i) = a(i,i)
  end do
  AM = (1/omega) * D - E
  AN = AM - a
  IAM = inv(n,AM)
  epsilon = 1
  do while ( epsilon .gt. (10.0**(-3)) )
    x = MATMUL(IAM, MATMUL(AN, x) + b)
    epsilon = SQRT(SUM((MATMUL(a, x) - b)**2, DIM=1))
  end do


end subroutine Relaxation

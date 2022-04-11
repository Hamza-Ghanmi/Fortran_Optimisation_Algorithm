subroutine gdcs(n, a, b, x)

  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(INOUT) :: x(n)
  REAL :: alpha, beta,d(n), g1(n), g2(n) ! alpha depend on spectral value "in ]0, 2/rho(a)"
  REAL :: epsilon = 1
! Gradient conjugated step Method
d = b - MATMUL(a,x)
g2 = -d
do while ( epsilon .gt. (10.0**(-9)) )
  g1 = g2
  d = b - MATMUL(a,x)
  alpha = DOT_PRODUCT(g1, g1) / DOT_PRODUCT(MATMUL(a, d), d)
  x = x + alpha * d
  g2 = MATMUL(a, x) - b
  beta = DOT_PRODUCT(g2, g2) / DOT_PRODUCT(g1, g1)
  d = -g2 + beta * d
  epsilon = SQRT(SUM(d**2, DIM=1))
end do

end subroutine gdcs

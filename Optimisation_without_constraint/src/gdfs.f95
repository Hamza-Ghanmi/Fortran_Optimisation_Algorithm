subroutine gdfs(n, a, b, x)

  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(INOUT) :: x(n)
  REAL :: alpha = 0.1, d(n) ! alpha depend on spectral value "in ]0, 2/rho(a)"
  REAL :: epsilon = 1

! Gradient with fixed step Method

do while ( epsilon .gt. (10.0**(-9)) )
  d = b - MATMUL(a,x)
  x = x + alpha * d
  epsilon = SQRT(SUM(d**2, DIM=1))
end do

end subroutine gdfs

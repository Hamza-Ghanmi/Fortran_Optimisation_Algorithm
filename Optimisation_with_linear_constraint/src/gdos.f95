subroutine gdos(n, a, b, x)

  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(INOUT) :: x(n)
  REAL :: alpha, d(n), ad ! alpha depend on spectral value "in ]0, 2/rho(a)"
  REAL :: epsilon = 1

! Gradient with optimal step Method

do while ( epsilon .gt. (10.0**(-9)) )
  d = b - MATMUL(a,x)
  alpha = 0
  ad = DOT_PRODUCT(MATMUL(a, d), d)
  if (ad.NE.0.0) then
    alpha = DOT_PRODUCT(d, d) / ad
  end if
  x = x + alpha * d
  epsilon = SQRT(SUM(d**2, DIM=1))
end do

end subroutine gdos

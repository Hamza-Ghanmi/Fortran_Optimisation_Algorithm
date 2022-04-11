subroutine Cholesky(n, a, b, C)

  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(out) :: C(n,n)
  INTEGER :: i, j

! Factorisation Cholesky
C(1,1) = SQRT(a(1,1))
do j = 2, n
  C(j,1) = a(j,1) / C(1,1)
end do
do i = 2, n-1
  C(i,i) = SQRT(a(i,i) - sum(C(i,1:i-1) ** 2, dim = 1))
  do j = i+1, n
    C(j,i) = (a(j,i) - sum(C(i,1:i-1) * C(j,1:i-1), dim = 1)) / C(i,i)
  end do
end do
C(n,n) = SQRT(a(n,n) - sum(C(n,1:n-1) ** 2, dim = 1))

end subroutine Cholesky

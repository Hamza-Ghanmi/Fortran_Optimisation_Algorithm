subroutine LU(n, a, b, L, U)

  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(out) :: L(n,n), U(n,n)
  INTEGER :: i, j, k

! LU Factorisation
do i = 1, n
  L(i,i) = 1
end do

do k = 1, n
  do j = k, n
    U(k,j) = a(k,j) - sum(L(k,1:k-1) * U(1:k-1,j), dim = 1)
  end do
  do i = k+1, n
    L(i,k) = (a(i,k) - sum(L(i,1:k-1) * U(1:k-1,k), dim = 1)) / U(k,k)
  end do
end do



end subroutine LU

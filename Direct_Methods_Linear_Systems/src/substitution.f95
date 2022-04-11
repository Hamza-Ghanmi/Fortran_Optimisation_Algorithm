subroutine back_subs(n, a, b, x)

  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(INOUT) :: x(n)
  INTEGER :: k

! x-elements using back substitution
  x(n) = b(n)/a(n,n)
  do k = n-1, 1, -1
    x(k) = (b(k)-sum(a(k,k+1:n) * x(k+1:n), dim = 1))/a(k,k)
  end do

end subroutine back_subs

subroutine forward_subs(n, a, b, x)

  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a(n,n), b(n)
  REAL, INTENT(INOUT) :: x(n)
  INTEGER :: k

! x-elements using back substitution
  x(1) = b(1)/a(1,1)
  do k = 2, n
    x(k) = (b(k)-sum(a(k,1:k-1) * x(1:k-1), dim = 1))/a(k,k)
  end do

end subroutine forward_subs

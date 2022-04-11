subroutine Gauss_elim(n, a, b)

  implicit none

  INTEGER, INTENT(IN) :: n
  REAL, INTENT(INOUT) :: a(n,n), b(n)
  Real :: pivot_line(n+1)
  INTEGER :: i, j, k, index_pivot
  REAL :: p

! Gaussian elimination algorithm
  do k = 1, n-1
    index_pivot = MAXLOC(ABS(a(k:n,k)), DIM = 1) + k - 1

!    if ( a(k,k).eq.0 ) then
!      j = k
!      do while ( (a(j,k).eq.0 ).and.(j.lt.n) )
!        j = j + 1
!        index_pivot = j
!      end do

    pivot_line(1:n) = a(index_pivot,:)
    pivot_line(n+1) = b(index_pivot)
    a(index_pivot,1:n) = a(k,:)
    b(index_pivot) = b(k)
    a(k,:) = pivot_line(1:n)
    b(k) = pivot_line(n+1)
!    end if
    do i = k+1, n
      p = a(i,k)/a(k,k)
      a(i, k:n) = a(i,k:n) - p * a(k,k:n)
      b(i) = b(i) - p * b(k)
    end do
  end do

end subroutine Gauss_elim

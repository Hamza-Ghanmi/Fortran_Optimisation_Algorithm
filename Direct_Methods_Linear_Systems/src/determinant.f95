module determinant

  implicit none

contains

  recursive function det(n, a) result(d)

    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: a(n,n)
    INTEGER :: i
    REAL,ALLOCATABLE ,DIMENSION(:,:) :: aij
    REAL :: d

    d = 0

    if ( n.eq.1 ) then
      d = a(1,1)
    elseif ( n.gt.1 ) then
      do i = 1, n
        ALLOCATE(aij(n-1,n-1))
        aij(1:i-1,:) = a(1:i-1,2:)
        aij(i:n-1,:) = a(i+1:n,2:)
        d = d + (-1)**(i+1)*a(i,1)*det(n-1, aij)
        DEALLOCATE(aij)
      end do
    end if

  end function det

end module determinant

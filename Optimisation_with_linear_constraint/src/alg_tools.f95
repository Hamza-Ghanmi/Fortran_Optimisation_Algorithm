module alg_tools

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

  function inv(n, a) result(ia)

    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: a(n,n)
    REAL :: ia(n,n), com(n-1, n-1)
    integer :: i, j

    do i = 1, n
      do j = 1, n
        com(1:i-1,1:j-1) = a(1:i-1,1:j-1)
        com(i:n-1,1:j-1) = a(i+1:n,1:j-1)
        com(1:i-1,j:n-1) = a(1:i-1,j+1:n)
        com(i:n-1,j:n-1) = a(i+1:n,j+1:n)
        ia(i,j) = (-1) ** (i+j) * det(n-1,com)
      end do
    end do
    ia = TRANSPOSE(ia) / det(n, a)

  end function inv


end module alg_tools

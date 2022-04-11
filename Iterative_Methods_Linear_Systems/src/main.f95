program main

  implicit none

  INTEGER :: i, n = 4
  REAL, ALLOCATABLE, DIMENSION(:) :: xj, xr, b
  REAL, ALLOCATABLE, DIMENSION(:,:) :: a

  ALLOCATE(a(n,n))
  ALLOCATE(xj(n))
  ALLOCATE(xr(n))
  ALLOCATE(b(n))

  a = TRANSPOSE(reshape((/10.0, -1.0, 2.0, 0.0, -1.0, 11.0, -1.0, 3.0, 2.0, -1.0, 10.0, -1.0, 0.0, 3.0, -1.0, 8.0/), shape(a)))
  b = (/6.0, 25.0, -11.0, 15.0/)

  Call Jacobi(n,a,b,xj)
  Call Relaxation(n,a,b,xr)



  do i = 1, n
    print *, a(i,:), b(i), xj(i), xr(i)
  end do

  DEALLOCATE(a)
  DEALLOCATE(xj)
  DEALLOCATE(xr)
  DEALLOCATE(b)


end program main

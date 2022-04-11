program main

  implicit none

  INTEGER :: i, n = 4
  REAL, ALLOCATABLE, DIMENSION(:) :: xf, xo, xc, b
  REAL, ALLOCATABLE, DIMENSION(:,:) :: a

  ALLOCATE(a(n,n))
  ALLOCATE(xf(n))
  ALLOCATE(xo(n))
  ALLOCATE(xc(n))
  ALLOCATE(b(n))

  a = TRANSPOSE(reshape((/10.0, -1.0, 2.0, 0.0, -1.0, 11.0, -1.0, 3.0, 2.0, -1.0, 10.0, -1.0, 0.0, 3.0, -1.0, 8.0/), shape(a)))
  b = (/6.0, 25.0, -11.0, 15.0/)

  Call gdfs(n,a,b,xf)
  Call gdos(n,a,b,xo)
  Call gdcs(n,a,b,xc)


  do i = 1, n
    print *, a(i,:), b(i), xf(i), xo(i), xc(i)
  end do

  DEALLOCATE(a)
  DEALLOCATE(xf)
  DEALLOCATE(xo)
  DEALLOCATE(xc)
  DEALLOCATE(b)


end program main

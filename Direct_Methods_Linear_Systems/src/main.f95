program main

  implicit none

  INTEGER :: i, n = 4
  REAL, ALLOCATABLE, DIMENSION(:) :: x, y, b, bg, xc, xg, xlu, xcl, ycl
  REAL, ALLOCATABLE, DIMENSION(:,:) :: a, ag, U, L, C

  ALLOCATE(a(n,n))
  ALLOCATE(ag(n,n))
  ALLOCATE(x(n))
  ALLOCATE(xc(n))
  ALLOCATE(xg(n))
  ALLOCATE(xlu(n))
  ALLOCATE(xcl(n))
  ALLOCATE(ycl(n))
  ALLOCATE(y(n))
  ALLOCATE(b(n))
  ALLOCATE(bg(n))
  ALLOCATE(U(n,n))
  ALLOCATE(L(n,n))
  ALLOCATE(C(n,n))
  a = TRANSPOSE(reshape((/10.0, -1.0, 2.0, 0.0, -1.0, 11.0, -1.0, 3.0, 2.0, -1.0, 10.0, -1.0, 0.0, 3.0, -1.0, 8.0/), shape(a)))
  b = (/6.0, 25.0, -11.0, 15.0/)
  CALL LU(n, a, b, L, U)
  CALL forward_subs(n, L, b, y)
  CALL back_subs(n, U, y, xlu)
  CALL cramer(n, a, b, xc)
  ag = a
  bg = b
  CALL Gauss_elim(n, ag, bg)
  CALL back_subs(n, ag, bg, xg)
  Call Cholesky(n, a, b, C)
  do i = 1, n
    print *, C(i,:)
  end do
  Call forward_subs(n, C, b, ycl)
  CALL back_subs(n, TRANSPOSE(C), ycl, xcl)
  do i = 1, n
    print *, a(i,:), b(i), xc(i), xg(i), xlu(i), xcl(i)
  end do

  DEALLOCATE(a)
  DEALLOCATE(ag)
  DEALLOCATE(x)
  DEALLOCATE(xc)
  DEALLOCATE(xg)
  DEALLOCATE(xlu)
  DEALLOCATE(y)
  DEALLOCATE(ycl)
  DEALLOCATE(xcl)
  DEALLOCATE(b)
  DEALLOCATE(bg)
  DEALLOCATE(U)
  DEALLOCATE(L)
  DEALLOCATE(C)

end program main

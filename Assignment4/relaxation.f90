module relaxation

!implicit none
integer, parameter:: DOUBLE = kind(1.0d0)

contains 

subroutine jacobi1d(n , niter, omega, f, v, resid, error, info) !jacobi iteration

  integer(DOUBLE), intent(in):: n, niter 
  integer(DOUBLE), intent(out):: info
  real(DOUBLE), intent(in):: omega, f(n-1)
  real(DOUBLE), intent(inout):: v(n-1), resid(n-1), error(niter)

  real(DOUBLE):: dx, temp(n-1)
  integer(DOUBLE):: it, i

  !initializing the values
  dx = 1.0/n
  info = 1
  it = 1

  do i =1, niter
    error(i) = 0
  end do

  !checking for proper values
  if(omega>1 .or. omega<0) then
    info = -1
  endif

  if(n<2) then
    info = -1
  endif

  if(niter<1) then
    info = -1
  endif
  
  !start the jacobi iterations
  if(info>0) then
    do while(it<=niter)
      temp(1) = 0.5*(v(2) + f(1)*dx*dx)
      temp(n-1) = 0.5*(v(n-2)+f(n-1)*dx*dx)
      temp(2:n-2) = 0.5*(f(2:n-2)*dx*dx + v(1:n-3) + v(3:n-1))

      v(1:n-1) = (1-omega)*v(1:n-1) + omega*temp(1:n-1)


      !calculating residuals
      resid(1) = f(1) - (2*v(1) - v(2))/(dx*dx)
      resid(n-1) = f(n-1) - (2*v(n-1)-v(n-2))/(dx*dx)
      do i = 2, n-2
        resid(i) = f(i) - (2*v(i) - v(i-1) - v(i+1))/(dx*dx)
      end do

      !calculating error
      temp(1) = 0.5*(resid(2) + f(1)*dx*dx)
      temp(n-1) = 0.5*(resid(n-2)+f(n-1)*dx*dx)
      temp(2:n-2) = 0.5*(f(2:n-2)*dx*dx + resid(1:n-3) + resid(3:n-1))

      error(it) = 0
      do i = 1, n-1
        error(it) = error(it) + temp(i)*temp(i)
      end do
      error(it) = sqrt(error(it)/niter) !average root mean squared error
      it = it + 1
    end do
  endif
  return
end subroutine jacobi1d

subroutine gs1d(n, niter, f, v, resid, error, info)

  integer(DOUBLE), intent(in):: n, niter
  integer(DOUBLE), intent(out):: info
  real(DOUBLE), intent(in):: f(n-1)
  real(DOUBLE), intent(inout):: v(n-1), resid(n-1), error(niter)

  real(DOUBLE):: dx, temp(n-1)
  integer(DOUBLE):: it , i

  !initializing values
  dx = 1.0/n
  info = 1
  it = 1

  do i =1, niter
    error(i) = 0
  end do

  !checking for proper values
  if(n<2) then
    info = -1
  endif

  if(niter<1) then
    info = -1
  endif
  
  !starting red-black gauss sidel iterations
  if(info>0) then
    do while(it<=niter)
      v(1) = 0.5*(v(2) + f(1)*dx*dx)
      v(n-1) = 0.5*(v(n-2)+f(n-1)*dx*dx)
      v(3:n-2:2) = 0.5*(f(3:n-2:2)*dx*dx + v(2:n-3:2) + v(4:n-1:2))
      v(2:n-2:2) = 0.5*(f(2:n-2:2)*dx*dx + v(1:n-3:2) + v(3:n-1:2))


      !calculating residuals
      resid(1) = f(1) - (2*v(1) - v(2))/(dx*dx)
      resid(n-1) = f(n-1) - (2*v(n-1)-v(n-2))/(dx*dx)
      resid(2:n-2) = f(2:n-2) - (2*v(2:n-2) - v(1:n-3) - v(3:n-1))/(dx*dx)

      !calculating error
      temp(1) = 0.5*(resid(2) + f(1)*dx*dx)
      temp(n-1) = 0.5*(resid(n-2)+f(n-1)*dx*dx)
      temp(2:n-2) = 0.5*(f(2:n-2)*dx*dx + resid(1:n-3) + resid(3:n-1))

      error(it) = 0
      do i = 1, n-1
        error(it) = error(it) + temp(i)*temp(i)
      end do
      error(it) = sqrt(error(it)/n)
      it = it + 1
    end do
  end if
  return
end subroutine gs1d

end module relaxation

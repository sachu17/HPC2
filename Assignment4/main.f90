PROGRAM main

USE relaxation

!n grid size, niter number of iterations
!info flag for error
!omega constraint value
!v solution vector
!resid residual values
!error average error
!resultTable table containing error for different initial guesses

integer(DOUBLE), parameter:: n=64, niter=100
real(DOUBLE), parameter:: pi = 3.1416
integer(DOUBLE):: info, i
real(DOUBLE):: omega, f(n-1)

real(DOUBLE):: v(n-1), resid(n-1), error(niter)
real(DOUBLE):: resultTable(niter, 4)

!initializing the values
omega = 2.0/3.0

!v = sin(k*pi/n)
do i =1,n-1
  f(i) = 0
  v(i) = sin(i*pi/n)
end do
call jacobi1d(n, niter, omega, f, v, resid, error, info)
do i = 1, niter
  resultTable(i, 1) = i
  resultTable(i, 2) = error(i)
end do

!v = sin(3k*pi/n)
do i =1,n-1
  f(i) = 0
  v(i) = sin(3*i*pi/n)
end do
call jacobi1d(n, niter, omega, f, v, resid, error, info)
do i = 1, niter
  resultTable(i, 1) = i
  resultTable(i, 3) = error(i)
end do

!v = sin(6k*pi/n)
do i =1,n-1
  f(i) = 0
  v(i) = sin(6*i*pi/n)
end do
call jacobi1d(n, niter, omega, f, v, resid, error, info)
do i = 1, niter
  resultTable(i, 1) = i
  resultTable(i, 4) = error(i)
end do

!writing values to a file
open(unit = 6, file='resultsjacobi.txt')
do i =1, niter
  write(6, '(i5, 3es12.3)') i, resultTable(i,2), resultTable(i,3), resultTable(i,4)
end do
close(6)

!gauss sidel red black method

!v = sin(k*pi/n)
do i =1,n-1
  f(i) = 0
  v(i) = sin(i*pi/n)
end do
call gs1d(n, niter, f, v, resid, error, info)
do i = 1, niter
  resultTable(i, 1) = i
  resultTable(i, 2) = error(i)
end do

!v = sin(3k*pi/n)
do i =1,n-1
  f(i) = 0
  v(i) = sin(3*i*pi/n)
end do
call gs1d(n, niter, f, v, resid, error, info)
do i = 1, niter
  resultTable(i, 1) = i
  resultTable(i, 3) = error(i)
end do

!v = sin(6k*pi/n)
do i =1,n-1
  f(i) = 0
  v(i) = sin(6*i*pi/n)
end do
call gs1d(n, niter, f, v, resid, error, info)
do i = 1, niter
  resultTable(i, 1) = i
  resultTable(i, 4) = error(i)
end do

open(unit = 6, file='resultsred-black.txt')
do i =1, niter
  write(6, '(i5, 3es12.3)') i, resultTable(i,2), resultTable(i, 3), resultTable(i, 4)
end do
close(6)

if(info<0) then
 print*, 'Info status negative, some of the input parameters maybe improper'
end if

end PROGRAM main
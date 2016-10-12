PROGRAM main
USE precision

implicit NONE
integer, PARAMETER:: N=4097
integer:: Bmin=0, Bplu=0
real(DOUBLE):: ep=2.0**(-12)


real(DOUBLE), PARAMETER:: lockout = 100
integer, PARAMETER:: maxiter = 100


!initializing the basin grid
real(DOUBLE), DIMENSION(N):: X, Y
real(DOUBLE):: tempx, tempy
integer(DOUBLE), DIMENSION(N, N):: infFlag, infFlagminep, infFlagpluep
integer(DOUBLE):: i, j
integer:: temp

do while (ep>(2.0**(-22)))
 Bmin = 0
 Bplu = 0
 do i=1,N
  X(i)=((i-1.0)/N)*6 - 3
  Y(i)=((i-1.0)/N)*6 - 3
 enddo 
 do i=1,N
  do j=1,N
   infFlag(i,j) = 0
   infFlagminep(i,j) = 0
   infFlagpluep(i,j) = 0
  enddo
 enddo

 do i=1,N
  do j=1,N
   tempx = X(i)
   tempy = Y(j)
   CALL henon(tempx, tempy, lockout, maxiter, temp)
   infFlag(i,j) = temp
  enddo
 enddo

 do i=1,N
  X(i) = X(i) + ep
  Y(i) = Y(i) + ep
 enddo

 do i=1,N
  do j=1,N
   tempx = X(i)
   tempy = Y(j)
   CALL henon(tempx, tempy, lockout, maxiter, temp)
   infFlagpluep(i,j) = temp
  enddo
 enddo

 do i=1,N
  X(i) = X(i) - 2*ep
  Y(i) = Y(i) - 2*ep
 enddo

 do i=1,N
  do j=1,N
   tempx = X(i)
   tempy = Y(j)
   CALL henon(tempx, tempy, lockout, maxiter, temp)
   infFlagminep(i,j) = temp
  enddo
 enddo


 open(unit=2, file='results.txt')
 write(2, *) , 'Epsilon', 'Bminus count', 'Bplus count'

 do i = 1,N
  do j = 1,N
   if(infFlag(i,j)/=infFlagminep(i,j)) then
    Bmin=Bmin+1
   endif
   if(infFlag(i,j)/=infFlagpluep(i,j)) then
    Bplu=Bplu+1
   endif
  enddo
 enddo
 print*, Bmin, Bplu
 write(2, *), ep, Bmin, Bplu
 ep = ep/2.0_DOUBLE
enddo

END PROGRAM

!  PROGRAMMER
!   Your name, sshivak8@asu.edu.
!              Sachin Shivakumar

module substitution
USE precision
implicit none
contains

!Forward substitution method
subroutine forward_sub(a, n, b, info)
integer, intent(in):: n !dimension of matrix
real(DOUBLE), intent(in):: a(n,n) !coefficient matrix
real(DOUBLE), intent(inout):: b(n) !target vector
integer, intent(out):: info !error handling flag
integer:: i,j !dummy variables

!checking for zero diagonal elements
do i=1, n
 if(a(i,i)==0) then
  info = 1
  exit
 else
  info = 0
 endif
enddo

!proceed with the forward substitution method
if(info==0) then
 do j=1,n
  do i=1,j-1
   b(j) = b(j) - a(i,j)*b(i)
  enddo
  b(j) = b(j)/a(j,j)
 enddo
else
  continue
endif
return
end subroutine forward_sub



!Backward substitution method
subroutine backward_sub(a, n, b, info)
integer, intent(in):: n !dimension of matrix
real(DOUBLE), intent(in):: a(n,n) !coefficient matrix
real(DOUBLE), intent(inout):: b(n) !target vector
integer, intent(out):: info !error handling flag
integer:: i,j !dummy variables

!checking for zero diagonal elements
do i=1, n
 if(a(i,i)==0) then
  info = 1
  exit
 else
  info = 0
 endif
enddo

!proceed with the backward substitution method
if(info==0) then
 do i=n,1,-1
  do j=n,i-1,-1
   b(i) = b(i) - a(i,j)*b(j)
  enddo
  b(i) = b(i)/a(i,i)
 enddo
else
 continue
endif
return 
end subroutine backward_sub
end module substitution


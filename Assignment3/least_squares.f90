!  PROGRAMMER
!   Your name, sshivak8@asu.edu.
!              Sachin Shivakumar

PROGRAM least_squares

USE precision
USE, intrinsic:: iso_fortran_env

implicit none

integer, parameter:: MAXDATA = 1000 !assuming maximum 1000 xy pairs
real(SINGLE):: x(MAXDATA), y(MAXDATA) !intialize x and y vectors
integer(SINGLE) :: m, k, xflag, yflag, count, negative_log, i !dummy variables
character(80):: arg, filename !taking command line arguments

xflag = 0 !flag to check whether to take log x 
yflag = 0 !flag to check whether to take log y 
negative_log = 0 !flag to handle log of negative number error
count = command_argument_count() !checking number of command line inputs

!get all command line arguments and act accordingly
do k=1, count
 call get_command_argument(k, arg)
 if(arg(1:2)=='-L') then !checking for filename or logarithmic flag
  if(arg(3:3)=='x') then 
   xflag = 1 !check log flag for x
  else
   if(arg(3:3)=='y') then !check log flag for y
    yflag = 1
   else
    xflag = -1 !invalid log command arguments
    yflag = -1
   endif
  endif
 else
  filename = arg !store filename
 endif
enddo

call readfile(filename, x, m, y) !open the file and store the values

if(xflag==1) then !skip transformation if the log of negative number is encountered
 do i=1, m
  if (x(i)<=0) then
   negative_log = 1
   exit
  endif
 enddo
endif

if (negative_log==0) then !transformation of x if negative flag is false
 do i =1, m
  x(i) = log(x(i))
 enddo
else
	print*, 'negative number found, cannot take log and proceeding with normal fit'
endif

if(yflag==1) then !skip transformation if the log of negative number is encountered
 do i=1, m
  if (y(i)<=0) then
   negative_log = 1
   exit
  endif
 enddo
endif

if (negative_log==0) then !transformation of y if negative flag is false
 do i =1, m
  y(i) = log(y(i))
 enddo
else
 print*, 'negative number found, cannot take log and proceeding with normal fit'
endif


call computells(x(m), m, y(m))

print*, y(1), y(2)

 
contains
 subroutine readfile(filename, x, m, y) !reading the file
  USE precision
  character(80), intent(in):: filename
  real(SINGLE), intent(inout):: x(m)
  real(SINGLE), intent(inout):: y(m)
  integer(SINGLE), intent(inout):: m
  integer(SINGLE):: MAXDATA

  MAXDATA = 1000

  open(unit=2, file = filename)
  do i=1, MAXDATA !counting number of the rows
   read(2, *, end=100)
  enddo
  100 m = MAXDATA - 1
  close(2)
  open(unit=2, file = filename)
  do i=1, m ! storing the values pairs in x and y
   read(2, *) x(i), y(i)
  enddo
  return
 end subroutine readfile
!!!!!
!!!!!
 subroutine computells(x, m, y)
  USE precision
  real(SINGLE), intent(in):: x(m)
  real(SINGLE), intent(in):: y(m)
  integer(SINGLE), intent(in):: m
  real(SINGLE):: X1(m,2)
  real(SINGLE):: wsize
  integer:: i

  do i=1,m
   X1(i,1) = x(m)
   X1(i,2) = 1
  enddo
  CALL sgels('N', m, 2, 1, X1, m, y, m , wsize,-1) !set worksize
  CALL sgels('N', m, 2, 1, X1, m, y, m , wsize) !solve for the values
  return
 end subroutine computells
!!!!!
!!!!!
end PROGRAM least_squares
!!!!!
!!!!!
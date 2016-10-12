    subroutine henon(x,y,lockout,maxiter,isin)
!  HENON - Sample implementation of a Henon subroutine that could be included
!  in your HENONDIM module.  This routine iterates a *single* initial condition
!  up to MAXITER times.  If |X| > LOCKOUT within MAXITER steps, then we
!  set ISIN to 1 to indicate that the initial condition lies in the basin
!  of attraction of infinity.  Otherwise, we set ISIN to 0.
!  Arguments:
!  X, Y :=: on entry, the initial condition; on return, the orbit after
!    MAXITER steps if ISIN is set to 0.
!  LOCKOUT : the lockout condition; we assume that the orbit diverges to
!   infinity if the absolute value of X exceeds LOCKOUT.
!  MAXITER : the maximum number of steps to iterate the orbit.
!  ISIN := set to 1 if the initial condition is deemed to be in the basin
!   of infinity and to 0 otherwise.
!  This routine assumes that the PRECISION module is in scope.
!
    USE precision
    real(DOUBLE), intent(inout):: x, y
    real(DOUBLE), intent(in):: lockout
    integer, intent(in):: maxiter
    integer, intent(out):: isin
    real(DOUBLE):: temp
!   
!  Local variables.  I have set A and B as local fixed parameters, but
!  you could make them module variables and include a method for changing
!  them from these defaults, particularly if you are intereted in the
!  bifurcation behavior as, say, A changes.  The _DOUBLE suffix is
!  essential to get parameters accurate to the full double precision.
!
    real(DOUBLE), parameter:: A=2.12_DOUBLE, B=-0.3_DOUBLE
    real(DOUBLE):: xnew
    integer:: j
    do j=1,maxiter
     xnew = A - x**2 + B*y
     y = x
     x = xnew
     temp = ABS(x)
     if(temp>lockout) exit
    enddo
!
    if(j <= maxiter) then  ! the point is diverging to infinity
     isin = 1
    else
     isin = 0
    endif
    return
    end subroutine henon
